#include "MoleculeUtils.h"
#include "Programs.h"
#include "TestUtils.h"

#include <cctype>
#include <set>

namespace ProgramsTests {
	namespace {
		std::set<int> ReadPositionRestraintAtoms(const GenericItpFile& file) {
			std::set<int> result;
			for (const auto& line : file.GetSection(position_restraints)) {
				const auto fields = StringUtils::SplitWords(line);
				if (fields.size() == 5 && std::isdigit(static_cast<unsigned char>(fields[0][0]))) result.insert(std::stoi(fields[0]));
			}
			return result;
		}

		std::set<int> ReadPositionRestraintAtoms(const fs::path& path) {
			return ReadPositionRestraintAtoms(GenericItpFile{ path });
		}

		void WriteTwoChainPdb(const fs::path& sourcePath, const fs::path& outputPath) {
			std::ifstream source(sourcePath);
			std::ofstream output(outputPath);
			if (!source || !output) throw std::runtime_error("Could not create multi-chain PDB test input");
			for (std::string line; std::getline(source, line);) {
				if (line.starts_with("COMPND") || line.starts_with("CRYST1")) {
					output << line << '\n';
				}
				else if (line.starts_with("ATOM  ") && line.size() >= 26) {
					const int residue = std::stoi(line.substr(22, 4));
					if (residue == 1 || residue == 2) {
						line[21] = residue == 1 ? 'A' : 'B';
						output << line << '\n';
					}
				}
			}
		}

	}

	LimaUnittestResult TestBuildMembrane(EnvMode envmode) {
		const fs::path workDir = HeavyTestsDir() / "etc";
		Lipids::Selection lipids{
			Lipids::Select{ "DPPC", workDir, 60.0 },
			Lipids::Select{ "DOPC", workDir, 40.0 }
		};
		GroFile grofile;
		grofile.box_size = Float3{ 12.f, 10.f, 8.f };
		grofile.title = "Membrane";
		TopologyFile topfile;
		topfile.SetSystem("Membrane");
		SimulationBuilder::CreateMembrane(grofile, topfile, lipids, 3.f);
		Programs::EnergyMinimize(grofile, topfile, true, workDir, envmode, true, 50.f);

		return LimaUnittestResult{ true, "Success", envmode == Full };
	}

	LimaUnittestResult TestToGmx_pdbfile(EnvMode envmode) {
		const fs::path directory = AutomatedTestsDir() / "pdb2gmx";
		const auto conversion = Programs::ToGmx(directory / "6lzm.pdb");
		LimaUnittestResult gResult = CompareGroFiles(conversion.grofile, GroFile{ directory / "conf_ref.gro" }, envmode);
		if (!gResult.success) 
			return LimaUnittestResult{ false, gResult.error_description, envmode == Full };
		
		LimaUnittestResult tResult = CompareTopologyFiles(conversion.topology, TopologyFile{ directory / "topol_ref.top" }, envmode);
		if (!tResult.success) 
			return LimaUnittestResult{ false, tResult.error_description, envmode == Full };

		ASSERT(conversion.positionRestraints.size() == 1
			&& ReadPositionRestraintAtoms(conversion.positionRestraints.front()) == ReadPositionRestraintAtoms(directory / "posre_ref.itp"),
			"Position-restraint atom sets differ");

		const auto spceConversion = Programs::ToGmx(directory / "6lzm.pdb", Programs::WaterModel::Spce);
		ASSERT(std::ranges::contains(spceConversion.topology.otherIncludes, fs::path{ "charmm27.ff/spce.itp" }),
			"Selected SPC/E water topology was not included");
		return LimaUnittestResult{ true, "Success", envmode == Full };
	}

	LimaUnittestResult TestToGmx_ciffile(EnvMode envmode) {
		const fs::path directory = AutomatedTestsDir() / "pdb2gmx";
		const auto conversion = Programs::ToGmx(directory / "7LZM.cif", Programs::WaterModel::Tip3p);
		LimaUnittestResult gResult = CompareGroFiles(
			conversion.grofile, GroFile{ directory / "conf_ref.gro" }, envmode, 0.75f, 0.05f, 0.075f);
		if (!gResult.success)
			return LimaUnittestResult{ false, gResult.error_description, envmode == Full };

		LimaUnittestResult tResult = CompareTopologyFiles(conversion.topology, TopologyFile{ directory / "topol_ref.top" }, envmode);
		if (!tResult.success)
			return LimaUnittestResult{ false, tResult.error_description, envmode == Full };

		ASSERT(conversion.positionRestraints.size() == 1
			&& ReadPositionRestraintAtoms(conversion.positionRestraints.front()) == ReadPositionRestraintAtoms(directory / "posre_ref.itp"),
			"CIF position-restraint atom sets differ");
		return LimaUnittestResult{ true, "Success", envmode == Full };
	}

	LimaUnittestResult TestToGmx_multichain(EnvMode envmode) {
		const fs::path directory = AutomatedTestsDir() / "pdb2gmx";
		const fs::path input = directory / "generated_multichain.pdb";
		TryDeleteFile(input);

		WriteTwoChainPdb(directory / "6lzm.pdb", input);
		auto conversion = Programs::ToGmx(input);
		ASSERT(conversion.positionRestraints.size() == 2
			&& !conversion.positionRestraints[0].GetSection(position_restraints).empty()
			&& !conversion.positionRestraints[1].GetSection(position_restraints).empty(),
			"Multi-chain position restraints were not generated correctly");

		ASSERT(conversion.topology.moleculetypes.size() == 2,
			"Multi-chain topology does not contain one molecule type per chain");
		ASSERT(conversion.topology.moleculetypes.contains("Protein_chain_A")
			&& conversion.topology.moleculetypes.contains("Protein_chain_B"),
			"Multi-chain topology is missing chain molecule names");

		const TopologyFile& parsedTopology = conversion.topology;
		GroFile parsedCoordinates = conversion.grofile;
		std::size_t topologyAtomCount = 0;
		for (const auto& molecule : parsedTopology.GetSystem().molecules) {
			topologyAtomCount += molecule.moleculetype->atoms.size();
		}
		ASSERT(parsedTopology.GetSystem().molecules.size() == 2
			&& topologyAtomCount == parsedCoordinates.atoms.size(),
			"Parsed multi-chain topology does not map one-to-one to GRO atoms");

		const auto& molecules = parsedTopology.GetSystem().molecules;
		const std::size_t secondMoleculeOffset = molecules.front().moleculetype->atoms.size();
		const auto& secondMolecule = *molecules[1].moleculetype;
		ASSERT(!secondMolecule.singlebonds.empty(), "Multi-chain test molecule contains no bond to fragment");
		const auto fragmentedBond = secondMolecule.singlebonds.front().ids;
		const Float3 untouchedFirstChainAtom = parsedCoordinates.atoms.front().position;
		parsedCoordinates.atoms[secondMoleculeOffset + fragmentedBond[1]].position.x += parsedCoordinates.box_size.x;

		MoleculeUtils::MakeMoleculeWholeAfterPBCFragmentation(parsedCoordinates, parsedTopology);
		ASSERT(parsedCoordinates.atoms.front().position == untouchedFirstChainAtom,
			"Making the second molecule whole modified the first molecule");
		const float reconstructedBondDx = std::abs(
			parsedCoordinates.atoms[secondMoleculeOffset + fragmentedBond[0]].position.x
			- parsedCoordinates.atoms[secondMoleculeOffset + fragmentedBond[1]].position.x);
		ASSERT(reconstructedBondDx < parsedCoordinates.box_size.x / 2.f,
			"System-level molecule reconstruction did not unwrap the fragmented chain");

		constexpr float padding = 0.5f;
		MoleculeUtils::FitMoleculeInBox(parsedCoordinates, padding);
		for (const auto& atom : parsedCoordinates.atoms) {
			ASSERT(atom.position.x >= padding - 1e-5f && atom.position.x <= parsedCoordinates.box_size.x - padding + 1e-5f
				&& atom.position.y >= padding - 1e-5f && atom.position.y <= parsedCoordinates.box_size.y - padding + 1e-5f
				&& atom.position.z >= padding - 1e-5f && atom.position.z <= parsedCoordinates.box_size.z - padding + 1e-5f,
				"Fitted render box does not encompass every molecule atom");
		}
		return LimaUnittestResult{ true, "Success", envmode == Full };
	}
}
