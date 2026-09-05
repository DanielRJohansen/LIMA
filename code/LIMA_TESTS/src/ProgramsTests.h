#include "../LIMA/include/BuildMembrane.h"
#include "MoleculeUtils.h"
#include "Programs.h"
#include "TestUtils.h"

#include <cctype>
#include <set>

namespace ProgramsTests {
	namespace {
		std::set<int> ReadPositionRestraintAtoms(const fs::path& path) {
			std::ifstream input(path);
			if (!input) throw std::runtime_error(std::format("Could not open position restraints {}", path.string()));
			std::set<int> result;
			for (std::string line; std::getline(input, line);) {
				if (const auto comment = line.find(';'); comment != std::string::npos) line.erase(comment);
				const auto fields = StringUtils::SplitWords(line);
				if (fields.size() == 5 && std::isdigit(static_cast<unsigned char>(fields[0][0]))) result.insert(std::stoi(fields[0]));
			}
			return result;
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

		std::size_t CountOccurrences(std::string_view text, std::string_view value) {
			std::size_t count = 0;
			for (std::size_t position = 0; (position = text.find(value, position)) != std::string_view::npos;
				position += value.size()) ++count;
			return count;
		}
	}

	LimaUnittestResult TestBuildMembrane(EnvMode envmode) {
		const fs::path workDir = HeavyTestsDir() / "etc";

		std::vector<std::string> args = {
			"lima", "buildMembrane",
			"-lipids", "DPPC", "60", "DOPC", "40",
			"-centerz", "3.0",
			"-boxsize", "12.0", "10.0","8.0",
			"-emtol", "50",
			"-working_dir", workDir.string().c_str(),
			"-d"
		};

		std::vector<char*> argv;
		for (auto& arg : args)
			argv.push_back(arg.data());

		int ret = buildMembrane(args.size(), argv.data());

		return LimaUnittestResult{ ret == 0, "Success", envmode == Full };
	}

	LimaUnittestResult TestPdb2Gmx_pdbfile(EnvMode envmode) {
		const fs::path directory = AutomatedTestsDir() / "pdb2gmx";
		const fs::path generatedGro = directory / "generated.gro";
		const fs::path generatedTop = directory / "generated.top";
		const fs::path generatedPosre = directory / "generated_posre.itp";
		TryDeleteFile(generatedGro);
		TryDeleteFile(generatedTop);
		TryDeleteFile(generatedPosre);

		Programs::pdb2gmx(directory / "6lzm.pdb", "generated");
		LimaUnittestResult gResult = CompareGroFiles(GroFile{ generatedGro }, GroFile{ directory / "conf_ref.gro" }, envmode);
		if (!gResult.success) 
			return LimaUnittestResult{ false, gResult.error_description, envmode == Full };
		
		LimaUnittestResult tResult = CompareTopologyFiles(generatedTop, directory / "topol_ref.top", envmode);	
		if (!tResult.success) 
			return LimaUnittestResult{ false, tResult.error_description, envmode == Full };

		ASSERT(ReadPositionRestraintAtoms(generatedPosre) == ReadPositionRestraintAtoms(directory / "posre_ref.itp"),
			"Position-restraint atom sets differ");

		Programs::pdb2gmx(directory / "6lzm.pdb", "generated_spce", Programs::WaterModel::Spce);
		std::ifstream spceTopology(directory / "generated_spce.top");
		const std::string spceContents(std::istreambuf_iterator<char>{ spceTopology }, {});
		ASSERT(spceContents.contains("#include \"charmm27.ff/spce.itp\""),
			"Selected SPC/E water topology was not included");
		return LimaUnittestResult{ true, "Success", envmode == Full };
	}

	LimaUnittestResult TestCif2Gmx_ciffile(EnvMode envmode) {
		const fs::path directory = AutomatedTestsDir() / "pdb2gmx";
		const fs::path generatedGro = directory / "generated_cif.gro";
		const fs::path generatedTop = directory / "generated_cif.top";
		const fs::path generatedPosre = directory / "generated_cif_posre.itp";
		TryDeleteFile(generatedGro);
		TryDeleteFile(generatedTop);
		TryDeleteFile(generatedPosre);

		const auto conversion = Programs::cif2gmx(
			directory / "7LZM.cif", "generated_cif", Programs::WaterModel::Tip3p, directory);
		ASSERT(conversion.gro == generatedGro && conversion.topology == generatedTop
			&& conversion.positionRestraints == generatedPosre, "cif2gmx returned incorrect output paths");
		LimaUnittestResult gResult = CompareGroFiles(
			GroFile{ generatedGro }, GroFile{ directory / "conf_ref.gro" }, envmode, 0.75f, 0.05f, 0.075f);
		if (!gResult.success)
			return LimaUnittestResult{ false, gResult.error_description, envmode == Full };

		LimaUnittestResult tResult = CompareTopologyFiles(generatedTop, directory / "topol_ref.top", envmode);
		if (!tResult.success)
			return LimaUnittestResult{ false, tResult.error_description, envmode == Full };

		ASSERT(ReadPositionRestraintAtoms(generatedPosre) == ReadPositionRestraintAtoms(directory / "posre_ref.itp"),
			"CIF position-restraint atom sets differ");
		return LimaUnittestResult{ true, "Success", envmode == Full };
	}

	LimaUnittestResult TestPdb2Gmx_multichain(EnvMode envmode) {
		const fs::path directory = AutomatedTestsDir() / "pdb2gmx";
		const fs::path input = directory / "generated_multichain.pdb";
		const fs::path gro = directory / "generated_multichain.gro";
		const fs::path top = directory / "generated_multichain.top";
		const fs::path primaryPosre = directory / "generated_multichain_posre.itp";
		const fs::path secondPosre = directory / "generated_multichain_posre_Protein_chain_B.itp";
		for (const auto& path : { input, gro, top, primaryPosre, secondPosre }) TryDeleteFile(path);

		WriteTwoChainPdb(directory / "6lzm.pdb", input);
		const auto conversion = Programs::pdb2gmx(input, "generated_multichain");
		ASSERT(conversion.additionalPositionRestraints.size() == 1
			&& conversion.additionalPositionRestraints.front() == secondPosre
			&& fs::is_regular_file(primaryPosre) && fs::is_regular_file(secondPosre),
			"Multi-chain position-restraint files were not generated correctly");

		std::ifstream topology(top);
		const std::string contents(std::istreambuf_iterator<char>{ topology }, {});
		ASSERT(CountOccurrences(contents, "[ moleculetype ]") == 2,
			"Multi-chain topology does not contain one molecule type per chain");
		ASSERT(contents.contains("Protein_chain_A") && contents.contains("Protein_chain_B"),
			"Multi-chain topology is missing chain molecule names");

		const TopologyFile parsedTopology{ top };
		GroFile parsedCoordinates{ gro };
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
