#include "Tests.h"
#include "MoleculeUtils.h"
#include "Programs.h"
#include "TestUtils.h"

#include <cctype>
#include <set>

namespace ProgramsTests {
	using namespace TestUtils;
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

		struct AtomsAndCharge { std::size_t atoms = 0; double charge = 0.0; };
		AtomsAndCharge CountAtomsAndCharge(const Programs::GmxConversionResult& conversion) {
			AtomsAndCharge result;
			for (const auto& molecule : conversion.topology.GetSystem().Instances())
				for (const auto& atom : molecule.moleculetype->atoms) {
					++result.atoms;
					result.charge += atom.charge.value_or(0.f);
				}
			return result;
		}

		// 6lzm with CYS 97's SG moved 0.2 nm from CYS 54's, so the two cysteines form a disulfide
		void WriteDisulfidePdb(const fs::path& sourcePath, const fs::path& outputPath) {
			std::ifstream source(sourcePath);
			std::ofstream output(outputPath);
			if (!source || !output) throw std::runtime_error("Could not create disulfide PDB test input");
			std::string sg54;
			std::vector<std::string> lines;
			for (std::string line; std::getline(source, line);) {
				if (line.starts_with("ATOM  ") && line.substr(12, 4) == " SG " && std::stoi(line.substr(22, 4)) == 54) sg54 = line;
				lines.push_back(line);
			}
			const double x = std::stod(sg54.substr(30, 8)) + 2.0, y = std::stod(sg54.substr(38, 8)), z = std::stod(sg54.substr(46, 8));
			for (auto& line : lines) {
				if (line.starts_with("ATOM  ") && line.substr(12, 4) == " SG " && std::stoi(line.substr(22, 4)) == 97) {
					char coordinates[32];
					std::snprintf(coordinates, sizeof(coordinates), "%8.3f%8.3f%8.3f", x, y, z);
					line.replace(30, 24, coordinates);
				}
				output << line << '\n';
			}
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
		auto job = Programs::MakeMembraneJob(workDir, Lipids::Composition{ { { "DPPC", 60.0 }, { "DOPC", 40.0 } }, workDir },
			Float3{ 12.f, 10.f, 8.f }, MembraneGeometry::Plane{ 3.f }, 0, SimParams::BasicEMSimParams(50.f), envmode);
		job.grofile->Get().title = "Membrane";
		job.topfile->Get().SetSystem("Membrane");
		Environment::Get().Submit(std::move(job)).Get();

		return LimaUnittestResult{ true, "Success", envmode == Full };
	}

	TestRoutine TestToGmx_pdbfile(Environment&, EnvMode envmode) {
		const fs::path directory = AutomatedTestsDir() / "pdb2gmx";
		const auto conversion = Programs::ToGmx(directory / "6lzm.pdb");
		LimaUnittestResult gResult = CompareGroFiles(conversion.grofile, GroFile{ directory / "conf_ref.gro" }, envmode);
		if (!gResult.success) 
			co_return LimaUnittestResult::Fail(gResult.error_description, envmode == Full);
		
		LimaUnittestResult tResult = CompareTopologyFiles(conversion.topology, TopologyFile{ directory / "topol_ref.top" }, envmode);
		if (!tResult.success) 
			co_return LimaUnittestResult::Fail(tResult.error_description, envmode == Full);

		TEST_ASSERT(conversion.positionRestraints.size() == 1
			&& ReadPositionRestraintAtoms(conversion.positionRestraints.front()) == ReadPositionRestraintAtoms(directory / "posre_ref.itp"),
			"Position-restraint atom sets differ");

		const auto spceConversion = Programs::ToGmx(directory / "6lzm.pdb", Programs::WaterModel::Spce);
		TEST_ASSERT(std::ranges::contains(spceConversion.topology.otherIncludes, fs::path{ "charmm27.ff/spce.itp" }, [](const std::string& include) { return fs::path{ include }; }),
			"Selected SPC/E water topology was not included");
		co_return LimaUnittestResult{ true, "Success", envmode == Full };
	}

	TestRoutine TestToGmx_ciffile(Environment&, EnvMode envmode) {
		const fs::path directory = AutomatedTestsDir() / "pdb2gmx";
		const auto conversion = Programs::ToGmx(directory / "7LZM.cif", Programs::WaterModel::Tip3p);
		LimaUnittestResult gResult = CompareGroFiles(
			conversion.grofile, GroFile{ directory / "conf_ref.gro" }, envmode, 0.75f, 0.05f, 0.075f);
		if (!gResult.success)
			co_return LimaUnittestResult::Fail(gResult.error_description, envmode == Full);

		LimaUnittestResult tResult = CompareTopologyFiles(conversion.topology, TopologyFile{ directory / "topol_ref.top" }, envmode);
		if (!tResult.success)
			co_return LimaUnittestResult::Fail(tResult.error_description, envmode == Full);

		TEST_ASSERT(conversion.positionRestraints.size() == 1
			&& ReadPositionRestraintAtoms(conversion.positionRestraints.front()) == ReadPositionRestraintAtoms(directory / "posre_ref.itp"),
			"CIF position-restraint atom sets differ");
		co_return LimaUnittestResult{ true, "Success", envmode == Full };
	}

	// The togmx options change exactly the atoms and charges they should, and disulfides are bonded
	TestRoutine TestToGmx_options(Environment&, EnvMode envmode) {
		const fs::path directory = AutomatedTestsDir() / "pdb2gmx";
		const fs::path input = directory / "6lzm.pdb";
		const auto convert = [&](auto configure) {
			Programs::ToGmxOptions options;
			configure(options);
			return Programs::ToGmx(input, options);
		};
		using O = Programs::ToGmxOptions;
		const auto base = CountAtomsAndCharge(Programs::ToGmx(input));
		const auto check = [&](std::string_view what, auto configure, int atomDelta, int chargeDelta) {
			const auto variant = CountAtomsAndCharge(convert(configure));
			const bool ok = static_cast<long long>(variant.atoms) - static_cast<long long>(base.atoms) == atomDelta
				&& std::abs(variant.charge - base.charge - chargeDelta) < 1e-3;
			return ok ? std::string{} : Lima::Format("{}: {} atoms, charge {:.3f}; expected {:+d} atoms, {:+d} charge vs {} / {:.3f}",
				what, variant.atoms, variant.charge, atomDelta, chargeDelta, base.atoms, base.charge);
		};
		for (const std::string& error : {
			check("--ter-n neutral", [](O& o) { o.nTerminus = O::NTerminus::Neutral; }, -1, -1),
			check("--ter-c neutral", [](O& o) { o.cTerminus = O::CTerminus::Neutral; }, +1, +1),
			check("--ter-c amide", [](O& o) { o.cTerminus = O::CTerminus::Amide; }, +2, +1),
			check("--ter-c methylamide", [](O& o) { o.cTerminus = O::CTerminus::Methylamide; }, +5, +1),
			check("--his hsp", [](O& o) { o.his = O::His::Hsp; }, +1, +1),
			check("--his hse", [](O& o) { o.his = O::His::Hse; }, 0, 0),
			check("--ignore-h", [](O& o) { o.ignoreHydrogens = true; }, 0, 0) }) {
			TEST_ASSERT(error.empty(), error);
		}

		const fs::path disulfideInput = directory / "generated_disulfide.pdb";
		WriteDisulfidePdb(input, disulfideInput);
		const auto bonded = Programs::ToGmx(disulfideInput);
		const auto& molecule = *bonded.topology.GetSystem().molecules.front().moleculetype;
		int sg54 = -1, sg97 = -1;
		for (const auto& atom : molecule.atoms) {
			if (atom.atomname == "SG" && atom.resnr == 54) sg54 = atom.id;
			if (atom.atomname == "SG" && atom.resnr == 97) sg97 = atom.id;
		}
		const bool hasBond = std::ranges::any_of(molecule.singlebonds, [&](const auto& bond) {
			return (bond.ids[0] == sg54 && bond.ids[1] == sg97) || (bond.ids[0] == sg97 && bond.ids[1] == sg54);
		});
		TEST_ASSERT(sg54 >= 0 && sg97 >= 0 && hasBond, "No SG-SG bond for the disulfide CYS 54 - CYS 97");
		TEST_ASSERT(CountAtomsAndCharge(bonded).atoms == base.atoms - 2, "The disulfide cysteines kept their HG hydrogens");
		TEST_ASSERT(std::ranges::any_of(bonded.notes, [](const std::string& note) { return note.starts_with("Disulfides: "); }),
			"The disulfide was not reported");
		const auto unbonded = convert([](O& o) { o.disulfides = false; });
		TEST_ASSERT(CountAtomsAndCharge(unbonded).atoms == base.atoms, "--no-disulfides changed the plain structure");
		TryDeleteFile(disulfideInput);

		TEST_ASSERT(Programs::ToGmx(directory / "7LZM.cif").notes.size() >= 1, "The CIF's HETATM waters were not reported");
		co_return LimaUnittestResult{ true, "Success", envmode == Full };
	}

	// AMBER99SB: Villin's coordinates, written as a PDB without hydrogens, give pdb2gmx's topology (ForceComparisons/VillinAmber).
	// The reference's charges are zeroed, so they are checked by the net charge instead
	TestRoutine TestToGmx_amber(Environment&, EnvMode envmode) {
		const fs::path referenceDir = AutomatedTestsDir() / "ForceComparisons" / "VillinAmber";
		const fs::path input = AutomatedTestsDir() / "pdb2gmx" / "generated_villin.pdb";
		{
			const GroFile gro{ referenceDir / "conf.gro" };
			std::ofstream pdb(input);
			for (std::size_t i = 0; i < gro.atoms.size(); ++i) {
				const auto& atom = gro.atoms[i];
				const std::string name{ atom.atomName.View() };
				char line[96];
				std::snprintf(line, sizeof(line), "ATOM  %5d %-4s %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00",
					static_cast<int>(i + 1), (name.size() < 4 ? " " + name : name).c_str(), std::string{ atom.residueName.View() }.c_str(), atom.residue_number,
					atom.position.x * 10.f, atom.position.y * 10.f, atom.position.z * 10.f);
				pdb << line << '\n';
			}
		}
		Programs::ToGmxOptions options;
		options.forceField = Programs::ToGmxOptions::ForceField::Amber99sb;
		options.ignoreHydrogens = true;
		const auto conversion = Programs::ToGmx(input, options);
		TryDeleteFile(input);

		const TopologyFile reference{ referenceDir / "topol.top" };
		const auto& expected = *reference.GetSystem().molecules.front().moleculetype;
		const auto& actual = *conversion.topology.GetSystem().molecules.front().moleculetype;
		TEST_ASSERT(actual.atoms.size() == expected.atoms.size(),
			Lima::Format("{} atoms, pdb2gmx {}", actual.atoms.size(), expected.atoms.size()));
		double netCharge = 0.0;
		for (std::size_t i = 0; i < actual.atoms.size(); ++i) {
			const auto& a = actual.atoms[i];
			const auto& e = expected.atoms[i];
			TEST_ASSERT(a.atomname == e.atomname && a.type == e.type && a.resnr == e.resnr,
				Lima::Format("Atom {}: {} {} {}, pdb2gmx {} {} {}", i + 1, a.resnr, a.atomname, a.type, e.resnr, e.atomname, e.type));
			netCharge += a.charge.value_or(0.f);
		}
		TEST_ASSERT(std::abs(netCharge - 2.0) < 1e-3, Lima::Format("Net charge {:.4f}, expected +2", netCharge));

		const auto sameIds = [](const auto& left, const auto& right) {
			auto ids = [](const auto& entries) {
				std::vector<std::vector<int>> result;
				for (const auto& entry : entries) {
					std::vector<int> v(entry.ids.begin(), entry.ids.end());
					if (v.back() < v.front()) std::ranges::reverse(v);
					result.push_back(std::move(v));
				}
				std::ranges::sort(result);
				return result;
			};
			return ids(left) == ids(right);
		};
		TEST_ASSERT(sameIds(actual.singlebonds, expected.singlebonds), "Bonds differ from pdb2gmx");
		TEST_ASSERT(sameIds(actual.pairbonds, expected.pairbonds), "Pairs differ from pdb2gmx");
		TEST_ASSERT(sameIds(actual.anglebonds, expected.anglebonds), "Angles differ from pdb2gmx");
		TEST_ASSERT(sameIds(actual.dihedralbonds, expected.dihedralbonds), "Proper dihedrals differ from pdb2gmx");
		TEST_ASSERT(sameIds(actual.improperdihedralbonds, expected.improperdihedralbonds), "Impropers differ from pdb2gmx");

		Programs::ToGmxOptions neutral = options;
		neutral.cTerminus = Programs::ToGmxOptions::CTerminus::Neutral;
		bool rejected = false;
		try { Programs::ToGmx(AutomatedTestsDir() / "pdb2gmx" / "6lzm.pdb", neutral); }
		catch (const std::runtime_error&) { rejected = true; }
		TEST_ASSERT(rejected, "AMBER99SB accepted a neutral C-terminus");
		co_return LimaUnittestResult{ true, OneLineMessage("{} atoms match pdb2gmx", actual.atoms.size()), envmode == Full };
	}

	TestRoutine TestToGmx_multichain(Environment&, EnvMode envmode) {
		const fs::path directory = AutomatedTestsDir() / "pdb2gmx";
		const fs::path input = directory / "generated_multichain.pdb";
		TryDeleteFile(input);

		WriteTwoChainPdb(directory / "6lzm.pdb", input);
		auto conversion = Programs::ToGmx(input);
		TEST_ASSERT(conversion.positionRestraints.size() == 2
			&& !conversion.positionRestraints[0].GetSection(position_restraints).empty()
			&& !conversion.positionRestraints[1].GetSection(position_restraints).empty(),
			"Multi-chain position restraints were not generated correctly");

		TEST_ASSERT(conversion.topology.moleculetypes.size() == 2,
			"Multi-chain topology does not contain one molecule type per chain");
		TEST_ASSERT(conversion.topology.moleculetypes.contains("Protein_chain_A")
			&& conversion.topology.moleculetypes.contains("Protein_chain_B"),
			"Multi-chain topology is missing chain molecule names");

		const TopologyFile& parsedTopology = conversion.topology;
		GroFile parsedCoordinates = conversion.grofile;
		std::size_t topologyAtomCount = 0;
		for (const auto& molecule : parsedTopology.GetSystem().Instances()) {
			topologyAtomCount += molecule.moleculetype->atoms.size();
		}
		TEST_ASSERT(parsedTopology.GetSystem().molecules.size() == 2 && parsedTopology.GetSystem().MoleculeCount() == 2
			&& topologyAtomCount == parsedCoordinates.atoms.size(),
			"Parsed multi-chain topology does not map one-to-one to GRO atoms");

		const auto& molecules = parsedTopology.GetSystem().molecules;
		const std::size_t secondMoleculeOffset = molecules.front().moleculetype->atoms.size();
		const auto& secondMolecule = *molecules[1].moleculetype;
		TEST_ASSERT(!secondMolecule.singlebonds.empty(), "Multi-chain test molecule contains no bond to fragment");
		const auto fragmentedBond = secondMolecule.singlebonds.front().ids;
		const Float3 untouchedFirstChainAtom = parsedCoordinates.atoms.front().position;
		parsedCoordinates.atoms[secondMoleculeOffset + fragmentedBond[1]].position.x += parsedCoordinates.box_size.x;

		MoleculeUtils::MakeMoleculeWholeAfterPBCFragmentation(parsedCoordinates, parsedTopology);
		TEST_ASSERT(parsedCoordinates.atoms.front().position == untouchedFirstChainAtom,
			"Making the second molecule whole modified the first molecule");
		const float reconstructedBondDx = std::abs(
			parsedCoordinates.atoms[secondMoleculeOffset + fragmentedBond[0]].position.x
			- parsedCoordinates.atoms[secondMoleculeOffset + fragmentedBond[1]].position.x);
		TEST_ASSERT(reconstructedBondDx < parsedCoordinates.box_size.x / 2.f,
			"System-level molecule reconstruction did not unwrap the fragmented chain");

		constexpr float padding = 0.5f;
		MoleculeUtils::FitMoleculeInBox(parsedCoordinates, padding);
		for (const auto& atom : parsedCoordinates.atoms) {
			TEST_ASSERT(atom.position.x >= padding - 1e-5f && atom.position.x <= parsedCoordinates.box_size.x - padding + 1e-5f
				&& atom.position.y >= padding - 1e-5f && atom.position.y <= parsedCoordinates.box_size.y - padding + 1e-5f
				&& atom.position.z >= padding - 1e-5f && atom.position.z <= parsedCoordinates.box_size.z - padding + 1e-5f,
				"Fitted render box does not encompass every molecule atom");
		}
		co_return LimaUnittestResult{ true, "Success", envmode == Full };
	}
}
