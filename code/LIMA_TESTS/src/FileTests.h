#pragma once

#include "TestUtils.h"
#include "MoleculeGraph.h"


namespace FileTests {
	using namespace TestUtils;
	namespace fs = std::filesystem;

	LimaUnittestResult TestFilesAreCachedAsBinaries(EnvMode envmode) {
		const fs::path workDir = HeavyTestsDir() / "filetests";

		// Remove any file in workdir ending in .itp.bin
		for (const auto& entry : fs::directory_iterator(workDir / "molecule")) {
			auto a = entry.path().extension();
			if (entry.path().extension() == ".bin") {
				fs::remove(entry.path());
			}
		}


		// Check grofiles first 
		{
			TimeIt time1(".gro read", envmode==Full);
			GroFile grofile{ workDir / "molecule/em.gro" };
			time1.stop();
			if (!fs::exists(workDir / "molecule/em.gro.bin")) {
				return LimaUnittestResult{ false , "GroFile did not make a bin cached file", envmode == Full };
			}

			TimeIt time2(".bin read", envmode==Full);
			GroFile grofileLoadedFromCache{ workDir / "molecule/em.gro" };
			time2.stop();
			if (!grofileLoadedFromCache.readFromCache) {
				return LimaUnittestResult{ false , "GroFile was not read from cached file", envmode == Full };
			}


			if (grofile.title != grofileLoadedFromCache.title
				|| grofile.box_size != grofileLoadedFromCache.box_size
				|| grofile.atoms.size() != grofileLoadedFromCache.atoms.size()
				|| grofile.atoms.back().position != grofileLoadedFromCache.atoms.back().position
				) {
				return LimaUnittestResult{ false , "Grofile did not match cached grofile", envmode == Full };
			}
		}


		// Check topol files now
		{
			TimeIt time1(".top read", envmode == Full);
			TopologyFile topolfile{ workDir / "molecule/topol.top" };
			time1.stop();
			if (!fs::exists(workDir / "molecule/topol.top.bin")) {
				return LimaUnittestResult{ false , "ParsedTopolFile did not make a bin cached file", envmode == Full };
			}

			TimeIt time2(".bin read", envmode == Full);
			TopologyFile topolfileLoadedFromCache{ workDir / "molecule/topol.top" };
			time2.stop();
	/*		if (!topolfileLoadedFromCache.readFromCache) {
				return LimaUnittestResult{ false , "ParsedTopolFile was not read from cached file", envmode == Full };
			}*/

			// TODONOW Fix
		/*	if (topolfile.title != topolfileLoadedFromCache.title
				|| topolfile.molecules.entries.size() != topolfileLoadedFromCache.molecules.entries.size()
				|| topolfile.molecules.entries.back().includeTopologyFile->atoms.entries.size() != topolfileLoadedFromCache.molecules.entries.back().includeTopologyFile->atoms.entries.size()
				|| topolfile.molecules.entries.back().includeTopologyFile->atoms.entries.back().atomname != topolfileLoadedFromCache.molecules.entries.back().includeTopologyFile->atoms.entries.back().atomname
				|| topolfileLoadedFromCache.molecules.entries[0].includeTopologyFile->atoms.entries.back().resnr != 250
				) {
				return LimaUnittestResult{ false , "Topolfile did not match cached topolfile", envmode == Full };
			}*/
		}


		return LimaUnittestResult{ true , "No error", envmode == Full };
	}


	namespace {
		void WriteTextFile(const fs::path& path, std::string_view contents) {
			std::ofstream file(path);
			if (!file) throw std::runtime_error("Failed to write " + path.string());
			file << contents;
		}

		std::vector<std::string> AtomNames(const TopologyFile& topology, const std::string& moleculetype) {
			std::vector<std::string> names;
			for (const auto& atom : topology.moleculetypes.at(moleculetype)->atoms)
				names.push_back(atom.atomname);
			return names;
		}
	}

	// #defines behave like in the C preprocessor, which GROMACS uses: a #define in an included file applies to all lines after
	// that #include, also in the parent file. The topology parser scans files in parallel, so this verifies it resolves in order
	TestRoutine TestTopologyPreprocessor(Environment&, EnvMode envmode) {
		const fs::path dir = fs::temp_directory_path() / "lima_topology_preprocessor_test";
		fs::remove_all(dir);
		fs::create_directories(dir);

		// b.itp defines ruleB, which must exclude c.itp in the parent. It also tests #else, nested conditionals,
		// multiple moleculetypes in 1 file, and bonds to missing atoms
		WriteTextFile(dir / "b.itp", R"(#define ruleB
[ moleculetype ]
MolB 3

[ atoms ]
1 CT 1 RES C1 1 0.0 12.0
#ifdef ruleB
2 CT 1 RES C2 2 0.0 12.0
#else
2 CT 1 RES ELSE_WRONG 2 0.0 12.0
#endif
#ifdef ruleB
  #ifndef ruleB
3 CT 1 RES NESTED_WRONG 3 0.0 12.0
  #else
3 CT 1 RES C3 3 0.0 12.0
  #endif
#endif
#ifdef NOT_DEFINED
  #ifdef ruleB
4 CT 1 RES NESTED_WRONG 4 0.0 12.0
  #endif
#endif

[ bonds ]
1 2 1
2 3 1
3 99 1 ; atom 99 does not exist, so this bond is discarded

[ moleculetype ]
MolB2 3

[ atoms ]
1 CT 1 RES D1 1 0.0 12.0
)");
		WriteTextFile(dir / "c.itp", R"([ moleculetype ]
MolC 3

[ atoms ]
1 CT 1 RES C1 1 0.0 12.0
)");
		// d.itp undefines ruleB for the lines after it
		WriteTextFile(dir / "d.itp", R"(#undef ruleB
)");
		WriteTextFile(dir / "e.itp", R"([ moleculetype ]
MolE 3

[ atoms ]
1 CT 1 RES E1 1 0.0 12.0
)");

		WriteTextFile(dir / "topol.top", R"(; Preprocessor test

#include "b.itp"
#ifndef ruleB
#include "c.itp"
#endif
#ifdef NOT_DEFINED
#include "does_not_exist.itp"
#endif
#include "d.itp"
#ifndef ruleB
#include "e.itp"
#endif

[ system ]
Test

[ molecules ]
MolB 2
MolB2 1
MolE 1
)");

		const TopologyFile topology{ dir / "topol.top" };

		ASSERT(topology.moleculetypes.contains("MolB") && topology.moleculetypes.contains("MolB2"), "Moleculetypes from b.itp are missing");
		ASSERT(!topology.moleculetypes.contains("MolC"), "c.itp was included, even though b.itp defined ruleB before the #ifndef");
		ASSERT(topology.moleculetypes.contains("MolE"), "e.itp was not included, even though d.itp undefined ruleB");
		ASSERT((AtomNames(topology, "MolB") == std::vector<std::string>{ "C1", "C2", "C3" }), "MolB atoms do not match the active #ifdef branches");
		ASSERT((AtomNames(topology, "MolB2") == std::vector<std::string>{ "D1" }), "Second moleculetype in b.itp did not get its own atoms");
		ASSERT(topology.moleculetypes.at("MolB")->singlebonds.size() == 2, "Bond to a missing atom was not discarded");
		ASSERT(topology.GetSystem().MoleculeCount() == 4 && topology.GetSystem().molecules.size() == 3, "Wrong molecule count");
		ASSERT(std::ranges::distance(topology.GetAllElements<TopologyFile::AtomsEntry>()) == 3 * 2 + 1 + 1, "Wrong total number of atoms");

		fs::remove_all(dir);
		co_return LimaUnittestResult{ true, "No error", envmode == Full };
	}
}

