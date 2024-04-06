#pragma once

#include "TestUtils.h"
#include "MoleculeGraph.h"


namespace FileTests {
	using namespace TestUtils;
	namespace fs = std::filesystem;

	// This is not ready to be a test, for now i just need it to create some lipids files
	LimaUnittestResult TestFilesAreCachedAsBinaries(EnvMode envmode) {
		const fs::path work_dir = simulations_dir + "/filetests";


		if (fs::exists(work_dir / "molecule/em.gro.bin"))
			fs::remove(work_dir / "molecule/em.gro.bin");
		if (fs::exists(work_dir / "molecule/topol.top.bin"))
			fs::remove(work_dir / "molecule/topol.top.bin");


		// Check grofiles first 
		{
			ParsedGroFile grofile{ work_dir / "molecule/em.gro" };
			if (!fs::exists(work_dir / "molecule/em.gro.bin")) {
				return LimaUnittestResult{ LimaUnittestResult::FAIL , "ParsedGroFile did not make a bin cached file", envmode == Full };
			}

			ParsedGroFile grofileLoadedFromCache{ work_dir / "molecule/em.gro" };
			if (!grofileLoadedFromCache.readFromCache) {
				return LimaUnittestResult{ LimaUnittestResult::FAIL , "ParsedGroFile was not read from cached file", envmode == Full };
			}


			if (grofile.title != grofileLoadedFromCache.title
				|| grofile.box_size != grofileLoadedFromCache.box_size
				|| grofile.atoms.size() != grofileLoadedFromCache.atoms.size()
				|| grofile.atoms.back().position != grofileLoadedFromCache.atoms.back().position
				) {
				return LimaUnittestResult{ LimaUnittestResult::FAIL , "Grofile did not match cached grofile", envmode == Full };
			}
		}


		//// Check topol files now
		//{
		//	ParsedTopologyFile topolfile{ work_dir / "molecule/topol.top" };
		//	if (!fs::exists(work_dir / "molecule/topol.top.bin")) {
		//		return LimaUnittestResult{ LimaUnittestResult::FAIL , "ParsedTopolFile did not make a bin cached file", envmode == Full };
		//	}

		//	ParsedTopologyFile topolfileLoadedFromCache{ work_dir / "molecule/topol.top" };
		//	if (!topolfileLoadedFromCache.readFromCache) {
		//		return LimaUnittestResult{ LimaUnittestResult::FAIL , "ParsedTopolFile was not read from cached file", envmode == Full };
		//	}

		//	if (topolfile.title != topolfileLoadedFromCache.title
		//		|| topolfile.molecules.size() != topolfileLoadedFromCache.molecules.size()
		//		|| topolfile.molecules.back().name != topolfileLoadedFromCache.molecules.back().name
		//		|| topolfile.molecules.back().atoms.size() != topolfileLoadedFromCache.molecules.back().atoms.size()
		//		|| topolfile.molecules.back().atoms.back().name != topolfileLoadedFromCache.molecules.back().atoms.back().name
		//		) {
		//		return LimaUnittestResult{ LimaUnittestResult::FAIL , "Topolfile did not match cached topolfile", envmode == Full };
		//	}
		//}


		return LimaUnittestResult{ LimaUnittestResult::SUCCESS , "No error", envmode == Full };
	}
}

