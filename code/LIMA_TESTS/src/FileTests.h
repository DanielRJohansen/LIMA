#pragma once

#include "TestUtils.h"
#include "MoleculeGraph.h"


namespace FileTests {
	using namespace TestUtils;
	namespace fs = std::filesystem;

	LimaUnittestResult TestFilesAreCachedAsBinaries(EnvMode envmode) {
		const fs::path work_dir = simulations_dir + "/filetests";

		// Remove any file in workdir ending in .itp.bin
		for (const auto& entry : fs::directory_iterator(work_dir / "molecule")) {
			auto a = entry.path().extension();
			if (entry.path().extension() == ".bin") {
				fs::remove(entry.path());
			}
		}


		// Check grofiles first 
		{
			TimeIt time1(envmode, ".gro read");
			ParsedGroFile grofile{ work_dir / "molecule/em.gro" };
			time1.stop();
			if (!fs::exists(work_dir / "molecule/em.gro.bin")) {
				return LimaUnittestResult{ LimaUnittestResult::FAIL , "ParsedGroFile did not make a bin cached file", envmode == Full };
			}

			TimeIt time2(envmode, ".bin read");
			ParsedGroFile grofileLoadedFromCache{ work_dir / "molecule/em.gro" };
			time2.stop();
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


		// Check topol files now
		{
			TimeIt time1(envmode, ".top read");
			ParsedTopologyFile topolfile{ work_dir / "molecule/topol.top" };
			time1.stop();
			if (!fs::exists(work_dir / "molecule/topol.top.bin")) {
				return LimaUnittestResult{ LimaUnittestResult::FAIL , "ParsedTopolFile did not make a bin cached file", envmode == Full };
			}

			TimeIt time2(envmode, ".bin read");
			ParsedTopologyFile topolfileLoadedFromCache{ work_dir / "molecule/topol.top" };
			time2.stop();
			if (!topolfileLoadedFromCache.readFromCache) {
				return LimaUnittestResult{ LimaUnittestResult::FAIL , "ParsedTopolFile was not read from cached file", envmode == Full };
			}

			// TODONOW Fix
		/*	if (topolfile.title != topolfileLoadedFromCache.title
				|| topolfile.molecules.entries.size() != topolfileLoadedFromCache.molecules.entries.size()
				|| topolfile.molecules.entries.back().includeTopologyFile->atoms.entries.size() != topolfileLoadedFromCache.molecules.entries.back().includeTopologyFile->atoms.entries.size()
				|| topolfile.molecules.entries.back().includeTopologyFile->atoms.entries.back().atomname != topolfileLoadedFromCache.molecules.entries.back().includeTopologyFile->atoms.entries.back().atomname
				|| topolfileLoadedFromCache.molecules.entries[0].includeTopologyFile->atoms.entries.back().resnr != 250
				) {
				return LimaUnittestResult{ LimaUnittestResult::FAIL , "Topolfile did not match cached topolfile", envmode == Full };
			}*/
		}


		return LimaUnittestResult{ LimaUnittestResult::SUCCESS , "No error", envmode == Full };
	}
}

