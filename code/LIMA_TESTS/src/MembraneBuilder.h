#pragma once

#include <filesystem>
#include <fstream>
#include <algorithm>

#include "TestUtils.h"

namespace TestMembraneBuilder {
	using namespace TestUtils;
	namespace fs = std::filesystem;





	static LimaUnittestResult testBuildmembraneSmall(EnvMode envmode, bool do_em)
	{
		const fs::path work_dir = simulations_dir + "/BuildMembraneSmall";
		Environment env{ work_dir.string(), envmode, false};

		env.CreateSimulation(7.f);
		env.createMembrane(do_em);

		const fs::path mol_dir = work_dir / "molecule";


		std::vector<std::array<std::string, 2>> files = { {"monolayer.gro", "monolayer_reference.gro"}, {"monolayer.top", "monolayer_reference.top"} };

		if (!do_em) {
			// These files are altered by the em, and thus the comparison cannot be made
			files.push_back({ "membrane.gro", "membrane_reference.gro" });
			files.push_back({ "membrane.top", "membrane_reference.top" });				
		}		

		for (const auto& pair : files) {
			const string error = compareFilesBitwise(mol_dir / pair[0], mol_dir / pair[1]);
			if (error != "") {
				return LimaUnittestResult{ LimaUnittestResult::FAIL , error, envmode == Full };
			}
		}

		return LimaUnittestResult{ LimaUnittestResult::SUCCESS , "No error", envmode == Full};
	}
}