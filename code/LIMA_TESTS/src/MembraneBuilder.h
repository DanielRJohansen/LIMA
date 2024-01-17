#pragma once

#include <filesystem>
#include <fstream>
#include <algorithm>

#include "TestUtils.h"

namespace TestMembraneBuilder {
	using namespace TestUtils;
	namespace fs = std::filesystem;

	string compareFiles(const std::filesystem::path& path1, const std::filesystem::path& path2) {
		// Open the files
		std::ifstream file1(path1, std::ifstream::ate);
		std::ifstream file2(path2, std::ifstream::ate);

		// Check if both files are open
		if (!file1.is_open() || !file2.is_open()) {
			return std::format("Failed to open either or both files \n\t\t{} \n\t\t{}", path1.string(), path2.string());
		}

		// Validate the files. If they are not even 50 bytes long, something is surely wrong
		if ( file1.tellg() < 50 || file2.tellg() < 50) {
			return std::format("Expected files to be atleast 50 bytes long \n\t\t{} \n\t\t{}", path1.string(), path2.string());
		}
		//// Compare file sizes
		//file1.seekg(0, std::ifstream::end);
		//file2.seekg(0, std::ifstream::end);
		//if (file1.tellg() != file2.tellg()) {
		//	return "Files are of different length";
		//}
		// 
		
		// Move ptr back to beginning of file
		file1.seekg(0, std::ifstream::beg);
		file2.seekg(0, std::ifstream::beg);

		// Compare the contents
		if (!std::equal(std::istreambuf_iterator<char>(file1.rdbuf()), std::istreambuf_iterator<char>(), std::istreambuf_iterator<char>(file2.rdbuf()))) {
			return std::format("Files did not match bit for bit \n\t\t{} \n\t\t{}", path1.string(), path2.string());
		};
		return "";
	}





	static LimaUnittestResult testBuildmembraneSmall(EnvMode envmode, bool do_em)
	{
		const fs::path work_dir = simulations_dir + "/BuildMembraneSmall";
		Environment env{ work_dir.string(), envmode, false};

		env.CreateSimulation(7.f);
		env.createMembrane(do_em);

		const fs::path mol_dir = work_dir / "molecule";
		const string err1 = compareFiles(mol_dir / "membrane.gro", mol_dir / "membrane_reference.gro");
		if (err1 != "") {
			return LimaUnittestResult{ LimaUnittestResult::FAIL , err1, envmode == Full };
		}
		const string err2 = compareFiles(mol_dir / "membrane.top", mol_dir / "membrane_reference.top");
		if (err2 != "") {
			return LimaUnittestResult{ LimaUnittestResult::FAIL , err2, envmode == Full };
		}

		return LimaUnittestResult{ LimaUnittestResult::SUCCESS , "No error", envmode == Full};

		//// Do em
		//env->run(true);
		////Analyzer::findAndDumpPiecewiseEnergies(*env->getSimPtr(), env->getWorkdir());

		//// Do sim
		////InputSimParams simparams{ 100, 2000 };
		//const std::string work_folder = simulations_dir + folder_name + "/";
		//const std::string simpar_path = work_folder + "sim_params.txt";
		//const InputSimParams simparams = env->loadInputSimParams(simpar_path);
		//auto sim = env->getSim();
		//env->CreateSimulation(*sim, simparams);
		//env->run();
		////Analyzer::findAndDumpPiecewiseEnergies(*env->getSimPtr(), env->getWorkdir());

		//const auto analytics = env->getAnalyzedPackage();

		//if (envmode != Headless) {
		//	Analyzer::printEnergy(analytics);
		//	LIMA_Print::printMatlabVec("cv", std::vector<float>{ analytics->variance_coefficient});
		//	LIMA_Print::printMatlabVec("energy_gradients", std::vector<float>{ analytics->energy_gradient});
		//}

		//const auto result = evaluateTest({ analytics->variance_coefficient }, max_vc, { analytics->energy_gradient }, max_gradient);
		//const auto status = result.first == true ? LimaUnittestResult::SUCCESS : LimaUnittestResult::FAIL;

		//return LimaUnittestResult{ status, result.second, envmode == Full };



	}
}