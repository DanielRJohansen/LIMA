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


		std::vector<std::array<std::string, 2>> files = { {"monolayer.gro", "monolayer_reference.gro"}, {"monolayer.top", "monolayer_reference.top"} };

		if (!do_em) {
			// These files are altered by the em, and thus the comparison cannot be made
			files.push_back({ "membrane.gro", "membrane_reference.gro" });
			files.push_back({ "membrane.top", "membrane_reference.top" });				
		}		

		for (const auto& pair : files) {
			const string error = compareFiles(mol_dir / pair[0], mol_dir / pair[1]);
			if (error != "") {
				return LimaUnittestResult{ LimaUnittestResult::FAIL , error, envmode == Full };
			}
		}

		return LimaUnittestResult{ LimaUnittestResult::SUCCESS , "No error", envmode == Full};
	}
}