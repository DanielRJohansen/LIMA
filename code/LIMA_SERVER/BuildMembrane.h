#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

#include "CommandlineUtils.h"
namespace fs = std::filesystem;

struct BuildMembraneSetup{

	BuildMembraneSetup(int argc, char** argv) {
		work_dir = std::filesystem::current_path();

        int sizeValue = 0;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "-size") {
                boxsize = true;
                if (i + 1 < argc && isInteger(argv[i + 1])) {
                    sizeValue = std::stoi(argv[++i]);
                }
                else {
                    std::cerr << "Invalid -size argument. It must have exactly one integer value." << std::endl;
                }
            }
            else if (arg == "-lipids") {
                while (i + 2 < argc && argv[i + 1][0] != '-' && isInteger(argv[i + 2])) {
                    lipids.emplace_back(argv[++i], std::stoi(argv[++i]));
                }
                if ((i + 1 < argc && argv[i + 1][0] != '-') || (lipids.size() * 2 != argc - i - 1)) {
                    std::cerr << "Invalid -lipids argument. It must have a multiple of two values." << std::endl;
                }
            }
        }
	}

	EnvMode envmode = Full;

	fs::path work_dir;
    int boxsize = 0;
    std::vector<std::pair<std::string, int>> lipids;

private:
    bool isInteger(const std::string& s) {
        for (char c : s) {
            if (!isdigit(c)) return false;
        }
        return true;
    }
};


int buildMembrane(int argc, char** argv) {
	BuildMembraneSetup setup(argc, argv);

    auto a = setup.work_dir.string();
	Environment env{ a, setup.envmode, false};

	env.CreateSimulation(7.f);
	LipidsSelection lipidselection;
	for (const auto& lipid : setup.lipids) {
		lipidselection.emplace_back(LipidSelect{ lipid.first, lipid.second });
	}
	env.createMembrane(lipidselection, true);

	return 0;
}