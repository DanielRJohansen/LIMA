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
            const std::string arg = argv[i];

            if (arg == "-lipids") {
                // If we have atleasat 2 more args, and next arg is not a keyword, and second arg is an integer
                while (i + 2 < argc && argv[i + 1][0] != '-' && isInteger(argv[i + 2])) {
                    lipids.emplace_back(argv[i+1], std::stoi(argv[i+2]));
                    i+=2;
                }
                if ((i + 1 < argc && argv[i + 1][0] != '-') || (lipids.size() * 2 != argc - i - 1)) {
                    std::cerr << "Invalid -lipids argument. It must have a multiple of two values." << std::endl;
                }
            }
        }
	}

	EnvMode envmode = Full;

	fs::path work_dir;
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
    const SimParams params{ SimParams::defaultPath() };



	env.CreateSimulation(params.box_size);
	LipidsSelection lipidselection;
	for (const auto& lipid : setup.lipids) {
		lipidselection.emplace_back(LipidSelect{ lipid.first, lipid.second });
	}
	env.createMembrane(lipidselection, true);

	return 0;
}
