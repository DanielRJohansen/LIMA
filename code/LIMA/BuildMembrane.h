#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

#include "Programs.h"

namespace fs = std::filesystem;

struct BuildMembraneSetup{

	BuildMembraneSetup(int argc, char** argv) : work_dir(std::filesystem::current_path())
    {		
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];

            if (arg == "-lipids") {
                // If we have atleasat 2 more args, and next arg is not a keyword, and second arg is a float
                while (i + 2 < argc && argv[i + 1][0] != '-' && isFloatingpoint(argv[i + 2])) {
                    lipids.emplace_back(argv[i+1], std::stod(argv[i+2]));
                    i+=2;
                }
                if ((i + 1 < argc && argv[i + 1][0] != '-') || (lipids.size() * 2 != argc - i - 1)) {
                    std::cerr << "Invalid -lipids argument. It must have a multiple of two values." << std::endl;
                }
            }
            else if (arg== "-centerZ") {
				if (i + 1 < argc) {
					membraneCenterZ = std::stof(argv[i + 1]);
					i++;
				}
				else {
					std::cerr << "-centerZ expected an argument." << std::endl;
				}
			}
            else if (arg == "-boxsize") {
                if (i+1 < argc) {
					boxsize = std::stof(argv[i+1]);
					i++;
				}
				else {
					std::cerr << "boxsize expected an argument." << std::endl;
				}
            }
            else {
                throw std::runtime_error("Unknown argument: " + arg);
            }
        }
	}

	EnvMode envmode = Full;
	const fs::path work_dir;

    std::vector<std::pair<std::string, double>> lipids; // {name, percentage}
    float membraneCenterZ = 0.0f;
    float boxsize;

private:
    bool isFloatingpoint(const std::string& s) {
        std::istringstream iss(s);
        float f;
        iss >> std::noskipws >> f; // no skipws to check for trailing characters
        return iss.eof() && !iss.fail();
    }
};


int buildMembrane(int argc, char** argv) {
	BuildMembraneSetup setup(argc, argv);

    const SimParams params{ SimParams::defaultPath() };

	LipidsSelection lipidselection;
	for (const auto& lipid : setup.lipids) {
		lipidselection.emplace_back(LipidSelect{ lipid.first, setup.work_dir, lipid.second });
	}
    auto [grofile, topfile] = Programs::CreateMembrane(setup.work_dir, lipidselection, setup.membraneCenterZ, true, setup.envmode);

    grofile->printToFile(setup.work_dir / "membrane.gro");
    topfile->printToFile(setup.work_dir / "membrane.top");
   
	return 0;
}
