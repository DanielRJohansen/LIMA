#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include "CommandlineUtils.h"
#include "Programs.h"

namespace fs = std::filesystem;


struct BuildMembraneSetup{
	BuildMembraneSetup(int argc, char** argv) : work_dir(std::filesystem::current_path())
    {		
        for (int i = 1; i < argc; ++i) {
            std::string arg = CmdLineUtils::ToLowercase(argv[i]);

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
            else if (arg== "-centerz") {
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
            else if (arg == "-help" || arg == "-h") {
                std::cout << helpText;
                exit(0);
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
    float boxsize = -1.f;

private:
    bool isFloatingpoint(const std::string& s) {
        std::istringstream iss(s);
        float f;
        iss >> std::noskipws >> f; // no skipws to check for trailing characters
        return iss.eof() && !iss.fail();
    }
    const std::string helpText = R"(
Usage: buildMembrane [OPTIONS]

Description:
    This command builds a membrane structure with specified lipid components and optional membrane positioning and box size.
    By default we support Stockholm lipids 2020. To use your own lipids, simply have the lipid files in the working
    directory, and provide their names (without extension) similarly to the examples below.

Options:
    -lipids {name] [percentage] ...    
        Specifies the lipid types and their respective percentages in the membrane. 
        You must provide an even number of arguments: each lipid name followed by its percentage (as a floating-point value. 
        Example: -lipids DPPC 70.5 DOPC 29.5)

    -centerz [value]    
        Sets the Z-coordinate (nm) of the membrane's center. 
        Example: -centerz 3.0

    -boxsize [value]    
        Defines the size (nm) of the simulation box. Non-cubic boxes are not yet supported
        Example: -boxsize 10.0

    -help, -h    
        Displays this help text and exits.

Example:
    buildMembrane -lipids DPPC 60 DOPC 40 -centerz 3.0 -boxsize 10.0
    This command creates a membrane with 60% DPPC and 40% DOPC, centered at Z = 3.0, within a box of size 10.
)";
};


int buildMembrane(int argc, char** argv) {
	BuildMembraneSetup setup(argc, argv);

    const SimParams params{ SimParams::defaultPath() };

	Lipids::Selection lipidselection;
	for (const auto& lipid : setup.lipids) {
		lipidselection.emplace_back(Lipids::Select{ lipid.first, setup.work_dir, lipid.second });
	}
    auto [grofile, topfile] = Programs::CreateMembrane(setup.work_dir, lipidselection, setup.membraneCenterZ, true, setup.envmode);

    grofile->printToFile(setup.work_dir / "membrane.gro");
    topfile->printToFile(setup.work_dir / "membrane.top");
   
	return 0;
}
