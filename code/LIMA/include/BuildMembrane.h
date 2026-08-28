#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include "CommandlineUtils.h"
#include "Programs.h"
#include "argparser.h"

namespace fs = std::filesystem;


int buildMembrane(int argc, char** argv) {

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
    
    -emtol [value]
        Sets the force tolerance (kJ/mol/nm) for the energy minimization. Default is 100.
        Example: -emtol 100

    -display, -d    
        Flag to enable the display, rendering the simulation and displaying information such as temperature, step and more.
    -working_dir, -workdir, -wd [path]
        Specifies the directory where the membrane files will be created. 
        Default is the current directory.
        Example: -working_dir /path/to/working/dir
    -help, -h    
        Displays this help text and exits.

Example:
    buildMembrane -lipids DPPC 60 DOPC 40 -centerz 3.0 -boxsize 10.0
    This command creates a membrane with 60% DPPC and 40% DOPC, centered at Z = 3.0, within a box of size 10.
)";


    EnvMode envmode = ConsoleOnly;
    fs::path workDir{ std::filesystem::current_path() };

    std::vector<std::pair<std::string, double>> lipids; // {name, percentage}
    std::optional<float> membraneCenterZ = std::nullopt;
    Float3 boxsize{};
    float emtol = 100.f;    

    ArgParser argparser(helpText);

	argparser.AddOption({ "-lipids" }, true, 
        [&lipids](const std::vector<std::string>& args) {
		if (args.size() % 2 != 0) {
			throw std::runtime_error("Invalid -lipids argument. It must have a multiple of two values.");
		}
		for (size_t i = 0; i < args.size(); i += 2) {
            std::string lipidname = args[i];
            double lipidPercentage = 0;
            try {
                lipidPercentage = std::stod(args[i+1]);
            }
            catch (...) {
				throw std::runtime_error(std::format("Invalid lipid percentage: '{}'. Expected a floating-point value.", args[i]));
                exit(1);
            }
			lipids.emplace_back(lipidname, lipidPercentage);
		}
		}
    );

    argparser.AddOption({ "-centerz", "-c"}, false, membraneCenterZ);
	argparser.AddOption({ "-boxsize", "-b" }, true, boxsize, true);
	argparser.AddOption({ "-emtol", "-tolerance" }, false, emtol);
    argparser.AddOption({ "-working_dir", "-workdir", "-wd" }, false, workDir);
	argparser.AddFlag({ "-display", "-d" }, [&envmode]() { envmode = Full; });

    argparser.Parse(argc, argv);

	Lipids::Selection lipidselection;
	for (const auto& lipid : lipids) {
		lipidselection.emplace_back(Lipids::Select{ lipid.first, workDir, lipid.second });
	}
    
    GroFile grofile;
    grofile.box_size = boxsize;
    grofile.title = "Membrane";
    TopologyFile topfile;
    topfile.SetSystem("Membrane");
    SimulationBuilder::CreateMembrane(grofile, topfile, lipidselection, membraneCenterZ.value_or(boxsize.z/2.f));
    auto sim = Programs::EnergyMinimize(grofile, topfile, true, workDir, envmode, true, emtol);

    grofile.printToFile(workDir / "membrane.gro");
    topfile.printToFile(workDir / "membrane.top");
        
    auto [step, force] = *std::min_element(sim->maxForceBuffer.begin(), sim->maxForceBuffer.end(),
        [](const std::pair<int64_t, float>& a, const std::pair<int64_t, float>& b) {
            return a.second < b.second;
        }
    );
    std::cout << std::format("buildmembrane finished with a min max-force of {:.3f}\n", force);

	return 0;
}
