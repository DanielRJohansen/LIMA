#include <iostream>
#include "render.h"
#include "selftest.h"

#include "CommandlineUtils.h"


#include "BuildMembrane.h"
#include "mdrun.h"
#include "GetForcefieldParams.h"

namespace fs = std::filesystem;


static const std::string helpText = R"(
Usage: lima <program> [OPTIONS]

Description:
    LIMA is a suite of tools for molecular dynamics and membrane simulation tasks. 
    Each program in the suite has its own set of options and parameters.
	Input "lima <program> -help" for help on a specific program.

Programs:
    mdrun               Runs the molecular dynamics simulation.
    buildmembrane       Builds a membrane structure with specified lipid composition and box size. 
    makesimparams       Generates default simulation parameters and writes them to a file.
    selftest            Runs internal self-tests to validate functionality.
	render				Render a molecule in a 3D viewer.	

Options:
    -help, -h           Displays this help message and exits.

Examples:
    lima buildmembrane -help
        Displays help for the buildmembrane program.
        
    lima makesimparams
        Generates and outputs default simulation parameters to a file in the current directory
)";



int main(int argc, char** argv) 
{
	try {
		const std::string program = argv[1];

		if (program == "mdrun") { mdrun(argc, argv); }
		else if (program == "buildmembrane") { buildMembrane(argc, argv); }
		else if (program == "makesimparams") {SimParams params{}; params.dumpToFile();}
		else if (program == "-help" || program =="-h"||program == "help") { std::cout << helpText; }
		else if (program == "selftest") { SelfTest(); }
		else if (program == "render") { render(argc, argv); }
		//else if (program == "getforcefieldparams") { GetForcefieldParams(); }
		else {
			std::cout << "Unregcognized lima program: " << program<< "\n";
		}
	}
	catch (const std::runtime_error& ex) {
		std::cerr << "LIMA encountered an exception:\n\t " << ex.what() << std::endl;
		return 1;
	}
	catch (...) {
		std::cerr << "LIMA caught an unknown exception\n";
		return 1;
	}

	return 0;
}
