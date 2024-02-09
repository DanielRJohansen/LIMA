#include <iostream>
#include <format>
#include <filesystem>
#include <algorithm>
#include <cctype>
#include "Environment.h"

#include "MoleculeGraph.h"

#include "CommandlineUtils.h"


#include "BuildMembrane.h"
#include "mdrun.h"

namespace fs = std::filesystem;


int reorderMoleculeParticles(int argc, char** argv) {
	if (argc != 3)
		throw std::runtime_error(std::format("reordermoleculeparticles expected 2 arguments (.gro file & .top file), but got {}", argc - 1));

	const fs::path gro_path_in = argv[1];
	const fs::path top_path_in = argv[2];
	const fs::path gro_path_out = argv[3];
	const fs::path top_path_out = argv[4];

	LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(gro_path_in, top_path_in, gro_path_out, top_path_out);

	return 0;
}

// TODO: Expand this with the options for each program
void printHelp() {
	std::cout << "LIMA - the faster Molecular Dynamics Engine\n";
	std::cout << "Usage: lima <program> <args>\n";
	std::cout << "Available programs:\n";
	std::cout << "  mdrun\n";
	std::cout << "  buildmembrane\n";
	//std::cout << "  reordermoleculeparticles\n";
	std::cout << "  makesimparams\n";
}


int main(int argc, char** argv) 
{
	try {
		const std::string program = argv[1];

		if (program == "mdrun") { mdrun(argc, argv); }
		else if (program == "buildmembrane") { buildMembrane(argc, argv); }
		else if (program == "makesimparams") {SimParams params{}; params.dumpToFile();}
		else if (program == "-help" || program =="-h"||program == "help") { printHelp(); }
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
