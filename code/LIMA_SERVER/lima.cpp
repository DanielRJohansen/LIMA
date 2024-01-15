#include <iostream>
#include <format>
#include <filesystem>
#include <algorithm>
#include <cctype>
#include "Environment.h"

#include "MoleculeGraph.h"

#include "CommandlineUtils.h"

namespace fs = std::filesystem;

struct MembraneBuilderSetup {

	MembraneBuilderSetup(int argc, char** argv) {
		work_dir = std::filesystem::current_path().string();
		coordinate_file = work_dir + "/molecule/conf.gro";
		topol_file = work_dir + "/molecule/topol.top";
		simparams = work_dir + "/sim_params.txt";



		char* user_structure = getCmdOption(argv, argv + argc, "-coordinates");
		if (user_structure)
		{
			coordinate_file = work_dir + user_structure;
		}

		char* user_topol = getCmdOption(argv, argv + argc, "-topology");
		if (user_topol)
		{
			topol_file = work_dir + user_topol;
		}

		char* user_params = getCmdOption(argv, argv + argc, "-simparams");
		if (user_params)
		{
			simparams = work_dir + user_params;
		}

		char* user_envmode = getCmdOption(argv, argv + argc, "-envmode");
		if (user_envmode)
		{
			const std::string user_envmode_str(user_envmode);

			if (user_envmode_str == "full") envmode = Full;
			else if (user_envmode_str == "console") envmode = ConsoleOnly;
			else if (user_envmode_str == "headless") envmode = Headless;
			else {
				throw std::runtime_error(std::format("Got illegal envmode parameter {}", user_envmode_str).c_str());
			}
		}


	}

	EnvMode envmode;

	std::string work_dir;
	std::string coordinate_file;
	std::string topol_file;
	std::string simparams;
};


int membraneBuilder(int argc, char** argv) {

	MembraneBuilderSetup setup(argc, argv);

	Environment env{ setup.work_dir, setup.envmode, true };
	
	const SimParams ip = env.loadSimParams(setup.simparams);

	env.CreateSimulation(setup.coordinate_file, setup.topol_file, ip);

	return 0;
}

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


int main(int argc, char** argv) 
{
	try {
		std::string program = argv[1];
		std::transform(program.begin(), program.end(), program.begin(),
			[](unsigned char c) { return std::tolower(c); });

		if (program == "membranebuilder") { membraneBuilder(argc, argv); }
		else {
			std::printf("Unregcognized lima program");
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
