#include <iostream>
#include <format>
#include <filesystem>

#include <CommandlineUtils.h>

#include "Environment.h"

namespace fs = std::filesystem;
struct MdrunSetup {

	MdrunSetup() {
		work_dir = std::filesystem::current_path();
		structure = work_dir / "molecule/conf.gro";
		topol = work_dir / "molecule/topol.top";
		simpar = work_dir / "sim_params.txt";
	}


	EnvMode envmode = Full;

	fs::path work_dir;
	fs::path structure;
	fs::path topol;
	fs::path simpar;
};


MdrunSetup parseProgramArguments(int argc, char** argv) {
	MdrunSetup setup{};

	char* user_structure = CmdLineUtils::getCmdOption(argv, argv + argc, "-structure");
	if (user_structure)
	{
		setup.structure = setup.work_dir / user_structure;
	}

	char* user_topol = CmdLineUtils::getCmdOption(argv, argv + argc, "-topology");
	if (user_topol)
	{
		setup.topol = setup.work_dir / user_topol;
	}

	char* user_params = CmdLineUtils::getCmdOption(argv, argv + argc, "-simparams");
	if (user_params)
	{
		setup.simpar = setup.work_dir / user_params;
	}

	char* user_envmode = CmdLineUtils::getCmdOption(argv, argv + argc, "-envmode");
	if (user_envmode)
	{
		const std::string user_envmode_str(user_envmode);

		if (user_envmode_str == "full") setup.envmode = Full;
		else if (user_envmode_str == "console") setup.envmode = ConsoleOnly;
		else if (user_envmode_str == "headless") setup.envmode = Headless;
		else {
			throw std::runtime_error(std::format("Got illegal envmode parameter {}", user_envmode_str).c_str());
		}
	}

	return setup;
}





int mdrun(int argc, char** argv) 
{
	std::cout << "LIMA is preparing simulation in dir " << std::filesystem::current_path().string() << "\n";
	MdrunSetup setup = parseProgramArguments(argc, argv);
	auto env = std::make_unique<Environment>(setup.work_dir, setup.envmode, true);

	const SimParams ip(setup.simpar);		
	GroFile grofile{setup.structure};
	TopologyFile topfile{setup.topol};

	env->CreateSimulation(grofile, topfile, ip);
	env->run();

	return 0;
}
