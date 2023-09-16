#include <iostream>
#include <format>
#include <filesystem>

#include "Environment.h"


char* getCmdOption(char** begin, char** end, const std::string& option)
{
	char** itr = std::find(begin, end, option);
	if (itr != end && std::next(itr) != end) {
		return *std::next(itr);
	}
	return nullptr;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
	return std::find(begin, end, option) != end;
}


struct MdrunSetup {

	MdrunSetup() {
		work_dir = std::filesystem::current_path();
	}


	EnvMode envmode = Headless;

	std::string work_dir;
	std::string structure = work_dir + "molecule/conf.gro";
	std::string topol = work_dir + "molecule/topol.top";
	std::string simpar = work_dir + "sim_params.txt";
};


MdrunSetup parseProgramArguments(int argc, char** argv) {
	MdrunSetup setup{};

	char* user_structure = getCmdOption(argv, argv + argc, "-structure");
	if (user_structure)
	{
		setup.structure = setup.work_dir + user_structure;
	}

	char* user_topol = getCmdOption(argv, argv + argc, "-topology");
	if (user_topol)
	{
		setup.topol = setup.work_dir + user_topol;
	}

	char* user_params = getCmdOption(argv, argv + argc, "-simparams");
	if (user_params)
	{
		setup.simpar = setup.work_dir + user_params;
	}

	char* user_envmode = getCmdOption(argv, argv + argc, "-envmode");
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





int main(int argc, char** argv) 
{
	try {
		MdrunSetup setup = parseProgramArguments(argc, argv);

		printf(std::format("LIMA is preparing simulation in dir {}", setup.work_dir));

		auto env = std::make_unique<Environment>(setup.work_dir, setup.envmode);

		const InputSimParams ip = env->loadInputSimParams(setup.simpar);

		env->CreateSimulation(setup.structure, setup.topol, ip);

		env->run();
	}
	catch (const std::runtime_error& ex) {
		std::cerr << "LIMA encountered an exception:\n\t " << ex.what();
		return 1;
	}
	catch (...) {
		std::cerr << "LIMA caught an unknown exception";
		return 1;
	}

	return 0;
}