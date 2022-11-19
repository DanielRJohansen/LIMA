#include "MDCalibrations.h"
#include "Environment.h"

#include "LimaTypes.cuh"

#include <iostream>
#include <string>


bool basicBenchmark(Environment& env) {
	const string conf_path = "C:\\PROJECTS\\Quantom\\Simulation\\Molecule\\conf.gro";
	const string topol_path = "C:\\PROJECTS\\Quantom\\Simulation\\Molecule\\topol.top";
	const std::string work_folder = "C:\\PROJECTS\\Quantom\\Simulation\\T4Lysozyme\\";

	env.prepFF(conf_path, topol_path);

	env.CreateSimulation(conf_path, topol_path, work_folder);
	env.run();
	return true;
}

bool doPoolBenchmark(Environment& env) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Pool/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";

	env.loadSimParams(work_folder + "sim_params.txt");
	auto sim_params = env.getSimparamRef();
	sim_params->n_steps = 200;

	env.CreateSimulation(conf, topol, work_folder);	

	env.run();

	return true;
}




int main() {
	Environment env;

	//basicBenchmark(env);
	doPoolBenchmark(env);

	return 0;
}