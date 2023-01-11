#include "MDCalibrations.h"
#include "Environment.h"

#include "LimaTypes.cuh"

#include <iostream>
#include <string>
#include <algorithm>
#include "Printer.h"
#include "Utilities.h"

//bool coordPrecesionBenchmark() {
//	Float3 pos1{ 3.5, 4, 4 };
//	Float3 pos2{ 4, 4, 4 };
//
//	Coord c1{ pos1 };
//	Coord c2{ pos2 };
//
//	double dist_c = sqrt(c1.distSqAbs(&c2));
//	double dist_f = (pos1 - pos2).len();
//
//	printf("dist_c %.12f dist_f %.12f\n", dist_c, dist_f);
//
//	return true;
//}



bool basicBenchmark(Environment& env) {
	const string conf_path = "C:/PROJECTS/Quantom/Simulation/Molecule/conf.gro";
	const string topol_path = "C:/PROJECTS/Quantom/Simulation/Molecule/topol.top";
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/T4Lysozyme";


	env.CreateSimulation(conf_path, topol_path, work_folder);
	env.run();
	return true;
}

bool doPoolBenchmark(Environment& env) {

	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Pool/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";
	const float particle_mass = 12.011000f * 1e-3f;
	
	env.loadSimParams(work_folder + "sim_params.txt");
	const float dt = env.getSimparamRef()->dt;

	auto* sim_params = env.getSimparamRef();

	std::vector<float> particle_temps{ 400, 1200, 2400, 4800 };// , 1000, 2000, 5000, 10000
	std::vector<float> std_devs;

	for (auto temp : particle_temps) {
		auto vel = EngineUtils::calcSpeedOfParticle(particle_mass, temp) * 2.f; //  *2 as 1 particles stores the velocity of both at temperature temp. [m/s]
		int steps_for_full_interaction = 10000000 / static_cast<int>(vel);
		sim_params->n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
		env.CreateSimulation(conf, topol, work_folder);

		
		auto coordarray_prev_ptr = env.getCoordarrayPtr("prev");
		coordarray_prev_ptr[1].rel_positions[0] += (Float3(1, 0, 0) * vel) * dt;	// convert vel from fm/fs to Lm/fs

		env.run();

		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		std_devs.push_back(Analyzer::getStdDevNorm(analytics->total_energy));
	}

	LIMA_Print::printMatlabVec("temperature", particle_temps);
	LIMA_Print::printMatlabVec("std_devs", std_devs);

	return true;
}

bool doSpringBenchmark(Environment& env) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Spring/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";
	const float particle_mass = 12.011000f * 1e-3f;

	env.loadSimParams(work_folder + "sim_params.txt");
	float dt = env.getSimparamRef()->dt;

	auto* sim_params = env.getSimparamRef();

	std::vector<float> bond_len_errors{ 0, 0.01, 0.05, 0.1 }; //(r-r0) [nm]
	std::vector<float> std_devs;

	for (auto bond_len_error : bond_len_errors) {		
		env.CreateSimulation(conf, topol, work_folder);

		auto coordarray_ptr = env.getCoordarrayPtr("current");
		auto coordarray_prev_ptr = env.getCoordarrayPtr("prev");
		coordarray_ptr[0].rel_positions[1].x += static_cast<int32_t>(bond_len_error * NANO_TO_LIMA);
		coordarray_prev_ptr[0].rel_positions[1].x += static_cast<int32_t>(bond_len_error * NANO_TO_LIMA);

		env.run();

		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		std_devs.push_back(Analyzer::getStdDevNorm(analytics->total_energy));
	}

	LIMA_Print::printMatlabVec("bond_len_errors", bond_len_errors);
	LIMA_Print::printMatlabVec("std_devs", std_devs);

	return true;
}


int main() {



	/*coordPrecesionBenchmark();
	return 0;*/
	Environment env;

	//basicBenchmark(env);
	//doPoolBenchmark(env);
	doSpringBenchmark(env);
	return 0;
}