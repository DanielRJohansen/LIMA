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
	float dt = env.getSimparamRef()->dt;

	auto* sim_params = env.getSimparamRef();

	//float particle_temps[] = { 10.f, 100.f, 200.f, 273.f, 300.f, 500.f };	// Simulated temp of a single particle
	//std::vector<float> particle_temps{ 300, 600, 2000, 4000, 8000 };// , 1000, 2000, 5000, 10000
	std::vector<float> particle_temps{ 1000};// , 1000, 2000, 5000, 10000
	std::vector<float> std_devs;

	for (auto temp : particle_temps) {
		auto vel = EngineUtils::calcSpeedOfParticle(particle_mass, temp * 2.f); // temp * 2 as 1 particles stores the temp of both in this test. [m/s] or [fm/fs]
		int steps_for_full_interaction = 8000000 / static_cast<int>(vel);
		sim_params->n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100) * 50;
		env.CreateSimulation(conf, topol, work_folder);

		
		auto coordarray_prev_ptr = env.getCoordarrayPtr("prev");
		coordarray_prev_ptr[1].rel_positions[0] += Float3(1, 0, 0) * vel * dt * FEMTO_TO_LIMA;	// convert vel from fm/fs to Lm/fs

		env.run();

		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		std_devs.push_back(Analyzer::getStdDevNorm(analytics->total_energy));
	}

	LIMA_Print::printMatlabVec("temperature", particle_temps);
	LIMA_Print::printMatlabVec("std_devs", std_devs);

	return true;
}




int main() {



	/*coordPrecesionBenchmark();
	return 0;*/
	Environment env;

	//basicBenchmark(env);
	doPoolBenchmark(env);

	return 0;
}