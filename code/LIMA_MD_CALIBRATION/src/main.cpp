#include "MDCalibrations.h"
#include "Environment.h"

#include "LimaTypes.cuh"

#include <iostream>
#include <string>
#include <algorithm>


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
	const float particle_mass = 12.011000 * 1e-3f;
	
	env.loadSimParams(work_folder + "sim_params.txt");
	float dt = env.getSimparamRef()->dt;

	//float particle_temps[] = { 10.f, 100.f, 200.f, 273.f, 300.f, 500.f };	// Simulated temp of a single particle
	float particle_temps[] = { 200000 };
	for (auto temp : particle_temps) {
		env.CreateSimulation(conf, topol, work_folder);
		auto* box = env.getSim()->box;



		box->compounds[0].prev_positions[0] += Float3(-1, 0, 0) / NORMALIZER * EngineUtils::calcSpeedOfParticle(particle_mass, temp) * dt;
		env.run();

		Analyzer::printEnergy(env.getAnalyzedPackage());
		break;
	}



	return true;
}




int main() {

	/*coordPrecesionBenchmark();
	return 0;*/
	Environment env;

	basicBenchmark(env);
	//doPoolBenchmark(env);

	return 0;
}