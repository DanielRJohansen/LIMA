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
	//const string conf_path = "C:/PROJECTS/Quantom/Simulation/Molecule/conf.gro";
	//const string topol_path = "C:/PROJECTS/Quantom/Simulation/Molecule/topol.top";
	//const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/T4Lysozyme";
	//env.CreateSimulation(conf_path, topol_path, work_folder);
	//env.run();
	//return true;


	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Met/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";

	env.loadSimParams(work_folder + "sim_params.txt");
	env.CreateSimulation(conf, topol, work_folder);

	env.run();

	auto analytics = env.getAnalyzedPackage();
	Analyzer::printEnergy(analytics);

	std::vector<float> std_devs;
	std_devs.push_back(Analyzer::getStdDevNorm(analytics->total_energy));


	LIMA_Print::printMatlabVec("std_devs", std_devs);

	return true;
}

bool doProteinBenchmark(Environment& env) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/T4LysozymeNoSolvent/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";

	env.loadSimParams(work_folder + "sim_params.txt");
	env.CreateSimulation(conf, topol, work_folder);

	env.run();

	auto analytics = env.getAnalyzedPackage();
	Analyzer::printEnergy(analytics);

	std::vector<float> std_devs;
	std_devs.push_back(Analyzer::getStdDevNorm(analytics->total_energy));


	LIMA_Print::printMatlabVec("std_devs", std_devs);

	return true;
}

// Test assumes two carbons particles in conf
bool doPoolBenchmark(Environment& env, const string& foldername) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/" +foldername + "/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";
	const float particle_mass = 12.011000f / 1000.f;	// kg/mol
	
	env.loadSimParams(work_folder + "sim_params.txt");
	const float dt = env.getSimparamRef()->dt;

	auto* sim_params = env.getSimparamRef();

	//std::vector<float> particle_temps{ 400, 1200, 2400, 4800 };// , 1000, 2000, 5000, 10000
	std::vector<float> particle_temps{ 400, 800, 1200};
	std::vector<float> std_devs;

	for (auto temp : particle_temps) {
		//const float vel = EngineUtils::calcSpeedOfParticle(particle_mass, temp) * 2.f; //  *2 as 1 particles stores the velocity of both at temperature temp. [m/s]
		const float vel = EngineUtils::tempToVelocity(temp, particle_mass);	// [m/s] <=> [lm/ls]
		int steps_for_full_interaction = 8000000 / static_cast<int>(vel);
		sim_params->n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
		env.CreateSimulation(conf, topol, work_folder);

		
		auto coordarray_prev_ptr = env.getCoordarrayPtr("prev");
		coordarray_prev_ptr[0].rel_positions[0] += Coord{ (Float3(-1, 0, 0) * vel) * dt };	
		coordarray_prev_ptr[1].rel_positions[0] += Coord{ (Float3(1, 0, 0) * vel) * dt };


		env.run();

		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		std_devs.push_back(Analyzer::getStdDevNorm(analytics->total_energy));
	}

	LIMA_Print::printMatlabVec("temperature", particle_temps);
	LIMA_Print::printMatlabVec("std_devs", std_devs);

	return true;
}

// Test assumes two carbons particles in conf
bool doPoolCompSolBenchmark(Environment& env, const string& foldername) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/" + foldername + "/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";
	//const float particle_mass = 12.011000f / 1000.f;	// kg/mol

	env.loadSimParams(work_folder + "sim_params.txt");
	const float dt = env.getSimparamRef()->dt;

	auto* sim_params = env.getSimparamRef();

	//std::vector<float> particle_temps{ 400, 1200, 2400, 4800 };// , 1000, 2000, 5000, 10000
	std::vector<float> particle_temps{ 400, 800, 1200 };
	std::vector<float> std_devs;

	for (auto temp : particle_temps) {
		// Give the carbon a velocity
		{
			const float vel = EngineUtils::tempToVelocity(temp, 12.011000f / 1000.f);	// [m/s] <=> [lm/ls]
			int steps_for_full_interaction = 8000000 / static_cast<int>(vel);
			sim_params->n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
			env.CreateSimulation(conf, topol, work_folder);

			auto coordarray_prev_ptr = env.getCoordarrayPtr("prev");
			coordarray_prev_ptr[0].rel_positions[0] += Coord{ Float3(-1, 0, 0) * vel * dt };
		}
		
		// Give the solvent a velocty
		{
			// We need to also give the solvent velocity. But we don't know which block it is in....
			const float vel = EngineUtils::tempToVelocity(temp, SOLVENT_MASS);	// [m/s] <=> [lm/ls]
			auto solventblockgridprev = env.getAllSolventBlocksPrev();
			for (int i = 0; i < solventblockgridprev->blocks_total; i++) {
				auto block = solventblockgridprev->getBlockPtr(i);
				if (block->n_solvents == 1) {
					block->rel_pos[0] += Coord{ Float3{1, 0, 0} *vel * dt };
					int a = 0;
				}
			}
		}
		

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

	std::vector<float> bond_len_errors{ 0.01f, 0.02f }; //(r-r0) [nm]
	std::vector<float> std_devs;

	for (auto bond_len_error : bond_len_errors) {		
		env.CreateSimulation(conf, topol, work_folder);

		auto coordarray_ptr = env.getCoordarrayPtr("current");
		auto coordarray_prev_ptr = env.getCoordarrayPtr("prev");
		coordarray_ptr[0].rel_positions[0].x -= static_cast<int32_t>(bond_len_error * NANO_TO_LIMA);
		coordarray_prev_ptr[0].rel_positions[0].x -= static_cast<int32_t>(bond_len_error * NANO_TO_LIMA);

		env.run();

		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		std_devs.push_back(Analyzer::getStdDevNorm(analytics->total_energy));
	}

	LIMA_Print::printMatlabVec("bond_len_errors", bond_len_errors);
	LIMA_Print::printMatlabVec("std_devs", std_devs);

	return true;
}

bool doAngleBenchmark(Environment& env) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/AngleBenchmark/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";

	env.loadSimParams(work_folder + "sim_params.txt");

	auto* sim_params = env.getSimparamRef();

	const float relaxed_angle = 1.8849f; // [rad]
	std::vector<float> angle_errors{ 0.5f }; //(t-t0) [rad]
	std::vector<float> std_devs;

	for (auto angle_error : angle_errors) {
		env.CreateSimulation(conf, topol, work_folder);

		auto coordarray_ptr = env.getCoordarrayPtr("current");
		auto coordarray_prev_ptr = env.getCoordarrayPtr("prev");

		// First rotate particle #3 to the relaxed position + the error angle
		Float3 p3_pos = coordarray_ptr[0].rel_positions[2].toFloat3();
		p3_pos.rotateAroundOrigo(Float3{ 0.f, relaxed_angle + angle_error, 0.f });
		//coordarray_ptr[0].rel_positions[2] = p3_pos;		// Temp disabled, fix soon plz
		//coordarray_prev_ptr[0].rel_positions[2] = p3_pos;		// Temp disabled, fix soon plz

		// Now center all 3 particles
		for (auto i = 0; i < 3; i++) {
			coordarray_ptr[0].rel_positions[i] += Coord{ 350'000'000 };
			coordarray_prev_ptr[0].rel_positions[i] += Coord{ 350'000'000 };
		}

		env.run();

		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		std_devs.push_back(Analyzer::getStdDevNorm(analytics->total_energy));
	}

	LIMA_Print::printMatlabVec("bond_angle_errors", angle_errors);
	LIMA_Print::printMatlabVec("std_devs", std_devs);

	return true;
}

bool doBasicBenchmark(Environment& env, const string& folder_name) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/" + folder_name + "/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";

	env.loadSimParams(work_folder + "sim_params.txt");
	env.CreateSimulation(conf, topol, work_folder);

	env.run();

	auto analytics = env.getAnalyzedPackage();
	Analyzer::printEnergy(analytics);

	std::vector<float> std_devs;
	std_devs.push_back(Analyzer::getStdDevNorm(analytics->total_energy));


	LIMA_Print::printMatlabVec("std_devs", std_devs);

	return true;
}



#include "EngineUtils.cuh"
int main() {

	/*coordPrecesionBenchmark();
	return 0;*/
	Environment env;

	//basicBenchmark(env);

	//doProteinBenchmark(env);
	//doPoolBenchmark(env, "Pool");			// Two 1-particle molecules colliding
	//doPoolCompSolBenchmark(env, "PoolCompSol");	// One 1-particle molecule colliding with 1 solvent
	//doSpringBenchmark(env);
	//doAngleBenchmark(env);	// Doesn't work currently
	//doBasicBenchmark(env, "TorsionBenchmark");
	//doBasicBenchmark(env, "Met");
	doBasicBenchmark(env, "T4LysozymeNoSolvent");
	//doBasicBenchmark(env, "SolventBenchmark");
	//doBasicBenchmark(env, "T4Lysozyme");
	//doBasicBenchmark(env, "4ake");// TOO big, almost 20 nm long!
	return 0;
}