#pragma once

#include "LIMA_TESTS/src/TestUtils.h"
#include "LIMA_TESTS/src/MDStability.cuh"	//Gives access to the testGenericBox


//Test assumes two carbons particles in conf
bool doPoolBenchmark(Environment::Mode envmode, float max_dev=0.007) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Pool/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";
	Environment env{ work_folder, envmode };

	const float particle_mass = 12.011000f / 1000.f;	// kg/mol
	std::vector<float> particle_temps{ 400 };
	//std::vector<float> particle_temps{ 400, 800, 1200 };
	std::vector<float> std_devs;

	for (auto temp : particle_temps) {
		const float vel = EngineUtils::tempToVelocity(temp, particle_mass);	// [m/s] <=> [lm/ls]
		int steps_for_full_interaction = 2000000 / static_cast<int>(vel);
		//sim_params->n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
		InputSimParams ip{};
		ip.n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
		ip.n_steps = 50;
		env.CreateSimulation(conf, topol, ip);

		Box* box_host = env.getSimPtr()->box_host.get();
		CompoundCoords* coordarray_prev_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, CompoundCoords::firststep_prev, 0);
		coordarray_prev_ptr[0].rel_positions[0] += Coord{ (Float3(-1, 0, 0) * vel) * ip.dt };
		coordarray_prev_ptr[1].rel_positions[0] += Coord{ (Float3(1, 0, 0) * vel) * ip.dt };

		env.run();

		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		std_devs.push_back(Analyzer::getVarianceCoefficient(analytics->total_energy));
	}

	LIMA_Print::printMatlabVec("temperature", particle_temps);
	LIMA_Print::printMatlabVec("std_devs", std_devs);

	for (auto& stddev : std_devs) {
		if (stddev > max_dev) { return false; }
	}
	return true;
}

bool doPoolCompSolBenchmark(Environment::Mode envmode, float max_dev = 0.01) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/PoolCompSol/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";
	Environment env{ work_folder, envmode };
	const auto a = sizeof(SolventBlockGrid);
	auto ip = env.loadInputSimParams(work_folder + "sim_params.txt");
	const float dt = ip.dt;

	//std::vector<float> particle_temps{ 400, 1200, 2400, 4800 };// , 1000, 2000, 5000, 10000
	std::vector<float> particle_temps{ 400, 800, 1200 };
	std::vector<float> std_devs;

	for (auto temp : particle_temps) {
		// Give the carbon a velocity
		{
			const float particle_mass = 12.011000f / 1000.f;	// kg/mol
			const float vel = EngineUtils::tempToVelocity(temp, particle_mass);	// [m/s] <=> [lm/ls]
			const int steps_for_full_interaction = 2000000 / static_cast<int>(vel);
			InputSimParams ip{};
			ip.n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
			ip.n_steps = 20;
			env.CreateSimulation(conf, topol, ip);


			Box* box_host = env.getSimPtr()->box_host.get();
			CompoundCoords* coordarray_prev_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, CompoundCoords::firststep_prev, 0);
			coordarray_prev_ptr[0].rel_positions[0] += Coord{ (Float3(-1, 0, 0) * vel) * ip.dt };
		}

		// Give the solvent a velocty
		{
			// We need to also give the solvent velocity. But we don't know which block it is in....
			const float vel = EngineUtils::tempToVelocity(temp, SOLVENT_MASS);	// [m/s] <=> [lm/ls]
			auto solventblockgridprev = env.getSolventBlocksPrevRef();
			for (int i = 0; i < solventblockgridprev->blocks_total; i++) {
				auto block = solventblockgridprev->getBlockPtr(i);
				if (block->n_solvents == 1) {
					block->rel_pos[0] += Coord{ Float3{1, 0, 0} *vel * dt };
				}
			}
		}


		env.run();

		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		std_devs.push_back(Analyzer::getVarianceCoefficient(analytics->total_energy));
	}

	LIMA_Print::printMatlabVec("temperature", particle_temps);
	LIMA_Print::printMatlabVec("std_devs", std_devs);

	for (auto& stddev : std_devs) {
		if (stddev > max_dev) { return false; }
	}
	return true;
}

bool doSinglebondBenchmark(Environment::Mode envmode, float max_dev = 0.1) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Spring/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";
	const std::string simpar = work_folder + "sim_params.txt";
	Environment env{ work_folder, envmode };

	const float particle_mass = 12.011000f * 1e-3f;

	InputSimParams ip = env.loadInputSimParams(simpar);
	ip.n_steps = 100;


	std::vector<float> bond_len_errors{ 0.01f, 0.02f }; //(r-r0) [nm]
	std::vector<float> std_devs;

	for (auto bond_len_error : bond_len_errors) {
		env.CreateSimulation(conf, topol, ip);

		Box* box_host = env.getSimPtr()->box_host.get();
		CompoundCoords* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, 0, 0);
		CompoundCoords* coordarray_prev_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, CompoundCoords::firststep_prev, 0);

		coordarray_ptr[0].rel_positions[0].x -= static_cast<int32_t>(bond_len_error * NANO_TO_LIMA);
		coordarray_prev_ptr[0].rel_positions[0].x -= static_cast<int32_t>(bond_len_error * NANO_TO_LIMA);

		env.run();

		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		std_devs.push_back(Analyzer::getVarianceCoefficient(analytics->total_energy));
	}

	LIMA_Print::printMatlabVec("bond_len_errors", bond_len_errors);
	LIMA_Print::printMatlabVec("std_devs", std_devs);

	for (auto& stddev : std_devs) {
		if (stddev > max_dev) { 
			std::cout << std::format("Stddev of {} superceeded the max of {}", stddev, max_dev);
			return false; 
		}
	}
	return true;
}

// Benchmarks anglebonds + singlebonds (for stability)
bool doAnglebondBenchmark(Environment::Mode envmode, float max_dev = 0.01) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/AngleBenchmark/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";
	const std::string simpar = work_folder + "sim_params.txt";

	Environment env{ work_folder, envmode };
	auto ip =  env.loadInputSimParams(simpar);
	ip.n_steps = 100;

	const float relaxed_angle = 1.8849f; // [rad]
	std::vector<float> angle_errors{ 0.5f }; //(t-t0) [rad]
	std::vector<float> std_devs;

	for (auto angle_error : angle_errors) {
		env.CreateSimulation(conf, topol, ip);

		Box* box_host = env.getSimPtr()->box_host.get();
		CompoundCoords* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, 0, 0);
		CompoundCoords* coordarray_prev_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, CompoundCoords::firststep_prev, 0);

		// First rotate particle #3 to the relaxed position + the error angle
		Float3 p3_pos = coordarray_ptr[0].rel_positions[2].toFloat3();
		p3_pos.rotateAroundOrigo(Float3{ 0.f, relaxed_angle + angle_error, 0.f });
		
		coordarray_ptr[0].rel_positions[2] = Coord{ p3_pos };		// Temp disabled, fix soon plz
		coordarray_prev_ptr[0].rel_positions[2] = Coord{ p3_pos };		// Temp disabled, fix soon plz

		// Now center all 3 particles
		for (auto i = 0; i < 3; i++) {			
			coordarray_ptr->origo += NodeIndex{3, 3, 3};
			coordarray_prev_ptr->origo += NodeIndex{ 3, 3, 3 };
		}

		env.run();

		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		std_devs.push_back(Analyzer::getVarianceCoefficient(analytics->total_energy));
	}

	LIMA_Print::printMatlabVec("bond_angle_errors", angle_errors);
	LIMA_Print::printMatlabVec("std_devs", std_devs);

	for (auto& stddev : std_devs) {
		if (stddev > max_dev) { return false; }
	}
	return true;
}

bool doDihedralbondBenchmark(Environment::Mode envmode) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/TorsionBenchmark/";
	const std::string simpar = work_folder + "sim_params.txt";

	auto ip = Environment::loadInputSimParams(simpar);
	ip.n_steps = 100;

	return TestUtils::loadAndRunBasicSimulation("TorsionBenchmark", envmode, 0.01, ip);
}

bool doMethionineBenchmark(Environment::Mode envmode) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Met/";
	const std::string simpar = work_folder + "sim_params.txt";

	auto ip = Environment::loadInputSimParams(simpar);
	ip.n_steps = 10;

	return TestUtils::loadAndRunBasicSimulation("Met", envmode, 0.01, ip);
}

