#pragma once

#include "TestUtils.h"
#include "MDStability.cuh"	//Gives access to the testGenericBox


namespace ForceCorrectness {
	using namespace TestUtils;

	//Test assumes two carbons particles in conf
	LimaUnittestResult doPoolBenchmark(EnvMode envmode, float max_dev = 4e-5f) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Pool/";
		const std::string conf = work_folder + "molecule/conf.gro";
		const std::string topol = work_folder + "molecule/topol.top";
		Environment env{ work_folder, envmode };

		const float particle_mass = 12.011000f / 1000.f;	// kg/mol
		std::vector<float> particle_temps{ 400 };
		//std::vector<float> particle_temps{ 400, 800, 1200 };
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto temp : particle_temps) {
			const float vel = EngineUtils::tempToVelocity(temp, particle_mass);	// [m/s] <=> [lm/ls]
			int steps_for_full_interaction = 3000000 / static_cast<int>(vel);

			InputSimParams ip{};
			ip.n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
			env.CreateSimulation(conf, topol, ip);

			Box* box_host = env.getSimPtr()->box_host.get();
			box_host->compounds[0].vels_prev[0] = Float3(1, 0, 0) * vel;
			box_host->compounds[1].vels_prev[0] = Float3(-1, 0, 0) * vel;

			env.run();

			const auto analytics = env.getAnalyzedPackage();
			varcoffs.push_back(analytics->variance_coefficient);
			energy_gradients.push_back(analytics->energy_gradient);
			if (envmode != Headless) { Analyzer::printEnergy(analytics); }
		}

		if (envmode != Headless) {
			LIMA_Print::printMatlabVec("temperature", particle_temps);
			LIMA_Print::printMatlabVec("varcoffs", varcoffs);
			LIMA_Print::printMatlabVec("energy_gradients", energy_gradients);
		}

		const auto result = evaluateTest(varcoffs, max_dev, energy_gradients, 1e-7);
		const auto status = result.first == true ? LimaUnittestResult::SUCCESS : LimaUnittestResult::FAIL;

		return LimaUnittestResult{status, result.second, envmode == Full};
	}

	LimaUnittestResult doPoolCompSolBenchmark(EnvMode envmode, float max_dev = 1e-4) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/PoolCompSol/";
		const std::string conf = work_folder + "molecule/conf.gro";
		const std::string topol = work_folder + "molecule/topol.top";
		Environment env{ work_folder, envmode };
		auto ip = env.loadInputSimParams(work_folder + "sim_params.txt");
		const float dt = ip.dt;

		//std::vector<float> particle_temps{ 400, 1200, 2400, 4800 };// , 1000, 2000, 5000, 10000
		std::vector<float> particle_temps{ 400, 800, 1200 };
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto temp : particle_temps) {
			// Give the carbon a velocity
			{
				const float particle_mass = 12.011000f / 1000.f;	// kg/mol
				const float vel = EngineUtils::tempToVelocity(temp, particle_mass);	// [m/s] <=> [lm/ls]
				const int steps_for_full_interaction = 6000000 / static_cast<int>(vel);
				InputSimParams ip{};
				ip.n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
				env.CreateSimulation(conf, topol, ip);


				Box* box_host = env.getSimPtr()->box_host.get();
				box_host->compounds[0].vels_prev[0] = Float3(1, 0, 0) * vel;
			}

			// Give the solvent a velocty
			{
				const float vel = EngineUtils::tempToVelocity(temp, SOLVENT_MASS);	// [m/s] <=> [lm/ls]
				env.getSimPtr()->box_host->solvents[0].vel_prev = Float3{ -1, 0, 0 } *vel;
			}


			env.run();

			auto analytics = env.getAnalyzedPackage();
			if (envmode != Headless) {
				Analyzer::printEnergy(analytics);
			}
			
			varcoffs.push_back(analytics->variance_coefficient);
			energy_gradients.push_back(analytics->energy_gradient);
		}

		if (envmode != Headless) {
			LIMA_Print::printMatlabVec("temperature", particle_temps);
			LIMA_Print::printMatlabVec("varcoffs", varcoffs);
		}	

		const auto result = evaluateTest(varcoffs, max_dev, energy_gradients, 2e-7);
		const auto status = result.first == true ? LimaUnittestResult::SUCCESS : LimaUnittestResult::FAIL;

		return LimaUnittestResult{status, result.second, envmode == Full };
	}

	LimaUnittestResult doSinglebondBenchmark(EnvMode envmode, float max_dev = 0.0031) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Spring/";
		const std::string conf = work_folder + "molecule/conf.gro";
		const std::string topol = work_folder + "molecule/topol.top";
		const std::string simpar = work_folder + "sim_params.txt";
		Environment env{ work_folder, envmode };

		const float particle_mass = 12.011000f * 1e-3f;

		InputSimParams ip = env.loadInputSimParams(simpar);


		std::vector<float> bond_len_errors{ 0.02f }; //(r-r0) [nm]
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto bond_len_error : bond_len_errors) {
			env.CreateSimulation(conf, topol, ip);

			Box* box_host = env.getSimPtr()->box_host.get();
			CompoundCoords* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, 0, 0);

			coordarray_ptr[0].rel_positions[0].x -= static_cast<int32_t>(bond_len_error * NANO_TO_LIMA);

			env.run();

			const auto analytics = env.getAnalyzedPackage();
			varcoffs.push_back(analytics->variance_coefficient);
			energy_gradients.push_back(analytics->energy_gradient);

			if (envmode != Headless) {
				Analyzer::printEnergy(analytics);
			}
		}

		if (envmode != Headless) {
			LIMA_Print::printMatlabVec("bond_len_errors", bond_len_errors);
			LIMA_Print::printMatlabVec("varcoffs", varcoffs);
			LIMA_Print::printMatlabVec("energy_gradients", energy_gradients);
		}

		const auto result = evaluateTest(varcoffs, max_dev, energy_gradients, 1.3e-7);
		const auto status = result.first == true ? LimaUnittestResult::SUCCESS : LimaUnittestResult::FAIL;

		return LimaUnittestResult{status, result.second, envmode == Full };
	}

	// Benchmarks anglebonds + singlebonds (for stability)
	LimaUnittestResult doAnglebondBenchmark(EnvMode envmode, float max_dev = 0.0007) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/AngleBenchmark/";
		const std::string conf = work_folder + "molecule/conf.gro";
		const std::string topol = work_folder + "molecule/topol.top";
		const std::string simpar = work_folder + "sim_params.txt";

		Environment env{ work_folder, envmode };
		auto ip = env.loadInputSimParams(simpar);

		const float relaxed_angle = 1.8849f; // [rad]
		std::vector<float> angle_errors{ 0.5f, 0.7f }; //(t-t0) [rad]
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto angle_error : angle_errors) {
			env.CreateSimulation(conf, topol, ip);

			Box* box_host = env.getSimPtr()->box_host.get();
			CompoundCoords* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, 0, 0);

			// First rotate particle #3 to the relaxed position + the error angle
			Float3 p3_pos = coordarray_ptr[0].rel_positions[2].toFloat3();
			p3_pos.rotateAroundOrigo(Float3{ 0.f, relaxed_angle + angle_error, 0.f });

			coordarray_ptr[0].rel_positions[2] = Coord{ p3_pos };		// Temp disabled, fix soon plz

			// Now center all 3 particles
			for (auto i = 0; i < 3; i++) {
				coordarray_ptr->origo += NodeIndex{ 3, 3, 3 };
			}

			env.run();

			const auto analytics = env.getAnalyzedPackage();
			varcoffs.push_back(analytics->variance_coefficient);
			energy_gradients.push_back(analytics->energy_gradient);

			if (envmode != Headless) {
				Analyzer::printEnergy(analytics);
			}
		}

		if (envmode != Headless) {
			LIMA_Print::printMatlabVec("bond_angle_errors", angle_errors);
			LIMA_Print::printMatlabVec("varcoffs", varcoffs);
			LIMA_Print::printMatlabVec("energy_gradients", energy_gradients);
		}

		const auto result = evaluateTest(varcoffs, max_dev, energy_gradients, 1e-7);
		const auto status = result.first == true ? LimaUnittestResult::SUCCESS : LimaUnittestResult::FAIL;

		return LimaUnittestResult{status, result.second, envmode == Full };
	}

	LimaUnittestResult doDihedralbondBenchmark(EnvMode envmode) {
		return TestUtils::loadAndRunBasicSimulation("TorsionBenchmark", envmode, 0.0006f, 2.7e-7);
	}

	LimaUnittestResult doImproperDihedralBenchmark(EnvMode envmode) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/ImproperDihedral/";
		const std::string conf = work_folder + "molecule/conf.gro";
		const std::string topol = work_folder + "molecule/topol.top";
		const std::string simpar = work_folder + "sim_params.txt";

		Environment env{ work_folder, envmode };
		auto ip = env.loadInputSimParams(simpar);

		std::vector<float> angle_errors{ 0.4f, -0.4f, 1.f }; //(t-t0) [rad]
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto angle_error : angle_errors) {
			env.CreateSimulation(conf, topol, ip);

			Box* box_host = env.getSimPtr()->box_host.get();
			CompoundCoords* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, 0, 0);

			auto atom_ids = box_host->compounds[0].impropers[0].atom_indexes;

			Float3 i = coordarray_ptr[0].rel_positions[atom_ids[0]].toFloat3();
			Float3 j = coordarray_ptr[0].rel_positions[atom_ids[1]].toFloat3();
			Float3 k = coordarray_ptr[0].rel_positions[atom_ids[2]].toFloat3();
			Float3 l = coordarray_ptr[0].rel_positions[atom_ids[3]].toFloat3();
			
			// Move i to origo
			j -= i;
			k -= i;
			l -= i;
			i -= i;	// Do this one last



			const Float3 plane_normal = (j-i).cross(k-i).norm();
			const Float3 l_vec = (l-i).norm();

			const Float3 rotatevec = (plane_normal.cross(l_vec)).norm();

			const Float3 l_point = l / l.len();
			const Float3 l_rotated = Float3::rodriguesRotatation(l_point, rotatevec, angle_error);


			Float3 l_diff = (l_rotated - l_point) *l.len();

			coordarray_ptr[0].rel_positions[atom_ids[3]] += Coord{ l_diff};

			env.run();

			const auto analytics = env.getAnalyzedPackage();
			varcoffs.push_back(analytics->variance_coefficient);
			energy_gradients.push_back(analytics->energy_gradient);

			if (envmode != Headless) {
				Analyzer::printEnergy(analytics);
			}
		}

		if (envmode != Headless) {
			LIMA_Print::printMatlabVec("bond_angle_errors", angle_errors);
			LIMA_Print::printMatlabVec("varcoffs", varcoffs);
			LIMA_Print::printMatlabVec("energy_gradients", energy_gradients);
		}

		const auto result = evaluateTest(varcoffs, 0.005, energy_gradients, 6e-5);
		const auto status = result.first == true ? LimaUnittestResult::SUCCESS : LimaUnittestResult::FAIL;

		return LimaUnittestResult{status, result.second, envmode == Full };
	}

	LimaUnittestResult doMethionineBenchmark(EnvMode envmode) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Met/";
		const std::string simpar = work_folder + "sim_params.txt";

		return TestUtils::loadAndRunBasicSimulation("Met", envmode, 4.1e-4, 9e-7);
	}

	LimaUnittestResult doPhenylalanineBenchmark(EnvMode envmode) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Phe/";
		const std::string simpar = work_folder + "sim_params.txt";

		return TestUtils::loadAndRunBasicSimulation("Phe", envmode, 4e-4f, 8e-8f);
	}

}



namespace StressTesting {
	bool doPool50x(EnvMode envmode) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Pool/";
		const std::string simpar = work_folder + "sim_params.txt";

		auto ip = Environment::loadInputSimParams(simpar);
		ip.n_steps = 100;

		auto func = [&]() {
			TestUtils::loadAndRunBasicSimulation("Pool", envmode, 0.0001f, 1e-7, ip, false);
		};
		TestUtils::stressTest(func, 50);
		return true;
	}
}
