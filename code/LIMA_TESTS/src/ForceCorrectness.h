#pragma once

#include "LIMA_TESTS/src/TestUtils.h"
#include "LIMA_TESTS/src/MDStability.cuh"	//Gives access to the testGenericBox


namespace ForceCorrectness {
	using namespace TestUtils;

	//Test assumes two carbons particles in conf
	LimaUnittest doPoolBenchmark(EnvMode envmode, float max_dev = 0.0001) {
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
		const auto status = result.first == true ? LimaUnittest::SUCCESS : LimaUnittest::FAIL;

		return LimaUnittest{ "doPoolBenchmark", status, result.second, envmode == Full };
	}

	LimaUnittest doPoolCompSolBenchmark(EnvMode envmode, float max_dev = 0.0012) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/PoolCompSol/";
		const std::string conf = work_folder + "molecule/conf.gro";
		const std::string topol = work_folder + "molecule/topol.top";
		Environment env{ work_folder, envmode };
		const auto a = sizeof(SolventBlockGrid);
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
				const int steps_for_full_interaction = 3000000 / static_cast<int>(vel);
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
		const auto status = result.first == true ? LimaUnittest::SUCCESS : LimaUnittest::FAIL;

		return LimaUnittest{ "doPoolCompSolBenchmark", status, result.second, envmode == Full };
	}

	LimaUnittest doSinglebondBenchmark(EnvMode envmode, float max_dev = 0.0031) {
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
			CompoundCoords* coordarray_prev_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, CompoundCoords::firststep_prev, 0);

			coordarray_ptr[0].rel_positions[0].x -= static_cast<int32_t>(bond_len_error * NANO_TO_LIMA);
			coordarray_prev_ptr[0].rel_positions[0].x -= static_cast<int32_t>(bond_len_error * NANO_TO_LIMA);

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

		const auto result = evaluateTest(varcoffs, max_dev, energy_gradients, 1e-7);
		const auto status = result.first == true ? LimaUnittest::SUCCESS : LimaUnittest::FAIL;

		return LimaUnittest{ "doSinglebondBenchmark", status, result.second, envmode == Full };
	}

	// Benchmarks anglebonds + singlebonds (for stability)
	LimaUnittest doAnglebondBenchmark(EnvMode envmode, float max_dev = 0.0007) {
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
			CompoundCoords* coordarray_prev_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, CompoundCoords::firststep_prev, 0);

			// First rotate particle #3 to the relaxed position + the error angle
			Float3 p3_pos = coordarray_ptr[0].rel_positions[2].toFloat3();
			p3_pos.rotateAroundOrigo(Float3{ 0.f, relaxed_angle + angle_error, 0.f });

			coordarray_ptr[0].rel_positions[2] = Coord{ p3_pos };		// Temp disabled, fix soon plz
			coordarray_prev_ptr[0].rel_positions[2] = Coord{ p3_pos };		// Temp disabled, fix soon plz

			// Now center all 3 particles
			for (auto i = 0; i < 3; i++) {
				coordarray_ptr->origo += NodeIndex{ 3, 3, 3 };
				coordarray_prev_ptr->origo += NodeIndex{ 3, 3, 3 };
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
		const auto status = result.first == true ? LimaUnittest::SUCCESS : LimaUnittest::FAIL;

		return LimaUnittest{ "doAnglebondBenchmark", status, result.second, envmode == Full };
	}

	LimaUnittest doDihedralbondBenchmark(EnvMode envmode) {
		return TestUtils::loadAndRunBasicSimulation("TorsionBenchmark", envmode, 0.0006f);
	}

	LimaUnittest doImproperDihedralBenchmark(EnvMode envmode) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/ImproperDihedral/";
		const std::string conf = work_folder + "molecule/conf.gro";
		const std::string topol = work_folder + "molecule/topol.top";
		const std::string simpar = work_folder + "sim_params.txt";

		Environment env{ work_folder, envmode };
		auto ip = env.loadInputSimParams(simpar);

		//std::vector<float> angle_errors{ 0.4f, 1.f }; //(t-t0) [rad]
		std::vector<float> angle_errors{ 1.f }; //(t-t0) [rad]
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto angle_error : angle_errors) {
			env.CreateSimulation(conf, topol, ip);

			Box* box_host = env.getSimPtr()->box_host.get();
			CompoundCoords* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, 0, 0);
			CompoundCoords* coordarray_prev_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box_host->coordarray_circular_queue, CompoundCoords::firststep_prev, 0);

			auto l1 = (coordarray_ptr[0].rel_positions[1] - coordarray_ptr[0].rel_positions[2]).toFloat3().len();

			auto atom_ids = box_host->compounds[0].impropers[0].atom_indexes;
			std::array<Float3, 4> pos;
			for (int i = 0; i < 4; i++) {
				pos[i] = coordarray_ptr[0].rel_positions[atom_ids[i]].toFloat3();
			}

			Float3 i = pos[1];
			Float3 j = pos[0];
			Float3 k = pos[3];
			Float3 l = pos[2];

			// Move i to origo
			j -= i;
			k -= i;
			l -= i;
			i -= i;	// Do this one last!





			Float3 plane_normal = (j-i).cross(k-i).norm();
			Float3 l_vec = (l-i).norm();

			//float angle = Float3::getAngle(plane_normal)

			Float3 rotatevec = (plane_normal.cross(l_vec)).norm();

			Float3 l_point = l / 100000.f;
			//Float3 l_rotated = l_point.rotateAroundVector(Float3{ 0.f,0.f,angle_error }, rotatevec);
			Float3 l_rotated = Float3::rodriguesRotatation(l_point, rotatevec, angle_error);

			Float3 l_diff = l_rotated - l_point;

			coordarray_ptr[0].rel_positions[3] += Coord{ l_diff * 100000.f };		// Temp disabled, fix soon plz
			coordarray_prev_ptr[0].rel_positions[3] += Coord{ l_diff * 100000.f };		// Temp disabled, fix soon plz

			auto l2 = (coordarray_ptr[0].rel_positions[1] - coordarray_ptr[0].rel_positions[2]).toFloat3().len();
			printf("%f\n", l2/l1);

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

		const auto result = evaluateTest(varcoffs, 0.0001, energy_gradients, 1e-5);
		const auto status = result.first == true ? LimaUnittest::SUCCESS : LimaUnittest::FAIL;

		return LimaUnittest{ "doImproperDihedralBenchmark", status, result.second, envmode == Full };
	}

	LimaUnittest doMethionineBenchmark(EnvMode envmode) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Met/";
		const std::string simpar = work_folder + "sim_params.txt";

		return TestUtils::loadAndRunBasicSimulation("Met", envmode, 0.00015f);
	}

	LimaUnittest doPhenylalanineBenchmark(EnvMode envmode) {
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Phe/";
		const std::string simpar = work_folder + "sim_params.txt";

		return TestUtils::loadAndRunBasicSimulation("Phe", envmode, 0.0002f);
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