#pragma once

#include "TestUtils.h"
#include "MDStability.cuh"	//Gives access to the testGenericBox
#include "PhysicsUtils.cuh"
#include "format"


namespace ForceCorrectness {
	using namespace TestUtils;

	//Test assumes two carbons particles in conf
	LimaUnittestResult doPoolBenchmark(EnvMode envmode, float target_vc = 6.45e-5) {
		const fs::path work_folder = simulations_dir / "Pool/";
		Environment env{ work_folder, envmode, false };

		const float particle_mass = 12.011000f / 1000.f;	// kg/mol
		std::vector<float> particle_temps{ 400 };
		//std::vector<float> particle_temps{ 400, 800, 1200 };
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto temp : particle_temps) {
			const float vel = PhysicsUtils::tempToVelocity(temp, particle_mass);	// [m/s] <=> [lm/ls]
			int steps_for_full_interaction = 3000000 / static_cast<int>(vel);

			SimParams params{};
			params.n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
			GroFile grofile{ work_folder / "molecule/conf.gro" };
			TopologyFile topfile{ work_folder / "molecule/topol.top" };
			env.CreateSimulation(grofile, topfile, params);

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

		const auto result = evaluateTest(varcoffs, target_vc, energy_gradients, 2e-7);

		return LimaUnittestResult{ result.first, result.second, envmode == Full};
	}

	LimaUnittestResult doPoolCompSolBenchmark(EnvMode envmode, float max_vc = 9.e-5) {
		const fs::path work_folder = simulations_dir / "PoolCompSol/";
		Environment env{ work_folder, envmode, false };
		SimParams params{ work_folder / "sim_params.txt"};
		const float dt = params.dt;

		//std::vector<float> particle_temps{ 400, 1200, 2400, 4800 };// , 1000, 2000, 5000, 10000
		std::vector<float> particle_temps{ 400, 1200 };
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto temp : particle_temps) {
			// Give the carbon a velocity
			{
				const float particle_mass = 12.011000f / 1000.f;	// kg/mol
				const float vel = PhysicsUtils::tempToVelocity(temp, particle_mass);	// [m/s] <=> [lm/ls]
				const int steps_for_full_interaction = 6000000 / static_cast<int>(vel);

				params.n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
				GroFile grofile{ work_folder / "molecule/conf.gro" };
				TopologyFile topfile{ work_folder / "molecule/topol.top" };
				env.CreateSimulation(grofile, topfile, params);


				Box* box_host = env.getSimPtr()->box_host.get();
				box_host->compounds[0].vels_prev[0] = Float3(1, 0, 0) * vel;
			}

			// Give the solvent a velocty
			{
				const float vel = PhysicsUtils::tempToVelocity(temp, SOLVENT_MASS);	// [m/s] <=> [lm/ls]
				env.getSimPtr()->box_host->solvents[0].vel_prev = Float3{ -1, 0, 0 } * vel;
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

		const auto result = evaluateTest(varcoffs, max_vc, energy_gradients, 2e-7);

		return LimaUnittestResult{ result.first, result.second, envmode == Full };
	}



	LimaUnittestResult SinglebondForceAndPotentialSanityCheck(EnvMode envmode) {
		const fs::path work_folder = simulations_dir / "Singlebond/";
		Environment env{ work_folder, envmode, false };

		SimParams params{ work_folder / "sim_params.txt" };
		params.n_steps = 1;
		params.data_logging_interval = 1;
		const float bondlenErrorNM = 0.02f; //(r-r0) [nm]
		const float bondlenErrorLM = bondlenErrorNM * NANO_TO_LIMA;

		GroFile grofile{ work_folder / "molecule/conf.gro" };
		TopologyFile topfile{ work_folder / "molecule/topol.top" };
		env.CreateSimulation(grofile, topfile, params);

		Box& box_host = *env.getSimPtr()->box_host.get();
		CompoundCoords* coordarray_ptr = box_host.compoundcoordsCircularQueue->getCoordarrayRef(0, 0);
		coordarray_ptr[0].rel_positions[1].x = coordarray_ptr[0].rel_positions[0].x + static_cast<int32_t>(bondlenErrorLM + box_host.compounds[0].singlebonds[0].params.b0);


		// Now figure the expected force and potential
		const double kB = box_host.compounds[0].singlebonds[0].params.kb / 2.; // [J/(mol lm^2)]
		const Float3 dir = (coordarray_ptr[0].rel_positions[1].toFloat3() - coordarray_ptr[0].rel_positions[0].toFloat3()).norm();
		const Float3 expectedForce = dir * 2.f * kB * bondlenErrorLM;			// [1/lima N/mol)]
		const float expectedPotential = kB * bondlenErrorLM * bondlenErrorLM;	// [J/mol]


		env.run();

		const auto sim = env.getSim();
		// Fetch the potE from a buffer. Remember the potE is split between the 2 particles, so we need to sum them here
		const float actualPotE = sim->potE_buffer->getCompoundparticleDatapointAtIndex(0, 0, 0) + sim->potE_buffer->getCompoundparticleDatapointAtIndex(0, 1, 0);
		const Float3 actualForce = sim->box_host->compounds[0].forces_prev[0];


		const float forceError = (actualForce - expectedForce).len() / expectedForce.len();
		ASSERT(forceError < 0.0001f, std::format("Expected force: {:.2e} {:.2e} {:.2e} Actual force: {:.2e} {:.2e} {:.2e} Error: {:.2e}", 
			expectedForce.x, expectedForce.y, expectedForce.z, actualForce.x, actualForce.y, actualForce.z, forceError));

		const float potEError = std::abs(actualPotE - expectedPotential) / expectedPotential;
		ASSERT(potEError < 0.0001f, std::format("Expected potential: {:.2e} Actual potential: {:.2e} Error: {:.2f}", expectedPotential, actualPotE, potEError));

		return LimaUnittestResult{ true, "Success", envmode == Full };
	}


	// Test that a singlebond oscillates at the correct frequency
	LimaUnittestResult SinglebondOscillationTest(EnvMode envmode) {
		const fs::path work_folder = simulations_dir / "Singlebond/";
		const fs::path conf = work_folder / "molecule/conf.gro";
		const fs::path topol = work_folder / "molecule/topol.top";
		const fs::path simpar = work_folder / "sim_params.txt";
		Environment env{ work_folder, envmode, false };

		const float particle_mass = 12.011000f * 1e-3f;

		SimParams params{ simpar };
		// Set time to 1000 [fs]
		params.dt = 100.f;
		params.n_steps = 1000; 
		params.data_logging_interval = 1;
		const float timeElapsed = params.dt * params.n_steps * LIMA_TO_FEMTO; // [fs]
		float bond_len_error = 0.04f ; //(r-r0) [nm]
			

		GroFile grofile{ conf };
		TopologyFile topfile{ topol };
		env.CreateSimulation(grofile, topfile, params);

		Box& box_host = *env.getSimPtr()->box_host.get();
		CompoundCoords* coordarray_ptr = box_host.compoundcoordsCircularQueue->getCoordarrayRef(0, 0);
		coordarray_ptr[0].rel_positions[1].x -= static_cast<int32_t>(bond_len_error * NANO_TO_LIMA + box_host.compounds[0].singlebonds[0].params.b0);




		// Now figure out how fast the bond should oscillate
		const auto atomtypeA = env.getSimPtr()->forcefield.particle_parameters[box_host.compounds[0].atom_types[0]];
		const auto atomtypeB = env.getSimPtr()->forcefield.particle_parameters[box_host.compounds[0].atom_types[1]];

		const double reducedMass = atomtypeA.mass * atomtypeB.mass / (atomtypeA.mass + atomtypeB.mass); // [kg/mol]
		const double kB = box_host.compounds[0].singlebonds[0].params.kb / LIMA / LIMA; // [J/(mol m^2)]

		const double expectedFrequency = sqrt(kB / reducedMass) / (2.f * PI) * FEMTO;	// [1/fs]


		env.run();

		const auto sim = env.getSim();

		std::vector<float> bondlenOverTime(params.n_steps);	// [nm]
		for (int i = 0; i < params.n_steps; i++) {
			bondlenOverTime[i] = (sim->traj_buffer->getCompoundparticleDatapointAtIndex(0, 0, i) - sim->traj_buffer->getCompoundparticleDatapointAtIndex(0, 1, i)).len();
		}

		const int nOscillations = SimAnalysis::CountOscillations(bondlenOverTime);
		const float actualFrequency = static_cast<float>(nOscillations) / timeElapsed; // [1/fs]

		const float error = std::abs(actualFrequency - expectedFrequency) / expectedFrequency;
		const float errorThreshold = 1e-2;

		return LimaUnittestResult{ error < errorThreshold ? true : false, 
			std::format("Expected freq: {:.2e} [1/fs], Actual: {:.2e} [1/fs], Error: {:.2e}", expectedFrequency, actualFrequency, error),
			envmode == Full };
	}









	LimaUnittestResult doSinglebondBenchmark(EnvMode envmode, float max_dev = 0.00746) {
		const fs::path work_folder = simulations_dir / "Singlebond/";
		Environment env{ work_folder, envmode, false };

		const float particle_mass = 12.011000f * 1e-3f;

		SimParams params{ work_folder / "sim_params.txt" };
		params.data_logging_interval = 1;
		params.n_steps = 5000;
		//params.dt = 50.f;
		std::vector<float> bond_len_errors{ 0.02f }; //(r-r0) [nm]
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		//const float bondEquilibrium = 0.149; // [nm]
		const float bondEquilibrium = 0.1335; // [nm]

		for (auto bond_len_error : bond_len_errors) {
			GroFile grofile{ work_folder / "molecule/conf.gro" };
			TopologyFile topfile{ work_folder / "molecule/topol.top" };
			env.CreateSimulation(grofile, topfile, params);

			Box* box_host = env.getSimPtr()->box_host.get();
			CompoundCoords* coordarray_ptr = box_host->compoundcoordsCircularQueue->getCoordarrayRef(0, 0);

			coordarray_ptr[0].rel_positions[1].x += static_cast<int32_t>((bondEquilibrium + bond_len_error) * NANO_TO_LIMA);

			//box_host->compounds[0].singlebonds[0].b0 = bondEquilibrium * NANO_TO_LIMA;
			//box_host->compounds[0].singlebonds[0].kb = 502080 * KILO / (NANO_TO_LIMA * NANO_TO_LIMA);

			env.run();

			const auto analytics = env.getAnalyzedPackage();
			varcoffs.push_back(analytics->variance_coefficient);
			energy_gradients.push_back(analytics->energy_gradient);

			if (envmode != Headless) {
				Analyzer::printEnergy(analytics);
			}

			//LIMA_Print::plotEnergies(analytics->pot_energy, analytics->kin_energy, analytics->total_energy);
			//LIMA_Print::printPythonVec("potE", analytics->pot_energy);
		}

		if (envmode != Headless) {
			LIMA_Print::printMatlabVec("bond_len_errors", bond_len_errors);
			LIMA_Print::printMatlabVec("varcoffs", varcoffs);
			LIMA_Print::printMatlabVec("energy_gradients", energy_gradients);
		}

		const auto result = evaluateTest(varcoffs, max_dev, energy_gradients, 2.e-7);

		return LimaUnittestResult{ result.first, result.second, envmode == Full };
	}

	// Benchmarks anglebonds + singlebonds (for stability)
	LimaUnittestResult doAnglebondBenchmark(EnvMode envmode, float max_vc = 4.7e-3) {
		const fs::path work_folder = simulations_dir / "Anglebond/";

		Environment env{ work_folder, envmode, false};
		SimParams params{ work_folder / "sim_params.txt" };

		const float relaxed_angle = 1.8849f; // [rad]
		std::vector<float> angle_errors{ 0.5f, 0.7f }; //(t-t0) [rad]
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto angle_error : angle_errors) {
			GroFile grofile{ work_folder / "molecule/conf.gro" };
			TopologyFile topfile{ work_folder / "molecule/topol.top" };
			env.CreateSimulation(grofile, topfile, params);

			Box* box_host = env.getSimPtr()->box_host.get();
			CompoundCoords* coordarray_ptr = box_host->compoundcoordsCircularQueue->getCoordarrayRef(0, 0);

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

		const auto result = evaluateTest(varcoffs, max_vc, energy_gradients, 1e-7);

		return LimaUnittestResult{ result.first, result.second, envmode == Full };
	}

	LimaUnittestResult doDihedralbondBenchmark(EnvMode envmode) {
		return TestUtils::loadAndRunBasicSimulation("Dihedralbond", envmode, 5.68e-4, 2.9e-7);
	}

	LimaUnittestResult doImproperDihedralBenchmark(EnvMode envmode, float max_vc=9.7e-3, float max_eg=6.037) {
		const fs::path work_folder = simulations_dir / "Improperbond/";

		Environment env{ work_folder, envmode, false };
		SimParams params{ work_folder / "sim_params.txt" };
		//params.n_steps = 5000;
		//params.data_logging_interval = 1;
		//params.dt = 50.f;
		std::vector<float> angle_errors{ 0.4f, -0.4f, 1.f }; //(t-t0) [rad]
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto angle_error : angle_errors) {
			GroFile grofile{ work_folder / "molecule/conf.gro" };
			TopologyFile topfile{ work_folder / "molecule/topol.top" };
			env.CreateSimulation(grofile, topfile, params);


		/*	env.getSimPtr()->forcefield.particle_parameters[env.getSimPtr()->box_host->compounds[0].atom_types[0]].mass = 100000000.0f;
			env.getSimPtr()->forcefield.particle_parameters[env.getSimPtr()->box_host->compounds[0].atom_types[1]].mass = 100000000.0f;
			env.getSimPtr()->forcefield.particle_parameters[env.getSimPtr()->box_host->compounds[0].atom_types[2]].mass = 100000000.0f;*/

			Box* box_host = env.getSimPtr()->box_host.get();
			CompoundCoords* coordarray_ptr = box_host->compoundcoordsCircularQueue->getCoordarrayRef(0, 0);

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
				//LIMA_Print::plotEnergies(env.getAnalyzedPackage()->pot_energy, env.getAnalyzedPackage()->kin_energy, env.getAnalyzedPackage()->total_energy);

			}
		}

		if (envmode != Headless) {
			//LIMA_Print::printMatlabVec("bond_angle_errors", angle_errors);
			//LIMA_Print::printMatlabVec("varcoffs", varcoffs);
			//LIMA_Print::printMatlabVec("energy_gradients", energy_gradients);
			//LIMA_Print::plotEnergies(env.getAnalyzedPackage()->pot_energy, env.getAnalyzedPackage()->kin_energy, env.getAnalyzedPackage()->total_energy);
		}

		const auto result = evaluateTest(varcoffs, max_vc, energy_gradients, max_eg);

		return LimaUnittestResult{ result.first, result.second, envmode == Full };
	}



}



namespace StressTesting {
	bool doPool50x(EnvMode envmode) {
		const fs::path work_folder = "C:/PROJECTS/Quantom/Simulation/Pool/";

		SimParams params{ work_folder / "sim_params.txt" };
		params.n_steps = 100;

		auto func = [&]() {
			TestUtils::loadAndRunBasicSimulation("Pool", envmode, 0.0001f, 1e-7, params);
		};
		TestUtils::stressTest(func, 50);
		return true;
	}
}


namespace VerletintegrationTesting {
	using namespace TestUtils;

	// Apply a constant force on a particle, and check that the particles achieves the expected kinetic energy
	LimaUnittestResult TestIntegration(EnvMode envmode) {
		const fs::path work_folder = simulations_dir / "Pool/";

		Environment env{ work_folder, envmode, false };

		SimParams params{};
		params.n_steps = 1000;
		params.enable_electrostatics = true;
		params.data_logging_interval = 1;
		params.snf_select = HorizontalChargeField;

		const double timeElapsed = params.dt * static_cast<double>(params.n_steps) * LIMA; // [s]

		GroFile grofile{ work_folder / "molecule/conf.gro" };
		TopologyFile topfile{ work_folder / "molecule/topol.top" };
		grofile.atoms.pop_back();
		topfile.GetLocalAtoms().pop_back();

		env.CreateSimulation(grofile, topfile, params);
		const float electricFieldStrength = .5f ; // [V/nm]
		env.getSimPtr()->box_host->uniformElectricField = UniformElectricField{ {-1, 0, 0 }, electricFieldStrength };



		const float particleCharge = static_cast<float>(env.getSimPtr()->box_host->compounds[0].atom_charges[0]) * KILO; // [C/mol]
		const float particleMass = env.getSimPtr()->forcefield.particle_parameters[1].mass; // [kg/mol]

		const float expectedVelocity = particleCharge * electricFieldStrength / NANO * timeElapsed / particleMass; // [m/s]
		const float expectedKinE = PhysicsUtils::calcKineticEnergy(expectedVelocity, particleMass); // [J/mol]

		env.run();



		const float actualKineticEnergy = env.getAnalyzedPackage()->kin_energy.back();

		const float error = std::abs(actualKineticEnergy - expectedKinE) / expectedKinE;

		ASSERT(error < 0.01f, std::format("Expected KE: {:.2e} Actual KE: {:.2e}", expectedKinE, actualKineticEnergy));


		return LimaUnittestResult{ true, "", envmode == Full};
	}

}