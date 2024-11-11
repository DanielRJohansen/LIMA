#pragma once

#include "TestUtils.h"
#include "MDStability.cuh"	//Gives access to the testGenericBox
#include "PhysicsUtils.cuh"
#include "format"


namespace ForceCorrectness {
	using namespace TestUtils;

	//Test assumes two carbons particles in conf
	LimaUnittestResult doPoolBenchmark(EnvMode envmode, float target_vc = 1.e-4) {
		const fs::path work_folder = simulations_dir / "Pool/";
		Environment env{ work_folder, envmode};

		const float particle_mass = 12.011000f / 1000.f;	// kg/mol
		std::vector<float> particle_temps{ 400 };
		//std::vector<float> particle_temps{ 400, 800, 1200 };
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto temp : particle_temps) {
			const float vel = PhysicsUtils::tempToVelocity(temp, particle_mass);	// [m/s] <=> [lm/ls]
			int64_t steps_for_full_interaction = 3000000 / static_cast<int>(vel);

			SimParams params{};
			params.enable_electrostatics = false;
			params.n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
			GroFile grofile{ work_folder / "molecule/conf.gro" };
			TopologyFile topfile{ work_folder / "molecule/topol.top" };
			env.CreateSimulation(grofile, topfile, params);

			Box* box_host = env.getSimPtr()->box_host.get();
			box_host->compoundInterimStates[0].vels_prev[0] = Float3(1, 0, 0) * vel;
			box_host->compoundInterimStates[1].vels_prev[0] = Float3(-1, 0, 0) * vel;

			env.run();

			const auto analytics = env.getAnalyzedPackage();
			varcoffs.push_back(analytics.variance_coefficient);
			energy_gradients.push_back(analytics.energy_gradient);
			if (envmode != Headless) { analytics.Print(); }
		}

		if (envmode != Headless) {
			LIMA_Print::printMatlabVec("temperature", particle_temps);
			LIMA_Print::printMatlabVec("varcoffs", varcoffs);
			LIMA_Print::printMatlabVec("energy_gradients", energy_gradients);
		}

		const auto result = evaluateTest(varcoffs, target_vc, energy_gradients, 2e-7);

		return LimaUnittestResult{ result.first, result.second, envmode == Full};
	}

	LimaUnittestResult doPoolCompSolBenchmark(EnvMode envmode, float max_vc = 1.66e-4) {
		const fs::path work_folder = simulations_dir / "PoolCompSol/";
		Environment env{ work_folder, envmode};
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
				const int64_t steps_for_full_interaction = 6000000 / static_cast<int>(vel);

				params.n_steps = LIMA_UTILS::roundUp(steps_for_full_interaction, 100);
				GroFile grofile{ work_folder / "molecule/conf.gro" };
				TopologyFile topfile{ work_folder / "molecule/topol.top" };
				env.CreateSimulation(grofile, topfile, params);


				Box* box_host = env.getSimPtr()->box_host.get();
				box_host->compoundInterimStates[0].vels_prev[0] = Float3(1, 0, 0) * vel;
			}

			// Give the solvent a velocty
			{	
				const float solventMass = env.getSimPtr()->forcefieldTinymol.types[env.getSimPtr()->box_host->tinyMols[0].tinymolTypeIndex].mass;
				const float vel = PhysicsUtils::tempToVelocity(temp, solventMass);	// [m/s] <=> [lm/ls]
				env.getSimPtr()->box_host->tinyMols[0].vel_prev = Float3{ -1, 0, 0 } * vel;
			}


			env.run();

			auto analytics = env.getAnalyzedPackage();
			if (envmode != Headless) {
				analytics.Print();
			}
			
			varcoffs.push_back(analytics.variance_coefficient);
			energy_gradients.push_back(analytics.energy_gradient);
		}

		if (envmode != Headless) {
			LIMA_Print::printMatlabVec("temperature", particle_temps);
			LIMA_Print::printMatlabVec("varcoffs", varcoffs);
		}	

		const auto result = evaluateTest(varcoffs, max_vc, energy_gradients, 3e-7);

		return LimaUnittestResult{ result.first, result.second, envmode == Full };
	}



	LimaUnittestResult SinglebondForceAndPotentialSanityCheck(EnvMode envmode) {
		const fs::path work_folder = simulations_dir / "Singlebond/";
		Environment env{ work_folder, envmode};

		SimParams params{ work_folder / "sim_params.txt" };
		params.n_steps = 1;
		params.data_logging_interval = 1;
		const float bondlenErrorNM = 0.02f; //(r-r0) [nm]
		const float bondlenErrorLM = bondlenErrorNM * NANO_TO_LIMA;

		GroFile grofile{ work_folder / "molecule/conf.gro" };
		TopologyFile topfile{ work_folder / "molecule/topol.top" };

		grofile.atoms[1].position = grofile.atoms[0].position + Float3{ bondlenErrorNM, 0.f, 0.f }; // Just so we dont get an 0 dist error

		env.CreateSimulation(grofile, topfile, params);

		Box& box_host = *env.getSimPtr()->box_host.get();
		//CompoundCoords* coordarray_ptr = box_host.compoundcoordsCircularQueue->getCoordarrayRef(0, 0);
		CompoundCoords* coordarray_ptr = &box_host.compoundCoordsBuffer[0];
		coordarray_ptr[0].rel_positions[1].x = coordarray_ptr[0].rel_positions[0].x + static_cast<int32_t>(bondlenErrorLM + box_host.compounds[0].singlebonds[0].params.b0);

		// Now figure the expected force and potential
		const double kB = box_host.compounds[0].singlebonds[0].params.kb / 2.; // [J/(mol lm^2)]
		const Float3 dir{ 1,0,0 };
		const Float3 expectedForce = dir * 2.f * kB * bondlenErrorLM;			// [1/lima N/mol)]
		const float expectedPotential = kB * bondlenErrorLM * bondlenErrorLM;	// [J/mol]


		env.run();

		const auto sim = env.getSim();
		// Fetch the potE from a buffer. Remember the potE is split between the 2 particles, so we need to sum them here
		const float actualPotE = sim->potE_buffer->getCompoundparticleDatapointAtIndex(0, 0, 0) + sim->potE_buffer->getCompoundparticleDatapointAtIndex(0, 1, 0);
		const Float3 actualForce = sim->box_host->compoundInterimStates[0].forces_prev[0];


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
		Environment env{ work_folder, envmode};

		const float particle_mass = 12.011000f * 1e-3f;

		SimParams params{ simpar };
		// Set time to 1000 [fs]
		params.dt = 100.f;
		params.n_steps = 1000; 
		params.data_logging_interval = 1;
		const float timeElapsed = params.dt * params.n_steps * LIMA_TO_FEMTO; // [fs]
		float bond_len_error = 0.04f ; //(r-r0) [nm]
			

		GroFile grofile{ conf };
		grofile.atoms[0].position.x += 0.1f; // Just so we dont get an 0 dist error
		TopologyFile topfile{ topol };
		env.CreateSimulation(grofile, topfile, params);

		Box& box_host = *env.getSimPtr()->box_host.get();
		//CompoundCoords* coordarray_ptr = box_host.compoundcoordsCircularQueue->getCoordarrayRef(0, 0);
		CompoundCoords* coordarray_ptr = &box_host.compoundCoordsBuffer[0];
		coordarray_ptr[0].rel_positions[1].x = coordarray_ptr[0].rel_positions[0].x - static_cast<int32_t>(bond_len_error * NANO_TO_LIMA + box_host.compounds[0].singlebonds[0].params.b0);




		// Now figure out how fast the bond should oscillate
		const float massA = box_host.compounds[0].atomMasses[0];
		const float massB = box_host.compounds[0].atomMasses[1];
		const double reducedMass = massA * massB / (massA + massB); // [kg/mol]
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
			std::format("Expected freq: {:.2e} [1/fs], Actual: {:.2e} [1/fs]", expectedFrequency, actualFrequency, error),
			envmode == Full };
	}









	LimaUnittestResult doSinglebondBenchmark(EnvMode envmode, float max_dev = 0.00746) {
		const fs::path work_folder = simulations_dir / "Singlebond/";
		Environment env{ work_folder, envmode};

		SimParams params{ work_folder / "sim_params.txt" };
		params.data_logging_interval = 1;
		params.n_steps = 5000;
		std::vector<float> bond_len_errors{ 0.02f }; //(r-r0) [nm]
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		//const float bondEquilibrium = 0.149; // [nm]
		const float bondEquilibrium = 0.1335; // [nm]

		for (auto bond_len_error : bond_len_errors) {
			GroFile grofile{ work_folder / "molecule/conf.gro" };
			grofile.atoms[1].position.x += bondEquilibrium + bond_len_error;

			TopologyFile topfile{ work_folder / "molecule/topol.top" };
			env.CreateSimulation(grofile, topfile, params);

			Box* box_host = env.getSimPtr()->box_host.get();

			env.run();

			const auto analytics = env.getAnalyzedPackage();
			varcoffs.push_back(analytics.variance_coefficient);
			energy_gradients.push_back(analytics.energy_gradient);

			if (envmode != Headless) {
				analytics.Print();
			}

			//LIMA_Print::plotEnergies(analytics.pot_energy, analytics.kin_energy, analytics.total_energy);
			//LIMA_Print::printPythonVec("potE", analytics.pot_energy);
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

		Environment env{ work_folder, envmode};
		SimParams params{ work_folder / "sim_params.txt" };

		const float relaxed_angle = 1.8849f; // [rad]
		std::vector<float> angle_errors{ 0.5f, 0.7f }; //(t-t0) [rad]
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto angle_error : angle_errors) {
			GroFile grofile{ work_folder / "molecule/conf.gro" };
			Float3& p3Pos = grofile.atoms[2].position;
			p3Pos.rotateAroundOrigo(Float3{ 0.f, relaxed_angle + angle_error, 0.f });
			for (auto& atom : grofile.atoms)
				atom.position += Float3{ 3.f };

			TopologyFile topfile{ work_folder / "molecule/topol.top" };
			env.CreateSimulation(grofile, topfile, params);

			env.run();

			const auto analytics = env.getAnalyzedPackage();
			varcoffs.push_back(analytics.variance_coefficient);
			energy_gradients.push_back(analytics.energy_gradient);

			if (envmode != Headless) {
				analytics.Print();
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

		Environment env{ work_folder, envmode};
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



			auto atom_ids = topfile.GetMoleculeType().improperdihedralbonds[0].ids;			
			Float3 i = grofile.atoms[atom_ids[0]].position;
			Float3 j = grofile.atoms[atom_ids[1]].position;
			Float3 k = grofile.atoms[atom_ids[2]].position;
			Float3 l = grofile.atoms[atom_ids[3]].position;

			// Move i to origo
			j -= i;
			k -= i;
			l -= i;
			i -= i;	// Do this one last

			const Float3 plane_normal = (j - i).cross(k - i).norm();
			const Float3 l_vec = (l - i).norm();

			const Float3 rotatevec = (plane_normal.cross(l_vec)).norm();

			const Float3 l_point = l / l.len();
			const Float3 l_rotated = Float3::rodriguesRotatation(l_point, rotatevec, angle_error);


			Float3 l_diff = (l_rotated - l_point) * l.len();

			//coordarray_ptr[0].rel_positions[atom_ids[3]] += Coord{ l_diff };
			grofile.atoms[atom_ids[3]].position += l_diff;


			env.CreateSimulation(grofile, topfile, params);

			Box* box_host = env.getSimPtr()->box_host.get();
			CompoundCoords* coordarray_ptr = &box_host->compoundCoordsBuffer[0];

			env.run();

			const auto analytics = env.getAnalyzedPackage();
			varcoffs.push_back(analytics.variance_coefficient);
			energy_gradients.push_back(analytics.energy_gradient);

			if (envmode != Headless) {
				analytics.Print();
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




	void TestForces1To1(EnvMode envmode) {
		const fs::path workDir = simulations_dir / "T4Lysozyme";
		GroFile grofile{ workDir / "molecule" / "conf.gro" };
		TopologyFile topfile{ workDir / "molecule" / "topol.top" };
		SimParams params{ workDir / "sim_params.txt" };
		Environment env{ workDir, envmode };
		params.n_steps = 11;
		params.data_logging_interval = 1;
		env.CreateSimulation(grofile, topfile, params);
		env.run();

		const ParticleDataBuffer<Float3>* forcebuffer = env.getSimPtr()->forceBuffer.get();
		std::vector<Float3> forces(forcebuffer->GetBufferAtStep(0), forcebuffer->GetBufferAtStep(0) + forcebuffer->n_particles_upperbound);
		//FileUtils::WriteVectorToBinaryFile(workDir/ "forces.bin", forces);

		const std::vector<Float3> forcesRef = FileUtils::ReadBinaryFileIntoVector<Float3>(workDir / "forces.bin");

		std::vector<float> errors(forces.size()); // Pre-allocate the vector
		std::transform(forces.begin(), forces.end(), forcesRef.begin(), errors.begin(),
			[](const Float3& a, const Float3& b) { return (a - b).len(); });



		FileUtils::WriteVectorToBinaryFile(workDir / "errors.bin", errors);

		std::string command = (FileUtils::GetLimaDir() / "dev/PyTools/pdf.py").string() + " \"" + (workDir / "errors.bin").string() + "\"";
		std::system(command.c_str());
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

		Environment env{ work_folder, envmode};

		SimParams params{};
		params.n_steps = 1000;
		params.enable_electrostatics = true;
		params.data_logging_interval = 1;
		params.snf_select = HorizontalChargeField;

		const double timeElapsed = params.dt * static_cast<double>(params.n_steps) * LIMA; // [s]

		GroFile grofile{ work_folder / "molecule/conf.gro" };
		TopologyFile topfile{ work_folder / "molecule/topol.top" };
		grofile.atoms.pop_back();
		topfile.GetMoleculeType().atoms.pop_back();

		env.CreateSimulation(grofile, topfile, params);
		const float electricFieldStrength = .5f ; // [V/nm]
		env.getSimPtr()->box_host->uniformElectricField = UniformElectricField{ Float3{-1.f, 0.f, 0.f }, electricFieldStrength };



		const float particleCharge = static_cast<float>(env.getSimPtr()->box_host->compounds[0].atom_charges[0]) * KILO; // [C/mol]
		const float particleMass = env.getSimPtr()->box_host->compounds[0].atomMasses[0]; // [kg/mol]

		const float expectedVelocity = particleCharge * electricFieldStrength / NANO * timeElapsed / particleMass; // [m/s]
		const float expectedKinE = PhysicsUtils::calcKineticEnergy(expectedVelocity, particleMass); // [J/mol]

		env.run();



		const float actualKineticEnergy = env.getAnalyzedPackage().kin_energy.back();

		const float error = std::abs(actualKineticEnergy - expectedKinE) / expectedKinE;

		ASSERT(error < 0.01f, std::format("Expected KE: {:.2e} Actual KE: {:.2e}", expectedKinE, actualKineticEnergy));


		return LimaUnittestResult{ true, "", envmode == Full};
	}

}
