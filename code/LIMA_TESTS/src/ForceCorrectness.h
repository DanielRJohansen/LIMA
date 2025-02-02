#pragma once

#include "TestUtils.h"
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
			const float vel = PhysicsUtils::tempToVelocity(temp, particle_mass);	// [m/s] <=> [nm/ns]
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
		params.data_logging_interval = 1;

		//std::vector<float> particle_temps{ 400, 1200, 2400, 4800 };// , 1000, 2000, 5000, 10000
		std::vector<float> particle_temps{ 400, 1200 };
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		for (auto temp : particle_temps) {
			// Give the carbon a velocity
			{
				const float particle_mass = 12.011000f / 1000.f;	// kg/mol
				const float vel = PhysicsUtils::tempToVelocity(temp, particle_mass);	// [m/s] <=> [nm/ns]
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
				float solventMass = 0;
				for (int i = 0; i < env.getSimPtr()->box_host->boxparams.nTinymolParticles; i++) {
					solventMass += env.getSimPtr()->forcefieldTinymol.types[env.getSimPtr()->box_host->tinyMolParticlesState[i].tinymolTypeIndex].mass;
				}
				const float vel = PhysicsUtils::tempToVelocity(temp, solventMass);	// [m/s] <=> [nm/ns]
				for (int i = 0; i < env.getSimPtr()->box_host->boxparams.nTinymolParticles; i++) {
					env.getSimPtr()->box_host->tinyMolParticlesState[i].vel_prev = Float3{ -vel, 0.f, 0.f };
				}
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

		GroFile grofile{ work_folder / "molecule/conf.gro" };
		TopologyFile topfile{ work_folder / "molecule/topol.top" };

		grofile.atoms[1].position = grofile.atoms[0].position + Float3{ bondlenErrorNM, 0.f, 0.f }; // Just so we dont get an 0 dist error as we load the simulation

		env.CreateSimulation(grofile, topfile, params);

		Box& box_host = *env.getSimPtr()->box_host.get();

		const SingleBond::Parameters bondparams = box_host.bondgroups[0].singlebonds[0].params;

		// Now we have the bond params, set the actual test position
		CompoundCoords* coordarray_ptr = &box_host.compoundCoordsBuffer[0];
		coordarray_ptr->rel_positions[1].x = coordarray_ptr->rel_positions[0].x + Coord{ Float3{bondlenErrorNM + bondparams.b0, 0.f, 0.f} }.x;

		// Now figure the expected force and potential
		const double kB = bondparams.kb / 2.;									// [J/mol/nm^2]
		const Float3 dir{ 1,0,0 };
		const Float3 expectedForce = dir * 2.f * kB * bondlenErrorNM;			// [J/mol/nm)]
		const float expectedPotential = kB * bondlenErrorNM * bondlenErrorNM;	// [J/mol]


		env.run();
		LIMA_UTILS::genericErrorCheck("Error during test");

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
		params.dt = 1.f * FEMTO_TO_NANO;
		params.n_steps = 1000; 
		params.data_logging_interval = 1;
		const float timeElapsed = params.dt * params.n_steps * NANO_TO_FEMTO; // [fs]
		float bond_len_error = 0.04f ; //(r-r0) [nm]
			

		GroFile grofile{ conf };
		grofile.atoms[0].position.x += 0.1f; // Just so we dont get an 0 dist error
		TopologyFile topfile{ topol };
		env.CreateSimulation(grofile, topfile, params);

		Box& box_host = *env.getSimPtr()->box_host.get();

		const SingleBond::Parameters bondparams = box_host.bondgroups[0].singlebonds[0].params;

		CompoundCoords* coordarray_ptr = &box_host.compoundCoordsBuffer[0];
		coordarray_ptr[0].rel_positions[1].x = coordarray_ptr[0].rel_positions[0].x - Coord{ Float3{bond_len_error + bondparams.b0, 0.f, 0.f } }.x;




		// Now figure out how fast the bond should oscillate
		const float massA = box_host.compounds[0].atomMasses[0];
		const float massB = box_host.compounds[0].atomMasses[1];
		const double reducedMass = massA * massB / (massA + massB); // [kg/mol]
		const double kB = bondparams.kb / NANO / NANO; // [J/(mol m^2)]

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

	LimaUnittestResult UreyBradleyForceAndPotentialSanityCheck(EnvMode envmode) {
		const fs::path work_folder = simulations_dir / "Anglebond/";
		Environment env{ work_folder, envmode };

		SimParams params{ work_folder / "sim_params.txt" };
		params.n_steps = 1;
		params.data_logging_interval = 1;

		// Deviation from equilibrium values
		
		//const float ubDistErrorNM = 0.02f;   // [nm]

		GroFile grofile{ work_folder / "molecule/conf.gro" };
		TopologyFile topfile{ work_folder / "molecule/topol.top" };

		// Initialize positions to avoid zero distance errors (temporary setup)
		grofile.atoms[0].position = Float3{ -0.13f, 0.0f, 0.0f };
		grofile.atoms[1].position = Float3{ 0.f, 0.0f, 0.0f };
		grofile.atoms[2].position = Float3{ 0.13f, 0.0f, 0.0f };

		env.CreateSimulation(grofile, topfile, params);

		Box& box_host = *env.getSimPtr()->box_host.get();
		
		//box_host.bondgroups[0].anglebonds[0].params.kTheta = 0.f;
		//box_host.bondgroups[0].anglebonds[0].params.kUB = 0.f;
		//box_host.bondgroups[0].nSinglebonds = 0; // Shouldn't be necessary..

		// First equilibrilize the singlebond
		{
			const SingleBond::Parameters bondParams = box_host.bondgroups[0].singlebonds[0].params;
			CompoundCoords* coordarray_ptr = &box_host.compoundCoordsBuffer[0];
			coordarray_ptr->rel_positions[0] = coordarray_ptr->rel_positions[1] + Coord{ Float3{ bondParams.b0, 0.0f, 0.0f } };
			coordarray_ptr->rel_positions[2] = coordarray_ptr->rel_positions[1] + Coord{ Float3{ bondParams.b0, 0.0f, 0.0f } };
		}

		// Now set the angle error
		const float angleErrorRad = 0.1f;    // [radians]
		const AngleUreyBradleyBond::Parameters angleparams = box_host.bondgroups[0].anglebonds[0].params;
		{
			const Float3 p2Pos = box_host.compoundCoordsBuffer[0].rel_positions[2].ToRelpos();
			const Float3 p2Rotated = Float3::rodriguesRotatation(p2Pos, Float3{ 0.f, 1.f, 0.f }, -(angleparams.theta0 + angleErrorRad));
			box_host.compoundCoordsBuffer[0].rel_positions[2] = Coord{ p2Rotated };
		}
		

		// Now calculate expected forces and potential energy

		const Float3 p0 = box_host.compoundCoordsBuffer[0].rel_positions[0].ToRelpos();
		const Float3 p1 = box_host.compoundCoordsBuffer[0].rel_positions[1].ToRelpos();
		const Float3 p2 = box_host.compoundCoordsBuffer[0].rel_positions[2].ToRelpos();

		//ASSERT((p1 - p0).len() == box_host.bondgroups[0].singlebonds[0].params.b0, std::format("Singlebondlen not as expected {}/{}", (p1 - p0).len(), box_host.bondgroups[0].singlebonds[0].params.b0));
		//ASSERT((p2 - p1).len() == box_host.bondgroups[0].singlebonds[0].params.b0, std::format("Anglebondlen not as expected {}/{}", (p2 - p1).len(), box_host.bondgroups[0].singlebonds[0].params.b0));

		// Angular component
		const float potAngle = angleparams.kTheta * angleErrorRad * angleErrorRad * 0.5f;		// Energy [J/mol]			
		const float torque = angleparams.kTheta * (angleErrorRad);				// Torque [J/(mol*rad)]
		// We know p0 will point directly up
		const Float3 forceAngle = Float3{ 0.f,0.f,1.f } *(torque / (p0 - p1).len());



		// Urey-Bradley component
		const float error = (p0-p2).len() - angleparams.ub0;
		ASSERT(error > 0.001, "UB error too small for test to be meaningful");
		const Float3 forceUB = (p0-p2).norm() * -angleparams.kUB * error;
		const float potUB = angleparams.kUB * error * error * 0.5f;


		// Total expected force on the middle atom (atom 1)
		const Float3 expectedForce = forceAngle + forceUB;
		const float expectedPotential = potAngle + potUB;

		env.run();
		LIMA_UTILS::genericErrorCheck("Error during test");

		const auto sim = env.getSim();

		// Fetch the potential energy from the buffer, summing over all three atoms
		const float actualPotE =
			sim->potE_buffer->getCompoundparticleDatapointAtIndex(0, 0, 0) +
			sim->potE_buffer->getCompoundparticleDatapointAtIndex(0, 1, 0) +
			sim->potE_buffer->getCompoundparticleDatapointAtIndex(0, 2, 0);

		// Fetch the actual force on the middle atom (atom 1)
		const Float3 actualForce = sim->box_host->compoundInterimStates[0].forces_prev[0];

		// Validate force and potential
		const float forceError = (actualForce - expectedForce).len() / expectedForce.len();
		ASSERT(forceError < 0.0001f, std::format(
			"Expected force: {:.2e} {:.2e} {:.2e} Actual force: {:.2e} {:.2e} {:.2e} Error: {:.2e}",
			expectedForce.x, expectedForce.y, expectedForce.z,
			actualForce.x, actualForce.y, actualForce.z, forceError));

		const float potEError = std::abs(actualPotE - expectedPotential) / expectedPotential;
		ASSERT(potEError < 0.0001f, std::format(
			"Expected potential: {:.2e} Actual potential: {:.2e} Error: {:.2e}",
			expectedPotential, actualPotE, potEError));

		return LimaUnittestResult{ true, "Success", envmode == Full };
	}


	LimaUnittestResult PairbondForceAndPotentialSanityCheck(EnvMode envmode) {
		const fs::path work_folder = simulations_dir / "Pairbond/";
		Environment env{ work_folder, envmode };

		SimParams params{ work_folder / "sim_params.txt" };
		params.n_steps = 1;
		params.data_logging_interval = 1;

		GroFile grofile{ work_folder / "molecule/conf.gro" };
		TopologyFile topfile{ work_folder / "molecule/topol.top" };

		env.CreateSimulation(grofile, topfile, params);

		Box& box_host = *env.getSimPtr()->box_host.get();

		ASSERT(box_host.bondgroups[0].nPairbonds == 1, std::format("Expected sim to contain 1 pairbond, found {}", box_host.bondgroups[0].nPairbonds));

		// Neutralize singlebonds
		box_host.bondgroups[0].nSinglebonds = 0;
		ASSERT(box_host.bondgroups[0].nAnglebonds == 0, "Expected 0 anglebonds");
		box_host.bondgroups[0].nDihedralbonds = 0;
		ASSERT(box_host.bondgroups[0].nImproperdihedralbonds == 0, "Expected 0 improperdihedralbonds");


		// Now calculate expected forces and potential energy
		const int p0 = box_host.bondgroups[0].particles[box_host.bondgroups[0].pairbonds[0].atom_indexes[0]].localIdInCompound;
		const int p1 = box_host.bondgroups[0].particles[box_host.bondgroups[0].pairbonds[0].atom_indexes[1]].localIdInCompound;

		const Float3 pos0 = box_host.compoundCoordsBuffer[0].rel_positions[p0].ToRelpos();
		const Float3 pos1 = box_host.compoundCoordsBuffer[0].rel_positions[p1].ToRelpos();
		const Float3 diff = pos1 - pos0;

		const PairBond::Parameters bondparams = box_host.bondgroups[0].pairbonds[0].params;
		const float s = std::powf(bondparams.sigma / (diff.len()), 6.f);
		float force_scalar = 24.f * bondparams.epsilon * s / diff.lenSquared() * (1.f - 2.f * s);
		const Float3 expectedForce = diff * force_scalar;
		const float expectedPotential = 4.f * bondparams.epsilon * s * (s - 1.f) * 0.5f;

		env.run();
		LIMA_UTILS::genericErrorCheck("Error during test");

		const auto sim = env.getSim();

		// Fetch the potential energy from the buffer, summing over all three atoms
		const float actualPotE = sim->potE_buffer->getCompoundparticleDatapointAtIndex(0, p0, 0);

		// Fetch the actual force on the middle atom (atom 1)
		const Float3 actualForce = sim->box_host->compoundInterimStates[0].forces_prev[p0];

		ASSERT(expectedForce.len() > 1.f, "Force too small for test to be meaningful");
		ASSERT(expectedPotential > 1.f, "Potential energy too small for test to be meaningful");

		// Validate force and potential
		const float forceError = (actualForce - expectedForce).len() / expectedForce.len();
		ASSERT(forceError < 0.0001f, std::format(
			"Expected force: {:.2e} {:.2e} {:.2e} Actual force: {:.2e} {:.2e} {:.2e} Error: {:.2e}",
			expectedForce.x, expectedForce.y, expectedForce.z,
			actualForce.x, actualForce.y, actualForce.z, forceError));

		const float potEError = std::abs(actualPotE - expectedPotential) / expectedPotential;
		ASSERT(potEError < 0.0001f, std::format(
			"Expected potential: {:.2e} Actual potential: {:.2e} Error: {:.2e}",
			expectedPotential, actualPotE, potEError));

		return LimaUnittestResult{ true, "Success", envmode == Full };
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

		TestUtils::CompareForces1To1(workDir, env, false);
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

		const double timeElapsed = params.dt * static_cast<double>(params.n_steps) * NANO; // [s]

		GroFile grofile{ work_folder / "molecule/conf.gro" };
		TopologyFile topfile{ work_folder / "molecule/topol.top" };
		grofile.atoms.pop_back();
		topfile.GetMoleculeType().atoms.pop_back();

		env.CreateSimulation(grofile, topfile, params);
		const float electricFieldStrength = .5f ; // [V/nm]
		env.getSimPtr()->box_host->uniformElectricField = UniformElectricField{ Float3{1.f, 0.f, 0.f }, electricFieldStrength };



		const float particleCharge = static_cast<float>(env.getSimPtr()->box_host->compounds[0].atom_charges[0]) * KILO; // [C/mol]
		const float particleMass = env.getSimPtr()->box_host->compounds[0].atomMasses[0]; // [kg/mol]

		const float expectedVelocity = particleCharge * electricFieldStrength / NANO * timeElapsed / particleMass; // [m/s]
		const float expectedKinE = PhysicsUtils::calcKineticEnergy(expectedVelocity, particleMass); // [J/mol]

		env.run();


		const float actualKineticEnergy = env.getAnalyzedPackage().kin_energy.back();
		const float error = std::abs(actualKineticEnergy - expectedKinE) / expectedKinE;
		ASSERT(error < 0.01f, std::format("Expected KE: {:.2e} Actual KE: {:.2e}", expectedKinE, actualKineticEnergy));

		{
			//const Float3 pos0 = env.getSimPtr()->traj_buffer->GetMostRecentCompoundparticleDatapoint(0, 0, 0);
			//const Float3 pos1 = env.getSimPtr()->traj_buffer->GetMostRecentCompoundparticleDatapoint(0, 0, params.n_steps - 1);
			//const float distanceTraveled = (pos1 - pos0).len();						// [nm]
			//const float expectedDistance = expectedVelocity * (timeElapsed / NANO); // [nm]
			//const float error = std::abs(distanceTraveled - expectedDistance) / expectedDistance;
			//ASSERT(error < 0.01f, std::format("Expected distance: {:.2e} Actual distance: {:.2e}", expectedDistance, distanceTraveled));
		}


		return LimaUnittestResult{ true, "", envmode == Full};
	}

}
