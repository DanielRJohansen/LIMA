#pragma once

#include "TestUtils.h"
#include "Environment.h"
#include "LimaTypes.cuh"
#include "LimaPositionSystem.cuh"
#include "Statistics.h"
#include "PhysicsUtils.cuh"
#include "PlotUtils.h"

#include <iostream>
#include <string>
#include <map>

namespace ElectrostaticsTests {
	using namespace TestUtils;


	static LimaUnittestResult CoulombForceSanityCheck(EnvMode envmode) {
		const float calcedForce = PhysicsUtils::CalcCoulumbForce(1.f*elementaryChargeToKiloCoulombPerMole, 1.f*elementaryChargeToKiloCoulombPerMole, Float3{ 1.f, 0.f, 0.f }).len(); // [1/l N / mol]
		const float expectedForce = 2.307078e-10 * AVOGADROSNUMBER * NANO;  // [J/mol/nm] https://www.omnicalculator.com/physics/coulombs-law

		ASSERT(std::abs(calcedForce - expectedForce) / expectedForce < 0.0001f, std::format("Expected {:.2e} Actual {:.2e}", expectedForce, calcedForce));
		// TODO: add potE to this also
		return LimaUnittestResult{ true, "Success", envmode == Full};
	}

	//static ForceEnergy CalcImmediateMirrorForceEnergy(const Float3& diff, const float chargeProduct, Float3 boxSize) {
	//	ForceEnergy forceEnergy{};
	//	for (int zOffset = -1; zOffset <= 1; zOffset += 1) {
	//		for (int yOffset = -1; yOffset <= 1; yOffset += 1) {
	//			for (int xOffset = -1; xOffset <= 1; xOffset += 1) {
	//				if (xOffset == 0 && yOffset == 0 && zOffset == 0)
	//					continue;

	//				const Float3 diffMirror = diff + Float3{ xOffset, yOffset, zOffset } *boxSize;
	//				forceEnergy.force += PhysicsUtils::CalcCoulumbForce(chargeProduct, 1, diffMirror);
	//				forceEnergy.potE += PhysicsUtils::CalcCoulumbPotential(chargeProduct, 1, diffMirror.len()) * 0.5f;
	//			}
	//		}
	//	}
	//		return forceEnergy;
	//}




	Float3 GetPositionOfParticleRelativeToSelfUsingTheWierdLogicOfTheKernel(const Float3& posOtherAbs, const NodeIndex nodeindexSelf) {
		// We cant use hyperdist, since the engine will hyperpos the block, and not individual particles
		//BoundaryConditionPublic::applyHyperposNM(posSelf, posOther, sim->simparams_host.box_size, PBC);
		// Instead we do this bullshit. Figure out the relative nodeindex compared to the self nodeindex
		// 
		// The final Coulomb force is calculated using the position from the second-to-last step, thus -2 not -1
		const NodeIndex nodeindexOther = LIMAPOSITIONSYSTEM::PositionToNodeIndexNM(posOtherAbs);
		const Float3 posOtherRel = posOtherAbs - LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(nodeindexOther);

		NodeIndex nodeindexOtherHyper = nodeindexOther;
		BoundaryConditionPublic::applyHyperpos(nodeindexSelf, nodeindexOtherHyper, 3, PBC);
		const NodeIndex nodeindexOfOtherRelativeToSelf = nodeindexOtherHyper - nodeindexSelf;


		if (nodeindexOfOtherRelativeToSelf.largestMagnitudeElement() > 1)
			throw std::runtime_error("Expected to get a relative nodeindex that was immediately adjacent to nodeindexSelf");		

		const Float3 posOtherRelativeToSelf = (posOtherRel + LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(nodeindexOfOtherRelativeToSelf));
		return posOtherRelativeToSelf;
	}


	LimaUnittestResult TestAttractiveParticlesInteractingWithESandLJ(EnvMode envmode) {
		const fs::path work_folder = simulations_dir / "Pool/";
		Environment env{ work_folder, envmode};


		const int nSteps = 1000;



		SimParams params{};
		params.n_steps = nSteps;
		params.enable_electrostatics = true;
		params.data_logging_interval = 1;
		params.cutoff_nm = 2.f;
		GroFile grofile{ work_folder / "molecule/conf.gro" };
		grofile.box_size = Float3{ 3.f };
		grofile.atoms[0].position = Float3{ 1.f, 1.5f, 1.5f };
		grofile.atoms[1].position = Float3{ 2.f, 1.5f, 1.5f };
		TopologyFile topfile{ work_folder / "molecule/topol.top" };


		env.CreateSimulation(grofile, topfile, params);

		//env.getSimPtr()->forcefield.particle_parameters[1].epsilon = 0.f;

		// Make the particles attractive
		env.getSimPtr()->box_host->compounds[1].atom_charges[0] = -env.getSimPtr()->box_host->compounds[0].atom_charges[0];

		env.run();


		//LIMA_Print::printPythonVec("potE", env.getAnalyzedPackage()->pot_energy);

		const float actualVC = env.getAnalyzedPackage().variance_coefficient;
		const float maxVC = 1e-3;
		ASSERT(actualVC < maxVC, std::format("VC {:.3e} / {:.3e}", actualVC, maxVC));

		return LimaUnittestResult{ true, "", envmode == Full };
	}

	static void MakeChargeParticlesSim(const std::string& dirName, const float boxLen, const AtomsSelection& atomsSelection, float particlesPerNm3) {
		Environment env(simulations_dir / dirName, EnvMode::Headless);

		env.createSimulationFiles(boxLen);

		MDFiles::SimulationFilesCollection simfiles(env.getWorkdir());
		for (const auto& atom : atomsSelection) {
			auto moltype = std::make_shared<TopologyFile::Moleculetype>( atom.atomtype.atomname, 3 );
			moltype->atoms.push_back(atom.atomtype);
			simfiles.topfile->moleculetypes.insert({ atom.atomtype.atomname, moltype });
		}
		simfiles.topfile->SetSystem("ElectroStatic Field Test");
		SimulationBuilder::DistributeParticlesInBox(*simfiles.grofile, *simfiles.topfile, atomsSelection, 0.24f, particlesPerNm3);

		// Overwrite the forcefield
		simfiles.topfile->forcefieldInclude.emplace("lima_custom_forcefield.itp");
		simfiles.topfile->forcefieldInclude->contents = std::move(GenericItpFile(FileUtils::GetLimaDir() / "resources" / "forcefields" / "lima_custom_forcefield.itp"));

		simfiles.grofile->title = "ElectroStatic Field Test";
		simfiles.topfile->title = "ElectroStatic Field Test";
		simfiles.grofile->printToFile();
		simfiles.topfile->printToFile();
	}

	static LimaUnittestResult TestChargedParticlesVelocityInUniformElectricField(EnvMode envmode) {
		MakeChargeParticlesSim("ElectrostaticField", 7.f, 
			AtomsSelection{
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt2", 0, "lxx", "lx1", 0, -1.f, 10.f}, 15},
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt2", 0, "lxx", "lx2", 0, -.5f, 10.f}, 15},
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt2", 0, "lxx", "lx3", 0, -0.f, 10.f}, 40},
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt2", 0, "lxx", "lx4", 0, 0.5f, 10.f}, 15},
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt2", 0, "lxx", "lx5", 0, 1.f, 10.f},  15}
			}, 
			5.f
			);

		SimParams simparams;
		simparams.dt = 0.2f * FEMTO_TO_NANO;
		simparams.coloring_method = ColoringMethod::Charge;
		simparams.data_logging_interval = 1;
		simparams.snf_select = HorizontalChargeField;
		auto env = basicSetup("ElectrostaticField", { simparams }, envmode);

		env->getSimPtr()->box_host->uniformElectricField = UniformElectricField{ Float3{-1.f, 0.f, 0.f }, 12.f};

		env->run();	

		if (envmode == Full)
			TestUtils::CompareForces1To1(simulations_dir / "ElectrostaticField", *env, false);


		auto sim = env->getSim();


		std::map<float, std::vector<float>> velDistributions;

		// Go through each particle in each compound, and assert that their velocities are as we expect in this horizontal electric field
		for (int cid = 0; cid < sim->box_host->boxparams.n_compounds; cid++) {
			const auto& compound = sim->box_host->compounds[cid];
			const auto& compoundInterimState = sim->box_host->compoundInterimStates[cid];

			for (int pid = 0; pid < compound.n_particles; pid++) {
				const float charge = static_cast<float>(compound.atom_charges[pid]);
				const float velHorizontal = compoundInterimState.vels_prev[pid].x;

				velDistributions[charge].push_back(velHorizontal);
			}
		}

		if (envmode == Full) {
			for (const auto& pair : velDistributions) {
				const int charge = pair.first;
				const auto& velocities = pair.second;
				std::cout << "Charge: " << charge << " | Mean Velocity: " << Statistics::Mean(velocities) << " | Standard Deviation: " << Statistics::StdDev(velocities) << std::endl;
			}
		}

		std::vector<float> x;
		std::vector<float> y;
		for (const auto& pair : velDistributions) {
			x.insert(x.end(), pair.second.size(), pair.first);
			y.insert(y.end(), pair.second.begin(), pair.second.end());
		}

		const auto [slope, intercept] = Statistics::linearFit(x, y);

		if (slope >= 0.f) {
			std::string errorMsg = std::format("Slope of velocity distribution should be negative, but got {:.4f} ", slope);
			return LimaUnittestResult{ false, errorMsg, envmode == Full };
		}
		if (std::abs(intercept) > 50.f) {
			std::string errorMsg = std::format("Intercept of velocity distribution should be close to 0, but got {:.2f}",intercept);
			return LimaUnittestResult{ false, errorMsg, envmode == Full };
		}

		const float r2 = Statistics::calculateR2(x, y, slope, intercept);
		ASSERT(!std::isnan(r2), "R2 value is nan");
		if (r2 < 0.5f) {
			//std::string errorMsg = "R2 value " + std::to_string(r2) + " of velocity distribution should be close to 1";
			std::string errorMsg = std::format("R2 value {:.2f} of velocity distribution should be close to 1", r2);
			return LimaUnittestResult{ false, errorMsg, envmode == Full };
		}

		return LimaUnittestResult{ true, std::format("R2 Value: {:.2f}", r2), envmode == Full};
	}

	static LimaUnittestResult TestElectrostaticsManyParticles(EnvMode envmode) {
		MakeChargeParticlesSim("ShortrangeElectrostaticsCompoundOnly", 5.f,
			AtomsSelection{
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt1", 0, "lxx", "lxx", 0, 1.f, 12.011}, 100}, // by naming the residue lxx we let these particles be full molecules, instead of tinymols, is that ideal? Does it matter? 
			},
			5.f // TODO: If we set this density to 32 as it should be, the result diverge too much. I should look into that later. And do a similar stresstest for a simple LJ system
		);

		const int nSteps = 1000;

		SimParams simparams;
		simparams.n_steps = nSteps;
		simparams.dt = 1.f * FEMTO_TO_NANO;
		simparams.coloring_method = ColoringMethod::Charge;
		simparams.data_logging_interval = 1;
		simparams.stepsPerNlistupdate = 1;
		simparams.enable_electrostatics = true;
		simparams.cutoff_nm = 2.f;
		auto env = basicSetup("ShortrangeElectrostaticsCompoundOnly", { simparams }, envmode);

		env->run();

		auto sim = env->getSim();

		//LIMA_Print::plotEnergies(env->getAnalyzedPackage()->pot_energy, env->getAnalyzedPackage()->kin_energy, env->getAnalyzedPackage()->total_energy);

		
		//ASSERT(sim->boxparams_host.boxSize == BoxGrid::blocksizeNM * 3, "This test assumes entire BoxGrid is in Shortrange range");

		// First check that the potential energy is calculated as we would expect if we do it the simple way
		float maxForceError = 0.f;
		for (int cidSelf = 0; cidSelf < sim->box_host->boxparams.n_compounds; cidSelf++) {
			double potESum{};
			Float3 forceSum{};

			const Compound& compoundSelf = sim->box_host->compounds[cidSelf];
			const CompoundInterimState& compoundInterimSelf = sim->box_host->compoundInterimStates[cidSelf];
			const float chargeSelf = compoundSelf.atom_charges[0];

			// The final Coulomb force is calculated using the position from the second-to-last step, thus -2 not -1
			//const Float3 posSelfAbs = sim->traj_buffer->GetMostRecentCompoundparticleDatapoint(cidSelf, 0, simparams.n_steps - 2);
			const Float3 posSelfAbs = sim->traj_buffer->getCompoundparticleDatapointAtIndex(cidSelf, 0, simparams.n_steps - 2);	// We MUST have the position at this index, to get accurate forces
			const NodeIndex nodeindexSelf = LIMAPOSITIONSYSTEM::PositionToNodeIndexNM(posSelfAbs);
			const Float3 posSelfRel = posSelfAbs - LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(nodeindexSelf);

			for (int cidOther = 0; cidOther < sim->box_host->boxparams.n_compounds; cidOther++) {
				if (cidSelf == cidOther)
					continue;
				
				const auto& compoundOther = sim->box_host->compounds[cidOther];
				const float chargeOther = compoundOther.atom_charges[0];
				const Float3 posOtherAbs = sim->traj_buffer->getCompoundparticleDatapointAtIndex(cidOther, 0, simparams.n_steps - 2);
				const Float3 posOtherRelativeToSelf = GetPositionOfParticleRelativeToSelfUsingTheWierdLogicOfTheKernel(posOtherAbs, nodeindexSelf);

				const Float3 diff = posSelfRel - posOtherRelativeToSelf;

				potESum += PhysicsUtils::CalcCoulumbPotential(chargeSelf, chargeOther, diff.len()) * 0.5f;
				forceSum += PhysicsUtils::CalcCoulumbForce(chargeSelf, chargeOther, diff);
			}
			
			// Need a expected error because in the test we do true hyperdist, but in sim we do no hyperdist
			// The error arises because a particle is moved 1 boxlen, not when it is correct for hyperPos, but when it moves into the next node in the boxgrid
			// Thus this error arises only when the box is so small that a particle go directly from nodes such as (-1, 0 0) to (1,0,0)
			//const float potEError = std::abs(compoundInterimSelf.sumPotentialenergy(0) - potESum) / potESum;
			const float forceError = std::abs((compoundInterimSelf.forces_prev[0] - forceSum).len()) / forceSum.len();
			maxForceError = std::max(maxForceError, forceError);

			//ASSERT(potEError < 1e-4, std::format("Actual PotE {:.7e} Expected potE: {:.7e} Error {:.7e}", compoundSelf.potE_interim[0], potESum, potEError));
			//ASSERT(forceError < 1e-4, std::format("Actual Force {:.7e} Expected force {:.7e} Error {:.7e}", compoundSelf.forces_interim[0].len(), forceSum.len(), forceError));
		}

		// Now do the normal VC check
		const float targetVarCoeff = 8.66e-3f;
		auto analytics = SimAnalysis::analyzeEnergy(sim.get());


		ASSERT(analytics.variance_coefficient < targetVarCoeff, std::format("VC {:.3e} / {:.3e}", analytics.variance_coefficient, targetVarCoeff));

		return LimaUnittestResult{ 
			true, 
			std::format("VC {:.3e} / {:.3e} Max F error {:.3e}", analytics.variance_coefficient, targetVarCoeff, maxForceError),
			envmode == Full };
	}


	LimaUnittestResult TestLongrangeEsNoLJTwoParticles(EnvMode envmode) {
		const fs::path work_folder = simulations_dir / "Pool/";
		Environment env{ work_folder, envmode};

		struct TestSetup {
			Float3 p0, p1;
			std::string name{};
			Float3 mirrorDir{ -1,0, 0 };
		};

		TopologyFile topfile{ work_folder / "molecule/topol.top" };
		GroFile grofile{ work_folder / "molecule/conf.gro" };
		grofile.box_size = Float3{ 30.f };

		// First check with 2 particles exactly on the nodeindices, such that the longrange approximation is perfect
		std::vector<TestSetup> testSetups{
			 {Float3{ 6.025f, 10.025f, 10.f }, Float3{ 8.025f, 10.025f, 10.025f }, "x-dim only 2 nm"}
			,{Float3{ 4.f, 10.f, 10.f }, Float3{ 10.07f, 10.025f, 10.025f },  "x-dim only 8 nm"}
			,{Float3{ 4.025f, 11.025f, 9.025f }, Float3{ 10.025f, 10.025f, 10.025f }, "3d force"}
			,{Float3{ 9.6f, 10.f, 10.f }, Float3{ 10.f, 10.f, 10.f}, "Within SR range"}
			,{Float3{ .1f, 10.f, 10.f }, Float3{ 5.6f, 10.f, 10.f }, "Close to boundary"}
			,{Float3{ .1f, 10.f, 10.f }, Float3{ grofile.box_size.x - 2.f, 10.f, 10.f }, "With hyperpos closest", {0.f,0.f, 0.f}}
		};





		SimParams params{};
		params.n_steps = 1;
		params.data_logging_interval = 1;
		const float c0 = -1.f * elementaryChargeToKiloCoulombPerMole;
		const float c1 = 1.f * elementaryChargeToKiloCoulombPerMole;

		for (int testIndex = 0; testIndex < testSetups.size(); testIndex++) {
			const auto setup = testSetups[testIndex];

			grofile.atoms[0].position = setup.p0;
			grofile.atoms[1].position = setup.p1;

			env.CreateSimulation(grofile, topfile, params);
			env.getSimPtr()->box_host->compounds[0].atom_charges[0] = c0;
			env.getSimPtr()->box_host->compounds[1].atom_charges[0] = c1;



			/*std::vector<ForceEnergy> forceEnergy(2);
			PMEtest::computePME({ setup.p0, setup.p1 }, { c0, c1 }, grofile.box_size.x, 2.5f, forceEnergy);*/




			Float3 hyperposOther = grofile.atoms[1].position;
			BoundaryConditionPublic::applyHyperposNM(grofile.atoms[0].position, hyperposOther, grofile.box_size.x, PBC);
			
			const Float3 diff = grofile.atoms[0].position - hyperposOther;

			const Float3 diffFromMirror = setup.p0 - (setup.p1 + (Float3(grofile.box_size) * setup.mirrorDir));
			const Float3 mirrorForce = PhysicsUtils::CalcCoulumbForce(c0, c1, diffFromMirror);
			const float mirrorPotential = PhysicsUtils::CalcCoulumbPotential(c0, c1, diffFromMirror.len()) * 0.5f;


			const float expectedPotential = PhysicsUtils::CalcCoulumbPotential(c0, c1, diff.len()) * 0.5f + mirrorPotential;
			const Float3 expectedForce = PhysicsUtils::CalcCoulumbForce(c0, c1, diff) + mirrorForce;



			env.run();
			const auto sim = env.getSim();

			const Float3 actualForce = sim->forceBuffer->getCompoundparticleDatapointAtIndex(0, 0, 0);
			const float actualPotential = sim->potE_buffer->getCompoundparticleDatapointAtIndex(0, 0, 0);
			const float potEError = std::abs((actualPotential - expectedPotential) / expectedPotential);
			const float forceError = (actualForce - expectedForce).len() / expectedForce.len();

			
			//ghostforce.print('G');
			if (envmode == Full) {
				printf("Exp Force %.3e %.3e %.3e\n", expectedForce.x, expectedForce.y, expectedForce.z);
				//printf("CPU Force %.3e %.3e %.3e\n", forceEnergy[0].force.x, forceEnergy[0].force.y, forceEnergy[0].force.z);
				printf("GPU force %.3e %.3e %.3e\n", actualForce.x, actualForce.y, actualForce.z);
				printf("Mir force %.3e %.3e %.3e\n", mirrorForce.x, mirrorForce.y, mirrorForce.z);
				printf("Expected pot %.3e GPU %.3e mirror %.3e\n", expectedPotential, actualPotential, mirrorPotential);
			}

			ASSERT(forceError < 0.07f, std::format("{}\n\tActual Force {:.3e} {:.3e} {:.3e} Expected force {:.3e} {:.3e} {:.3e} Error {:.3f}", setup.name, actualForce.x, actualForce.y, actualForce.z, expectedForce.x, expectedForce.y, expectedForce.z, forceError));
			// Potential is hopeless to match realspace and kspace
			ASSERT(potEError < 3.f, std::format("{}\n\tActual PotE {:.5e} Expected potE: {:.5e} Error {:.3}", setup.name, actualPotential, expectedPotential, potEError));

			const Float3 actualForceP1 = sim->forceBuffer->getCompoundparticleDatapointAtIndex(1, 0, 0);
			ASSERT((actualForce + actualForceP1).len() / actualForce.len() < 0.001f,
				std::format("{}\n\tExpected forces to be equal and opposite. P0 {:.3e} {:.3e} {:.3e} P1 {:.3e} {:.3e} {:.3e}", setup.name,
					actualForce.x, actualForce.y, actualForce.z, actualForceP1.x, actualForceP1.y, actualForceP1.z));			
		}

		return LimaUnittestResult{ true, "Success", envmode == Full };
	}

	LimaUnittestResult PlotPmePotAsFactorOfDistance(EnvMode envmode) {
		const fs::path work_folder = simulations_dir / "Pool/";
		Environment env{ work_folder, envmode };

		TopologyFile topfile{ work_folder / "molecule/topol.top" };
		GroFile grofile{ work_folder / "molecule/conf.gro" };
		grofile.box_size = Float3{ 30.f };

		SimParams params{};
		params.n_steps = 1;
		params.data_logging_interval = 1;

		const float c0 = -1.f * elementaryChargeToKiloCoulombPerMole;
		const float c1 = 1.f * elementaryChargeToKiloCoulombPerMole;

		std::vector<Float3> actualForce, expectedForce;
		std::vector<float> actualPot, expectedPot;
		std::vector<float> distances;

		for (float dist = 15.f; dist > 0.4; dist -= 1.13f) {
			grofile.atoms[0].position = Float3{ 15.025f, 10.025f, 10.025f };
			grofile.atoms[1].position = grofile.atoms[0].position - Float3{ dist, 0.f, 0.f };

			env.CreateSimulation(grofile, topfile, params);
			env.getSimPtr()->box_host->compounds[0].atom_charges[0] = c0;
			env.getSimPtr()->box_host->compounds[1].atom_charges[0] = c1;


			Float3 hyperposOther = grofile.atoms[1].position;
			BoundaryConditionPublic::applyHyperposNM(grofile.atoms[0].position, hyperposOther, grofile.box_size.x, PBC);
			const Float3 diff = grofile.atoms[0].position - hyperposOther;

//			Float3 ghostforce = PhysicsUtils::CalcCoulumbForce(c0, c1, Float3{ diff.x - grofile.box_size.x , 0.f, 0.f });
			//ForceEnergy mirrorFE = CalcImmediateMirrorForceEnergy(diff, c0 * c1, grofile.box_size);
			Float3 mirrorForce = PhysicsUtils::CalcCoulumbForce(c0, c1, Float3{ diff.x - grofile.box_size.x , 0.f, 0.f });
			float mirrorPotential = PhysicsUtils::CalcCoulumbPotential(c0, c1, (Float3{ diff.x - grofile.box_size.x , 0.f, 0.f }).len()) * 0.5f;
			Float3 force = PhysicsUtils::CalcCoulumbForce(c0, c1, diff);
			float pot = PhysicsUtils::CalcCoulumbPotential(c0, c1, diff.len()) * 0.5f;

			expectedPot.push_back(pot);
			expectedForce.push_back(force + mirrorForce);

			env.run();
			const auto sim = env.getSim();			

			actualPot.push_back(sim->potE_buffer->getCompoundparticleDatapointAtIndex(0, 0, 0));	
			actualForce.push_back(sim->forceBuffer->getCompoundparticleDatapointAtIndex(0, 0, 0));

			distances.push_back(dist);
		}

		//PlotUtils::PlotData({ actual, expected }, {"Actual Pot", "Expected Pot"}, distances);
		std::vector<float> potError(actualPot.size());
		std::vector<float> forceError(actualPot.size());
		for (int i = 0; i < potError.size(); i++) {
			potError[i] = std::abs((actualPot[i] - expectedPot[i]) / expectedPot[i]);
			forceError[i] = std::abs((actualForce[i] - expectedForce[i]).len() / expectedForce[i].len());
		}

		//PlotUtils::PlotData({ actual, expected, errors }, { "Actual Force mag", "Expected Force mag", "Error"}, distances);
		PlotUtils::PlotData({ potError, forceError}, { "PotE Error", "Force Error"}, distances);
		return LimaUnittestResult{ true, "", envmode == Full };
	}

	LimaUnittestResult TestConsistentEnergyWhenGoingFromLresToSres(EnvMode envmode) {
		const fs::path work_folder = simulations_dir / "Pool/";
		Environment env{ work_folder, envmode };


		TopologyFile topfile{ work_folder / "molecule/topol.top" };
		GroFile grofile{ work_folder / "molecule/conf.gro" };
		grofile.box_size = Float3{ 20.f };

		SimParams params{};
		params.n_steps = 500;
		params.data_logging_interval = 1;
		params.dt = 1.f * FEMTO_TO_NANO;

		grofile.atoms[0].position = Float3{ 7.0f, 10.f, 10.f };
		grofile.atoms[1].position = Float3{ 9.f, 10.f, 10.f };

		const float c0 = 1.f * elementaryChargeToKiloCoulombPerMole;
		const float c1 = -c0;

		env.CreateSimulation(grofile, topfile, params);
		env.getSimPtr()->box_host->compounds[0].atom_charges[0] = c0;
		env.getSimPtr()->box_host->compounds[1].atom_charges[0] = c1;

		//env.getSimPtr()->box_host->compoundInterimStates[0].vels_prev[0] = Float3{ 5000, 0, 0 };

		Float3 hyperposOther = grofile.atoms[1].position;
		BoundaryConditionPublic::applyHyperposNM(grofile.atoms[0].position, hyperposOther, grofile.box_size.x, PBC);
		const Float3 diff = grofile.atoms[0].position - hyperposOther;
		const float expectedPotential = PhysicsUtils::CalcCoulumbPotential(c0, c1, diff.len()) * 0.5f;
		const Float3 expectedForce = PhysicsUtils::CalcCoulumbForce(c0, c1, diff);



		env.run();

		// First go trough the traj data and find the step where the particles are less than 0.5 nm apart
		int step = -1;
		for (int i = 0; i < params.n_steps; i++)
		{
			const Float3 pos0 = env.getSimPtr()->traj_buffer->getCompoundparticleDatapointAtIndex(0, 0, i);
			const Float3 pos1 = env.getSimPtr()->traj_buffer->getCompoundparticleDatapointAtIndex(1, 0, i);
			const float dist = (pos0 - pos1).len();
			if (dist < 0.4f) {
				step = i;
				break;
			}
		}

		// Only use the energies up untill the step found above
		auto anal = env.getAnalyzedPackage();
		//LIMA_Print::plotEnergies(std::span(anal.pot_energy).subspan(0, step), std::span(anal.kin_energy).subspan(0, step), std::span(anal.total_energy).subspan(0, step));

		return LimaUnittestResult{ anal.variance_coefficient < 1e-3f, "", envmode == Full };
	}



	// Create many pos charged Ions as compounds. Set all LJ to 0. Compute exact SR and LR interactions between all particles. Run simulation 1 step, and compare the errors
	LimaUnittestResult TestLongrangeEsNoLJManyParticles(EnvMode envmode) {
		const float boxlen = 20.f;
		const float chargeExtern = 1.f;
		const float charge = chargeExtern * elementaryChargeToKiloCoulombPerMole;
		MakeChargeParticlesSim("ShortrangeElectrostaticsCompoundOnly", boxlen,
			AtomsSelection{
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt1", 0, "lxx", "lxx", 0, chargeExtern, 12.011}, 100}, // by naming the residue lxx we let these particles be full molecules, instead of tinymols, is that ideal? Does it matter? 
			},
			1.f
			);

		const fs::path work_folder = simulations_dir / "ShortrangeElectrostaticsCompoundOnly/";
		Environment env{ work_folder, envmode };

		
		
		SimParams params{};
		params.n_steps = 2;
		params.data_logging_interval = 1;
		GroFile grofile{ work_folder / "molecule/conf.gro" };

		TopologyFile topfile{ work_folder / "molecule/topol.top" };
		//topfile.GetSystemMutable().molecules.resize(3);

		env.CreateSimulation(grofile, topfile, params);
		env.getSimPtr()->forcefield.particle_parameters[0].epsilonSqrt = 0.f; // There is only 1 particle type
		env.run();


		// Now compute all expected forces and potentials
		std::vector<float> expectedPotentials(grofile.atoms.size());
		std::vector<Float3> expectedForces(grofile.atoms.size());

		for (int i = 0; i < grofile.atoms.size(); i++) {
			const auto& atom = grofile.atoms[i];
			float potential{};
			Float3 force{};
			for (int j = 0; j < grofile.atoms.size(); j++) {
				if (i == j)
					continue;

				Float3 hyperposOther = grofile.atoms[j].position;
				BoundaryConditionPublic::applyHyperposNM(atom.position, hyperposOther, boxlen, PBC);
				const Float3 diff = atom.position - hyperposOther;
				const float expectedPotential = PhysicsUtils::CalcCoulumbPotential(charge, charge, diff.len()) * 0.5f;
				const Float3 expectedForce = PhysicsUtils::CalcCoulumbForce(charge, charge, diff);

				potential += expectedPotential;
				force += expectedForce;
			}
			expectedPotentials[i] = potential;
			expectedForces[i] = force;
		}

		const auto sim = env.getSim();
		const Float3 actualForce = sim->forceBuffer->getCompoundparticleDatapointAtIndex(0, 0, 0);

		std::vector<float> potErrors(grofile.atoms.size());
		std::vector<float> forceErrors(grofile.atoms.size());

		for (int i = 0; i < grofile.atoms.size(); i++) {
			const float potEError = std::abs(sim->potE_buffer->getCompoundparticleDatapointAtIndex(i, 0, 0) - expectedPotentials[i]) / expectedPotentials[i];
			Float3 actualForce = sim->forceBuffer->getCompoundparticleDatapointAtIndex(i, 0, 0);
			Float3 expectedForce = expectedForces[i];
			Float3 position = grofile.atoms[i].position;
			const float forceError = (sim->forceBuffer->getCompoundparticleDatapointAtIndex(i, 0, 0) - expectedForces[i]).len() / expectedForces[i].len();

			if (expectedForces[i].len() < 50'000.f) // [J/mol/nm
				continue; // Force is quite small, hard to be relative accurate here

			if (forceError > 4.f)
				int a = 0;

			potErrors[i] = potEError;
			forceErrors[i] = forceError;			
		}

		/*const float maxPotError = *std::max_element(potErrors.begin(), potErrors.end());
		const float meanPotError = Statistics::Mean(potErrors);
		ASSERT(meanPotError < 5e-2, std::format("Mean PotE Error {:.3e}", meanPotError));
		ASSERT(maxPotError < 1, std::format("Max PotE Error {:.3e}", maxPotError));*/
		
		const float maxForceError = *std::max_element(forceErrors.begin(), forceErrors.end());
		const float meanForceError = Statistics::Mean(forceErrors);
		ASSERT(meanForceError < 0.18f, std::format("Mean Force Error {:.3f}", meanForceError));
		//ASSERT(maxForceError < 0.8f, std::format("Max Force Error {:.3e}", maxForceError));

		return LimaUnittestResult{ true, "", envmode == Full };
	}
}

