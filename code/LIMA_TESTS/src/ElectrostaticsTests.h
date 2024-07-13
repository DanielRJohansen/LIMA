#pragma once

#include "TestUtils.h"
#include "Environment.h"
#include "Printer.h"
#include "Utilities.h"
#include "LimaTypes.cuh"
//#include "EngineUtils.cuh"
#include "LimaPositionSystem.cuh"
#include "Statistics.h"

#include "EngineCore.h"
#include "PhysicsUtils.cuh"
#include <unordered_map>

#include <iostream>
#include <string>
#include <algorithm>


namespace ElectrostaticsTests {
	using namespace TestUtils;


	static LimaUnittestResult CoulombForceSanityCheck(EnvMode envmode) {
		const float calcedForce = PhysicsUtils::CalcCoulumbForce(1.f*elementaryChargeToKiloCoulombPerMole, 1.f*elementaryChargeToKiloCoulombPerMole, Float3{ 1.f, 0.f, 0.f }).len(); // [1/l N / mol]
		const float expectedForce = 2.307078e-10 * AVOGADROSNUMBER / UNIT_TO_LIMA;  // [1/l N / mol] https://www.omnicalculator.com/physics/coulombs-law

		ASSERT(std::abs(calcedForce - expectedForce) / expectedForce < 0.0001f, std::format("Expected {:.2f} Actual {:.2f}", expectedForce, calcedForce));
		return LimaUnittestResult{ LimaUnittestResult::SUCCESS, "", envmode == Full};
	}






	LimaUnittestResult doPoolBenchmarkES(EnvMode envmode) {
		const std::string work_folder = simulations_dir + "Pool/";
		const std::string conf = work_folder + "molecule/conf.gro";
		const std::string topol = work_folder + "molecule/topol.top";
		Environment env{ work_folder, envmode, false };

		const float particle_mass = 12.011000f / 1000.f;	// kg/mol
		float temperature =  400.f;
		std::vector<float> varcoffs;
		std::vector<float> energy_gradients;

		SimParams params{};
		params.n_steps = 1000;
		params.enable_electrostatics = true;
		params.data_logging_interval = 1;
		GroFile grofile{ conf };
		grofile.box_size = Float3{ 3.f };
		grofile.atoms[0].position = Float3{ 1.f, 1.5f, 1.5f };
		grofile.atoms[1].position = Float3{ 1.3f, 1.5f, 1.5f };
		TopologyFile topfile{ topol };
			
		env.CreateSimulation(grofile, topfile, params);

		Box* box_host = env.getSimPtr()->box_host.get();

		// Disable LJ force
		ASSERT(box_host->compounds[0].atom_types[0] == 1, "Expected atom type 1");
		ASSERT(box_host->compounds[1].atom_types[0] == 1, "Expected atom type 1");
		env.getSimPtr()->forcefield.particle_parameters[1].epsilon = 0.f;
	

		env.run();

		const auto analytics = env.getAnalyzedPackage();
		varcoffs.push_back(analytics->variance_coefficient);
		energy_gradients.push_back(analytics->energy_gradient);
		if (envmode != Headless) { Analyzer::printEnergy(analytics); }
		

		if (envmode != Headless) {
			LIMA_Print::printPythonVec("kinE", analytics->kin_energy);
			LIMA_Print::printPythonVec("potE", analytics->pot_energy);
			LIMA_Print::printPythonVec("totE", analytics->total_energy);
		}



		const auto result = evaluateTest(varcoffs, 4.f, energy_gradients, 1e-7);
		const auto status = result.first == true ? LimaUnittestResult::SUCCESS : LimaUnittestResult::FAIL;

		return LimaUnittestResult{ status, result.second, envmode == Full };
	}





	/// <summary>
	/// 
	/// </summary>
	/// <param name="dirName"></param>
	/// <param name="boxLen"></param>
	/// <param name="atomType">For determining LJ params, but also coloring</param>
	static void MakeChargeParticlesSim(const std::string& dirName, const float boxLen, const AtomsSelection& atomsSelection) {
		Environment env(simulations_dir + dirName, EnvMode::Headless, false);

		env.createSimulationFiles(boxLen);

		auto a = env.getWorkdir();
		MDFiles::SimulationFilesCollection simfiles(env.getWorkdir());
		SimulationBuilder::DistributeParticlesInBox(*simfiles.grofile, *simfiles.topfile, atomsSelection, 0.15f, 5.f);

		simfiles.grofile->title = "ElectroStatic Field Test";
		simfiles.topfile->title = "ElectroStatic Field Test";
		simfiles.grofile->printToFile();
		simfiles.topfile->printToFile();
	}

	static LimaUnittestResult TestChargedParticlesVelocityInUniformElectricField(EnvMode envmode) {
		MakeChargeParticlesSim("ElectrostaticField", 7.f, 
			AtomsSelection{
				{TopologyFile::AtomsEntry{";residue_X", 0, "C", 0, "XXX", "lxx", 0, -1.f, 10.f}, 15},
				{TopologyFile::AtomsEntry{";residue_X", 0, "C", 0, "XXX", "lxx", 0, -.5f, 10.f}, 15},
				{TopologyFile::AtomsEntry{";residue_X", 0, "C", 0, "XXX", "lxx", 0, -0.f, 10.f}, 40},
				{TopologyFile::AtomsEntry{";residue_X", 0, "C", 0, "XXX", "lxx", 0, 0.5f, 10.f}, 15},
				{TopologyFile::AtomsEntry{";residue_X", 0, "C", 0, "XXX", "lxx", 0, 1.f, 10.f},  15}
			}
			);


		SimParams simparams{ 5000, 20, true, PBC };
		simparams.coloring_method = ColoringMethod::Charge;
		simparams.data_logging_interval = 20;
		simparams.snf_select = HorizontalChargeField;
		auto env = basicSetup("ElectrostaticField", { simparams }, envmode);

		env->getSimPtr()->box_host->uniformElectricField = UniformElectricField{ {-1, 0, 0 }, 15.f *KILO};
	
		env->run();	

		auto sim = env->getSim();


		std::unordered_map<int, std::vector<float>> velDistributions;

		// Go through each particle in each compound, and assert that their velocities are as we expect in this horizontal electric field
		for (int cid = 0; cid < sim->box_host->boxparams.n_compounds; cid++) {
			const auto& compound = sim->box_host->compounds[cid];

			for (int pid = 0; pid < compound.n_particles; pid++) {

				const int charge = static_cast<int>(compound.atom_charges[pid]);
				const float velHorizontal = compound.vels_prev[pid].x;

				velDistributions[charge].push_back(velHorizontal);
			}
		}

		if (envmode == Full) {
			for (const auto& pair : velDistributions) {
				const int charge = pair.first;
				const auto& velocities = pair.second;
				std::cout << "Charge: " << charge << " | Mean Velocity: " << getMean(velocities) << " | Standard Deviation: " << getStdDev(velocities) << std::endl;
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
			return LimaUnittestResult{ LimaUnittestResult::FAIL, errorMsg, envmode == Full };
		}
		if (std::abs(intercept) > 50.f) {
			std::string errorMsg = std::format("Intercept of velocity distribution should be close to 0, but got {:.2f}",intercept);
			return LimaUnittestResult{ LimaUnittestResult::FAIL, errorMsg, envmode == Full };
		}

		const float r2 = Statistics::calculateR2(x, y, slope, intercept);

		if (r2 < 0.8) {
			//std::string errorMsg = "R2 value " + std::to_string(r2) + " of velocity distribution should be close to 1";
			std::string errorMsg = std::format("R2 value {:.2f} of velocity distribution should be close to 1", r2);
			return LimaUnittestResult{ LimaUnittestResult::FAIL, errorMsg, envmode == Full };
		}

		return LimaUnittestResult{ LimaUnittestResult::SUCCESS, std::format("R2 Value: {:.2f}", r2), envmode == Full};
	}

	static LimaUnittestResult TestShortrangeElectrostaticsCompoundsOnly(EnvMode envmode) {
		MakeChargeParticlesSim("ShortrangeElectrostaticsCompoundOnly", 3.f,
			AtomsSelection{
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt1", 0, "XXX", "lxx", 0, 1.f, 10.f}, 100},				
			}
		);

		const int nSteps = 2000;

		SimParams simparams{ nSteps, 20, false, PBC };
		simparams.coloring_method = ColoringMethod::Charge;
		simparams.data_logging_interval = 5;
		simparams.enable_electrostatics = true;
		auto env = basicSetup("ShortrangeElectrostaticsCompoundOnly", { simparams }, envmode);

		env->run();

		auto sim = env->getSim();

		
		ASSERT(sim->boxparams_host.boxSize == BoxGrid::blocksizeNM * 3, "This test assumes entire BoxGrid is in Shortrange range");

		// First check that the potential energy is calculated as we would expect if we do it the simple way
		float maxForceError = 0.f;
		for (int cidSelf = 0; cidSelf < sim->box_host->boxparams.n_compounds; cidSelf++) {
			float potESum = 0.f;
			Float3 forceSum = 0.f;

			const Compound& compoundSelf = sim->box_host->compounds[cidSelf];
			const int chargeSelf = static_cast<int>(compoundSelf.atom_charges[0]);
			const Float3 posSelf = sim->traj_buffer->GetMostRecentCompoundparticleDatapoint(cidSelf, 0, nSteps-1);

			for (int cidOther = 0; cidOther < sim->box_host->boxparams.n_compounds; cidOther++) {
				if (cidSelf == cidOther)
					continue;

				const auto& compoundOther = sim->box_host->compounds[cidOther];
				const int chargeOther = static_cast<int>(compoundOther.atom_charges[0]);
				Float3 posOther = sim->traj_buffer->GetMostRecentCompoundparticleDatapoint(cidOther, 0, nSteps-1);
				BoundaryConditionPublic::applyHyperposNM(posSelf, posOther, sim->simparams_host.box_size, PBC);
				
				const Float3 diff = posSelf - posOther;

				potESum += PhysicsUtils::CalcCoulumbPotential(chargeSelf, chargeOther, diff.len()) * 0.5f;
				forceSum += PhysicsUtils::CalcCoulumbForce(chargeSelf, chargeOther, diff);
			}
			
			// Need a expected error because in the test we do true hyperdist, but in sim we do no hyperdist
			// The error arises because a particle is moved 1 boxlen, not when it is correct for hyperPos, but when it moves into the next node in the boxgrid
			// Thus this error arises only when the box is so small that a particle go directly from nodes such as (-1, 0 0) to (1,0,0)
			const float potEError = std::abs(compoundSelf.potE_interim[0] - potESum) / potESum;
			const float forceError = std::abs((compoundSelf.forces_interim[0] - forceSum).len()) / forceSum.len();
			maxForceError = std::max(maxForceError, forceError);

			ASSERT(potEError < 1.f, std::format("Actual PotE {:.2e} Expected potE: {:.2e}", compoundSelf.potE_interim[0], potESum));
			//ASSERT(forceError < 0.5f, std::format("Actual Force {:.2e} Expected force {:.2e}", compoundSelf.forces_interim[0].len(), forceSum.len()));
		}

		// Now do the normal VC check
		const float targetVarCoeff = 1.5e-3f;
		auto analytics = env->getAnalyzedPackage();
		const auto result = evaluateTest({ analytics->variance_coefficient }, targetVarCoeff, { analytics->energy_gradient }, 2e-5);
		if (result.first==false)
			return LimaUnittestResult{ LimaUnittestResult::FAIL, result.second, envmode == Full };

		return LimaUnittestResult{ 
			LimaUnittestResult::SUCCESS, 
			std::format("VC {:.3e} / {:.3e} Max F error {:.2f}", analytics->variance_coefficient, targetVarCoeff, maxForceError),
			envmode == Full };
	}
}

