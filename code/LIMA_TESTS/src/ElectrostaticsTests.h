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

	/// <summary>
	/// 
	/// </summary>
	/// <param name="dirName"></param>
	/// <param name="boxLen"></param>
	/// <param name="atomType">For determining LJ params, but also coloring</param>
	static void MakeChargeParticlesSim(const std::string& dirName, const float boxLen, const std::string& atomType) {
		Environment env(simulations_dir + dirName, EnvMode::Headless, false);

		env.createSimulationFiles(boxLen);

		// Create 5 custom particles
		AtomsSelection atomsSelection{
			{TopologyFile::AtomsEntry{";residue_X", 0, atomType, 0, "XXX", "lxx", 0, -1.f, 10.f}, 15},
			{TopologyFile::AtomsEntry{";residue_X", 0, atomType, 0, "XXX", "lxx", 0, -.5f, 10.f}, 15},
			{TopologyFile::AtomsEntry{";residue_X", 0, atomType, 0, "XXX", "lxx", 0, -0.f, 10.f}, 40},
			{TopologyFile::AtomsEntry{";residue_X", 0, atomType, 0, "XXX", "lxx", 0, 0.5f, 10.f}, 15},
			{TopologyFile::AtomsEntry{";residue_X", 0, atomType, 0, "XXX", "lxx", 0, 1.f, 10.f},  15}
		};
		auto a = env.getWorkdir();
		MDFiles::SimulationFilesCollection simfiles(env.getWorkdir());
		SimulationBuilder::DistributeParticlesInBox(*simfiles.grofile, *simfiles.topfile, atomsSelection, 0.15f, 5.f);

		simfiles.grofile->title = "ElectroStatic Field Test";
		simfiles.topfile->title = "ElectroStatic Field Test";
		simfiles.grofile->printToFile();
		simfiles.topfile->printToFile();
	}

	static LimaUnittestResult TestChargedParticlesVelocityInUniformElectricField(EnvMode envmode) {
		MakeChargeParticlesSim("ElectrostaticField", 7.f, "C");


		SimParams simparams{ 5000, 20, true, PBC };
		simparams.coloring_method = ColoringMethod::Charge;
		simparams.data_logging_interval = 20;
		simparams.snf_select = HorizontalChargeField;
		auto env = basicSetup("ElectrostaticField", { simparams }, envmode);

		env->getSimPtr()->box_host->uniformElectricField = UniformElectricField{ {-1, 0, 0 }, 15.f };
	
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
		MakeChargeParticlesSim("ShortrangeElectrostaticsCompoundOnly", 3.f, "lt1");

		const int nSteps = 1000;

		SimParams simparams{ nSteps, 20, true, PBC };
		simparams.coloring_method = ColoringMethod::Charge;
		simparams.data_logging_interval = 20;
		auto env = basicSetup("ShortrangeElectrostaticsCompoundOnly", { simparams }, envmode);

		env->run();

		auto sim = env->getSim();

		
		ASSERT(sim->boxparams_host.boxSize == BoxGrid::blocksizeNM * 3, "This test assumes entire BoxGrid is in Shortrange range");

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

				potESum += PhysicsUtils::CalcCoulumbPotential(chargeSelf, chargeOther, diff.len());
				forceSum += PhysicsUtils::CalcCoulumbForce(chargeSelf, chargeOther, diff);
			}

			forceSum = forceSum / 1'000'000'000.f;
			potESum = potESum / 1'000'000'000.f;
			
			const float potEError = std::abs(potESum - compoundSelf.potE_interim[0]);
			const float forceError = (forceSum - compoundSelf.forces_interim[0].x).len();
			maxForceError = std::max(maxForceError, forceError);

			const float potEErrorThreshold = 1.f;
			ASSERT(potEError < potEErrorThreshold, std::format("PotE Error {:.2f} above threshold {:.2f}. True potE: {:.2f}", potEError, potEErrorThreshold, potESum));
			const float forceErrorThreshold = 1.f;
			ASSERT(forceError < forceErrorThreshold, std::format("Force Error {:.2f} above threshold {:.2f}", forceError, forceErrorThreshold));
		}

		return LimaUnittestResult{ LimaUnittestResult::SUCCESS, std::format("Max force error: {:.2f}", maxForceError), envmode == Full };
	}
}

