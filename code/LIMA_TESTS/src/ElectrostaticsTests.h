#pragma once

#include "TestUtils.h"
#include "Environment.h"
#include "Printer.h"
#include "Utilities.h"
#include "LimaTypes.cuh"
//#include "EngineUtils.cuh"
#include "LimaPositionSystem.cuh"
#include "Statistics.h"

#include <unordered_map>

#include <iostream>
#include <string>
#include <algorithm>


namespace ElectrostaticsTests {
	using namespace TestUtils;

	static void MakeChargeParticlesSim() {
		Environment env("C:/Users/Daniel/git_repo/LIMA_data/ElectrostaticField", EnvMode::Headless, false);

		env.createSimulationFiles(7.f);

		// Create 5 custom particles
		AtomsSelection atomsSelection{
			{ParsedTopologyFile::AtomsEntry{";residue_X", 0, "C", 0, "lxx", "lxx", 0, -1.f, 10.f}, 15},
			{ParsedTopologyFile::AtomsEntry{";residue_X", 0, "C", 0, "lxx", "lxx", 0, -.5f, 10.f}, 15},
			{ParsedTopologyFile::AtomsEntry{";residue_X", 0, "C", 0, "lxx", "lxx", 0, -0.f, 10.f}, 40},
			{ParsedTopologyFile::AtomsEntry{";residue_X", 0, "C", 0, "lxx", "lxx", 0, 0.5f, 10.f}, 15},
			{ParsedTopologyFile::AtomsEntry{";residue_X", 0, "C", 0, "lxx", "lxx", 0, 1.f, 10.f},  15}
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
		MakeChargeParticlesSim();


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
}

