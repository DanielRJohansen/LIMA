#pragma once

#include "TestUtils.h"
#include "Environment.h"
#include "Printer.h"
#include "Utilities.h"
#include "LimaTypes.cuh"
#include "EngineUtils.cuh"



#include <iostream>
#include <string>
#include <algorithm>


namespace TestMDStability {
	using namespace TestUtils;

	static LimaUnittestResult loadAndEMAndRunBasicSimulation(const string& folder_name, EnvMode envmode, float max_vc = 0.05, float max_gradient=1e-5) {
		SimParams emparams{ 2000, 20, true, PBC };
		auto env = basicSetup(folder_name, { emparams }, envmode);

		// Do em
		env->run(true);
		//Analyzer::findAndDumpPiecewiseEnergies(*env->getSimPtr(), env->getWorkdir());

		// Do sim
		//InputSimParams simparams{ 100, 2000 };
		const std::string work_folder = simulations_dir + folder_name;
		const std::string simpar_path = work_folder + "/sim_params.txt";
		const SimParams simparams = env->loadSimParams(simpar_path);
		auto sim = env->getSim();
		env->CreateSimulation(*sim, simparams);
		env->run();
		//Analyzer::findAndDumpPiecewiseEnergies(*env->getSimPtr(), env->getWorkdir());

		const auto analytics = env->getAnalyzedPackage();
		
		if (envmode != Headless) {
			Analyzer::printEnergy(analytics);
			LIMA_Print::printMatlabVec("cv", std::vector<float>{ analytics->variance_coefficient});
			LIMA_Print::printMatlabVec("energy_gradients", std::vector<float>{ analytics->energy_gradient});
		}		

		const auto result = evaluateTest({ analytics->variance_coefficient }, max_vc, { analytics->energy_gradient }, max_gradient);
		const auto status = result.first == true ? LimaUnittestResult::SUCCESS : LimaUnittestResult::FAIL;

		return LimaUnittestResult{status, result.second, envmode == Full };
	}

	LimaUnittestResult doEightResiduesNoSolvent(EnvMode envmode) {
		const std::string name = "8ResNoSol";
		const std::string work_folder = simulations_dir + name + "/";
		const std::string simpar = work_folder + "sim_params.txt";

		return loadAndRunBasicSimulation(name, envmode, 2.46e-4, 8.6e-7);
	}

	static bool doMoleculeTranslationTest(std::string foldername) {
		//auto env = TestUtils::basicSetup(foldername, SimulationParams{100, 10000});

		//const float vel = EngineUtils::tempToVelocity(300, 0.012f);	// [m/s] <=> [lm/ls]
		//const float dt = env.getSimparamRef()->dt;

		//const auto dir = Float3{ -1.f, -1.f, -1.f }.norm();	// Inv since we move the pos_prev

		//auto& coordarray_prev_ptr = env.getCoordarrayRef("prev");

		//// Too lazy to figure out which coords are vacant.. So i give same velocity to all possible
		//// particles, even inactive ones. If everything works, it should be no problem
		//for (auto& coords : coordarray_prev_ptr) {
		//	for (auto& pos : coords.rel_positions) {
		//		pos += Coord{ dir * vel * dt };
		//	}
		//}

		//env.run();

		//return TestUtils::verifyStability(env, 0.08);
		return true;
	}









	// Tests proposed by GPT 3.5, lol
	/*
		Atom Shift Analysis Test: This test verifies that the positions of atoms are accurately updated according to the integration scheme being used.

		Velocity Correlation Test: This test checks the correctness of the velocity distribution, ensuring that it follows the expected Maxwell-Boltzmann distribution.

		Force Field Validation Test: This test assesses the accuracy of the interatomic force field being used, comparing predicted structural and energetic properties to experimental data.

		Temperature Equilibration Test: This test confirms that the system reaches a stable temperature and that the temperature remains constant over time.

		Trajectory Analysis Test: This test analyzes the trajectory of the system to verify that it is physically reasonable and corresponds to expected behavior.

		Pressure Control Test: This test ensures that the system is correctly maintained at the desired pressure and that fluctuations are within expected limits.

		Energy Conservation Test: This test verifies that the total energy of the system is conserved over time, which is essential for accurate and stable simulations.

		Structural Stability Test: This test assesses the stability of the system by measuring fluctuations in bond lengths, angles, and other structural parameters.

		Hydrogen Bond Analysis Test: This test analyzes the presence and behavior of hydrogen bonds within the system, which can be important for understanding molecular interactions.

		Diffusion Coefficient Test: This test measures the rate of diffusion of particles in the system, which is related to properties such as viscosity and can provide insights into the behavior of the system.

		Radial Distribution Function Test: This test calculates the radial distribution function of the system, which provides information on the distribution of particles around each other.

		Potential Energy Surface Test: This test examines the potential energy surface of the system to ensure that it is consistent with expected behavior and accurately reflects the interatomic interactions.

		Solvation Free Energy Test: This test calculates the solvation free energy of a solute molecule in a solvent, which can provide insights into the energetics of molecular interactions.

		Conformational Sampling Test: This test examines the ability of the simulation to explore the conformational space of the system, which is important for understanding the behavior of flexible molecules.

		Viscosity Calculation Test: This test measures the viscosity of the system, which is related to diffusion and can be important for understanding the behavior of complex fluids.
	*/
}

