#pragma once

#include "TestUtils.h"

namespace TestMembraneBuilder {
	using namespace TestUtils;


	static LimaUnittestResult testBuildmembraneSmall(EnvMode envmode)
	{
		Environment env{ simulations_dir + "/MembraneSmall", envmode, false};

		env.CreateSimulation(7.f);
		env.createMembrane();


		return LimaUnittestResult{ LimaUnittestResult::SUCCESS , "no err", envmode == Full };

		//// Do em
		//env->run(true);
		////Analyzer::findAndDumpPiecewiseEnergies(*env->getSimPtr(), env->getWorkdir());

		//// Do sim
		////InputSimParams simparams{ 100, 2000 };
		//const std::string work_folder = simulations_dir + folder_name + "/";
		//const std::string simpar_path = work_folder + "sim_params.txt";
		//const InputSimParams simparams = env->loadInputSimParams(simpar_path);
		//auto sim = env->getSim();
		//env->CreateSimulation(*sim, simparams);
		//env->run();
		////Analyzer::findAndDumpPiecewiseEnergies(*env->getSimPtr(), env->getWorkdir());

		//const auto analytics = env->getAnalyzedPackage();

		//if (envmode != Headless) {
		//	Analyzer::printEnergy(analytics);
		//	LIMA_Print::printMatlabVec("cv", std::vector<float>{ analytics->variance_coefficient});
		//	LIMA_Print::printMatlabVec("energy_gradients", std::vector<float>{ analytics->energy_gradient});
		//}

		//const auto result = evaluateTest({ analytics->variance_coefficient }, max_vc, { analytics->energy_gradient }, max_gradient);
		//const auto status = result.first == true ? LimaUnittestResult::SUCCESS : LimaUnittestResult::FAIL;

		//return LimaUnittestResult{ status, result.second, envmode == Full };



	}
}