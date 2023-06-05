#pragma once

#include "LIMA_MD/include/Environment.h"
#include "LIMA_BASE/include/Printer.h"
#include "LIMA_BASE/include/Utilities.h"
#include "LIMA_BASE/include/LimaTypes.cuh"
#include "LIMA_ENGINE/include/EngineUtils.cuh"


#include <iostream>
#include <string>
#include <algorithm>

#include <iostream>
#include <optional>



namespace TestUtils {

	// Creates a simulation from the folder which should contain a molecule with conf and topol
	// Returns an environment where solvents and compound can still be modified, and nothing (i hope) have
	// yet been moved to device. I should find a way to enforce this...
	static std::unique_ptr<Environment> basicSetup(const std::string& foldername, LAL::optional<InputSimParams> simparams, Environment::Mode envmode) {
		
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/" + foldername + "/";
		const std::string conf = work_folder + "molecule/conf.gro";
		const std::string topol = work_folder + "molecule/topol.top";
		const std::string simpar = work_folder + "sim_params.txt";

		auto env = std::make_unique<Environment>(work_folder, envmode );

		const auto ip = simparams
			? simparams.value()
			: env->loadInputSimParams(simpar);

		env->CreateSimulation(conf, topol, ip);

		return std::move(env);
	}

	static bool loadAndRunBasicSimulation(
		const string& folder_name,
		Environment::Mode envmode,
		float max_dev = 0.05,
		LAL::optional<InputSimParams> ip = {},
		bool em_variant = false
	)
	{
		auto env = TestUtils::basicSetup(folder_name, ip, envmode);
		env->run(em_variant);

		auto analytics = env->getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		float std_dev = Analyzer::getVarianceCoefficient(analytics->total_energy);

		LIMA_Print::printMatlabVec("std_devs", std::vector<float>{ std_dev });

		return std_dev < max_dev;
	}

	static bool verifyStability(Environment& env, float max_dev) {
		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		float std_dev = Analyzer::getVarianceCoefficient(analytics->total_energy);

		LIMA_Print::printMatlabVec("std_dev", std::vector<float>{std_dev});

		return (std_dev < max_dev);
	}

	void stressTest(std::function<void()> func, size_t reps) {
		for (size_t i = 0; i < reps; i++) {
			func();
		}
	}
}


