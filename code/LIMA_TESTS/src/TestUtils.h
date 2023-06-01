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

//#include <cuda/std/detail/libcxx/include/optional>


namespace TestUtils {

	// Creates a simulation from the folder which should contain a molecule with conf and topol
	// Returns an environment where solvents and compound can still be modified, and nothing (i hope) have
	// yet been moved to device. I should find a way to enforce this...
	Environment basicSetup(const std::string& foldername, LAL::optional<InputSimParams> simparams) {
		
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/" + foldername + "/";
		const std::string conf = work_folder + "molecule/conf.gro";
		const std::string topol = work_folder + "molecule/topol.top";

		Environment env{work_folder};

		auto ip = simparams
			? simparams.value()
			: env.loadInputSimParams(work_folder + "sim_params.txt");

		env.CreateSimulation(conf, topol, ip);

		return env;
	}

	bool verifyStability(Environment& env, float max_dev) {
		auto analytics = env.getAnalyzedPackage();
		Analyzer::printEnergy(analytics);
		float std_dev = Analyzer::getVarianceCoefficient(analytics->total_energy);

		LIMA_Print::printMatlabVec("std_dev", std::vector<float>{std_dev});

		return (std_dev < max_dev);
	}
}


