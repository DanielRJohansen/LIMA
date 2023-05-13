#pragma once

#include "LIMA_TESTS/include/TestUtils.h"

#include "LIMA_MD/include/Environment.h"
#include "LIMA_BASE/include/Printer.h"
#include "LIMA_BASE/include/Utilities.h"
#include "LIMA_BASE/include/LimaTypes.cuh"
#include "LIMA_ENGINE/include/EngineUtils.cuh"



#include "cuda_runtime.h"

#include "LIMA_ENGINE/include/engine.cuh"	//?=?
//#include "LIMA_ENGINE/src/engine.cu"

#include <iostream>
#include <string>
#include <algorithm>

#include <iostream>



bool loadAndRunBasicSimulation(const string& folder_name, float max_dev = 0.05) {
	auto env = TestUtils::basicSetup(folder_name, {});
	env.run();

	auto analytics = env.getAnalyzedPackage();
	Analyzer::printEnergy(analytics);
	float std_dev = Analyzer::getVarianceCoefficient(analytics->total_energy);

	LIMA_Print::printMatlabVec("std_devs", std::vector<float>{ std_dev });

	return std_dev < max_dev;
}

bool testNearestSolventSolventAfterEnergyMinimizationIsDecreasing() {
	//SimulationParams params{ .dt = 5, .n_steps = 1000 };

	//auto env = TestUtils::basicSetup("SolventBenchmark", { params });
	//env.prepareForRun();	// Lock down simulation

	//LimaLogger logger{ LimaLogger::LogMode::compact, "nearestSolSolTest", env.getWorkdir() };



	//std::vector<float> closestSolvent;
	//auto solventgrid = env.getCurrentSolventblockGrid();
	//DEBUGUTILS::findAllNearestSolventSolvent(solventgrid.get(), env.getSim()->box->n_solvents, closestSolvent);
	//logger.printToFile("solsoldist_prerun.bin", closestSolvent);

	//env.run();

	//solventgrid = env.getCurrentSolventblockGrid();
	//DEBUGUTILS::findAllNearestSolventSolvent(solventgrid.get(), env.getSim()->box->n_solvents, closestSolvent);
	//logger.printToFile("solsoldist_postrun.bin", closestSolvent);

	return 1;
}