#pragma once

#include <iostream>
#include <string>
#include <algorithm>

#include "LIMA_TESTS/src/TestUtils.h"
#include "LIMA_MD/include/Environment.h"

bool testNearestSolventSolventAfterEnergyMinimizationIsDecreasing() {
	//SimulationParams params{ .dt = 5, .n_steps = 1000 };
	SimulationParams params{ 5, 1000 };

	auto env = TestUtils::basicSetup("SolventBenchmark", { params });
	env.prepareForRun();	// Lock down simulation

	LimaLogger logger{ LimaLogger::LogMode::compact, "nearestSolSolTest", env.getWorkdir() };



	std::vector<float> closestSolvent;
	auto solventgrid = env.getCurrentSolventblockGrid();
	DEBUGUTILS::findAllNearestSolventSolvent(solventgrid.get(), env.getSim()->box->n_solvents, closestSolvent);
	logger.printToFile("solsoldist_prerun.bin", closestSolvent);
	const float minradius_before = *std::min_element(closestSolvent.begin(), closestSolvent.end());

	env.run(true);

	solventgrid = env.getCurrentSolventblockGrid();
	DEBUGUTILS::findAllNearestSolventSolvent(solventgrid.get(), env.getSim()->box->n_solvents, closestSolvent);
	logger.printToFile("solsoldist_postrun.bin", closestSolvent);
	const float minradius_after = *std::min_element(closestSolvent.begin(), closestSolvent.end());

	float mindist = 0.21;
	const bool success = minradius_after > mindist;

	if (!success) {
		printf("Test failed. Min radius for any solvent went from %.3f [nm] to %.3f [nm] during energy minimization. Min dist allowed was %.3f [nm]\n",
			minradius_before, minradius_after, mindist);
	}

	return success;
}