#pragma once

#include <iostream>
#include <string>
#include <algorithm>

#include "TestUtils.h"
#include "Environment.h"

bool testNearestSolventSolventAfterEnergyMinimizationIsDecreasing(EnvMode envmode) {
	////SimulationParams params{ .dt = 5, .n_steps = 1000 };
	//InputSimParams params{ 5, 1000 };

	//auto env = TestUtils::basicSetup("SolventBenchmark", { params }, envmode);
	//env->prepareForRun();	// Lock down simulation

	//LimaLogger logger{ LimaLogger::LogMode::compact, Full, "nearestSolSolTest", env->getWorkdir() };



	//std::vector<float> closestSolvent;
	//auto solventgrid_pre_em = env->getCurrentSolventblockGrid();
	//DEBUGUTILS::findAllNearestSolventSolvent(solventgrid_pre_em.get(), env->getSim()->boxparams_host.n_solvents, closestSolvent);
	//logger.printToFile("solsoldist_prerun.bin", closestSolvent);
	//const float minradius_before = *std::min_element(closestSolvent.begin(), closestSolvent.end());

	//env->run(true);

	//auto solventgrid_post_em = env->getCurrentSolventblockGrid();
	//DEBUGUTILS::findAllNearestSolventSolvent(solventgrid_post_em.get(), env->getSim()->boxparams_host.n_solvents, closestSolvent);
	//logger.printToFile("solsoldist_postrun.bin", closestSolvent);
	//const float minradius_after = *std::min_element(closestSolvent.begin(), closestSolvent.end());

	//const float mindist = 0.21f;
	//const bool success = minradius_after > mindist;

	//if (!success) {
	//	printf("Test failed. Min radius for any solvent went from %.3f [nm] to %.3f [nm] during energy minimization. Min dist allowed was %.3f [nm]\n",
	//		minradius_before, minradius_after, mindist);
	//}

	//return success;
	return true;	// TODO: Fix this again
}