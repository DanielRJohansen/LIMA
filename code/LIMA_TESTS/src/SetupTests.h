#include "TestUtils.h"

using namespace TestUtils;

LimaUnittestResult TestBoxIsSavedCorrectlyBetweenSimulations(EnvMode mode) {
	const fs::path workDir = simulations_dir / "T4Lysozyme";

	Environment env{ workDir , mode, false };

	SimParams simparams;
	simparams.n_steps = 10;
	simparams.data_logging_interval = 10;
	env.CreateSimulation(workDir / "conf.gro", workDir / "topol.top", simparams);
	env.run();
	auto sim = env.getSim();

	return LimaUnittestResult{ true, "Success", mode == Full };
}
