#include "TestUtils.h"

using namespace TestUtils;

LimaUnittestResult TestBoxIsSavedCorrectlyBetweenSimulations(EnvMode envmode) {
	//const fs::path workDir = simulations_dir / "pool";
	const fs::path workDir = simulations_dir / "T4Lysozyme";

	Environment env{ workDir , envmode, false };

	SimParams simparams;
	simparams.n_steps = 100;
	simparams.dt = 100.f;
	simparams.data_logging_interval = 1;

	env.CreateSimulation(workDir / "molecule/conf.gro", workDir / "molecule/topol.top", simparams);
	//env.getSimPtr()->box_host->compounds[0].vels_prev[0] = Float3(1, 0, 0) * 2000.f;
	env.run();
	auto sim1 = env.getSim();


	simparams.dt = 0.f;
	simparams.n_steps = 1;
	env.CreateSimulation(*sim1, simparams);
	env.run();
	auto sim2 = env.getSim();
	for (int cid = 0; cid < sim2->box_host->boxparams.n_compounds; cid++) {
		for (int pid = 0; pid < sim2->box_host->compounds[cid].n_particles; pid++) {
			Float3 pos1 = sim1->traj_buffer->GetMostRecentCompoundparticleDatapoint(cid, pid, 100-1);

			Float3 pos2 = sim2->traj_buffer->GetMostRecentCompoundparticleDatapoint(cid, pid, 1-1);

			ASSERT(pos1 == pos2, "Position of compound " + std::to_string(cid) + " particle " + std::to_string(pid) + " is not the same between simulations");
		}
	}



	return LimaUnittestResult{ true, "Success", envmode == Full };
}
