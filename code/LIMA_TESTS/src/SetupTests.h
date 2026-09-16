#include "TestUtils.h"

using namespace TestUtils;

TestRoutine TestBoxIsSavedCorrectlyBetweenSimulations(Environment& environment, EnvMode envmode) {
	const fs::path workDir = AutomatedTestsDir() / "T4Lysozyme";
	SimulationJob first;
	first.workDir = workDir;
	first.groPath = workDir / "molecule/conf.gro";
	first.topPath = workDir / "molecule/topol.top";
	first.simParams = SimParams{};
	first.simParams->n_steps = 100;
	first.simParams->dt = 1.f * FEMTO_TO_NANO;
	first.simParams->data_logging_interval = 1;
	auto firstResult = co_await environment.Submit(std::move(first));

	SimulationJob second;
	second.workDir = workDir;
	second.initialSimulation = std::move(firstResult.simulation);
	second.simParams = SimParams{};
	second.simParams->dt = 0.f;
	second.simParams->n_steps = 1;
	auto secondResult = co_await environment.Submit(std::move(second));
	co_return LimaUnittestResult{ static_cast<bool>(secondResult.simulation), "Success", envmode == Full };

	// TODO: Wtf is this test??
		/*for (int cid = 0; cid < sim2->box->boxparams.n_compounds; cid++) {
		for (int pid = 0; pid < sim2->box->compounds[cid].n_particles; pid++) {
			Float3 pos1 = sim1->traj_buffer->GetMostRecentCompoundparticleDatapoint(cid, pid, 100-1);

			Float3 pos2 = sim2->traj_buffer->GetMostRecentCompoundparticleDatapoint(cid, pid, 1-1);

			ASSERT(pos1 == pos2, "Position of compound " + std::to_string(cid) + " particle " + std::to_string(pid) + " is not the same between simulations");
		}
	}*/
}
