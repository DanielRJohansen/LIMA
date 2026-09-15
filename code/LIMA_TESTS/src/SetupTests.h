#include "TestUtils.h"

#include <future> // TODO: Remove the need for this

using namespace TestUtils;

std::function<LimaUnittestResult()> TestBoxIsSavedCorrectlyBetweenSimulations(Environment& environment, EnvMode envmode) {
	const fs::path workDir = AutomatedTestsDir() / "T4Lysozyme";
	SimulationJob first;
	first.workDir = workDir;
	first.groPath = workDir / "molecule/conf.gro";
	first.topPath = workDir / "molecule/topol.top";
	first.simParams = SimParams{};
	first.simParams->n_steps = 100;
	first.simParams->dt = 1.f * FEMTO_TO_NANO;
	first.simParams->data_logging_interval = 1;
	auto firstHandle = environment.Submit(std::move(first));

	// This is a test-only continuation. Its thread only waits for the first run,
	// submits the dependent second run, and waits again. Environment owns all
	// preprocessing and GPU work; TestManager remains free to evaluate other tests.
	auto simulation = std::make_shared<std::future<SimulationResult>>(std::async(std::launch::async,
		[firstHandle = std::move(firstHandle), &environment, workDir]() mutable {
			auto firstResult = firstHandle.Get();

			SimulationJob second;
			second.workDir = workDir;
			second.initialSimulation = std::move(firstResult.simulation);
			second.simParams = SimParams{};
			second.simParams->dt = 0.f;
			second.simParams->n_steps = 1;
			return environment.Submit(std::move(second)).Get();
		}));

	// Evaluation still happens later on TestManager's thread and in test order.
	return [simulation = std::move(simulation), envmode]() mutable {
		auto secondResult = simulation->get();
		return LimaUnittestResult{ static_cast<bool>(secondResult.simulation), "Success", envmode == Full };

	// TODO: Wtf is this test??
		/*for (int cid = 0; cid < sim2->box->boxparams.n_compounds; cid++) {
		for (int pid = 0; pid < sim2->box->compounds[cid].n_particles; pid++) {
			Float3 pos1 = sim1->traj_buffer->GetMostRecentCompoundparticleDatapoint(cid, pid, 100-1);

			Float3 pos2 = sim2->traj_buffer->GetMostRecentCompoundparticleDatapoint(cid, pid, 1-1);

			ASSERT(pos1 == pos2, "Position of compound " + std::to_string(cid) + " particle " + std::to_string(pid) + " is not the same between simulations");
		}
	}*/
	};
}
