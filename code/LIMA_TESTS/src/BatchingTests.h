#pragma once

#include "TestUtils.h"

#include <array>

namespace BatchingTests {
	using namespace TestUtils;

	static SimulationJob MakeT4Job(EnvMode envmode, bool mustRunAlone) {
		const fs::path workDir = AutomatedTestsDir() / "T4Lysozyme";
		SimulationJob job;
		job.workDir = workDir;
		job.groPath = workDir / "molecule/conf.gro";
		job.topPath = workDir / "molecule/topol.top";
		job.simParamsPath = workDir / "sim_params.txt";
		job.mode = envmode;
		job.mustRunAlone = mustRunAlone;
		job.preprocess = [](GroFile&, TopologyFile&, SimParams& params) {
			params.n_steps = 10;
			params.data_logging_interval = 2; // must NOT be 0
		};
		return job;
	}

	static bool HasIdenticalCoordinates(const Simulation& lhs, const Simulation& rhs) {
		const auto& lhsPclusters = lhs.box->persistentClusters;
		const auto& rhsPclusters = rhs.box->persistentClusters;
		if (lhsPclusters.size() != rhsPclusters.size())
			return false;

		for (size_t pcId = 0; pcId < lhsPclusters.size(); ++pcId) {
			for (int particleId = 0; particleId < PersistentCluster::maxParticles; ++particleId) {
				if (lhsPclusters[pcId].pqd[particleId].position != rhsPclusters[pcId].pqd[particleId].position)
					return false;
			}
		}
		return true;
	}

	// Enable this once Environment executes a real shared-kernel batch. It deliberately fails
	// against the current single-simulation scheduler rather than accepting sequential execution.
	static TestRoutine TestFourT4BatchMatchReference(Environment& environment, EnvMode envmode) {
		auto reference = co_await environment.Submit(MakeT4Job(envmode, true));
		if (!reference.simulation)
			co_return LimaUnittestResult{ false, "Reference T4 simulation did not complete", envmode == Full };

		std::array<SimulationHandle, 4> handles;
		for (auto& handle : handles)
			handle = environment.Submit(MakeT4Job(envmode, false));

		std::optional<int> batchId;
		for (size_t index = 0; index < handles.size(); ++index) {
			auto result = co_await std::move(handles[index]);
			if (!result.simulation || result.execution.batchSize != 4
				|| (batchId && result.execution.batchId != *batchId)
				|| !HasIdenticalCoordinates(*reference.simulation, *result.simulation)) {
				co_return LimaUnittestResult{ false, std::format("Batched T4 run {} did not execute as the expected batch", index), envmode == Full };
			}
			batchId = result.execution.batchId;
		}

		co_return LimaUnittestResult{ true, "Four T4 simulations shared one batch and matched the isolated reference", envmode == Full };
	}
}
