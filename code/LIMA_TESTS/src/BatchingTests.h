#pragma once

#include "TestUtils.h"
#include "EngineBatchTests.h"

#include <array>
#include <latch>

namespace BatchingTests {
	using namespace TestUtils;

	static SimulationJob MakeT4Job(EnvMode envmode, bool mustRunAlone) {
		const fs::path workDir = AutomatedTestsDir() / "T4Lysozyme";
		SimulationJob job;
		job.workDir = workDir;
		job.groPath = workDir / "molecule/conf.gro";
		job.topPath = workDir / "molecule/topol.top";
		job.simParamsPath = workDir / "sim_params.txt";
		job.mode = Headless;// envmode;
		job.mustRunAlone = mustRunAlone;
		job.preprocess = [](GroFile&, TopologyFile&, SimParams& params) {
			params.n_steps = 4000;
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

	static TestRoutine TestFourT4BatchMatchReference(Environment& environment, EnvMode envmode) {
		auto referenceHandle = environment.Submit(MakeT4Job(envmode, true));

		std::array<SimulationHandle, 4> handles;
		for (auto& handle : handles)
			handle = environment.Submit(MakeT4Job(envmode, false));

		// Everything has now been submitted. Only now do we start awaiting.
		auto reference = co_await std::move(referenceHandle);

		if (!reference.simulation)
			co_return LimaUnittestResult{ false, "Reference T4 simulation did not complete", envmode == Full };

		const double soloSimulatedNs = reference.simulation->getStep() * reference.simulation->simParams.dt;
		const double soloNsPerDay = soloSimulatedNs / reference.engineTime.count() * 86400.;

		std::optional<int> batchId;
		double batchSimulatedNs = 0.;
		double batchWallTime = 0.;

		for (size_t index = 0; index < handles.size(); ++index) {
			auto result = co_await std::move(handles[index]);

			if (!result.simulation || result.execution.batchSize != 4
				|| (batchId && result.execution.batchId != *batchId)
				|| !HasIdenticalCoordinates(*reference.simulation, *result.simulation)) {
				co_return LimaUnittestResult{
					false,
					std::format("Batched T4 run {} did not execute as the expected batch", index),
					envmode == Full
				};
			}

			batchId = result.execution.batchId;
			batchSimulatedNs += result.simulation->getStep() * result.simulation->simParams.dt;
			batchWallTime = (std::max)(batchWallTime, result.engineTime.count());
		}

		const double batchNsPerDay = batchSimulatedNs / batchWallTime * 86400.;
		const double throughputRatio = batchNsPerDay / soloNsPerDay;
		const bool throughputImproved = std::isfinite(throughputRatio) && throughputRatio >= 2.;
		const std::string result = std::format(
			"Batch perf.: {:.2f}x ({:.2f} vs {:.2f} [ns/day])",
			throughputRatio, batchNsPerDay, soloNsPerDay);

		co_return LimaUnittestResult{
			throughputImproved,
			throughputImproved ? result : result + "; expected at least 2x",
			envmode == Full
		};
	}

	// Require shared-kernel execution, not merely equivalent sequential results.
	//static TestRoutine TestFourT4BatchMatchReference(Environment& environment, EnvMode envmode) {
	//	auto reference = environment.Submit(MakeT4Job(envmode, true));

	//	std::array<SimulationHandle, 4> handles;
	//	for (auto& handle : handles)
	//		handle = environment.Submit(MakeT4Job(envmode, false));

	//	co_await reference;
	//	for (auto& handle : handles)
	//		co_await handle;

	//	std::optional<int> batchId;
	//	double batchSimulatedNs = 0.;
	//	double batchWallTime = 0.;
	//	for (size_t index = 0; index < handles.size(); ++index) {
	//		auto result = std::move(handles[index]);
	//		if (!result.simulation || result.execution.batchSize != 4
	//			|| (batchId && result.execution.batchId != *batchId)
	//			|| !HasIdenticalCoordinates(*reference.simulation, *result.simulation)) {
	//			LimaUnittestResult{ false, std::format("Batched T4 run {} did not execute as the expected batch", index), envmode == Full };
	//		}
	//		batchId = result.execution.batchId;
	//		batchSimulatedNs += result.simulation->getStep() * result.simulation->simParams.dt;
	//		batchWallTime = (std::max)(batchWallTime, result.engineTime.count());
	//	}


	//
	//

	//	if (!reference.simulation)
	//		LimaUnittestResult{ false, "Reference T4 simulation did not complete", envmode == Full };
	//	const double soloSimulatedNs = reference.simulation->getStep() * reference.simulation->simParams.dt;
	//	const double soloNsPerDay = soloSimulatedNs / reference.engineTime.count() * 86400.;

	//	const double batchNsPerDay = batchSimulatedNs / batchWallTime * 86400.;
	//	const double throughputRatio = batchNsPerDay / soloNsPerDay;
	//	const bool throughputImproved = std::isfinite(throughputRatio) && throughputRatio >= 2.;
	//	const std::string result = std::format(
	//		"4 T4 batch matched ref.; throughput {:.2f}x ({:.2f} vs {:.2f} [ns/day])",
	//		throughputRatio, batchNsPerDay, soloNsPerDay);
	//	LimaUnittestResult{ throughputImproved, throughputImproved ? result : result + "; expected at least 2x", envmode == Full };
	//}

	inline SimulationJob MakeSmallJob(int particles, int steps, float dt, bool electrostatics) {
		auto job = MakeT4Job(Headless, false);
		job.configureSimulation = [=](Simulation& simulation) {
			auto small = EngineBatchTests::MakeSimulation(particles, steps, dt, electrostatics);
			simulation.box = std::move(small->box);
			simulation.simParams = small->simParams;
		};
		return job;
	}

	inline std::vector<SimulationHandle> SubmitTogether(Environment& environment, std::vector<SimulationJob> jobs) {
		// Ensure all submissions are visible before preprocessing finishes, without
		// relying on sleeps or the relative speed of the scheduler and test thread.
		auto submitted = std::make_shared<std::latch>(1);
		auto configure = std::move(jobs.front().configureSimulation);
		jobs.front().configureSimulation = [submitted, configure](Simulation& simulation) {
			submitted->wait();
			if (configure) configure(simulation);
		};
		std::vector<SimulationHandle> handles;
		try {
			for (auto& job : jobs) handles.push_back(environment.Submit(std::move(job)));
		}
		catch (...) { submitted->count_down(); throw; }
		submitted->count_down();
		return handles;
	}

	inline void RunSchedulerTests(Environment& environment) {
		using EngineBatchTests::Require;
		for (bool electrostatics : {false, true}) {
			std::vector<SimulationJob> jobs;
			std::vector<SimulationResult> references;
			for (int i = 0; i < 4; ++i) {
				auto reference = MakeSmallJob(5 + i * 4, i * 5, (i + 1) * 0.00001f, electrostatics);
				reference.mustRunAlone = true;
				references.push_back(environment.Submit(std::move(reference)).Get());
				jobs.push_back(MakeSmallJob(5 + i * 4, i * 5, (i + 1) * 0.00001f, electrostatics));
			}
			auto handles = SubmitTogether(environment, std::move(jobs));
			int batchId = -1;
			double previousEngineTime = 0.;
			for (size_t i = 0; i < handles.size(); ++i) {
				auto result = handles[i].Get();
				if (i == 0) batchId = result.execution.batchId;
				Require(result.execution.batchSize == 4 && result.execution.batchId == batchId, "Environment did not form a four-member batch");
				Require(result.simulation->finished, "Environment did not mark the member finished");
				EngineBatchTests::Compare(*references[i].simulation, *result.simulation);
				Require(result.engineTime.count() >= previousEngineTime, "Member timing did not track its retirement");
				previousEngineTime = result.engineTime.count();
			}
		}
		{
			std::vector<SimulationJob> jobs;
			for (int i = 0; i < 6; ++i) {
				auto job = MakeSmallJob(5, 3, 0.00001f, false);
				auto configure = std::move(job.configureSimulation);
				job.configureSimulation = [configure, i](Simulation& simulation) {
					configure(simulation);
					simulation.simParams.ref_t = i % 2 ? 310.f : 300.f;
				};
				jobs.push_back(std::move(job));
			}
			auto handles = SubmitTogether(environment, std::move(jobs));
			std::vector<SimulationResult> results;
			for (auto& handle : handles) results.push_back(handle.Get());
			for (size_t i = 0; i < results.size(); ++i) {
				Require(results[i].simulation->getStep() == 3, "Bounded scheduler failed to drain jobs");
				for (size_t j = 0; j < i; ++j)
					if ((i % 2) != (j % 2))
						Require(results[i].execution.batchId != results[j].execution.batchId, "Scheduler ignored configured parameter incompatibility");
			}
			Require(results[0].execution.batchSize == 3 && results[0].execution.batchId == results[2].execution.batchId
				&& results[2].execution.batchId == results[4].execution.batchId
				&& results[1].execution.batchId == results[3].execution.batchId
				&& results[3].execution.batchId == results[5].execution.batchId,
				"Scheduler started before queued compatible simulations finished preparing");
		}
		{
			std::vector<SimulationJob> jobs;
			for (int i = 0; i < 4; ++i) {
				auto job = MakeSmallJob(5, 3, 0.00001f, false);
				job.mustRunAlone = i == 1;
				job.run = i != 2;
				jobs.push_back(std::move(job));
			}
			auto handles = SubmitTogether(environment, std::move(jobs));
			for (size_t i = 0; i < handles.size(); ++i) {
				auto result = handles[i].Get();
				Require(result.execution.batchSize == 1, "Standalone/preparation-only job was batched");
				Require(result.simulation->finished == (i != 2), "Preparation-only job executed");
				Require(result.simulation->getStep() == (i == 2 ? 0 : 3), "Standalone job step count changed");
			}
		}
		{
			std::vector<SimulationJob> jobs;
			for (int i = 0; i < 4; ++i) jobs.push_back(MakeSmallJob(5, 3, 0.00001f, false));
			jobs[1].configureSimulation = [](Simulation&) { throw std::runtime_error("Expected preprocessing failure"); };
			jobs[2].postprocess = [](SimulationResult&) { throw std::runtime_error("Expected postprocessing failure"); };
			auto handles = SubmitTogether(environment, std::move(jobs));
			for (size_t i = 0; i < handles.size(); ++i) {
				bool failed = false;
				try { auto result = handles[i].Get(); Require(result.simulation->finished, "Healthy job did not finish"); }
				catch (const std::runtime_error&) { failed = true; }
				Require(failed == (i == 1 || i == 2), "Callback failure affected the wrong job");
			}
		}
		const auto t4 = TestFourT4BatchMatchReference(environment, Headless).RunToCompletion();
		Require(t4.success, t4.error_description.c_str());
		std::cout << "Environment batching tests passed\n";
	}
}
