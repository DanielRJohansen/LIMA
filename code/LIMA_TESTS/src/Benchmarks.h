#include "Programs.h"
#include "TestUtils.h"
#include "TimeIt.h"

namespace Benchmarks {

	using namespace TestUtils;
	namespace fs = std::filesystem;
	constexpr int automatedTestRuns = 3;

	template<typename Duration>
	struct PerformanceBounds {
		Duration min;
		Duration max;
	};

	const fs::path TestsDir() {
		return HeavyTestsDir();
	}

	static TestRoutine ToGmxLargeCif(Environment&, EnvMode envmode) {
		const fs::path input = TestsDir() / "fileconversions" / "3J3Q.cif";
		if (!fs::is_regular_file(input))
			co_return LimaUnittestResult{ false, "Missing ToGmx benchmark input: " + input.string(), envmode == Full };
		TimeIt timer;
		const auto conversion = Programs::ToGmx(input);
		const auto elapsed = timer.stop();
		const std::chrono::seconds allowedTime{ 10 };

		if (conversion.grofile.atoms.empty())
			co_return LimaUnittestResult{ false, "ToGmx benchmark produced no atoms", envmode == Full };
		co_return LimaUnittestResult{ elapsed < allowedTime,
			std::format("3J3Q.cif elapsed: {:.3f} allowed: {:.3f}",
				std::chrono::duration<double>(elapsed).count(), std::chrono::duration<double>(allowedTime).count()),
			envmode == Full };
	}

	static TestRoutine Bench(Environment& environment, EnvMode envmode, fs::path workDir,
		fs::path groPath, fs::path topPath, fs::path simParamsPath,
		PerformanceBounds<std::chrono::microseconds> allowedTimePerStep, int nSteps, int nRuns)
	{
		std::vector<SimulationHandle> handles;
		handles.reserve(nRuns);
		for (int run = 0; run < nRuns; run++) {
			SimulationJob job;
			job.workDir = workDir;
			job.groPath = groPath;
			job.topPath = topPath;
			job.simParamsPath = simParamsPath;
			job.mode = EnvMode::Headless;
			job.preprocess = [nSteps](GroFile&, TopologyFile&, SimParams& params) {
				params.data_logging_interval = 20;
				params.enable_electrostatics = true;
				params.n_steps = nSteps;
			};

			handles.push_back(environment.Submit(std::move(job)));
		}
		std::vector<std::chrono::microseconds> timesPerStep;
		timesPerStep.reserve(nRuns);
		for (auto& handle : handles) {
			auto completed = co_await std::move(handle);
			if (!completed.simulation || completed.simulation->getStep() != completed.simulation->simParams.n_steps)
				co_return LimaUnittestResult{ false, "Simulation did not run fully", envmode != Headless };
			timesPerStep.push_back(std::chrono::duration_cast<std::chrono::microseconds>(completed.engineTime / nSteps));
		}

		const auto [fastest, slowest] = std::minmax_element(timesPerStep.begin(), timesPerStep.end());
		const bool withinBounds = *fastest >= allowedTimePerStep.min && *slowest <= allowedTimePerStep.max;
		co_return LimaUnittestResult{ withinBounds,
			std::format("({:.3f}-{:.3f}) / ({:.3f}-{:.3f}) [ms/step]",
				fastest->count() / 1000., slowest->count() / 1000., allowedTimePerStep.min.count() / 1000.,
				allowedTimePerStep.max.count() / 1000.), envmode != Headless };
	}

	static TestRoutine STMV(Environment& environment, EnvMode envmode, int nSteps, int nRuns = 1) {
		const fs::path workDir = TestsDir() / "benchmarking/stmv";
		return Bench(environment, envmode, workDir, workDir / "conf.gro", workDir / "topol.top",
			workDir / "sim_params.txt", { std::chrono::microseconds{ 12000 }, std::chrono::microseconds{ 14000 } },
			nSteps, nRuns);
	}

	static TestRoutine T4(Environment& environment, EnvMode envmode, int nSteps = 500, int nRuns = 1) {
		const fs::path workDir = TestsDir() / "benchmarking/t4";
		return Bench(environment, envmode, workDir, workDir / "conf.gro", workDir / "topol.top",
			workDir / "../sim_params.txt", { std::chrono::microseconds{ 200 }, std::chrono::microseconds{ 300 } },
			nSteps, nRuns);
	}
}
