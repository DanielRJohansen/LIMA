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

	// Profiles loading a large system: gro/top parsing, then BoxImage + Box building inside the Environment
	static TestRoutine Load3J3Q(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "3j3q";
		if (!fs::is_regular_file(workDir / "conf.gro") || !fs::is_regular_file(workDir / "topol.top"))
			co_return LimaUnittestResult{ false, "Missing 3j3q load benchmark input: " + workDir.string(), envmode == Full };


		TimeIt totalTimer;
		SimulationJob job;
		job.workDir = workDir;

		TimeIt groTimer;
		job.grofile.emplace(workDir / "conf.gro");
		const std::chrono::duration<double> groTime = groTimer.stop();

		TimeIt topTimer;
		job.topfile.emplace(workDir / "topol.top");
		const std::chrono::duration<double> topTime = topTimer.stop();

		job.simParams.emplace();
		job.mode = EnvMode::Headless;
		job.mustRunAlone = true;
		job.run = false;		
				
		const auto& top = *job.topfile;
		const size_t nTopAtoms = std::ranges::distance(top.GetAllElements<TopologyFile::AtomsEntry>());

		if (nTopAtoms != job.grofile->atoms.size())
			co_return LimaUnittestResult{ false, std::format("3j3q atom count mismatch between gro({}) and top({})", job.grofile->atoms.size(), nTopAtoms), envmode == Full};

		auto completed = co_await environment.Submit(std::move(job));		
		if (!completed.simulation)
			co_return LimaUnittestResult{ false, "3j3q load benchmark produced no simulation", envmode == Full };
		
		const std::chrono::duration<double> allowedFileTime{ 3.5 };
		const std::chrono::duration<double> allowedBuildTime{ 7.0 };

		auto buildtime = completed.environmentTime;
		auto success = (groTime + topTime) < allowedFileTime && buildtime < allowedBuildTime;

		co_return LimaUnittestResult{ success,
			std::format("files: ({:.2f}+{:.2f})/{:.2f}   build: {:.2f}/{:.2f} [s]",
				groTime.count(), topTime.count(), allowedFileTime.count(), buildtime.count(), allowedBuildTime.count()),
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
			job.grofile.emplace(groPath);
			job.topfile.emplace(topPath);
			job.simParams.emplace(simParamsPath);
			job.mode = EnvMode::Headless;
			job.mustRunAlone = true;
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
			workDir / "sim_params.txt", { std::chrono::microseconds{ 11500 }, std::chrono::microseconds{ 14000 } },
			nSteps, nRuns);
	}

	static TestRoutine T4(Environment& environment, EnvMode envmode, int nSteps = 500, int nRuns = 1) {
		const fs::path workDir = TestsDir() / "benchmarking/t4";
		return Bench(environment, envmode, workDir, workDir / "conf.gro", workDir / "topol.top",
			workDir / "../sim_params.txt", { std::chrono::microseconds{ 180 }, std::chrono::microseconds{ 300 } },
			nSteps, nRuns);
	}
}
