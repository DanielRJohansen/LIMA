#include "Tests.h"
#include "Format.h"
#include "Programs.h"
#include "TestUtils.h"
#include "TimeIt.h"

namespace Benchmarks {

	using namespace TestUtils;

	const fs::path TestsDir() {
		return HeavyTestsDir();
	}

	TestRoutine ToGmxLargeCif(Environment&, EnvMode envmode) {
		const fs::path input = TestsDir() / "fileconversions" / "3J3Q.cif";
		if (!fs::is_regular_file(input))
			co_return LimaUnittestResult::Fail("Missing ToGmx benchmark input: " + input.string(), envmode == Full);
		TimeIt timer;
		const auto conversion = Programs::ToGmx(input);
		const auto elapsed = timer.stop();
		const std::chrono::seconds allowedTime{ 12 };	// See TODO.md, was 10

		if (conversion.grofile.atoms.empty())
			co_return LimaUnittestResult::Fail("ToGmx benchmark produced no atoms", envmode == Full);
		co_return LimaUnittestResult{ elapsed < allowedTime,
			OneLineMessage("3J3Q.cif elapsed: {:.3f} allowed: {:.3f}",
				std::chrono::duration<double>(elapsed).count(), std::chrono::duration<double>(allowedTime).count()),
			envmode == Full };
	}

	// Profiles loading a large system: gro/top parsing, then BoxImage + Box building inside the Environment
	TestRoutine Load3J3Q(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "3j3q";
		if (!fs::is_regular_file(workDir / "conf.gro") || !fs::is_regular_file(workDir / "topol.top"))
			co_return LimaUnittestResult::Fail("Missing 3j3q load benchmark input: " + workDir.string(), envmode == Full);


		TimeIt totalTimer;
		SimulationJob job;
		job.workDir = workDir;

		// Parsed here rather than by the Environment, to time the files separately from the build
		TimeIt groTimer;
		const GroFile& gro = job.grofile.emplace(workDir / "conf.gro").Load();
		const std::chrono::duration<double> groTime = groTimer.stop();

		TimeIt topTimer;
		const TopologyFile& top = job.topfile.emplace(workDir / "topol.top").Load();
		const std::chrono::duration<double> topTime = topTimer.stop();

		job.simParams.emplace();
		job.mode = EnvMode::Headless;
		job.mustRunAlone = true;
		job.run = false;		
				
		const size_t nTopAtoms = std::ranges::distance(top.GetAllElements<TopologyFile::AtomsEntry>());

		if (nTopAtoms != gro.atoms.size())
			co_return LimaUnittestResult::Fail(Lima::Format("3j3q atom count mismatch between gro({}) and top({})", gro.atoms.size(), nTopAtoms), envmode == Full);

		auto completed = co_await environment.Submit(std::move(job));		
		if (!completed.simulation)
			co_return LimaUnittestResult::Fail("3j3q load benchmark produced no simulation", envmode == Full);
		
		const std::chrono::duration<double> allowedFileTime{ 3.5 };
		const std::chrono::duration<double> allowedBuildTime{ 10.0 };	// See TODO.md, was 7.0

		auto buildtime = completed.environmentTime;
		auto success = (groTime + topTime) < allowedFileTime && buildtime < allowedBuildTime;

		co_return LimaUnittestResult{ success,
			OneLineMessage("files {:.1f}+{:.1f}/{:.1f}s build {:.1f}/{:.1f}s",
				groTime.count(), topTime.count(), allowedFileTime.count(), buildtime.count(), allowedBuildTime.count()),
			envmode == Full };
	}

	// Profiles preparing a small system end to end inside the Environment: gro/top parsing, BoxImage + Box building,
	// Engine construction and the first step
	TestRoutine LoadT4(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "benchmarking/t4";

		SimulationJob job;
		job.workDir = workDir;
		job.grofile.emplace(workDir / "conf.gro");
		job.topfile.emplace(workDir / "topol.top");
		job.simParams.emplace(workDir / "../sim_params.txt");
		job.mode = EnvMode::Headless;
		job.mustRunAlone = true;
		job.preprocess = [](GroFile&, TopologyFile&, SimParams& params) { params.n_steps = 1; };

		auto completed = co_await environment.Submit(std::move(job));
		if (!completed.simulation || completed.simulation->getStep() != 1)
			co_return LimaUnittestResult::Fail("T4 load benchmark did not complete its step", envmode == Full);

		// environmentTime excludes time spent queued behind other jobs
		const std::chrono::duration<double> setupTime = completed.environmentTime;
		const std::chrono::duration<double> allowedTime{ 0.5 };

		co_return LimaUnittestResult{ setupTime < allowedTime,
			OneLineMessage("setup {:.0f} (1st step {:.1f}) / {:.0f} [ms]",
				setupTime.count() * 1000., completed.engineTime.count() * 1000., allowedTime.count() * 1000.),
			envmode == Full };
	}

	TestRoutine Bench(Environment& environment, EnvMode envmode, fs::path workDir,
		fs::path groPath, fs::path topPath, fs::path simParamsPath,
		PerformanceBounds<std::chrono::microseconds> allowedTimePerStep, int nSteps, int nRuns, int warmupSteps, bool printPerformance, bool barostat,
		bool constrained) {
		const auto MakeJob = [&](int steps) {
			SimulationJob job;
			job.workDir = workDir;
			job.grofile.emplace(groPath);
			job.topfile.emplace(topPath);
			job.simParams.emplace(simParamsPath);
			job.mode = EnvMode::Headless;
			job.mustRunAlone = true;
			job.preprocess = [steps, barostat, constrained](GroFile&, TopologyFile&, SimParams& params) {
				params.data_logging_interval = 20;
				params.enable_electrostatics = true;
				params.n_steps = steps;
				params.apply_barostat = barostat;
				if (constrained)
					params.constraints = ConstraintsSelect::HBonds;
			};
			return job;
		};

		// The GPU drops to idle clocks within ~1 s without work, and small systems load it too lightly
		// to ramp back up quickly. An untimed run first brings it back to full clocks.
		if (warmupSteps > 0) {
			auto warmup = co_await environment.Submit(MakeJob(warmupSteps));
			if (!warmup.simulation || warmup.simulation->getStep() != warmupSteps)
				co_return LimaUnittestResult::Fail("Warmup simulation did not run fully", envmode != Headless);
		}

		std::vector<std::chrono::microseconds> timesPerStep;
		timesPerStep.reserve(nRuns);
		// Await each run before submitting the next, so a run's preparation never overlaps the previous run's timed steps
		for (int run = 0; run < nRuns; run++) {
			auto completed = co_await environment.Submit(MakeJob(nSteps));
			if (!completed.simulation || completed.simulation->getStep() != completed.simulation->simParams.n_steps)
				co_return LimaUnittestResult::Fail("Simulation did not run fully", envmode != Headless);
			timesPerStep.push_back(std::chrono::duration_cast<std::chrono::microseconds>(completed.engineTime / nSteps));
			if (printPerformance)
				std::cout << StringUtils::FormatPerformance(std::chrono::duration<double>(completed.engineTime).count(),
					static_cast<double>(nSteps) * completed.simulation->simParams.dt);
			if (printPerformance && barostat && !completed.simulation->pressure_buffer.empty())
				std::cout << Lima::Format("Pressure: {:.1f} bar, box: {} nm\n", completed.simulation->pressure_buffer.back(),
					completed.simulation->box->boxparams.BoxSizeFloat().toString());
		}

		const auto [fastest, slowest] = std::minmax_element(timesPerStep.begin(), timesPerStep.end());
		const bool withinBounds = *fastest >= allowedTimePerStep.min && *slowest <= allowedTimePerStep.max;
		co_return LimaUnittestResult{ withinBounds,
			OneLineMessage("({:.3f}-{:.3f}) / ({:.3f}-{:.3f}) [ms/step]",
				fastest->count() / 1000., slowest->count() / 1000., allowedTimePerStep.min.count() / 1000.,
				allowedTimePerStep.max.count() / 1000.), envmode != Headless };
	}

	TestRoutine STMV(Environment& environment, EnvMode envmode, int nSteps, int nRuns, bool printPerformance, bool barostat, bool constrained) {
		const fs::path workDir = TestsDir() / "benchmarking/stmv";
		return Bench(environment, envmode, workDir, workDir / "conf.gro", workDir / "topol.top",
			workDir / "sim_params.txt", { std::chrono::microseconds{ 2700 }, std::chrono::microseconds{ 3000 } },
			nSteps, nRuns, 0, printPerformance, barostat, constrained);
	}

	TestRoutine T4(Environment& environment, EnvMode envmode, int nSteps, int nRuns, bool printPerformance) {
		const fs::path workDir = TestsDir() / "benchmarking/t4";
		return Bench(environment, envmode, workDir, workDir / "conf.gro", workDir / "topol.top",
			workDir / "../sim_params.txt", { std::chrono::microseconds{ 45 }, std::chrono::microseconds{ 200 } },
			nSteps, nRuns, 2000, printPerformance);
	}
}
