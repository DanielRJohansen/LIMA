#include "Tests.h"

#include <array>
#include <cstring>
#include <limits>

static int ProfileT4(int batchSize, int steps) {
	if (batchSize < 1 || steps < 1)
		throw std::invalid_argument("Invalid T4 profiling arguments");
	std::vector<SimulationHandle> handles;
	handles.reserve(batchSize);
	for (int i = 0; i < batchSize; ++i) {
		SimulationJob job;
		job.workDir = TestUtils::AutomatedTestsDir() / "T4Lysozyme";
		job.mode = Headless;
		job.profileCuda = true;
		job.preprocess = [steps](GroFile&, TopologyFile&, SimParams& params) {
			params.nsteps = steps;
		};
		handles.push_back(Environment::Get().Submit(std::move(job)));
	}
	for (auto& handle : handles)
		handle.Get();
	return 0;
}

// Runs 1, 2, 4, .. copies of T4 lysozyme (34k atoms, 7 nm box, with PME) as one batch each, up to the Environment's
// max batch size, and prints how the combined throughput grows with the batch size
static int ProfileBatchScaling(int steps) {
	if (steps < 1)
		throw std::invalid_argument("Invalid batch scaling arguments");
	Environment& environment = Environment::Get();
	const auto MakeJob = [](int nSteps) {
		SimulationJob job;
		job.workDir = TestUtils::AutomatedTestsDir() / "T4Lysozyme";
		job.mode = Headless;
		job.preprocess = [nSteps](GroFile&, TopologyFile&, SimParams& params) {
			params.nsteps = nSteps;
			params.enable_electrostatics = true;
			params.nstenergy = 50;
		};
		return job;
	};

	// Warm up CUDA and the file caches, so the first measurement is not penalized
	environment.Submit(MakeJob(100)).Get();

	std::cout << Lima::Format("{:>6} {:>12} {:>12} {:>14} {:>16} {:>9}\n",
		"batch", "engine [s]", "wall [s]", "ns/day per sim", "ns/day combined", "speedup");
	double singleThroughput = 0.0;
	// Submits n jobs at once and prints their throughput if they ran as one batch. Returns the size of the first batch
	const auto Measure = [&](int n) {
		const auto start = std::chrono::steady_clock::now();
		std::vector<SimulationHandle> handles;
		for (int i = 0; i < n; ++i)
			handles.push_back(environment.Submit(MakeJob(steps)));
		std::vector<SimulationResult> results;
		for (auto& handle : handles)
			results.push_back(handle.Get());
		const double wallTime = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

		const int firstBatch = results.front().execution.batchSize;
		if (firstBatch != n)
			return firstBatch;
		double engineTime = 0.0;	// [s] the slowest simulation of the batch
		for (const auto& result : results) {
			if (!result.simulation || result.simulation->getStep() != steps)
				throw std::runtime_error("A batch scaling simulation did not run fully");
			engineTime = (std::max)(engineTime, std::chrono::duration<double>(result.engineTime).count());
		}
		const double simulatedNs = static_cast<double>(steps) * results.front().simulation->simParams.dt;
		const double perSimulation = simulatedNs / engineTime * 86400.0;	// [ns/day]
		const double combined = perSimulation * n;
		if (n == 1) singleThroughput = combined;
		std::cout << Lima::Format("{:>6} {:>12.2f} {:>12.2f} {:>14.1f} {:>16.1f} {:>8.2f}x\n",
			n, engineTime, wallTime, perSimulation, combined, combined / singleThroughput);
		return n;
	};

	// Double the batch until the Environment's max batch size splits the jobs, then measure a batch of exactly that size
	int previous = 0;
	for (int n = 1;; n *= 2) {
		const int batchSize = Measure(n);
		if (batchSize < n) {
			if (batchSize > previous)
				Measure(batchSize);
			std::cout << Lima::Format("Max batch size for this system: {}\n", batchSize);
			break;
		}
		previous = n;
	}
	return 0;
}

// FNV-1a over the bits of a run's temperatures, one per thermostat coupling. Equal between two builds only if their dynamics are
// bitwise identical, so a performance change that should not change the numerics can be checked on the benchmark systems
static uint64_t TemperatureChecksum(const std::vector<float>& temperatures) {
	uint64_t hash = 14695981039346656037ull;
	for (const float temperature : temperatures) {
		uint32_t bits = 0;
		std::memcpy(&bits, &temperature, sizeof(bits));
		for (int byte = 0; byte < 4; byte++)
			hash = (hash ^ ((bits >> (8 * byte)) & 0xFF)) * 1099511628211ull;
	}
	return hash;
}

// Runs a system of the GROMACS benchmark set (LIMA_data/salad_gromacs_bench) with its own lima.mdp, changing only nsteps.
// warmupSteps > 0 first runs an untimed simulation, so the GPU is back at full clocks. Returns the fastest run's time per step [ms],
// and the last run's TemperatureChecksum in checksum
static double ProfileSalad(const std::string& name, int steps, int nRuns, int warmupSteps, uint64_t* checksum = nullptr) {
	if (steps < 1 || nRuns < 1)
		throw std::invalid_argument("Invalid salad profiling arguments");
	const std::filesystem::path dir = TestUtils::HeavyTestsDir() / "salad_gromacs_bench" / name;
	if (!std::filesystem::exists(dir / "lima.mdp"))
		throw std::invalid_argument(Lima::Format("Unknown salad system {} (no {})", name, (dir / "lima.mdp").string()));
	const auto MakeJob = [&dir](int nSteps) {
		SimulationJob job;
		job.workDir = dir;
		job.grofile.emplace(dir / "conf.gro");
		job.topfile.emplace(dir / "topol.top");
		job.simParams.emplace(dir / "lima.mdp");
		job.mode = Headless;
		job.mustRunAlone = true;
		job.preprocess = [nSteps](GroFile&, TopologyFile&, SimParams& params) { params.nsteps = nSteps; };
		return job;
	};

	Environment& environment = Environment::Get();
	if (warmupSteps > 0)
		environment.Submit(MakeJob(warmupSteps)).Get();
	double fastest = std::numeric_limits<double>::max();
	// Await each run before submitting the next, so a run's preparation never overlaps the previous run's timed steps
	for (int run = 0; run < nRuns; run++) {
		const SimulationResult result = environment.Submit(MakeJob(steps)).Get();
		if (!result.simulation || result.simulation->getStep() != steps)
			throw std::runtime_error(Lima::Format("{} did not run fully", name));
		const double msPerStep = std::chrono::duration<double, std::milli>(result.engineTime).count() / steps;
		fastest = (std::min)(fastest, msPerStep);
		if (checksum)
			*checksum = TemperatureChecksum(result.simulation->temperature_buffer);
	}
	return fastest;
}

// The performance suite: a small (villin, 5k atoms, launch-bound), a medium (adh_dodec, 191k) and a large (stmv, 1.07M,
// throughput-bound) system, all with the same physics (PME, h-bond constraints, v-rescale, 2 fs). Step counts give
// each run ~1.5 s. Prints the fastest of nRuns runs per system
static int ProfileSuite(int nRuns) {
	struct System { std::string name; int steps; };
	const std::array<System, 3> systems{ { { "villin", 20000 }, { "adh_dodec", 3000 }, { "stmv", 500 } } };
	std::cout << Lima::Format("{:>10} {:>10} {:>10}\n", "system", "ms/step", "ns/day");
	for (const auto& system : systems) {
		const double msPerStep = ProfileSalad(system.name, system.steps, nRuns, system.steps / 4);
		const double dtPs = 0.002;	// lima.mdp of every salad system
		std::cout << Lima::Format("{:>10} {:>10.4f} {:>10.1f}\n", system.name, msPerStep, dtPs * 86400.0 / msPerStep);
	}
	return 0;
}

int main(int argc, char** argv) {
	if (argc > 1 && std::string_view(argv[1]) == "--suite") {
		try {
			if (argc > 3) throw std::invalid_argument("Usage: limaprofile --suite [runs=3]");
			return ProfileSuite(argc > 2 ? std::stoi(argv[2]) : 3);
		}
		catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
	}
	// A single run without warmup, eg. for ncu. The checksum tells whether two builds' dynamics are bitwise identical
	if (argc > 1 && std::string_view(argv[1]) == "--salad") {
		try {
			if (argc != 4) throw std::invalid_argument("Usage: limaprofile --salad villin|rnase_cubic|adh_dodec|stmv steps");
			uint64_t checksum = 0;
			const double msPerStep = ProfileSalad(argv[2], std::stoi(argv[3]), 1, 0, &checksum);
			std::cout << Lima::Format("{:.4f} ms/step, temperature checksum {:016x}\n", msPerStep, checksum);
			return 0;
		}
		catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
	}
	if (argc > 1 && std::string_view(argv[1]) == "--batch-scaling") {
		try {
			if (argc > 3) throw std::invalid_argument("Usage: limaprofile --batch-scaling [steps=2000]");
			return ProfileBatchScaling(argc > 2 ? std::stoi(argv[2]) : 2000);
		}
		catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
	}
	if (argc > 1 && std::string_view(argv[1]) == "--batch-t4") {
		try {
			if (argc != 4) throw std::invalid_argument("Usage: limaprofile --batch-t4 batchSize steps");
			return ProfileT4(std::stoi(argv[2]), std::stoi(argv[3]));
		}
		catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
	}
	if (argc > 1 && std::string_view(argv[1]) == "--t4") {
		if (argc != 4) { std::cerr << "Usage: limaprofile --t4 steps runs\n"; return 1; }
		auto result = Benchmarks::T4(Environment::Get(), EnvMode::Headless, std::stoi(argv[2]), std::stoi(argv[3]), true).RunToCompletion();
		result.printStatus();
		return result.success ? 0 : 1;
	}
	int stmvSteps = 1000;
	bool barostat = false, constrained = false;
	if (argc > 1 && (std::string_view(argv[1]) == "--stmv" || std::string_view(argv[1]) == "--stmv-barostat" || std::string_view(argv[1]) == "--stmv-constrained")) {
		if (argc != 3) { std::cerr << "Usage: limaprofile --stmv[-barostat|-constrained] steps\n"; return 1; }
		barostat = std::string_view(argv[1]) == "--stmv-barostat";
		constrained = std::string_view(argv[1]) == "--stmv-constrained";
		stmvSteps = std::stoi(argv[2]);
	}
	auto result = Benchmarks::STMV(Environment::Get(), EnvMode::Headless, stmvSteps, 1, true, barostat, constrained).RunToCompletion();
	result.printStatus();
	return result.success ? 0 : 1;
}
