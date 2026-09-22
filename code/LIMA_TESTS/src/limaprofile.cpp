#include "Benchmarks.h"
#include "Engine.cuh"
#include <cuda_profiler_api.h>

static int ProfileT4(int batchSize, int steps, int loggingInterval, int nlistInterval) {
	if (batchSize < 1 || steps < 1 || loggingInterval < 0 || nlistInterval < 1)
		throw std::invalid_argument("Invalid T4 profiling arguments");
	std::vector<std::unique_ptr<Simulation>> simulations;
	std::vector<Simulation*> members;
	for (int i = 0; i < batchSize; ++i) {
		SimulationJob job;
		job.workDir = TestUtils::AutomatedTestsDir() / "T4Lysozyme";
		job.mode = Headless;
		job.run = false;
		job.preprocess = [=](GroFile&, TopologyFile&, SimParams& params) {
			params.n_steps = steps;
			params.data_logging_interval = loggingInterval;
			params.stepsPerNlistupdate = nlistInterval;
		};
		auto result = Environment::Get().Submit(std::move(job)).Get();
		result.simulation->PrepareDataBuffers();
		members.push_back(result.simulation.get());
		simulations.push_back(std::move(result.simulation));
	}
	Engine engine(members);
	cudaDeviceSynchronize();
	cudaProfilerStart();
	while (!engine.IsFinished()) engine.step();
	cudaDeviceSynchronize();
	cudaProfilerStop();
	return 0;
}

static int ProfileBaselineT4() {
	SimulationJob job;
	job.workDir = TestUtils::AutomatedTestsDir() / "T4Lysozyme";
	job.mode = Headless;
	job.mustRunAlone = true;
	job.preprocess = [](GroFile&, TopologyFile&, SimParams& params) {
		params.n_steps = 4;
	};
	Environment::Get().Submit(std::move(job)).Get();
	return 0;
}

int main(int argc, char** argv) {
	if (argc == 2 && std::string_view(argv[1]) == "--baseline") {
		try {
			return ProfileBaselineT4();
		}
		catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
	}
	if (argc > 1 && std::string_view(argv[1]) == "--batch-t4") {
		try {
			if (argc != 6) throw std::invalid_argument("Usage: limaprofile --batch-t4 batchSize steps loggingInterval nlistInterval");
			return ProfileT4(std::stoi(argv[2]), std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]));
		}
		catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
	}
	auto result = Benchmarks::STMV(Environment::Get(), EnvMode::Headless, 300).RunToCompletion();
	result.printStatus();
	return result.success ? 0 : 1;
}
