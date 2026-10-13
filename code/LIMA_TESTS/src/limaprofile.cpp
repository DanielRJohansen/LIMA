#include "Tests.h"
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
			params.n_steps = steps;
		};
		handles.push_back(Environment::Get().Submit(std::move(job)));
	}
	for (auto& handle : handles)
		handle.Get();
	return 0;
}

int main(int argc, char** argv) {
	if (argc > 1 && std::string_view(argv[1]) == "--batch-t4") {
		try {
			if (argc != 4) throw std::invalid_argument("Usage: limaprofile --batch-t4 batchSize steps");
			return ProfileT4(std::stoi(argv[2]), std::stoi(argv[3]));
		}
		catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
	}
	if (argc > 1 && std::string_view(argv[1]) == "--t4") {
		if (argc != 4) { std::cerr << "Usage: limaprofile --t4 steps runs\n"; return 1; }
		auto result = Benchmarks::T4(Environment::Get(), EnvMode::Headless, std::stoi(argv[2]), std::stoi(argv[3])).RunToCompletion();
		result.printStatus();
		return result.success ? 0 : 1;
	}
	int stmvSteps = 300;
	LongRangeElectrostaticsMethod longrangeMethod = LongRangeElectrostaticsMethod::PME;
	int fmmOrder = 6;
	int fmmLeafParticles = 64;
	if (argc > 1 && (std::string_view(argv[1]) == "--stmv" || std::string_view(argv[1]) == "--stmv-fmm")) {
		if (argc < 3 || argc > 5) { std::cerr << "Usage: limaprofile --stmv[[-fmm]] steps [fmm_order] [fmm_leaf_particles]\n"; return 1; }
		if (std::string_view(argv[1]) == "--stmv-fmm") longrangeMethod = LongRangeElectrostaticsMethod::FMM;
		stmvSteps = std::stoi(argv[2]);
		if (argc >= 4) fmmOrder = std::stoi(argv[3]);
		if (argc >= 5) fmmLeafParticles = std::stoi(argv[4]);
	}
	auto result = Benchmarks::STMV(Environment::Get(), EnvMode::Headless, stmvSteps, 1, longrangeMethod, fmmOrder, fmmLeafParticles).RunToCompletion();
	result.printStatus();
	return result.success ? 0 : 1;
}
