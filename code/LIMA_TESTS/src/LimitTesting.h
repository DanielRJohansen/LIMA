#pragma once

#include "Engine.cuh"
#include "TestUtils.h"

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string_view>
#include <thread>

#ifdef _WIN32
#include <Windows.h>
#include <crtdbg.h>
#else
#include <fcntl.h>
#include <signal.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char** environ;
#endif

// Reproduce a failure with: limatest --limit-case NAME
// Run just this suite with: limatest --limit-tests [NAME-SUBSTRING]
// For memory diagnostics, run an individual case under compute-sanitizer --tool memcheck.
// Each case guards an engine limit or contract that has broken before. Cases that still fail expose
// open engine defects; they remain failures until the engine handles the input safely.
// Never turn a crash, timeout, generic CUDA error, or nonfinite result into an expected pass.
namespace LimitTesting {

	enum class Kind { Transfer, Occupancy, ChargeBlock, Bootstrap, Pair, InvalidInput, Logging, Thermostat, Unwind };
	struct Case {
		std::string name;
		Kind kind;
		int value = 0;
	};

	inline std::vector<Case> Cases() {
		std::vector<Case> cases;
		// Capacities are probed at the limit and one past it
		for (int count : {8, 9}) cases.push_back({std::format("transfer-{}", count), Kind::Transfer, count});
		for (int count : {128, 129}) {
			cases.push_back({std::format("runtime-occupancy-{}", count), Kind::Occupancy, count});
			cases.push_back({std::format("bootstrap-occupancy-{}", count), Kind::Bootstrap, count});
		}
		for (int count : {384, 385}) cases.push_back({std::format("pme-entries-{}", count), Kind::ChargeBlock, count});
		for (int separation : {200, 10}) cases.push_back({std::format("pair-md-{}pm", separation), Kind::Pair, separation});
		cases.push_back({"invalid-input", Kind::InvalidInput});
		for (int steps : {0, 1, 15, 16, 31}) cases.push_back({std::format("logging-steps-{}", steps), Kind::Logging, steps});
		cases.push_back({"thermostat-logging-independent", Kind::Thermostat});
		cases.push_back({"cuda-failure-unwind", Kind::Unwind});
		return cases;
	}

	class Failure : public std::runtime_error {
	public:
		using std::runtime_error::runtime_error;
	};

	inline void Require(bool condition, const std::string& message) {
		if (!condition) throw Failure(message);
	}

	inline void Phase(const std::string& message) {
		std::cout << "PHASE " << message << std::endl;
	}

	inline void CheckCuda() {
		LIMA_UTILS::genericErrorCheck(cudaDeviceSynchronize());
		LIMA_UTILS::genericErrorCheck(cudaGetLastError());
	}

	inline bool Finite(const Float3& value) {
		return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
	}

	inline std::unique_ptr<Simulation> MakeSimulation(int particles, int steps = 3, int loggingInterval = 0) {
		SimParams params;
		params.n_steps = steps;
		params.dt = 0.1f * FEMTO_TO_NANO;
		params.enable_electrostatics = false;
		params.apply_thermostat = false;
		params.stepsPerNlistupdate = 1;
		params.steps_per_temperature_measurement = 1;
		params.data_logging_interval = loggingInterval;
		auto simulation = std::make_unique<Simulation>(params, std::make_unique<Box>(Float3{4.f}));
		auto& box = *simulation->box;
		box.boxparams.totalParticles = particles;
		box.boxparams.degreesOfFreedom = particles * 3;
		// One atom per persistent cluster makes occupancy independent of particle packing.
		box.persistentClusters.resize(particles);
		box.persistentClustersMetadata.resize(particles);
		box.pclusterInterimStates.resize(particles);
		box.particlesBondedToParticle.resize(particles);
		box.pclustersBondedToPcluster.resize(particles);
		for (int i = 0; i < particles; ++i) {
			box.persistentClusters[i].pqd[0] = PData{
				Float3{1.2f + 0.07f * (i % 5), 1.2f + 0.07f * ((i / 5) % 5), 1.2f + 0.07f * (i / 25)}, NBParams{0.f, 0.f, 0.f}};
			auto& meta = box.persistentClustersMetadata[i];
			meta.nParticles = 1;
			meta.mass[0] = 12.f;
			meta.particleIdsGlobal[0] = i;
		}
		simulation->PrepareDataBuffers();
		return simulation;
	}

	inline void VerifyState(Engine& engine, Simulation& simulation) {
		CheckCuda();
		engine.CopySimulationToHost();
		const auto clusters = engine.OffloadPclusterState().GetData();
		CheckCuda();
		for (size_t pc = 0; pc < simulation.box->persistentClusters.size(); ++pc) {
			const auto& state = simulation.box->pclusterInterimStates[pc];
			for (int lane = 0; lane < simulation.box->persistentClustersMetadata[pc].nParticles; ++lane) {
				Require(clusters[pc].pqd[lane].Valid(), "Particle disappeared");
				Require(Finite(clusters[pc].pqd[lane].position), std::format("Nonfinite position at step {}, cluster {}", simulation.getStep(), pc));
				Require(Finite(state.vels_prev[lane]) && Finite(state.forces_prev[lane]), std::format("Nonfinite velocity/force at step {}, cluster {}", simulation.getStep(), pc));
			}
		}
		const auto forceMagnitudes = engine.OffloadForcesMagnitudeBuffer().GetData();
		CheckCuda();
		for (float force : forceMagnitudes) Require(std::isfinite(force), "Nonfinite force magnitude");
	}

	inline void Run(Simulation& simulation) {
		Phase("construct engine");
		{
			Engine engine({&simulation});
			while (!engine.IsFinished()) {
				Phase(std::format("step {}", simulation.getStep() + 1));
				engine.step();
				VerifyState(engine, simulation);
			}
			Phase("finalize");
			engine.terminateSimulation();
			CheckCuda();
		}
		Phase("engine destroyed");
		CheckCuda();
	}

	inline void Pair(int separationPm) {
		auto simulation = MakeSimulation(2, 8, 1);
		auto& box = *simulation->box;
		for (auto& cluster : box.persistentClusters)
			cluster.pqd[0] = PData{Float3{1.5f, 1.5f, 1.5f}, NBParams{0.1f, 0.1f, 0.f}};
		box.persistentClusters[1].pqd[0].position.x += separationPm * 0.001f;
		const double initialDistance = (box.persistentClusters[1].pqd[0].position - box.persistentClusters[0].pqd[0].position).len();
		try {
			Run(*simulation);
			for (float energy : simulation->potE_buffer->GetBuffer()) Require(std::isfinite(energy), "Nonfinite logged potential energy");
			// Finite output alone misses overflow in the fixed-point MD accumulator. Compare the
			// first force and energy (evaluated at the initial coordinates) to the analytic pair.
			const double s = std::pow(0.2 / initialDistance, 6);
			const double expectedForce = -24. * 0.01 * s * (2. * s - 1.) / initialDistance;
			const double expectedEnergy = 4. * 0.01 * s * (s - 1.);
			const double actualForce = simulation->forceBuffer->GetDatapoint(0, 0, 0).x;
			const double actualEnergy = simulation->potE_buffer->GetDatapoint(0, 0, 0) + simulation->potE_buffer->GetDatapoint(1, 0, 0);
			std::cout << "Initial pair force expected=" << expectedForce << " actual=" << actualForce
				<< " energy expected=" << expectedEnergy << " actual=" << actualEnergy << std::endl;
			Require(std::abs(actualForce - expectedForce) <= (std::max)(1e-4, std::abs(expectedForce) * 1e-3), "Finite MD force disagrees with analytic Lennard-Jones force (possible accumulator overflow)");
			Require(std::abs(actualEnergy - expectedEnergy) <= (std::max)(1e-4, std::abs(expectedEnergy) * 1e-3), "Finite MD energy disagrees with analytic Lennard-Jones energy");
		}
		catch (const Failure&) { throw; }
		catch (const std::exception& error) {
			// Allow a future explicit, clean rejection of the close pair. Generic CUDA failures never qualify.
			if (separationPm == 200) throw;
			const std::string message = error.what();
			if (message.find("overlap") == std::string::npos && message.find("non-finite") == std::string::npos) throw;
			CheckCuda();
			std::cout << "Clean numerical rejection: " << message << std::endl;
		}
	}

	inline void InvalidInput() {
		const std::array<std::pair<const char*, void(*)(Simulation&)>, 8> inputs{{
			{"zero timestep", [](Simulation& s) { s.simParams.dt = 0.f; }},
			{"nan timestep", [](Simulation& s) { s.simParams.dt = std::numeric_limits<float>::quiet_NaN(); }},
			{"zero cutoff", [](Simulation& s) { s.simParams.cutoff_nm = 0.f; }},
			{"zero mass", [](Simulation& s) { s.box->persistentClustersMetadata[0].mass[0] = 0.f; }},
			{"nan position", [](Simulation& s) { s.box->persistentClusters[0].pqd[0].position.x = std::numeric_limits<float>::quiet_NaN(); }},
			{"zero nlist interval", [](Simulation& s) { s.simParams.stepsPerNlistupdate = 0; }},
			{"negative logging interval", [](Simulation& s) { s.simParams.data_logging_interval = -1; }},
			{"zero temperature interval", [](Simulation& s) { s.simParams.steps_per_temperature_measurement = 0; }},
		}};
		for (const auto& [name, mutate] : inputs) {
			// Mutate after buffer preparation to test the engine entry contract separately from builders.
			auto simulation = MakeSimulation(2);
			mutate(*simulation);
			bool rejected = false;
			try { Engine engine({simulation.get()}); }
			catch (const std::invalid_argument& error) {
				rejected = true;
				std::cout << "Rejected " << name << ": " << error.what() << std::endl;
			}
			Require(rejected, std::string("Engine accepted invalid input: ") + name);
			CheckCuda();
		}
		// A rejected job must not leave the process unable to run a subsequent valid job.
		auto healthy = MakeSimulation(2);
		Run(*healthy);
	}

	inline void ProbeCapacity(EngineLimitProbe probe, int count, int capacity) {
		try { Engine::TestLimit(probe, count); }
		catch (const std::length_error& error) {
			// A deliberate capacity rejection is acceptable above the limit. Fixture invariant
			// failures and generic CUDA errors are runtime_errors and must still fail the case.
			Require(count > capacity, "Engine rejected a within-capacity control");
			CheckCuda();
			std::cout << "Clean capacity rejection: " << error.what() << std::endl;
			auto healthy = MakeSimulation(2);
			Run(*healthy);
		}
	}

	inline void RunCase(const Case& test) {
		switch (test.kind) {
		case Kind::Transfer: ProbeCapacity(EngineLimitProbe::ClusterTransfer, test.value, 8); break;
		case Kind::Occupancy: ProbeCapacity(EngineLimitProbe::ClusterOccupancy, test.value, 128); break;
		case Kind::ChargeBlock: ProbeCapacity(EngineLimitProbe::ChargeBlock, test.value, 384); break;
		case Kind::Bootstrap: {
			auto simulation = MakeSimulation(test.value);
			if (test.value <= 128) Run(*simulation);
			else {
				bool rejected = false;
				try { Engine engine({simulation.get()}); }
				catch (const std::runtime_error& error) {
					rejected = std::string_view(error.what()).find("Too many pclusters in bootstrap grid bin") != std::string_view::npos;
					if (!rejected) throw;
				}
				Require(rejected, "Bootstrap did not reject excess cluster occupancy");
				CheckCuda();
				auto healthy = MakeSimulation(2);
				Run(*healthy);
			}
			break;
		}
		case Kind::Pair: Pair(test.value); break;
		case Kind::InvalidInput: InvalidInput(); break;
		case Kind::Logging: {
			auto simulation = MakeSimulation(3, test.value, 3);
			const auto initial = simulation->box->persistentClusters;
			Run(*simulation);
			Require(simulation->getStep() == test.value, "Incorrect final step");
			const int entries = (test.value + 2) / 3;
			for (int entry = 0; entry < entries; ++entry)
				for (int pc = 0; pc < initial.size(); ++pc)
					Require(simulation->traj_buffer->GetDatapoint(pc, 0, entry) == initial[pc].pqd[0].position, "Missing or incorrect trajectory frame across ring-buffer boundary");
			break;
		}
		case Kind::Thermostat: {
			// Logging intervals 0, 1 and 7 must not change the thermostat cadence or its dynamics
			constexpr int steps = 20, temperatureInterval = 4;
			std::array simulations{ MakeSimulation(2, steps, 0), MakeSimulation(2, steps, 1), MakeSimulation(2, steps, 7) };
			for (auto& simulation : simulations) {
				simulation->simParams.apply_thermostat = true;
				simulation->simParams.ref_t = 600.f;
				simulation->simParams.steps_per_temperature_measurement = temperatureInterval;
				for (auto& state : simulation->box->pclusterInterimStates) state.vels_prev[0] = Float3{0.01f, 0.02f, 0.f};
				Run(*simulation);
				Require(simulation->temperature_buffer.size() == steps / temperatureInterval,
					"Temperature measured " + std::to_string(simulation->temperature_buffer.size()) + " times, expected " + std::to_string(steps / temperatureInterval));
			}
			for (int i = 1; i < simulations.size(); ++i) {
				for (int pc = 0; pc < 2; ++pc) {
					const auto a = simulations[0]->box->pclusterInterimStates[pc].vels_prev[0];
					const auto b = simulations[i]->box->pclusterInterimStates[pc].vels_prev[0];
					Require((a - b).len() <= 1e-7f, "Changing trajectory logging changed thermostat dynamics");
				}
			}
			break;
		}
		case Kind::Unwind: {
			auto simulation = MakeSimulation(2);
			try {
				Engine engine({simulation.get()});
				CheckCuda();
				Phase("inject device trap; engine destruction must not terminate the process");
				Engine::TestLimit(EngineLimitProbe::DeviceFailure);
				const auto status = cudaDeviceSynchronize();
				Require(status != cudaSuccess, "Device fault injection did not fail");
				std::cout << "Injected CUDA error: " << cudaGetErrorString(status) << std::endl;
				throw Failure("original simulation failure");
			}
			catch (const Failure& error) {
				Require(std::string_view(error.what()) == "original simulation failure", "Cleanup replaced the original failure");
			}
			Phase("survived cleanup");
			break;
		}
		}
	}

	inline int RunChild(std::string_view name) {
#ifdef _WIN32
		SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
		_set_abort_behavior(0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
#endif
		std::cout << std::unitbuf;
		std::cerr << std::unitbuf;
		try {
			const auto cases = Cases();
			const auto found = std::ranges::find(cases, name, &Case::name);
			Require(found != cases.end(), "Unknown limit case: " + std::string(name));
			cudaDeviceProp properties{};
			LIMA_UTILS::genericErrorCheck(cudaGetDeviceProperties(&properties, 0));
			int driver = 0, runtime = 0;
			LIMA_UTILS::genericErrorCheck(cudaDriverGetVersion(&driver));
			LIMA_UTILS::genericErrorCheck(cudaRuntimeGetVersion(&runtime));
			std::cout << "CASE " << name << "\nGPU " << properties.name << "\nCUDA driver " << driver << " runtime " << runtime
				<< "\nBUILD " << __DATE__ << ' ' << __TIME__ << "\nINDEXING_CHECKS " << INDEXING_CHECKS << " FORCE_CHECKS " << FORCE_CHECKS << '\n';
			RunCase(*found);
			std::cout << "LIMIT_PASS\n";
			return 0;
		}
		catch (const std::exception& error) { std::cerr << "LIMIT_FAIL " << error.what() << '\n'; }
		catch (...) { std::cerr << "LIMIT_FAIL unknown exception\n"; }
		return 1;
	}

	inline std::filesystem::path ExecutablePath() {
#ifdef _WIN32
		std::wstring path(32768, L'\0');
		const DWORD size = GetModuleFileNameW(nullptr, path.data(), static_cast<DWORD>(path.size()));
		Require(size > 0 && size < path.size(), "Cannot locate limatest executable");
		path.resize(size);
		return path;
#else
		return std::filesystem::read_symlink("/proc/self/exe");
#endif
	}

	struct ProcessResult {
		uint64_t exitCode = 0;
		bool timedOut = false;
	};

	inline ProcessResult RunProcess(const Case& test, const std::filesystem::path& log) {
		constexpr auto timeout = std::chrono::seconds(120);
		const auto executable = ExecutablePath();
#ifdef _WIN32
		struct Handle {
			HANDLE value = INVALID_HANDLE_VALUE;
			~Handle() { if (value != INVALID_HANDLE_VALUE && value != nullptr) CloseHandle(value); }
		};
		SECURITY_ATTRIBUTES security{sizeof(security), nullptr, TRUE};
		Handle output{CreateFileW(log.c_str(), GENERIC_WRITE, FILE_SHARE_READ, &security, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr)};
		Handle input{CreateFileW(L"NUL", GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, &security, OPEN_EXISTING, 0, nullptr)};
		Require(output.value != INVALID_HANDLE_VALUE && input.value != INVALID_HANDLE_VALUE, "Cannot open child process output/input");
		STARTUPINFOW startup{sizeof(startup)};
		startup.dwFlags = STARTF_USESTDHANDLES;
		startup.hStdInput = input.value;
		startup.hStdOutput = startup.hStdError = output.value;
		PROCESS_INFORMATION process{};
		// Case names come from Cases(), contain no quotes, and are never interpreted by a shell.
		std::wstring command = L"\"" + executable.wstring() + L"\" --limit-case " + std::wstring(test.name.begin(), test.name.end());
		Require(CreateProcessW(executable.c_str(), command.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW, nullptr, nullptr, &startup, &process), "Cannot start limit test process: " + std::to_string(GetLastError()));
		Handle child{process.hProcess}, thread{process.hThread};
		const DWORD wait = WaitForSingleObject(child.value, static_cast<DWORD>(std::chrono::duration_cast<std::chrono::milliseconds>(timeout).count()));
		if (wait != WAIT_OBJECT_0) {
			TerminateProcess(child.value, 124);
			WaitForSingleObject(child.value, INFINITE);
			Require(wait == WAIT_TIMEOUT, "Failed waiting for limit test process");
			return {124, true};
		}
		DWORD code = 0;
		Require(GetExitCodeProcess(child.value, &code), "Cannot read limit test exit code");
		return {code, false};
#else
		posix_spawn_file_actions_t actions;
		Require(posix_spawn_file_actions_init(&actions) == 0, "Cannot initialize child process actions");
		const auto Destroy = [&actions](int*) { posix_spawn_file_actions_destroy(&actions); };
		int token = 0;
		std::unique_ptr<int, decltype(Destroy)> cleanup(&token, Destroy);
		Require(posix_spawn_file_actions_addopen(&actions, STDIN_FILENO, "/dev/null", O_RDONLY, 0) == 0
			&& posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, log.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600) == 0
			&& posix_spawn_file_actions_adddup2(&actions, STDOUT_FILENO, STDERR_FILENO) == 0, "Cannot redirect child output");
		std::string program = executable.string(), option = "--limit-case", name = test.name;
		char* args[]{program.data(), option.data(), name.data(), nullptr};
		pid_t pid;
		Require(posix_spawn(&pid, program.c_str(), &actions, nullptr, args, environ) == 0, "Cannot start limit test process");
		const auto deadline = std::chrono::steady_clock::now() + timeout;
		int status = 0;
		for (;;) {
			const pid_t result = waitpid(pid, &status, WNOHANG);
			if (result == pid) break;
			if (result < 0 && errno != EINTR) {
				kill(pid, SIGKILL);
				while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
				throw Failure("Cannot wait for limit test process");
			}
			if (std::chrono::steady_clock::now() >= deadline) {
				kill(pid, SIGKILL);
				while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
				return {124, true};
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
		return {static_cast<uint64_t>(WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status)), false};
#endif
	}

	inline TestUtils::TestRoutine TestIsolated(Case test, std::filesystem::path directory) {
		const auto log = directory / (test.name + ".log");
		const auto started = std::chrono::steady_clock::now();
		const auto result = RunProcess(test, log);
		const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
		std::ifstream input(log);
		const std::string output{std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
		const bool passed = !result.timedOut && result.exitCode == 0 && output.ends_with("LIMIT_PASS\n");
		std::ofstream report(directory / "results.csv", std::ios::app);
		Require(static_cast<bool>(report), "Cannot append limit test results");
		report << test.name << ',' << passed << ',' << result.timedOut << ',' << result.exitCode << ',' << seconds << ',' << log.filename().string() << '\n';
		const std::string summary = std::format("{}: {} (exit {}). Log: {}", test.name, passed ? "passed" : result.timedOut ? "timed out" : "FAILED", result.exitCode, log.string());
		if (!passed) std::cout << output.substr(output.size() > 2000 ? output.size() - 2000 : 0) << std::endl;
		co_return TestUtils::LimaUnittestResult{passed, summary, false};
	}

	inline void AddTests(TestUtils::LimaUnittestManager& manager, std::string_view filter = {}) {
		const auto directory = std::filesystem::current_path() / "build" / "limit-testing-results"
			/ std::to_string(std::chrono::system_clock::now().time_since_epoch().count());
		std::filesystem::create_directories(directory);
		std::ofstream report(directory / "results.csv");
		Require(static_cast<bool>(report), "Cannot create limit test results");
		report << "case,passed,timed_out,exit_code,seconds,log\n";
		report.close();
		std::cout << "Limit testing results: " << directory << std::endl;
		int selected = 0;
		for (const auto& test : Cases()) {
			if (test.name.find(filter) == std::string::npos) continue;
			++selected;
			manager.AddTest("LimitTesting: " + test.name, [test, directory] { return TestIsolated(test, directory); });
		}
		Require(selected > 0, "No limit tests match the filter");
	}
}
