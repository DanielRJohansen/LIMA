#pragma once

#include "Analyzer.h"
#include "Bodies.cuh"
#include "TimeIt.h"
#include "MDFiles.h"
#include "Trajectory.h"
#include "LiveEditCommands.h"

#include <memory>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <thread>

class Display;
struct BoxImage;
class Engine;
struct LiveEditData;
struct ScheduledSimulationState;

namespace fs = std::filesystem;

struct SimulationJob {
	SimulationJob() = default;
	SimulationJob(fs::path workDir, GroFile grofile, TopologyFile topfile, SimParams simParams, EnvMode mode)
		: workDir(std::move(workDir)), simParams(std::move(simParams)), grofile(std::move(grofile)),
		topfile(std::move(topfile)), mode(mode) {}

	fs::path workDir;
	fs::path groPath{ "molecule/conf.gro" };
	fs::path topPath{ "molecule/topol.top" };
	fs::path simParamsPath{ "sim_params.txt" };
	std::optional<SimParams> simParams;
	std::optional<GroFile> grofile;
	std::optional<TopologyFile> topfile;
	std::unique_ptr<Simulation> initialSimulation;
	EnvMode mode = EnvMode::Headless;
	std::function<void(SimParams&)> configureParams = [](SimParams&) {};
	std::function<void(GroFile&, TopologyFile&, SimParams&)> configureInput = [](GroFile&, TopologyFile&, SimParams&) {};
	std::function<void(Simulation&)> configure = [](Simulation&) {};
	bool run = true;
	bool analyze = false;
};

struct SimulationResult {
	std::unique_ptr<Simulation> simulation;
	std::optional<SimAnalysis::AnalyzedPackage> analysis;
	std::chrono::duration<double> engineTime{};
	std::vector<float> averageStepTimes;

	void WriteCoordinatesTo(GroFile& grofile, std::optional<int64_t> step = std::nullopt) const;
	Trajectory MakeTrajectory() const;
	void WriteTrajectoryAsUff(const fs::path& path) const;
};

class SimulationHandle {
public:
	SimulationHandle() = default;
	SimulationHandle(const SimulationHandle&) = default;
	SimulationHandle& operator=(const SimulationHandle&) = default;
	SimulationHandle(SimulationHandle&&) noexcept = default;
	SimulationHandle& operator=(SimulationHandle&&) noexcept = default;

	// Blocks the calling thread until the job finishes; does not consume the result.
	void Wait() const;
	// Checks completion without blocking or consuming the result.
	bool IsReady() const;
	// Blocks until completion, rethrows job errors, and transfers the result; may only be called once across all handle copies.
	SimulationResult Get();

private:
	explicit SimulationHandle(std::shared_ptr<ScheduledSimulationState> state);
	std::shared_ptr<ScheduledSimulationState> state;
	friend class Environment;
};


class Environment
{
public:
	static Environment& Get();

	// Queues a lightweight simulation description. All expensive construction and
	// GPU work is performed by Environment's bounded worker pipeline.
	[[nodiscard]] SimulationHandle Submit(SimulationJob job);

	////////////////// LIVE EDIT //////////////////

	/// <summary>
	/// A mode where the user can continously give inputs to the program
	/// </summary>
	std::tuple<GroFile, TopologyFile, SimParams> BeginLiveEdit(
		const fs::path& workDir, EnvMode mode, Float3 boxlen);
	void LiveEdit(GroFile& grofile, TopologyFile& topfile);
private:
	void InsertMolecule(LiveEditData*, GroFile& grofile, TopologyFile& topfile, LiveEdit::InsertMolecule& insertionCmd, SimParams simparams);
	void BuildMembrane(LiveEditData*, const LiveEdit::BuildMembrane& cmd, GroFile& grofile, TopologyFile& topfile);
	void HandleMoveMoleculeCommand(LiveEditData*, const LiveEdit::MoveMolecule& newMoveCommand);
	void UpdateForcemask(LiveEditData*, const LiveEdit::AddForcemaskToSelection&);
	void UpdateSelection(LiveEditData*, const LiveEdit::AtomSelected&);
	void UpdateSelection(LiveEditData*, const LiveEdit::SelectAtomsBasedOnQualifier&);
	void UpdateElasticPosition(LiveEditData*, const LiveEdit::ElasticPosition&);
	void EM(LiveEditData*);
	////////////////// ////////////////// ////////////////// 
public:



	void QueueLiveEditCommand(LiveEdit::Command command);

	// Development diagnostics for scheduler/worker utilization.
	void PrintDevPerformanceReport();

private:
	Environment();
	~Environment();

	struct QueuedSimulation {
		SimulationJob job;
		std::shared_ptr<ScheduledSimulationState> state;
	};

	struct PreparedSimulation {
		SimulationJob job;
		std::shared_ptr<ScheduledSimulationState> state;
		std::unique_ptr<Simulation> simulation;
	};

	void StartScheduling();
	void StopScheduling();
	void MainLoop();


	// Functions that are only run by their own dedicated worker thread
	void Preprocess(QueuedSimulation next);					// preprocessor thread	
	void RunPreparedSimulation(PreparedSimulation next);	// simulation thread
	//


	std::unique_ptr<Simulation> BuildSimulation(SimulationJob& job) const;
	void InitializeSimulation(
		const GroFile&, const TopologyFile&, const SimParams&, EnvMode mode, const fs::path& workDir);
	std::tuple<GroFile, TopologyFile, SimParams> CreateLiveEditSimulationFiles(
		Float3 boxlen, const fs::path& workDir);
	void UpdateLiveEditCoordinates(GroFile& grofile);
	std::chrono::duration<double> RunSimulation();

	struct SimulationSession {
		SimulationSession(std::unique_ptr<Simulation> simulation, EnvMode mode, const fs::path& workDir);
		~SimulationSession();

		std::unique_ptr<Simulation> simulation;
		std::unique_ptr<Engine> engine = nullptr;
		std::chrono::steady_clock::time_point time0;
		std::optional<TimeIt> simulationTimer;
		std::vector<float> avgStepTimes;
		std::optional<std::chrono::duration<double>> engineTime;
		std::deque<LiveEdit::Command> liveEditCommandsQueue;
		SimStatus simStatus{};
		bool forceWriteSimstatusToDisplay = false;
		int64_t stepAtLastRender = INT64_MIN;
		std::optional<SimAnalysis::AnalyzedPackage> analyzedPackage;
		EnvMode mode;
		fs::path workDir;
	};

	SimulationSession& Session();
	const SimulationSession& Session() const;
	void SetSimulation(
		std::unique_ptr<Simulation> simulation, EnvMode mode, const fs::path& workDir);
	void WriteTrajectoryAsUff(const fs::path& path) const;
	fs::path FixPath(const fs::path& path) const;
	const SimAnalysis::AnalyzedPackage& getAnalyzedPackage();
	
	void UpdateSimstatus(Engine& engine, bool printToConsole, bool alwaysUpdate/*Performance hit*/);

	// Returns false if display has been closed by user
	bool handleDisplay(Engine& engine, const BoxParams& boxparams, Display* const display, bool emVariant, bool stepwise);

	void sayHello();

	


	std::unique_ptr<Display> display = nullptr;
	std::unique_ptr<SimulationSession> simulationSession = nullptr;

	std::mutex schedulingMutex;
	std::condition_variable schedulerWakeup;
	std::deque<QueuedSimulation> pendingSimulations;
	static constexpr size_t maxPreparedSimulations = 3;
	std::deque<PreparedSimulation> preparedSimulations;
	bool preparingSimulation = false;
	bool runningSimulation = false;
	bool stopping = false;
	std::jthread coordinator;

	std::chrono::steady_clock::time_point timingStarted;
	std::chrono::duration<double> preprocessTime{};
	std::chrono::duration<double> simulationTime{};
};
