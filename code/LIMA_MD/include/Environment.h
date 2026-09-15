#pragma once

#include "Analyzer.cuh"
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

	void Wait() const;
	bool IsReady() const;
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
	bool prepareForRun();
	void WriteTrajectoryAsUff(const fs::path& path) const;
	fs::path FixPath(const fs::path& path) const;
	const SimAnalysis::AnalyzedPackage& getAnalyzedPackage();
	void verifySimulationParameters();			// Constants before doing anything
	void verifyBox();							// Checks wheter the box will break
	
	void UpdateSimstatus(bool printToConsole, bool alwaysUpdate/*Performance hit*/);

	// Returns false if display has been closed by user
	bool handleDisplay(const BoxParams& boxparams, Display* const display, bool emVariant, bool stepwise);

	void sayHello();

	


	std::unique_ptr<Display> display = nullptr;
	std::unique_ptr<SimulationSession> simulationSession = nullptr;

	std::mutex schedulingMutex;
	std::condition_variable schedulerWakeup;
	std::deque<QueuedSimulation> pendingSimulations;
	std::optional<PreparedSimulation> preparedSimulation;
	bool preparingSimulation = false;
	bool runningSimulation = false;
	bool stopping = false;
	std::jthread coordinator;
};
