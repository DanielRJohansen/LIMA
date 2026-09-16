#include <chrono>
#include <filesystem>
#include <string>
#include <optional>
#include <numeric>
#include "Environment.h"
#include "MDFiles.h"
#include "BoxImageBuilder.h"
#include "Display.h"
#include "BoxBuilder.cuh"
#include "Engine.cuh"
#include "UpgradeableFileFormat.h"
#include "SimulationBuilder.h"
#include "MoleculeUtils.h"

namespace lfs = FileUtils;
namespace fs = std::filesystem;

struct ScheduledSimulationState {
	mutable std::mutex mutex;
	std::condition_variable completed;
	std::optional<SimulationResult> result;
	std::exception_ptr error;
	bool consumed = false;

	void SetResult(SimulationResult value) {
		{
			const std::lock_guard lock(mutex);
			result.emplace(std::move(value));
		}
		completed.notify_all();
	}

	void SetError(std::exception_ptr value) {
		{
			const std::lock_guard lock(mutex);
			error = std::move(value);
		}
		completed.notify_all();
	}

	bool Ready() const {
		const std::lock_guard lock(mutex);
		return result.has_value() || error;
	}
};

SimulationHandle::SimulationHandle(std::shared_ptr<ScheduledSimulationState> state)
	: state(std::move(state)) {}

void SimulationHandle::Wait() const {
	if (!state)
		throw std::runtime_error("Cannot wait for an empty simulation handle");
	std::unique_lock lock(state->mutex);
	state->completed.wait(lock, [&] { return state->result.has_value() || state->error; });
}

bool SimulationHandle::IsReady() const {
	return state && state->Ready();
}

SimulationResult SimulationHandle::Get() {
	Wait();
	std::lock_guard lock(state->mutex);
	if (state->consumed)
		throw std::runtime_error("Simulation result has already been consumed");
	state->consumed = true;
	if (state->error)
		std::rethrow_exception(state->error);
	return std::move(*state->result);
}

void SimulationResult::WriteCoordinatesTo(GroFile& grofile, std::optional<int64_t>) const {
	if (!simulation)
		throw std::runtime_error("Cannot write coordinates without a simulation result");

	int particlesUpdated = 0;
	for (int clusterId = 0; clusterId < simulation->box->persistentClusters.size(); clusterId++) {
		const auto& metadata = simulation->box->persistentClustersMetadata[clusterId];
		const auto& cluster = simulation->box->persistentClusters[clusterId];
		for (int particleId = 0; particleId < PersistentCluster::maxParticles; particleId++) {
			const int globalId = metadata.particleIdsGlobal[particleId];
			if (globalId >= 0) {
				grofile.atoms[globalId].position = cluster.pqd[particleId].position;
				particlesUpdated++;
			}
		}
	}
	if (AllAtom && particlesUpdated != grofile.atoms.size())
		throw std::runtime_error(std::format(
			"Only {} out of {} particles were updated", particlesUpdated, grofile.atoms.size()));
}

Trajectory SimulationResult::MakeTrajectory() const {
	if (!simulation || !simulation->boxImage)
		throw std::runtime_error("Cannot write a trajectory without a simulation result and BoxImage");
	return Trajectory{ static_cast<int>(simulation->getStep()),
		static_cast<int>(simulation->boxImage->grofile.atoms.size()),
		simulation->boxImage->grofile.box_size, simulation->simParams.dt };
}

void SimulationResult::WriteTrajectoryAsUff(const fs::path& path) const {
	if (!simulation || !simulation->boxImage)
		throw std::runtime_error("Cannot write a trajectory without a simulation result and BoxImage");
	UpgradeableFileFormat file(path);
	file.WriteSection("numAtoms", std::vector{ static_cast<int>(simulation->boxImage->grofile.atoms.size()) });
	file.WriteSection("numFrames", std::vector{
		static_cast<int>(simulation->getStep() / simulation->simParams.data_logging_interval) });
	file.WriteSection("trajectory", simulation->traj_buffer->GetBuffer());
}

// ------------------------------------------------ Display Parameters ------------------------------------------ //
const int STEPS_PER_UPDATE = 100;
constexpr float MIN_STEP_TIME = 0.f;		// [ms] Set to 0 for full speed sim
// -------------------------------------------------------------------------------------------------------------- //

Environment::SimulationSession::SimulationSession(std::unique_ptr<Simulation> simulation, EnvMode mode, const fs::path& workDir)
	: simulation(std::move(simulation))
	, mode(mode)
	, workDir(workDir)
	{}

Environment::SimulationSession::~SimulationSession() = default;

Environment::SimulationSession& Environment::Session() {
	if (!simulationSession)
		throw std::runtime_error("Environment has no simulation session");
	return *simulationSession;
}

const Environment::SimulationSession& Environment::Session() const {
	if (!simulationSession)
		throw std::runtime_error("Environment has no simulation session");
	return *simulationSession;
}

void Environment::SetSimulation(
	std::unique_ptr<Simulation> simulation, EnvMode mode, const fs::path& workDir) {
	if (!simulationSession) {
		simulationSession = std::make_unique<SimulationSession>(std::move(simulation), mode, workDir);
		return;
	}

	simulationSession->simulation = std::move(simulation);
	simulationSession->mode = mode;
	simulationSession->workDir = workDir;
}

Environment::Environment() {
	StartScheduling();
}

Environment& Environment::Get() {
	static Environment environment;
	return environment;
}

Environment::~Environment() {
	StopScheduling();
}

void Environment::StartScheduling() {
	const std::lock_guard lock(schedulingMutex);
	stopping = false;
	coordinator = std::jthread([this] { MainLoop(); });
}

void Environment::StopScheduling() {
	{
		const std::lock_guard lock(schedulingMutex);
		stopping = true;
	}

	// MainLoop may be asleep with an empty queue. Wake it so it can observe
	// stopping; accepted jobs are drained before the thread exits.
	schedulerWakeup.notify_one();
	if (coordinator.joinable())
		coordinator.join();
}

SimulationHandle Environment::Submit(SimulationJob job) {
	auto state = std::make_shared<ScheduledSimulationState>();
	{
		const std::lock_guard lock(schedulingMutex);
		if (stopping)
			throw std::runtime_error("Cannot submit a simulation while Environment is stopping");
		pendingSimulations.push_back({ std::move(job), state });
	}
	schedulerWakeup.notify_one();
	return SimulationHandle{ std::move(state) };
}

void Environment::MainLoop() {
	timingStarted = std::chrono::steady_clock::now();

	std::jthread preprocessThread;
	std::jthread simulationThread;


	auto MayPreprocessNextJob = [this]() -> bool {
		bool canStartNewSim = (!runningSimulation && !preparedSimulations.empty());
		bool canPrepareSimulation = !preparingSimulation && preparedSimulations.size() < maxPreparedSimulations && !pendingSimulations.empty();
		bool finishedStopping = stopping && pendingSimulations.empty() && !preparingSimulation && preparedSimulations.empty() && !runningSimulation;

		return canStartNewSim || canPrepareSimulation || finishedStopping;
		};

	while (true) {
		std::optional<PreparedSimulation> simulationToRun;
		std::optional<QueuedSimulation> simulationToPreprocess;
		{
			std::unique_lock lock(schedulingMutex);
			schedulerWakeup.wait(lock, MayPreprocessNextJob);

			if (stopping && pendingSimulations.empty() && !preparingSimulation
				&& preparedSimulations.empty() && !runningSimulation)
				break;

			if (!runningSimulation && !preparedSimulations.empty()) {
				simulationToRun.emplace(std::move(preparedSimulations.front()));
				preparedSimulations.pop_front();
				runningSimulation = true;
			}
			if (!preparingSimulation && preparedSimulations.size() < maxPreparedSimulations
				&& !pendingSimulations.empty()) {
				simulationToPreprocess.emplace(std::move(pendingSimulations.front()));
				pendingSimulations.pop_front();
				preparingSimulation = true;
			}
		}

		if (simulationToRun) {
			simulationThread = std::jthread(
				[this, next = std::move(*simulationToRun)]() mutable { RunPreparedSimulation(std::move(next)); });
		}
		if (simulationToPreprocess) {
			preprocessThread = std::jthread(
				[this, next = std::move(*simulationToPreprocess)]() mutable { Preprocess(std::move(next)); });
		}
	}
}

void Environment::Preprocess(QueuedSimulation next) {
	const auto started = std::chrono::steady_clock::now();
	try {
		auto simulation = BuildSimulation(next.job);
		if (next.job.configure)
			next.job.configure(*simulation);
		if (next.job.run)
			simulation->PrepareDataBuffers();
		const std::lock_guard lock(schedulingMutex);
		preparedSimulations.emplace_back(PreparedSimulation{
			std::move(next.job), std::move(next.state), std::move(simulation) });
	}
	catch (...) {
		next.state->SetError(std::current_exception());
	}
	{
		const std::lock_guard lock(schedulingMutex);
		preprocessTime += std::chrono::steady_clock::now() - started;
		preparingSimulation = false;
	}
	schedulerWakeup.notify_one();
}

void Environment::RunPreparedSimulation(PreparedSimulation next) {
	const auto started = std::chrono::steady_clock::now();
	try {
		simulationSession = std::make_unique<SimulationSession>(
			std::move(next.simulation), next.job.mode, next.job.workDir);
		const auto elapsed = next.job.run ? RunSimulation() : std::chrono::duration<double>{};
		SimulationSession& session = Session();
		std::optional<SimAnalysis::AnalyzedPackage> analysis;
		if (next.job.analyze)
			analysis.emplace(getAnalyzedPackage());

		SimulationResult result{
			std::move(session.simulation), std::move(analysis), elapsed, std::move(session.avgStepTimes) };
		next.state->SetResult(std::move(result));
	}
	catch (...) {
		next.state->SetError(std::current_exception());
	}
	{
		const std::lock_guard lock(schedulingMutex);
		simulationTime += std::chrono::steady_clock::now() - started;
		runningSimulation = false;
	}
	schedulerWakeup.notify_one();
}

void Environment::PrintDevPerformanceReport() {
	std::chrono::duration<double> elapsed;
	std::chrono::duration<double> preprocessing;
	std::chrono::duration<double> simulation;
	{
		const std::lock_guard lock(schedulingMutex);
		elapsed = std::chrono::steady_clock::now() - timingStarted;
		preprocessing = preprocessTime;
		simulation = simulationTime;
	}

	const auto Percentage = [total = elapsed.count()](double seconds) {
		return total > 0. ? seconds / total * 100. : 0.;
	};
	const auto Bar = [](double percentage) {
		constexpr int width = 30;
		const int filled = std::clamp(static_cast<int>(std::round(percentage / 100. * width)), 0, width);
		return std::string(filled, '#') + std::string(width - filled, '-');
	};
	// The GPU is the serialized resource. Preprocessing deliberately overlaps it,
	// so idle is the wall time for which the simulation worker was not running.
	const double idleSeconds = (std::max)(0., elapsed.count() - simulation.count());

	std::printf("\n");
	std::printf("========================================================================\n");
	std::printf("                    LIMA ENVIRONMENT UTILIZATION\n");
	std::printf("                    %.2f seconds observed\n", elapsed.count());
	std::printf("------------------------------------------------------------------------\n");
	std::printf("  GPU simulation [%s] %6.2f%%  %8.2f s\n",
		Bar(Percentage(simulation.count())).c_str(), Percentage(simulation.count()), simulation.count());
	std::printf("  Preprocessing  [%s] %6.2f%%  %8.2f s\n",
		Bar(Percentage(preprocessing.count())).c_str(), Percentage(preprocessing.count()), preprocessing.count());
	std::printf("  GPU idle       [%s] %6.2f%%  %8.2f s\n",
		Bar(Percentage(idleSeconds)).c_str(), Percentage(idleSeconds), idleSeconds);
	std::printf("========================================================================\n");
	std::printf("  Preprocessing and GPU simulation can overlap.\n");
}

std::unique_ptr<Simulation> Environment::BuildSimulation(SimulationJob& job) const {
	const fs::path simParamsPath = job.simParamsPath.is_absolute() ? job.simParamsPath : job.workDir / job.simParamsPath;
	if (!job.simParams)
		job.simParams.emplace(simParamsPath);
	job.configureParams(*job.simParams);
	const SimParams& simParams = *job.simParams;

	if (job.initialSimulation) {
		auto simulation = std::make_unique<Simulation>(simParams);
		BoxBuilder::copyBoxState(*simulation, std::move(job.initialSimulation->box), job.initialSimulation->getStep());
		simulation->boxImage = std::move(job.initialSimulation->boxImage);
		return simulation;
	}

	const fs::path groPath = job.groPath.is_absolute() ? job.groPath : job.workDir / job.groPath;
	const fs::path topPath = job.topPath.is_absolute() ? job.topPath : job.workDir / job.topPath;
	GroFile grofile = job.grofile ? std::move(*job.grofile) : GroFile{ groPath };
	TopologyFile topolfile = job.topfile ? std::move(*job.topfile) : TopologyFile{ topPath };
	job.configureInput(grofile, topolfile, *job.simParams);

	auto boxImage = LIMA_MOLECULEBUILD::buildMolecules(
		grofile, topolfile, V1,
		std::make_unique<LimaLogger>(LimaLogger::normal, job.mode, "moleculebuilder", job.workDir),
		IGNORE_HYDROGEN, simParams);
	auto simulation = std::make_unique<Simulation>(simParams, BoxBuilder::BuildBox(simParams, *boxImage));
	simulation->boxImage = std::shared_ptr<BoxImage>(std::move(boxImage));
	return simulation;
}


void Environment::InitializeSimulation(
	const GroFile& grofile, const TopologyFile& topolfile, const SimParams& params,
	EnvMode mode, const fs::path& workDir)
{
	auto boxImage = LIMA_MOLECULEBUILD::buildMolecules(
		grofile,
		topolfile,
		V1,
		std::make_unique<LimaLogger>(LimaLogger::normal, mode, "moleculebuilder", workDir),
		IGNORE_HYDROGEN,
		params
		);

	auto simulation = std::make_unique<Simulation>(params, BoxBuilder::BuildBox(params, *boxImage));
	simulation->boxImage = std::shared_ptr<BoxImage>(std::move(boxImage));
	SetSimulation(std::move(simulation), mode, workDir);
	SimulationSession& session = Session();

	if (display) {
		display->Render(std::make_unique<Rendering::AtomRenderTask>(
			session.simulation->box->persistentClusters, session.simulation->box->persistentClustersMetadata,
			session.simulation->box->boxparams, session.simStatus, session.simulation->box->backboneChains
		));
	}
}

std::tuple<GroFile, TopologyFile, SimParams> Environment::CreateLiveEditSimulationFiles(
	Float3 boxlen, const fs::path& workDir) {
	GroFile grofile{};
	grofile.m_path = workDir / "conf.gro";
	grofile.box_size = Float3{ boxlen };
	grofile.printToFile();

	TopologyFile topfile{};
	topfile.SetSystem("MySystem");
	topfile.path = workDir / "topol.top";
	//topfile.forcefieldInclude = TopologyFile::ForcefieldInclude("charmm27.ff/forcefield.itp");
	topfile.printToFile();
	topfile = TopologyFile{workDir / "topol.top"}; // Reload the topfile to parse the ffinclude


	SimParams simparams{};
	simparams.DumpToFile(workDir / "sim_params.txt");

	return { grofile, topfile, simparams };
}

void Environment::UpdateLiveEditCoordinates(GroFile& grofile) {
	SimulationSession& session = Session();
	if (session.engine) {
		CudaBuffer<PersistentCluster>& deviceState = session.engine->OffloadPclusterState();
		session.simulation->box->persistentClusters = GenericCopyToHost(
			deviceState.Get(), session.simulation->box->persistentClusters.size());
	}
	SimulationResult view;
	view.simulation = std::move(session.simulation);
	view.WriteCoordinatesTo(grofile);
	session.simulation = std::move(view.simulation);
}

fs::path Environment::FixPath(const fs::path& path) const {
	if (path.is_absolute())
		return path;
	if (fs::exists(Session().workDir / path))
		return Session().workDir / path;
	if (fs::exists( "./" / path))
		return "./" / path;
	return path;
}

void Environment::sayHello() {
	static bool hasSaidHello = false;
	if (hasSaidHello)
		return;
	hasSaidHello = true;

	std::ifstream file(FileUtils::GetLimaDir() / "resources/logo/logo_ascii.txt");
	if (!file) {
		throw std::runtime_error("Failed to open logo file");
	}

	std::string file_contents((std::istreambuf_iterator<char>(file)),
		std::istreambuf_iterator<char>());

	std::cout << file_contents;
}

std::chrono::duration<double> Environment::RunSimulation() {
	SimulationSession& session = Session();
	auto& simulation = session.simulation;
	auto& simStatus = session.simStatus;
	auto& time0 = session.time0;
	auto& simulationTimer = session.simulationTimer;
	auto& engineTime = session.engineTime;
	if (!simulation)
		throw std::runtime_error("Cannot run without a simulation");
	if (simulation->finished)
		throw std::runtime_error("Cannot run a simulation that has already finished");
	const bool emVariant = simulation->simParams.em_variant;
	const bool stepwise = simulation->simParams.stepwise;

	session.avgStepTimes.reserve((simulation->simParams.n_steps + 1) / STEPS_PER_UPDATE);
	Engine engine(simulation.get(), simulation->simParams.bc_select);

	std::unique_ptr<Display> display = nullptr;

	if (session.mode == Full) {
		display = std::make_unique<Display>();
		display->WaitForDisplayReady();
		display->Render(std::make_unique<Rendering::AtomRenderTask>(
			simulation->box->persistentClusters, simulation->box->persistentClustersMetadata,
			simulation->box->boxparams, simStatus, simulation->box->backboneChains
		), stepwise);
	}

	simulationTimer.emplace(TimeIt{ "Simulation" });
	time0 = std::chrono::steady_clock::now();
    auto t0 = std::chrono::steady_clock::now();
	while (true) {

		if (!handleDisplay(engine, simulation->box->boxparams, display.get(), emVariant, stepwise)) {
			break;
		}

		auto stepStartTime = std::chrono::steady_clock::now();
		
		engine.step();

		UpdateSimstatus(engine, true, true);
		
		if (engine.runstatus.simulation_finished) {
			break;
		}

		// Deadspin to slow down rendering for visual debugging :)
		while ((double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - stepStartTime).count() < MIN_STEP_TIME) {}
	}
    auto t1 = std::chrono::steady_clock::now();
	simulationTimer->stop();

	// Transfers the remaining traj data and more
	engine.terminateSimulation();
	CudaBuffer<PersistentCluster>& deviceState = engine.OffloadPclusterState();
	simulation->box->persistentClusters = GenericCopyToHost(
		deviceState.Get(), simulation->box->persistentClusters.size());

	simulation->finished = true;

	engineTime = t1 - t0;
    return t1-t0;
}







void Environment::WriteTrajectoryAsUff(const fs::path& path) const {
	const auto& simulation = Session().simulation;
	const auto& boximage = simulation->boxImage;
	if (!boximage)
		throw std::runtime_error("Cannot write a trajectory without a BoxImage");
	const int nSteps = simulation->getStep();
	const int nAtoms = boximage->grofile.atoms.size();

	UpgradeableFileFormat file(path);

	file.WriteSection("numAtoms", std::vector{ nAtoms });
	file.WriteSection("numFrames", std::vector{nSteps / simulation->simParams.data_logging_interval});
	file.WriteSection("trajectory", simulation->traj_buffer->GetBuffer());
}

void Environment::UpdateSimstatus(Engine& engine, bool printToConsole, bool alwaysUpdate) {
	SimulationSession& session = Session();
	auto& simulation = session.simulation;
	auto& simStatus = session.simStatus;
	auto& time0 = session.time0;
	auto& avgStepTimes = session.avgStepTimes;
	auto& simulationTimer = session.simulationTimer;
	auto& forceWriteSimstatusToDisplay = session.forceWriteSimstatusToDisplay;
	if (!simulation) {
		return;
	}

	const int64_t step = simulation->getStep();
	if ((step % STEPS_PER_UPDATE == STEPS_PER_UPDATE-1) || forceWriteSimstatusToDisplay) {		
		forceWriteSimstatusToDisplay = false;
		auto duration = std::chrono::steady_clock::now() - time0;		
		const double duration_ms = std::chrono::duration_cast<std::chrono::microseconds>(duration).count() * 1e-3;
		const double avgSteptime = duration_ms / (double) STEPS_PER_UPDATE;

		if (printToConsole && session.mode == Full) {
			//// First clear the current line
			//printf("\r\033[K");
			// Move cursor to the beginning of the line and clear it
			printf("\033[1000D\033[K");

			printf("Step #%06llu", step);
			printf("\tAvg. time: %.2fms", avgSteptime);
		}

		time0 = std::chrono::steady_clock::now();
		avgStepTimes.emplace_back(avgSteptime);




		const int nStepsSinceLast = engine.runstatus.current_step - *simStatus.step;
		const double totalNsSimulated = nStepsSinceLast * simulation->simParams.dt; // [ns]
		const double wall_time_sec = duration_ms * 1e-3;
		const double ns_per_day = totalNsSimulated / (wall_time_sec / 86400.0);  // 86400 seconds in a day
		const double completionFraction = (double)step / (double)simulation->simParams.n_steps;
		const std::optional<std::chrono::duration<double>> expectedTimeToFinish = simulation->simParams.n_steps > 0 && simulationTimer.has_value()
			? std::optional<std::chrono::duration<double>> {simulationTimer->Elapsed()* (1. / completionFraction * (1.-completionFraction))}
			: std::nullopt;

		SimStatus newStatus{};
		newStatus.step = engine.runstatus.current_step;
		newStatus.avgStepTime = avgStepTimes.empty() ? 0.f : avgStepTimes.back();
		newStatus.expectedTimeToFinish = expectedTimeToFinish;
		if (simulation->simParams.em_variant) {
		}
		else {
			newStatus.simulationPerformance = ns_per_day;
		}

		simStatus = newStatus;
	}

	// "Free" updates
	if (simulation->simParams.em_variant) {
		simStatus.maxForce = engine.runstatus.greatestForce;
	}
	else {
		if (!std::isnan(engine.runstatus.current_temperature))
			simStatus.temperature = engine.runstatus.current_temperature;
	}
}



bool Environment::handleDisplay(Engine& engine, const BoxParams& boxparams, Display* const display, bool emVariant, bool stepwise) {
	SimulationSession& session = Session();
	auto& simStatus = session.simStatus;
	auto& step_at_last_render = session.stepAtLastRender;
	if (session.mode != Full) {
		return true;
	}

	auto displayException = display->displayThreadException;
	if (displayException) {
		std::rethrow_exception(displayException);
	}

	int64_t stepForMostRecentData = engine.runstatus.stepForMostRecentData;
	Float3* renderPositions = engine.runstatus.most_recent_positions;
	const std::string info = emVariant
		? std::format("Step {:d} MaxForce {:.02f}", static_cast<int>(engine.runstatus.current_step), static_cast<float>(engine.runstatus.greatestForce))
		: std::format("Step {:d} Temp {:.02f}", static_cast<int>(engine.runstatus.current_step), static_cast<float>(engine.runstatus.current_temperature));

	if (stepForMostRecentData > step_at_last_render) {
		display->Render(std::make_unique<Rendering::SimulationTaskUpdate>(
			renderPositions, nullptr, simStatus
		), stepwise);
		step_at_last_render = stepForMostRecentData;
		//engine->runstatus.most_recent_positions = nullptr;
	}

	return !display->DisplaySelfTerminated();
}

const SimAnalysis::AnalyzedPackage& Environment::getAnalyzedPackage()
{
	SimulationSession& session = Session();
	auto& simulation = session.simulation;
	auto& postsim_anal_package = session.analyzedPackage;
	if (simulation == nullptr)
		throw std::runtime_error("Env has no simulation");
	// TODO: make some check here that the simulation has finished
	if (!postsim_anal_package.has_value())
		postsim_anal_package = SimAnalysis::analyzeEnergy(simulation.get());
	return postsim_anal_package.value();
}

void Environment::QueueLiveEditCommand(LiveEdit::Command command) {
	Session().liveEditCommandsQueue.push_back(std::move(command));
}
