#include <chrono>
#include <filesystem>
#include <string>
#include <optional>
#include <numeric>
#include <iterator>
#include "Environment.h"
#include "MDFiles.h"
#include "BoxImageBuilder.h"
#include "Display.h"
#include "BoxBuilder.cuh"
#include "Engine.cuh"
#include "BatchCompatibility.h"
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
Environment::SimulationSession::SimulationSession(SimulationSession&&) noexcept = default;

Environment::SimulationSession& Environment::LiveEditSession() {
	if (!liveEditSession)
		throw std::runtime_error("Environment has no live-edit session");
	return *liveEditSession;
}

const Environment::SimulationSession& Environment::LiveEditSession() const {
	if (!liveEditSession)
		throw std::runtime_error("Environment has no live-edit session");
	return *liveEditSession;
}

void Environment::SetLiveEditSimulation(
	std::unique_ptr<Simulation> simulation, EnvMode mode, const fs::path& workDir) {
	if (!liveEditSession) {
		liveEditSession = std::make_unique<SimulationSession>(std::move(simulation), mode, workDir);
		return;
	}

	liveEditSession->simulation = std::move(simulation);
	liveEditSession->mode = mode;
	liveEditSession->workDir = workDir;
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
		++unpreparedSimulations;
	}
	schedulerWakeup.notify_one();
	return SimulationHandle{ std::move(state) };
}

bool Environment::MustRunAlone(const PreparedSimulation& next) {
	return next.job.mustRunAlone || !next.job.run /*|| next.job.mode == Full*/
		|| next.simulation->simParams.stepwise;
}

bool Environment::IsDrained() const {
	return stopping && unpreparedSimulations == 0 && !preparingSimulation
		&& preparedSimulations.empty() && !runningSimulation
		&& processedSimulations.empty() && !postprocessingSimulation;
}

bool Environment::CanPrepare() const {
	return !preparingSimulation && !pendingSimulations.empty()
		&& preparedSimulations.size() < maxPreparedSimulations;
}

bool Environment::CanPostprocess() const {
	return !postprocessingSimulation && !processedSimulations.empty();
}

bool Environment::CanStartBatch() const {
	if (runningSimulation || preparedSimulations.empty()
		|| processedSimulations.size() + maxBatchSize > maxProcessedSimulations)
		return false;
	if (MustRunAlone(preparedSimulations.front()) || GetReadyBatchSize() == maxBatchSize)
		return true;
	// All already-submitted jobs must finish preparation before a partial batch
	// starts. At capacity we must dispatch to free preparation space.
	return unpreparedSimulations == 0 || preparedSimulations.size() == maxPreparedSimulations;
}

size_t Environment::GetReadyBatchSize() const {
	if (preparedSimulations.empty() || MustRunAlone(preparedSimulations.front()))
		return preparedSimulations.empty() ? 0 : 1;
	size_t count = 1;
	for (auto it = std::next(preparedSimulations.begin()); it != preparedSimulations.end()
		&& count < maxBatchSize; ++it) {
		if (MustRunAlone(*it)) break;
		if (!EngineBatch::FindIncompatibility(*preparedSimulations.front().simulation, *it->simulation))
			++count;
	}
	return count;
}

std::vector<Environment::PreparedSimulation> Environment::TakeReadyBatch() {
	std::vector<PreparedSimulation> batch;
	batch.reserve(maxBatchSize);
	batch.push_back(std::move(preparedSimulations.front()));
	preparedSimulations.pop_front();
	if (MustRunAlone(batch.front())) return batch;
	for (auto it = preparedSimulations.begin(); it != preparedSimulations.end()
		&& batch.size() < maxBatchSize;) {
		if (MustRunAlone(*it)) break;
		if (EngineBatch::FindIncompatibility(*batch.front().simulation, *it->simulation)) {
			++it;
			continue;
		}
		batch.push_back(std::move(*it));
		it = preparedSimulations.erase(it);
	}
	return batch;
}

void Environment::MainLoop() {
	timingStarted = std::chrono::steady_clock::now();

	std::jthread preprocessThread;
	std::jthread simulationThread;
	std::jthread postprocessThread;

	while (true) {
		std::optional<QueuedSimulation> toPrepare;
		std::optional<ProcessedSimulation> toPostprocess;
		std::vector<PreparedSimulation> toRun;
		int batchId = 0;
		{
			std::unique_lock lock(schedulingMutex);
			schedulerWakeup.wait(lock, [this] {
				return IsDrained() || CanPrepare() || CanPostprocess() || CanStartBatch();
			});
			if (IsDrained()) break;

			// Reserve worker slots while holding the lock. A later dispatch decision
			// sees the newly reserved preparation and every submitted job still
			// counted by unpreparedSimulations.
			if (CanPrepare()) {
				toPrepare.emplace(std::move(pendingSimulations.front()));
				pendingSimulations.pop_front();
				preparingSimulation = true;
			}
			if (CanPostprocess()) {
				toPostprocess.emplace(std::move(processedSimulations.front()));
				processedSimulations.pop_front();
				postprocessingSimulation = true;
			}
			if (CanStartBatch()) {
				toRun = TakeReadyBatch();
				batchId = nextBatchId++;
				runningSimulation = true;
			}
		}

		if (toPrepare) {
			preprocessThread = std::jthread(
				[this, next = std::move(*toPrepare)]() mutable { Preprocess(std::move(next)); });
		}
		if (toPostprocess) {
			postprocessThread = std::jthread(
				[this, next = std::move(*toPostprocess)]() mutable { Postprocess(std::move(next)); });
		}
		if (!toRun.empty()) {
			simulationThread = std::jthread(
				[this, next = std::move(toRun), batchId]() mutable { RunPreparedSimulations(std::move(next), batchId); });
		}
	}
}

void Environment::Preprocess(QueuedSimulation next) {
	const auto started = std::chrono::steady_clock::now();
	try {
		auto simulation = BuildSimulation(next.job);
		if (next.job.configureSimulation)
			next.job.configureSimulation(*simulation);
		if (next.job.run)
			simulation->PrepareDataBuffers();
		const auto elapsed = std::chrono::steady_clock::now() - started;
		const std::lock_guard lock(schedulingMutex);
		preparedSimulations.emplace_back(PreparedSimulation{
			std::move(next.job), next.state, std::move(simulation), elapsed });
	}
	catch (...) {
		next.state->SetError(std::current_exception());
	}
	{
		const std::lock_guard lock(schedulingMutex);
		preprocessTime += std::chrono::steady_clock::now() - started;
		preparingSimulation = false;
		--unpreparedSimulations;
	}
	schedulerWakeup.notify_one();
}

void Environment::RunPreparedSimulations(std::vector<PreparedSimulation> next, int batchId) {
	const auto started = std::chrono::steady_clock::now();
	try {
		BatchSession batch;
		batch.sessions.reserve(next.size());
		for (auto& member : next)
			batch.sessions.emplace_back(std::move(member.simulation), member.job.mode, member.job.workDir);
		if (next.front().job.run)
			RunSimulation(batch);
		// RunSimulation destroys Engine before any of its nonowning simulation
		// pointers can be transferred to (and consumed by) postprocessing.
		const auto processTime = std::chrono::steady_clock::now() - started;
		const std::lock_guard lock(schedulingMutex);
		for (size_t i = 0; i < next.size(); ++i) {
			auto& session = batch.sessions[i];
			auto& member = next[i];
			SimulationResult result{
				std::move(session.simulation), std::nullopt, session.engineTime.value_or(std::chrono::duration<double>{}),
				member.preprocessingTime + processTime, std::move(session.avgStepTimes),
				SimulationExecutionInfo{ batchId, static_cast<int>(next.size()) } };
			processedSimulations.emplace_back(ProcessedSimulation{
				std::move(member.job), member.state, std::move(result) });
		}
	}
	catch (...) {
		for (auto& member : next)
			member.state->SetError(std::current_exception());
	}
	{
		const std::lock_guard lock(schedulingMutex);
		simulationTime += std::chrono::steady_clock::now() - started;
		runningSimulation = false;
	}
	schedulerWakeup.notify_one();
}

void Environment::Postprocess(ProcessedSimulation next) {
	const auto started = std::chrono::steady_clock::now();
	try {
		if (next.job.postprocess)
			next.job.postprocess(next.result);
		next.result.environmentTime += std::chrono::steady_clock::now() - started;
		next.state->SetResult(std::move(next.result));
	}
	catch (...) {
		next.state->SetError(std::current_exception());
	}
	{
		const std::lock_guard lock(schedulingMutex);
		postprocessTime += std::chrono::steady_clock::now() - started;
		postprocessingSimulation = false;
	}
	schedulerWakeup.notify_one();
}

void Environment::PrintDevPerformanceReport() {
	std::chrono::duration<double> elapsed;
	std::chrono::duration<double> preprocessing;
	std::chrono::duration<double> simulation;
	std::chrono::duration<double> postprocessing;
	{
		const std::lock_guard lock(schedulingMutex);
		elapsed = std::chrono::steady_clock::now() - timingStarted;
		preprocessing = preprocessTime;
		simulation = simulationTime;
		postprocessing = postprocessTime;
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
	std::printf("  Postprocessing [%s] %6.2f%%  %8.2f s\n",
		Bar(Percentage(postprocessing.count())).c_str(), Percentage(postprocessing.count()), postprocessing.count());
	std::printf("  GPU idle       [%s] %6.2f%%  %8.2f s\n",
		Bar(Percentage(idleSeconds)).c_str(), Percentage(idleSeconds), idleSeconds);
	std::printf("========================================================================\n");
	std::printf("  Preprocessing, GPU simulation, and postprocessing can overlap.\n");
}

std::unique_ptr<Simulation> Environment::BuildSimulation(SimulationJob& job) const {
	const fs::path simParamsPath = job.simParamsPath.is_absolute() ? job.simParamsPath : job.workDir / job.simParamsPath;
	if (!job.simParams)
		job.simParams.emplace(simParamsPath);

	if (job.initialSimulation) {
		auto simulation = std::make_unique<Simulation>(*job.simParams);
		BoxBuilder::copyBoxState(*simulation, std::move(job.initialSimulation->box), job.initialSimulation->getStep());
		simulation->boxImage = std::move(job.initialSimulation->boxImage);
		return simulation;
	}

	const fs::path groPath = job.groPath.is_absolute() ? job.groPath : job.workDir / job.groPath;
	const fs::path topPath = job.topPath.is_absolute() ? job.topPath : job.workDir / job.topPath;
	GroFile grofile = job.grofile ? std::move(*job.grofile) : GroFile{ groPath };
	TopologyFile topolfile = job.topfile ? std::move(*job.topfile) : TopologyFile{ topPath };
	if (job.preprocess)
		job.preprocess(grofile, topolfile, *job.simParams);
	const SimParams& simParams = *job.simParams;

	auto boxImage = LIMA_MOLECULEBUILD::buildMolecules(
		grofile, topolfile, V1,
		std::make_unique<LimaLogger>(LimaLogger::normal, job.mode, "moleculebuilder", job.workDir),
		IGNORE_HYDROGEN, simParams);
	auto simulation = std::make_unique<Simulation>(simParams, BoxBuilder::BuildBox(simParams, *boxImage));
	simulation->boxImage = std::shared_ptr<BoxImage>(std::move(boxImage));
	return simulation;
}


void Environment::InitializeLiveEditSimulation(
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
	SetLiveEditSimulation(std::move(simulation), mode, workDir);
	SimulationSession& session = LiveEditSession();

	if (display) {
		display->Submit(0, std::make_unique<Rendering::AtomRenderTask>(
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
	SimulationSession& session = LiveEditSession();
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

fs::path Environment::FixLiveEditPath(const fs::path& path) const {
	if (path.is_absolute())
		return path;
	if (fs::exists(LiveEditSession().workDir / path))
		return LiveEditSession().workDir / path;
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

std::chrono::duration<double> Environment::RunSimulation(BatchSession& batch) {
	std::vector<Simulation*> simPointers;
	for (auto& session : batch.sessions) {
		session.avgStepTimes.reserve((session.simulation->simParams.n_steps + 1) / STEPS_PER_UPDATE);
		simPointers.push_back(session.simulation.get());
	}
	Engine engine(simPointers);

	auto& controlSession = batch.sessions.front();
	auto& simulation = controlSession.simulation;
	const bool emVariant = simulation->simParams.em_variant;
	const bool stepwise = simulation->simParams.stepwise;
	std::unique_ptr<Display> display;
	if (controlSession.mode == Full) {
		display = std::make_unique<Display>();
		display->WaitForDisplayReady();
		display->Submit(0, std::make_unique<Rendering::AtomRenderTask>(
			simulation->box->persistentClusters, simulation->box->persistentClustersMetadata,
			simulation->box->boxparams, controlSession.simStatus, simulation->box->backboneChains
		), stepwise);
	}

	const auto started = std::chrono::steady_clock::now();
	for (size_t i = 0; i < batch.sessions.size(); ++i) {
		auto& session = batch.sessions[i];
		session.simulationTimer.emplace("Simulation");
		session.time0 = started;
		if (engine.GetRunStatus(i).simulation_finished) {
			session.engineTime = std::chrono::duration<double>{};
			session.simulationTimer->stop();
		}
	}
	while (!engine.IsFinished()) {
		if (!HandleDisplay(controlSession, engine, simulation->box->boxparams, display.get(), emVariant, stepwise))
			break;
		const auto stepStarted = std::chrono::steady_clock::now();
		engine.step();
		for (size_t i = 0; i < batch.sessions.size(); ++i) {
			auto& session = batch.sessions[i];
			if (session.engineTime) continue;
			UpdateSimstatus(session, engine, true, true, i);
			if (engine.GetRunStatus(i).simulation_finished) {
				session.engineTime = std::chrono::steady_clock::now() - started;
				session.simulationTimer->stop();
			}
		}
		// Optional rendering throttle for visual debugging.
		while (std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - stepStarted).count() < MIN_STEP_TIME) {}
	}
	const auto elapsed = std::chrono::steady_clock::now() - started;
	// Finalization copies final coordinates, integration state and remaining logs
	// independently for every member, including when the display is closed early.
	engine.terminateSimulation();
	for (auto& session : batch.sessions) {
		if (!session.engineTime) {
			session.engineTime = elapsed;
			session.simulationTimer->stop();
		}
		session.simulation->finished = true;
	}
	return elapsed;
}

void Environment::UpdateSimstatus(SimulationSession& session, Engine& engine, bool printToConsole, bool alwaysUpdate, size_t simulationId) {
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




		const int nStepsSinceLast = engine.GetRunStatus(simulationId).current_step - *simStatus.step;
		const double totalNsSimulated = nStepsSinceLast * simulation->simParams.dt; // [ns]
		const double wall_time_sec = duration_ms * 1e-3;
		const double ns_per_day = totalNsSimulated / (wall_time_sec / 86400.0);  // 86400 seconds in a day
		const double completionFraction = (double)step / (double)simulation->simParams.n_steps;
		const std::optional<std::chrono::duration<double>> expectedTimeToFinish = step > 0 && simulation->simParams.n_steps > 0 && simulationTimer.has_value()
			? std::optional<std::chrono::duration<double>> {simulationTimer->Elapsed()* (1. / completionFraction * (1.-completionFraction))}
			: std::nullopt;

		SimStatus newStatus{};
		newStatus.step = engine.GetRunStatus(simulationId).current_step;
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
		simStatus.maxForce = engine.GetRunStatus(simulationId).greatestForce;
	}
	else {
		if (!std::isnan(engine.GetRunStatus(simulationId).current_temperature))
			simStatus.temperature = engine.GetRunStatus(simulationId).current_temperature;
	}
}



bool Environment::HandleDisplay(SimulationSession& session, Engine& engine, const BoxParams& boxparams, Display* const display, bool emVariant, bool stepwise) {
	auto& simStatus = session.simStatus;
	auto& step_at_last_render = session.stepAtLastRender;
	if (session.mode != Full) {
		return true;
	}

	auto displayException = display->displayThreadException;
	if (displayException) {
		std::rethrow_exception(displayException);
	}

	int64_t stepForMostRecentData = engine.GetRunStatus().stepForMostRecentData;
	Float3* renderPositions = engine.GetRunStatus().most_recent_positions;
	const std::string info = emVariant
		? std::format("Step {:d} MaxForce {:.02f}", static_cast<int>(engine.GetRunStatus().current_step), static_cast<float>(engine.GetRunStatus().greatestForce))
		: std::format("Step {:d} Temp {:.02f}", static_cast<int>(engine.GetRunStatus().current_step), static_cast<float>(engine.GetRunStatus().current_temperature));

	if (stepForMostRecentData > step_at_last_render) {
		display->Submit(0, std::make_unique<Rendering::SimulationTaskUpdate>(
			renderPositions, nullptr, simStatus
		), stepwise);
		step_at_last_render = stepForMostRecentData;
		//engine->GetRunStatus().most_recent_positions = nullptr;
	}

	return !display->DisplaySelfTerminated();
}

void Environment::QueueLiveEditCommand(LiveEdit::Command command) {
	LiveEditSession().liveEditCommandsQueue.push_back(std::move(command));
}
