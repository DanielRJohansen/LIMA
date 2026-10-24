#include "Engine.cuh"
#include <cstring>
#include <thread>
#include <condition_variable>
#include <functional>
#include <mutex>

#include "BoundaryCondition.cuh"
#include "EngineBodies.cuh"
#include "EngineKernels.cuh"
#include "EnergyMinimization.cuh"
#include "LimaPositionSystem.cuh"
#include "PME.cuh"
#include "SimulationData.h"
#include "Thermostat.cuh"
#include "Statistics.h"
#include "Utilities.h"
#include "EngineHostside.h"
#include "SuperclusterStagingControl.cuh"
#include "TaskBuilderControl.cuh"
#include "BatchData.cuh"
#include "LimitTesting.cuh"
#include "RenderDataPipe.h"

#include <random>
#include <numeric>

EngineBatchData::~EngineBatchData() {
	pmeController.reset();
	boxState.FreeMembers();
	if (forceEnergyInterims) forceEnergyInterims->Free();
	if (superClustersControl) superClustersControl->Free();
	if (pclusterTransfermodule) pclusterTransfermodule->Free();
	if (superclusterStagingControl) superclusterStagingControl->Free();
}

Engine::Engine(const std::vector<Simulation*>& simulations, EngineRunMode mode,
	const std::vector<RenderDataPipe*>& inputRenderDataPipes)
	: mode(mode), batch(std::make_unique<EngineBatchData>()), renderDataPipes(inputRenderDataPipes)
{
	if (mode == EngineRunMode::Interactive && simulations.size() != 1)
		throw std::invalid_argument("Interactive engines require one simulation");
	if (mode == EngineRunMode::Interactive && simulations.front() && simulations.front()->simParams.apply_barostat)
		throw std::invalid_argument("Pressure coupling is not supported during live editing");
	if (!renderDataPipes.empty() && renderDataPipes.size() != simulations.size())
		throw std::invalid_argument("Render data pipe count does not match simulation count");
	EngineBatch::Pack(*batch, simulations);
	for (auto& simulation : batch->simulations)
		if (!simulation.simulation->temperature_buffer.empty())
			simulation.runstatus.current_temperature = simulation.simulation->temperature_buffer.back();
	if (renderDataPipes.empty())
		renderDataPipes.resize(simulations.size(), nullptr);
	if (mode == EngineRunMode::Interactive) {
		batch->simulations.front().device.active = true;
		batch->simulations.front().runstatus.simulation_finished = false;
	}
	verifyEngine();
#if defined(LIMA_HIP)
	// The kernels assume 32-wide warps. AMD's wave64 GPUs (MI-series, older GCN) would compute wrong results silently
	static const int warpSize = [] {
		int device = 0;
		cudaDeviceProp properties{};
		LIMA_UTILS::genericErrorCheck(cudaGetDevice(&device));
		LIMA_UTILS::genericErrorCheck(cudaGetDeviceProperties(&properties, device));
		return properties.warpSize;
	}();
	if (warpSize != 32)
		throw std::runtime_error("LIMA requires a GPU with 32-wide warps, this one has " + std::to_string(warpSize));
#endif
	try {
		if (UsesEnergyMinimization()) UploadEnergyMinimizationPreconditioner();
		for (auto& stream : cudaStreams) cudaStreamCreate(&stream);
		cudaStreamCreate(&pmeStream);
		// The SMs' L1/shared split follows the resident kernels. With the default, the NB and bonded blocks of a small system leave
		// too little shared memory for PME's kernels, which run concurrently on the critical path, so their blocks waited for the
		// SMs to drain (T4: 2.7 us before SpreadChargesKernel, a 4 us shorter step without). No cost on STMV's kernels (ncu)
		cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
		cudaStreamCreateWithFlags(&logCopyStream, cudaStreamNonBlocking); // Not synchronized with the legacy default stream, which would wait for the copy
		for (auto& event : streamJoinEvents) cudaEventCreateWithFlags(&event, cudaEventDisableTiming);
		cudaEventCreateWithFlags(&stepStartEvent, cudaEventDisableTiming);
		LIMA_UTILS::genericErrorCheck(cudaMalloc(&stepCounterDevice, sizeof(int64_t)));
		LIMA_UTILS::genericErrorCheck(cudaMemset(stepCounterDevice, 0xFF, sizeof(int64_t))); // -1, before step 0
		if (batch->params.data_logging_interval != 0)
			for (size_t i = 0; i < batch->simulations.size(); ++i)
				logDrains.push_back(std::make_unique<LogDrain>());
		for (size_t simulationId = 0; simulationId < renderDataPipes.size(); ++simulationId) {
			if (!renderDataPipes[simulationId]) continue;
			const auto count = batch->simulations[simulationId].device.pclusters.count * PersistentCluster::maxParticles;
			renderDataPipes[simulationId]->Initialize(count);
		}
		batch->dataBuffersDevice = std::make_unique<DatabuffersDeviceController>(batch->nPclusters, batch->params.data_logging_interval);
		batch->superClustersControl = std::make_unique<SuperClustersControl>(batch->nGridnodes, batch->nPclusters);
		batch->pclusterTransfermodule = std::make_unique<PClusterTransfermodule>(PClusterTransfermodule::Create(batch->nGridnodes));
		batch->thermostat = std::make_unique<Thermostat>(batch->nPclusters);
		{
			std::vector<VRescale> rngs;
			for (const auto& sim : batch->simulations)
				rngs.push_back(sim.thermostat);
			deviceThermostat = std::make_unique<DeviceThermostat>(rngs, std::max<int>(1, static_cast<int>(batch->params.steps_per_temperature_measurement.value)), batch->step);
			UploadDeviceThermostatSimulations();
		}
		BootstrapClustering(cudaStreams[0]);
		MakeSuperClusterTasksGPU(cudaStreams[0]);
		InitializePME();
		UploadNbChunks();
		UploadIntegrationSimulationData();
		bool retired = false;
		for (auto& sim : batch->simulations) {
			BootstrapTrajbufferWithCoords(sim);
			if (!sim.device.active) {
				FinalizeSimulation(sim);
				retired = true;
			}
		}
		if (retired && !IsFinished()) RebuildActiveBatch();
	}
	catch (...) {
		Synchronize();
		StopRenderDataPipes();
		batch.reset();
		for (auto stream : cudaStreams) if (stream) cudaStreamDestroy(stream);
		if (pmeStream) cudaStreamDestroy(pmeStream);
		if (logCopyStream) cudaStreamDestroy(logCopyStream);
		for (auto event : streamJoinEvents) if (event) cudaEventDestroy(event);
		if (stepStartEvent) cudaEventDestroy(stepStartEvent);
		cudaFree(stepCounterDevice);
		for (auto& virialGraphs : stepGraphs) for (auto& graphs : virialGraphs) for (auto graph : graphs) if (graph) cudaGraphExecDestroy(graph);
		throw;
	}
}

Engine::~Engine() {
	Synchronize();
	StopRenderDataPipes();
	batch.reset(); // Controllers reference streams, so destroy them first.
	for (auto stream : cudaStreams) cudaStreamDestroy(stream);
	cudaStreamDestroy(pmeStream);
	logDrains.clear(); // Joins the unpacking threads
	if (logCopyStream) cudaStreamDestroy(logCopyStream);
	for (auto event : streamJoinEvents) if (event) cudaEventDestroy(event);
	if (stepStartEvent) cudaEventDestroy(stepStartEvent);
	cudaFree(stepCounterDevice);
	for (auto& virialGraphs : stepGraphs) for (auto& graphs : virialGraphs) for (auto graph : graphs) if (graph) cudaGraphExecDestroy(graph);
	LIMA_UTILS::genericErrorCheckNoSync("Error during Engine destruction");
}

void Engine::Synchronize() {
	cudaStreamSynchronize(pmeStream);
	for (auto stream : cudaStreams) cudaStreamSynchronize(stream);
}

// Makes work subsequently queued on cudaStreams[0] wait for all work queued so far on the other streams, on the GPU, so the host doesn't block
void Engine::JoinStreamsIntoMainStream() {
	cudaEventRecord(streamJoinEvents[0], pmeStream);
	cudaStreamWaitEvent(cudaStreams[0], streamJoinEvents[0]);
	for (size_t i = 1; i < cudaStreams.size(); i++) {
		cudaEventRecord(streamJoinEvents[i], cudaStreams[i]);
		cudaStreamWaitEvent(cudaStreams[0], streamJoinEvents[i]);
	}
}

// Makes work subsequently queued on the other streams wait for all work queued so far on cudaStreams[0], on the GPU
void Engine::ForkStreamsFromMainStream() {
	cudaEventRecord(stepStartEvent, cudaStreams[0]);
	cudaStreamWaitEvent(pmeStream, stepStartEvent);
	for (size_t i = 1; i < cudaStreams.size(); i++)
		cudaStreamWaitEvent(cudaStreams[i], stepStartEvent);
}

// Graphs remove the WDDM submit batching that idles the GPU on small systems' short steps (T4 -15-18%, STMV
// neutral). Encoding a launch order (NB after PME's spreading, bonds after NB) was slower on both. On HIP steps are launched
// directly; it only needs to be correct
#if defined(LIMA_HIP)
void Engine::BeginStepCapture() {}
void Engine::LaunchCapturedStep(bool, bool, bool) {}
void Engine::LaunchStepGraph(cudaGraphExec_t, bool) {}
#else
// ThreadLocal, so the log drain thread's CUDA calls are not affected by the capture
void Engine::BeginStepCapture() {
	LIMA_UTILS::genericErrorCheck(cudaStreamBeginCapture(cudaStreams[0], cudaStreamCaptureModeThreadLocal));
}
void Engine::LaunchCapturedStep(bool logData, bool computeVirial, bool measureTemperature) {
	cudaGraph_t graph = nullptr;
	LIMA_UTILS::genericErrorCheck(cudaStreamEndCapture(cudaStreams[0], &graph));
	cudaGraphExec_t& exec = stepGraphs[logData][computeVirial][measureTemperature];
	if (exec) {
		cudaGraphExecUpdateResultInfo updateInfo{};
		if (cudaGraphExecUpdate(exec, graph, &updateInfo) != cudaSuccess) {
			cudaGetLastError(); // A changed topology is expected after a task rebuild, the graph is then instantiated anew
			cudaGraphExecDestroy(exec);
			exec = nullptr;
		}
	}
	if (!exec)
		LIMA_UTILS::genericErrorCheck(cudaGraphInstantiate(&exec, graph, 0));
	LaunchStepGraph(exec, true);
	cudaGraphDestroy(graph);
	if (!logData && !computeVirial && mode != EngineRunMode::Interactive)
		replayableStepGraph[measureTemperature] = true;
}
void Engine::LaunchStepGraph(cudaGraphExec_t exec, bool captured) {
	LIMA_UTILS::genericErrorCheck(cudaGraphLaunch(exec, cudaStreams[0]));
	// When the GPU has drained its queue, after a host wait or while the host captured a step, WDDM holds launches back until a later
	// submit, up to ~750 us. Not every launch: replays keep the queue filled, and then the extra submits cost more than they save
	if (captured || flushNextStepLaunch) {
		cudaStreamQuery(cudaStreams[0]);
		flushNextStepLaunch = false;
	}
}
#endif

const RunStatus& Engine::GetRunStatus(size_t simulationId) const {
	return batch->simulations.at(simulationId).runstatus;
}

void Engine::StopSimulation(size_t simulationId) {
	auto& sim = batch->simulations.at(simulationId);
	if (sim.finalized) return;
	FinalizeSimulation(sim);
	rebuildPending = true;
}

bool Engine::IsFinished() const {
	return std::none_of(batch->simulations.begin(), batch->simulations.end(), [](const auto& sim) { return sim.device.active; });
}

// From the simulation's own size, so a simulation computes the same forces batched and alone: its supercluster count estimated
// from its particles, against a large GPU's resident warps of NbNonlocalKernel. A constant, so the forces don't depend on the GPU.
// More than 8 gained nothing on T4 (2.6k particles, -18% at 8), and STMV-sized simulations get 1, as more cost them (+7% at 8)
void Engine::UploadNbChunks() {
	constexpr int targetWarps = 4096, maxChunks = 8;
	std::vector<int> chunks;
	for (const auto& sim : batch->simulations) {
		const int superclusters = std::max(1, sim.device.particles.count / SuperCluster::maxParticles);
		chunks.push_back(std::clamp(targetWarps / superclusters, 1, maxChunks));
	}
	batch->nbMaxChunks = *std::max_element(chunks.begin(), chunks.end());
	batch->nbSimulationChunks.Expand(chunks.size());
	cudaMemcpy(batch->nbSimulationChunks.Get(), chunks.data(), sizeof(int) * chunks.size(), cudaMemcpyHostToDevice);
	InvalidateStepGraph();
}

void Engine::InitializePME() {
	if (ENABLE_ES_LR && batch->params.enable_electrostatics) {
		batch->pmeController = std::make_unique<PME::Controller>(batch->simulations, batch->params.cutoff_nm, batch->params.fourier_spacing, pmeStream);
	}
}

namespace {
	template<typename T>
	void CopyParticleBuffer(const std::optional<CudaBuffer<T>>& source, std::optional<CudaBuffer<T>>& destination,
		const EngineBatchData& oldBatch, const EngineBatchData& newBatch) {
		if (!source) return;
		destination.emplace();
		destination->Expand(newBatch.nParticles);
		for (size_t id = 0; id < newBatch.simulations.size(); ++id) {
			const auto& oldRange = oldBatch.simulations[id].device.particles;
			const auto& newRange = newBatch.simulations[id].device.particles;
			if (!newBatch.simulations[id].device.active) continue;
			cudaMemcpy(destination->Get() + newRange.offset, source->Get() + oldRange.offset,
				sizeof(T) * newRange.count, cudaMemcpyDeviceToDevice);
		}
	}
}

bool Engine::UsesEnergyMinimization() const {
	return mode == EngineRunMode::Interactive || batch->params.em_variant;
}

void Engine::PrepareEnergyMinimization(Simulation& simulation) {
	if (!simulation.box) throw std::invalid_argument("Cannot prepare energy minimization without a box");
	// Engines use the default EM::Config
	auto preconditioner = EM::MakePreconditioner(*simulation.box, EM::Config{}.nonbondedStiffness);
	simulation.emInverseStiffness = std::move(preconditioner.inverseStiffness);
	simulation.emWholeMolecule = std::move(preconditioner.wholeMolecule);
}

// Each simulation's preconditioner is made once, from its coordinates when this engine first sees it, and reused when
// the batch is rebuilt, so neither the cost nor the result depends on when other members retire
void Engine::UploadEnergyMinimizationPreconditioner() {
	std::vector<const Box*> boxes;
	std::vector<EM::Preconditioner*> preconditioners;
	for (auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		if (sim.emPreconditioner.inverseStiffness.empty()) {
			sim.emPreconditioner.inverseStiffness = std::move(sim.simulation->emInverseStiffness);
			sim.emPreconditioner.wholeMolecule = std::move(sim.simulation->emWholeMolecule);
			sim.simulation->emInverseStiffness.clear();
			sim.simulation->emWholeMolecule.clear();
		}
		boxes.push_back(sim.simulation->box.get());
		preconditioners.push_back(&sim.emPreconditioner);
	}
	EM::MakeMissingPreconditioners(boxes, preconditioners, batch->emConfig.nonbondedStiffness);

	std::vector<float> inverseStiffness(size_t(batch->nPclusters) * PersistentCluster::maxParticles, 0.f);
	std::vector<uint8_t> wholeMolecule(batch->nPclusters, 0);
	for (const auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		const auto range = sim.device.pclusters;
		std::ranges::copy(sim.emPreconditioner.inverseStiffness, inverseStiffness.begin() + size_t(range.offset) * PersistentCluster::maxParticles);
		std::ranges::copy(sim.emPreconditioner.wholeMolecule, wholeMolecule.begin() + range.offset);
	}
	batch->emInverseStiffness.SetData(inverseStiffness);
	batch->emWholeMolecule.SetData(wholeMolecule);
}

void Engine::RebuildActiveBatch() {
	StoreIntegrationStates(cudaStreams[0]);
	Synchronize();
	std::vector<Simulation*> simulations;
	std::vector<bool> active;
	for (auto& sim : batch->simulations) {
		simulations.push_back(sim.simulation);
		active.push_back(sim.device.active);
		if (!sim.device.active) continue;
		OffloadLoggingData(sim);
		const auto range = sim.device.pclusters;
		auto& box = *sim.simulation->box;
		cudaMemcpy(box.persistentClusters.data(), batch->pClusterDevice.Get() + range.offset,
			sizeof(PersistentCluster) * range.count, cudaMemcpyDeviceToHost);
		cudaMemcpy(box.pclusterInterimStates.data(), batch->boxState.pclusterInterimStates + range.offset,
			sizeof(PersistentclusterInterimState) * range.count, cudaMemcpyDeviceToHost);
	}

	JoinLogDrains(); // The snapshots must be taken before the logging ring is freed
	auto rebuilt = std::make_unique<EngineBatchData>();
	EngineBatch::Pack(*rebuilt, simulations, &active);
	rebuilt->step = batch->step;
	rebuilt->listBuffer = batch->listBuffer;	// Pack sets it only at the start
	for (size_t id = 0; id < simulations.size(); ++id) {
		const auto& oldSim = batch->simulations[id];
		auto& newSim = rebuilt->simulations[id];
		newSim.runstatus = oldSim.runstatus;
		newSim.step = oldSim.step;
		newSim.nLogEntriesTransferred = oldSim.nLogEntriesTransferred;
		newSim.finalized = oldSim.finalized;
		newSim.finalForcesMagnitudeSquared = oldSim.finalForcesMagnitudeSquared;
		newSim.device.thermostatScalar = oldSim.device.thermostatScalar;
		newSim.thermostat = oldSim.thermostat;
		if (!active[id]) continue;
		newSim.emPreconditioner = std::move(batch->simulations[id].emPreconditioner);
		cudaMemcpy(rebuilt->emParticles.Get() + newSim.device.pclusters.offset * PersistentCluster::maxParticles,
			batch->emParticles.Get() + oldSim.device.pclusters.offset * PersistentCluster::maxParticles,
			sizeof(EM::ParticleState) * newSim.device.pclusters.count * PersistentCluster::maxParticles, cudaMemcpyDeviceToDevice);
		cudaMemcpy(rebuilt->emStates.Get() + id, batch->emStates.Get() + id, sizeof(EM::SimState), cudaMemcpyDeviceToDevice);
	}
	CopyParticleBuffer(batch->fixedParticleMovementBuffer, rebuilt->fixedParticleMovementBuffer, *batch, *rebuilt);
	CopyParticleBuffer(batch->fixedParticleRotationBuffer, rebuilt->fixedParticleRotationBuffer, *batch, *rebuilt);
	CopyParticleBuffer(batch->forceMaskBuffer, rebuilt->forceMaskBuffer, *batch, *rebuilt);
	CopyParticleBuffer(batch->elasticPositionsBuffer, rebuilt->elasticPositionsBuffer, *batch, *rebuilt);

	auto oldBatch = std::move(batch);
	batch = std::move(rebuilt);
	rebuildPending = false;
	if (UsesEnergyMinimization()) UploadEnergyMinimizationPreconditioner();
	batch->dataBuffersDevice = std::make_unique<DatabuffersDeviceController>(batch->nPclusters, batch->params.data_logging_interval);
	batch->superClustersControl = std::make_unique<SuperClustersControl>(batch->nGridnodes, batch->nPclusters);
	batch->pclusterTransfermodule = std::make_unique<PClusterTransfermodule>(PClusterTransfermodule::Create(batch->nGridnodes));
	batch->thermostat = std::make_unique<Thermostat>(batch->nPclusters);
	UploadDeviceThermostatSimulations();
	if (batch->step == 0) BootstrapClustering(cudaStreams[0]);
	else RunClustering(cudaStreams[0]);
	MakeSuperClusterTasksGPU(cudaStreams[0]);
	InitializePME();
	UploadNbChunks();
	UploadIntegrationSimulationData();
	Synchronize();
}

namespace {
	__global__ void PackRenderPositions(const PersistentCluster* source, Float3* destination,
		int pclusterOffset, int positionCount) {
		const int positionId = blockIdx.x * blockDim.x + threadIdx.x;
		if (positionId >= positionCount) return;
		const int pclusterId = positionId / PersistentCluster::maxParticles;
		const int lane = positionId % PersistentCluster::maxParticles;
		destination[positionId] = source[pclusterOffset + pclusterId].pqd[lane].position;
	}
}

void Engine::StopRenderDataPipes() {
	for (auto* pipe : renderDataPipes)
		if (pipe)
			pipe->Stop();
}

void Engine::step() { // step isnt an accurate name
	if (!IsFinished())
		FullStep(); // isnt accurate either..
}

int Engine::StepUntilHostWork() {
	if (IsFinished()) return 0;
	FullStep();
	int nSteps = 1;
	if (nSteps < maxStepsPerCall && IsPlainReplayStep()) {
		do {
			LaunchStepGraph(stepGraphs[0][0][DeviceThermostatStep()], false);
			AdvanceStepCounters();
			for (auto& sim : batch->simulations)
				if (sim.device.active) sim.runstatus.current_step = sim.step;
			nSteps++;
		} while (nSteps < maxStepsPerCall && IsPlainReplayStep());
		cudaStreamQuery(cudaStreams[0]); // Submits the queued steps at once, see LaunchStepGraph
		LIMA_UTILS::genericErrorCheckNoSync("Error after replayed steps");
	}
	return nSteps;
}

void Engine::AdvanceStepCounters() {
	++batch->step;
	for (auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		++sim.simulation->step;
		sim.step = sim.simulation->getStep();
	}
}

bool Engine::IsPlainReplayStep() const {
	if (!replayableStepGraph[DeviceThermostatStep()] || mode != EngineRunMode::Simulation || IsFinished())
		return false;
	const auto& params = batch->params;
	const int64_t step = batch->step;	// The step to run, the steps' host work below sees step + 1
	if (params.em_variant || !params.snf_select.empty() || batch->nSuperclusters == 0)
		return false;
	if (params.data_logging_interval != 0 && step % params.data_logging_interval == 0)
		return false;
	if (params.apply_barostat && (step + 1) % params.steps_per_pressure_measurement == 0)
		return false;
	if ((!UsesDeviceThermostat() && (step + 1) % params.steps_per_temperature_measurement == 0) || (step + 1) % params.stepsPerNlistupdate == 0)
		return false;
	if ((step + 1) % StepsPerRender == 0 && std::any_of(renderDataPipes.begin(), renderDataPipes.end(), [](auto* pipe) { return pipe != nullptr; }))
		return false;
	for (const auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		const auto& simParams = sim.simulation->simParams;
		if (DatabuffersDeviceController::IsBufferFull(sim.step + 1, simParams.data_logging_interval)
			|| sim.device.thermostatScalar != 1.f || sim.step + 1 >= simParams.n_steps
			|| sim.runstatus.critical_error_occured || sim.runstatus.simulation_finished)
			return false;
	}
	return true;
}

void Engine::FullStep() {
	LIMA_UTILS::genericErrorCheckNoSync("Error before step");
	if (mode == EngineRunMode::Interactive) {
		auto& sim = batch->simulations.front();
		// Live editing changes EM/MD and force selections between steps.
		const bool wasEnergyMinimizing = batch->params.em_variant;
		batch->params = sim.simulation->simParams;
		if (batch->params.apply_barostat)
			throw std::invalid_argument("Pressure coupling is not supported during live editing");
		if (batch->params.em_variant && !wasEnergyMinimizing)
			ResetEnergyMinimization();
		if (sim.device.dt != batch->params.dt) {
			sim.device.dt = batch->params.dt;
			UploadIntegrationSimulationData();
		}
	}
	deviceMaster();
	AdvanceStepCounters();
	PublishRenderData();
	const bool rebuilt = hostMaster();
	flushNextStepLaunch |= rebuilt;
	if (!rebuilt && !IsFinished() && batch->step % batch->params.stepsPerNlistupdate == 0) {
		batch->superClustersControl->Reset(batch->nGridnodes, cudaStreams[0]);
		RunClustering(cudaStreams[0]);
		MakeSuperClusterTasksGPU(cudaStreams[0]);
		flushNextStepLaunch = true;
	}
	LIMA_UTILS::genericErrorCheckNoSync("Error after step");
}

void Engine::PublishRenderData() {
	if (batch->step % StepsPerRender != 0)
		return;
	for (size_t simulationIndex = 0; simulationIndex < batch->simulations.size(); ++simulationIndex)
		PublishRenderData(simulationIndex);
}

void Engine::PublishRenderData(size_t simulationIndex) {
	auto* pipe = renderDataPipes[simulationIndex];
	auto& simulation = batch->simulations[simulationIndex];
	if (!pipe || !simulation.device.active) return;
	Float3* destination = pipe->TryBeginWrite();
	if (!destination) return;
	StoreIntegrationStates(cudaStreams[0]);
	const int count = simulation.device.pclusters.count * PersistentCluster::maxParticles;
	PackRenderPositions<<<(count + 255) / 256, 256, 0, cudaStreams[0]>>>(
		batch->pClusterDevice.Get(), destination, simulation.device.pclusters.offset, count);
	if (cudaPeekAtLastError() != cudaSuccess) {
		pipe->CancelWrite();
		LIMA_UTILS::genericErrorCheckNoSync("Could not pack render positions");
		return;
	}
	pipe->Publish(cudaStreams[0], simulation.step);
}

bool Engine::hostMaster() {
	bool retired = false;
	bool thermostatScalarsChanged = false;
	float barostatCoordinateScale = 1.f;
	// Temperature cadence is independent of logging, otherwise the logging interval changes the thermostat's dynamics
	const bool measureTemperature = !UsesDeviceThermostat() && batch->step % batch->params.steps_per_temperature_measurement == 0;
	const bool measurePressure = batch->params.apply_barostat
		&& batch->step % batch->params.steps_per_pressure_measurement == 0;
	if (measureTemperature || measurePressure) {
		flushNextStepLaunch = true; // The measurements wait for the GPU
		StoreIntegrationStates(cudaStreams[0]);
		batch->thermostat->ComputeKineticEnergy(batch->boxState.pclusterInterimStates, batch->pClusterMetaDevice.Get(),
			batch->nPclusters, cudaStreams[0]);
	}
	for (auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		const auto step = sim.step;
		const auto& params = sim.simulation->simParams;
		if (DatabuffersDeviceController::IsBufferFull(step, params.data_logging_interval))
			OffloadLoggingData(sim);
		if (UsesDeviceThermostat())
			HarvestTemperatures(sim);
		// v-rescale is one impulse on a coupling step, not a factor spread over the interval.
		if ((!measureTemperature || !params.apply_thermostat) && sim.device.thermostatScalar != 1.f) {
			sim.device.thermostatScalar = 1.f;
			thermostatScalarsChanged = true;
		}
		double kineticEnergy = 0.0;
		double virial = 0.0;
		if (measurePressure) {
			const int kineticCount = sim.device.pclusters.count * PersistentCluster::maxParticles;
			const int reciprocalCount = batch->pmeController ? batch->pmeController->VirialPartialCount() : 0;
			const auto sums = CubWrappers::SumPressure(batch->thermostat->KineticEnergyValues(sim.device.pclusters), kineticCount,
				batch->nbVirialPartials.Get(), batch->nSuperclusters * batch->nbMaxChunks,
				batch->nBondgroups > 0 ? batch->bondVirialPartials.Get() : nullptr, batch->nBondgroups,
				batch->hasConstraints ? batch->constraintVirialPartials.Get() : nullptr, batch->hasConstraints ? batch->nSuperclusters : 0,
				reciprocalCount > 0 ? batch->pmeVirialPartials.Get() : nullptr, reciprocalCount, cudaStreams[0]);
			kineticEnergy = sums.kineticEnergy;
			virial = sums.virial * Virial::fixedPointScaleInv;
		}
		else if (measureTemperature)
			kineticEnergy = batch->thermostat->KineticEnergy(sim.device.pclusters, cudaStreams[0]);
		if (measureTemperature) {
			const auto& boxparams = sim.simulation->box->boxparams;
			const float temperature = Thermostat::Temperature(boxparams, kineticEnergy);
			sim.simulation->temperature_buffer.push_back(temperature);
			sim.runstatus.current_temperature = temperature;
			if (params.apply_thermostat) {
				const int64_t interval = params.steps_per_temperature_measurement;
				const double tauInCouplings = params.tau_t * PICO_TO_NANO / (interval * double(params.dt));
				const double factor = sim.thermostat.Factor(temperature, params.ref_t, boxparams.degreesOfFreedom, tauInCouplings);
				const float scalar = static_cast<float>(factor);
				thermostatScalarsChanged |= sim.device.thermostatScalar != scalar;
				sim.device.thermostatScalar = scalar;
			}
		}
		if (measurePressure) {
			const auto& boxparams = sim.simulation->box->boxparams;
			const Float3 boxLength = boxparams.BoxSizeFloat();
			const double volume = double(boxLength.x) * boxLength.y * boxLength.z;
			constexpr double pressureUnit = 1e22 / AVOGADROSNUMBER; // (J/mol)/nm^3 to bar
			const double pressure = (2.0 * kineticEnergy + virial) / (3.0 * volume) * pressureUnit;
			sim.simulation->pressure_buffer.push_back(static_cast<float>(pressure));
			sim.runstatus.current_pressure = static_cast<float>(pressure);
			const double couplingTime = params.steps_per_pressure_measurement * double(params.dt);
			const double tau = params.tau_p * PICO_TO_NANO;
			const double volumeScale = 1.0 + params.compressibility * couplingTime / tau * (pressure - params.ref_p);
			if (!(volumeScale > 0.0) || !std::isfinite(volumeScale))
				throw std::runtime_error("Barostat produced a non-positive or non-finite volume scale");
			barostatCoordinateScale = static_cast<float>(std::cbrt(volumeScale));
		}
		HandleEarlyStoppingInEM(sim);
		sim.runstatus.current_step = step;
		if (step >= params.n_steps || sim.runstatus.critical_error_occured) sim.runstatus.simulation_finished = true;
		if (mode == EngineRunMode::Simulation && sim.runstatus.simulation_finished) {
			FinalizeSimulation(sim);
			retired = true;
		}
	}
	rebuildPending |= retired;
	if (rebuildPending && !IsFinished() && batch->step % batch->params.stepsPerNlistupdate == 0) {
		RebuildActiveBatch(); // Uploads the integration data itself
		return true;
	}
	if (thermostatScalarsChanged)
		UploadIntegrationSimulationData();
	if (measurePressure && !IsFinished())
		return ApplyBarostat(barostatCoordinateScale);
	return false;
}

bool Engine::ApplyBarostat(float coordinateScale) {
	const Float3 oldLength = batch->boxLength;
	const Float3 newLength = oldLength * coordinateScale;
	if (!(newLength.x > 2.f * (batch->params.cutoff_nm + batch->listBuffer)
		&& newLength.y > 2.f * (batch->params.cutoff_nm + batch->listBuffer)
		&& newLength.z > 2.f * (batch->params.cutoff_nm + batch->listBuffer)))
		throw std::runtime_error("Barostat shrank the box below twice the neighbor-list radius");
	const int slots = batch->nSuperclusters * SuperCluster::maxParticles;
	ScaleBoxKernel<<<(slots + 255) / 256, 256, 0, cudaStreams[0]>>>(
		batch->superClustersControl->scData, batch->nSuperclusters, coordinateScale);
	auto& boxparams = batch->simulations.front().simulation->box->boxparams;
	batch->positionsConstrained = false;
	const bool engineGridNeedsRebuild = GridSizeRequiresRebuild(batch->boxSize, BoxParams::GridSizeForLength(newLength));
	const bool dimensionsChanged = engineGridNeedsRebuild
		|| (batch->pmeController && batch->pmeController->RequiresRebuild(newLength));
	if (dimensionsChanged) {
		boxparams.SetBoxLength(newLength);
		RebuildActiveBatch();
		return true;
	}
	boxparams.boxLength = newLength;
	batch->boxLength = newLength;
	if (batch->pmeController)
		batch->pmeController->UpdateBoxLength(newLength);
	batch->superClustersControl->Reset(batch->nGridnodes, cudaStreams[0]);
	RunClustering(cudaStreams[0]);
	MakeSuperClusterTasksGPU(cudaStreams[0]);
	return true;
}

void Engine::UploadIntegrationSimulationData() {
	InvalidateStepGraph();
	std::vector<IntegrationSimulationData> simulationData;
	simulationData.reserve(batch->simulations.size());
	for (const auto& sim : batch->simulations) {
		simulationData.push_back({
			sim.device.particles, sim.device.pclusters, sim.device.logOffset, sim.device.dt, sim.device.thermostatScalar });
	}
	batch->integrationSimulationDataDevice.Expand(simulationData.size());
	// Queued on the stream used by the integration kernel, so it is ordered before the next integration
	cudaMemcpyAsync(batch->integrationSimulationDataDevice.Get(), simulationData.data(),
		sizeof(IntegrationSimulationData) * simulationData.size(), cudaMemcpyHostToDevice, cudaStreams[0]);
}

void Engine::FinalizeSimulation(EngineSimulationData& sim) {
	if (sim.finalized) return;
	StoreIntegrationStates(cudaStreams[0]);
	Synchronize();
	// Results are only valid if no kernel dropped entries at a capacity limit
	if (batch->pclusterTransfermodule) batch->pclusterTransfermodule->overflow.Check();
	if (batch->pmeController)
		if (const CapacityOverflow* overflow = batch->pmeController->Overflow()) overflow->Check();
	OffloadLoggingData(sim);
	JoinLogDrains();
	if (UsesDeviceThermostat()) {
		// The temperature after its last step is measured by the next step, which it does not run
		if (sim.device.active && deviceThermostat->MeasuresAt(batch->step) && batch->nSuperclusters > 0) {
			deviceThermostat->Queue(batch->integrationStates.Get(), batch->nSuperclusters * SuperCluster::maxParticles, batch->pClusterMetaDevice.Get(),
				batch->thermostat->KineticEnergyValues(), stepCounterDevice, cudaStreams[0], batch->step, static_cast<int>(&sim - batch->simulations.data()));
			Synchronize();
		}
		HarvestTemperatures(sim);
	}
	const auto range = sim.device.pclusters;
	cudaMemcpy(sim.simulation->box->pclusterInterimStates.data(), batch->boxState.pclusterInterimStates + range.offset,
		sizeof(PersistentclusterInterimState) * range.count, cudaMemcpyDeviceToHost);
	cudaMemcpy(sim.simulation->box->persistentClusters.data(), batch->pClusterDevice.Get() + range.offset,
		sizeof(PersistentCluster) * range.count, cudaMemcpyDeviceToHost);
	if (sim.device.particles.count > 0)
		sim.finalForcesMagnitudeSquared = GenericCopyToHost(batch->forcesMagnitudeSquareDevice.Get() + sim.device.particles.offset, sim.device.particles.count);
	sim.device.active = false;
	sim.runstatus.simulation_finished = true;
	sim.finalized = true;
	UploadDeviceThermostatSimulations(); // It no longer measures this simulation
}

bool Engine::UsesDeviceThermostat() const {
	return deviceThermostat && mode == EngineRunMode::Simulation && !batch->params.em_variant && !batch->params.apply_barostat;
}

bool Engine::DeviceThermostatStep() const {
	return UsesDeviceThermostat() && deviceThermostat->MeasuresAt(batch->step);
}

void Engine::UploadDeviceThermostatSimulations() {
	if (!deviceThermostat)
		return;
	const int interval = std::max<int>(1, static_cast<int>(batch->params.steps_per_temperature_measurement.value));
	std::vector<_Thermostat::DeviceThermostatSimulation> simulations;
	std::vector<BatchRange> pclusters;
	for (const auto& sim : batch->simulations) {
		const SimParams& params = sim.simulation->simParams;
		_Thermostat::DeviceThermostatSimulation simulation;
		simulation.measure = sim.device.active && !sim.finalized && sim.device.pclusters.count > 0;
		simulation.degreesOfFreedom = sim.simulation->box->boxparams.degreesOfFreedom;
		simulation.referenceTemperature = params.ref_t;
		simulation.tauInCouplings = params.tau_t * PICO_TO_NANO / (interval * double(params.dt));	// As hostMaster
		simulation.applyThermostat = params.apply_thermostat;
		simulations.push_back(simulation);
		pclusters.push_back(simulation.measure ? sim.device.pclusters : BatchRange{});
	}
	deviceThermostat->SetSimulations(simulations, pclusters);
}

void Engine::HarvestTemperatures(EngineSimulationData& sim) {
	if (const auto latest = deviceThermostat->Harvest(&sim - batch->simulations.data(), sim.simulation->temperature_buffer))
		sim.runstatus.current_temperature = *latest;
}

void Engine::terminateSimulation() {
	for (auto& sim : batch->simulations) FinalizeSimulation(sim);
	StopRenderDataPipes();
	LIMA_UTILS::genericErrorCheckNoSync("Error during TerminateSimulation");
}

struct Engine::LogDrain {
	struct Unpack { void* dst; size_t offset; size_t bytes; };
	char* device = nullptr;		// Snapshot of the simulation's part of the logging ring
	char* pinned = nullptr;
	size_t capacity = 0;
	cudaEvent_t snapshotted = nullptr;
	std::array<cudaEvent_t, 2> copied{};	// Of the last two slices issued
	cudaError_t error = cudaSuccess;

	// One thread for the drain's lifetime: a thread's first CUDA call initializes its runtime state, which on WDDM blocked
	// the main thread's launches for ~10 ms, so a new thread per drain stalled the steps (T4 every 100 steps)
	std::thread unpacker;
	std::mutex mutex;
	std::condition_variable changed;
	std::function<void()> job;	// Set by Start, cleared by the thread once done
	bool stop = false;

	LogDrain() {
		cudaEventCreateWithFlags(&snapshotted, cudaEventDisableTiming);
		for (auto& event : copied) cudaEventCreateWithFlags(&event, cudaEventDisableTiming);
		int deviceId = 0;
		cudaGetDevice(&deviceId);
		unpacker = std::thread([this, deviceId] {
			cudaSetDevice(deviceId); // Also initializes the thread's runtime state, before the steps run
			std::unique_lock lock(mutex);
			while (true) {
				changed.wait(lock, [this] { return stop || job; });
				if (!job) return;
				lock.unlock();
				job();
				lock.lock();
				job = nullptr;
				changed.notify_all();
			}
		});
	}
	~LogDrain() {
		{
			std::lock_guard lock(mutex);
			stop = true;
		}
		changed.notify_all();
		unpacker.join(); // Finishes a started job first
		cudaFree(device);
		cudaFreeHost(pinned);
		cudaEventDestroy(snapshotted);
		for (auto event : copied) cudaEventDestroy(event);
	}
	// Must be joined
	void Start(std::function<void()> work) {
		{
			std::lock_guard lock(mutex);
			job = std::move(work);
		}
		changed.notify_all();
	}
	void Join() {
		{
			std::unique_lock lock(mutex);
			changed.wait(lock, [this] { return !job; });
		}
		if (error != cudaSuccess) {
			const cudaError_t e = error;
			error = cudaSuccess;
			throw std::runtime_error(std::string("Error offloading logging data: ") + cudaGetErrorString(e));
		}
	}
	// Must be joined
	void Reserve(size_t bytes) {
		if (bytes <= capacity) return;
		cudaFree(device);
		cudaFreeHost(pinned);
		device = pinned = nullptr;
		capacity = 0;
		LIMA_UTILS::genericErrorCheck(cudaMalloc(&device, bytes));
		LIMA_UTILS::genericErrorCheck(cudaMallocHost(&pinned, bytes));
		capacity = bytes;
	}
};

void Engine::JoinLogDrains() {
	for (auto& drain : logDrains)
		if (drain) drain->Join();
}

void Engine::OffloadLoggingData(EngineSimulationData& sim) {
	const int interval = batch->params.data_logging_interval;
	if (interval == 0 || sim.step == 0) return;
	const size_t entries = (sim.step - 1) / interval + 1;
	const size_t count = entries - sim.nLogEntriesTransferred;
	if (count == 0) return;
	if (count > DatabuffersDeviceController::nStepsInBuffer) throw std::runtime_error("Logging buffer was not drained");

	const size_t simulationIndex = &sim - batch->simulations.data();
	if (logDrains.size() <= simulationIndex) logDrains.resize(simulationIndex + 1);
	if (!logDrains[simulationIndex]) logDrains[simulationIndex] = std::make_unique<LogDrain>();
	LogDrain& drain = *logDrains[simulationIndex];
	drain.Join(); // Its buffers hold the previous drain until then

	// A full ring or its final partial segment needs only four contiguous copies.
	const size_t width = sim.device.pclusters.count * PersistentCluster::maxParticles;
	struct Copy { void* dst; const void* src; size_t bytes; };
	std::vector<Copy> copies;
	for (size_t entry = sim.nLogEntriesTransferred; entry < entries;) {
		const size_t slot = entry % DatabuffersDeviceController::nStepsInBuffer;
		const size_t chunk = std::min(entries - entry, DatabuffersDeviceController::nStepsInBuffer - slot);
		const size_t src = sim.device.logOffset + slot * width;
		copies.push_back({ sim.simulation->traj_buffer->getBufferAtIndex(entry), batch->dataBuffersDevice->traj_buffer + src, sizeof(Float3) * width * chunk });
		copies.push_back({ sim.simulation->potE_buffer->getBufferAtIndex(entry), batch->dataBuffersDevice->potE_buffer + src, sizeof(float) * width * chunk });
		copies.push_back({ sim.simulation->vel_buffer->getBufferAtIndex(entry), batch->dataBuffersDevice->vel_buffer + src, sizeof(float) * width * chunk });
		copies.push_back({ sim.simulation->forceBuffer->getBufferAtIndex(entry), batch->dataBuffersDevice->forceBuffer + src, sizeof(Float3) * width * chunk });
		entry += chunk;
	}
	size_t totalBytes = 0;
	for (const Copy& copy : copies) totalBytes += copy.bytes;
	drain.Reserve(totalBytes);

	// The snapshot is ordered after the steps that logged the entries, and before those that overwrite them
	constexpr size_t sliceBytes = size_t(8) << 20;
	std::vector<LogDrain::Unpack> slices;
	size_t offset = 0;
	for (const Copy& copy : copies) {
		cudaMemcpyAsync(drain.device + offset, copy.src, copy.bytes, cudaMemcpyDeviceToDevice, cudaStreams[0]);
		for (size_t sliceOffset = 0; sliceOffset < copy.bytes; sliceOffset += sliceBytes)
			slices.push_back({ static_cast<char*>(copy.dst) + sliceOffset, offset + sliceOffset, std::min(sliceBytes, copy.bytes - sliceOffset) });
		offset += copy.bytes;
	}
	cudaEventRecord(drain.snapshotted, cudaStreams[0]);
	cudaStreamWaitEvent(logCopyStream, drain.snapshotted);
	LIMA_UTILS::genericErrorCheckNoSync("Error queueing the logging data offload");

	// The thread copies the snapshot to the host on its own stream, so the steps continue meanwhile. It keeps only two slices queued:
	// the steps' own readbacks share the copy engine, and would otherwise wait for the whole copy. Each slice is unpacked while the
	// next is copied
	drain.Start([&drain, slices = std::move(slices), stream = logCopyStream] {
		auto issue = [&](size_t i) {
			cudaMemcpyAsync(drain.pinned + slices[i].offset, drain.device + slices[i].offset, slices[i].bytes, cudaMemcpyDeviceToHost, stream);
			cudaEventRecord(drain.copied[i % 2], stream);
		};
		for (size_t i = 0; i < std::min(slices.size(), size_t(2)); i++)
			issue(i);
		for (size_t i = 0; i < slices.size(); i++) {
			drain.error = cudaEventSynchronize(drain.copied[i % 2]);
			if (drain.error != cudaSuccess) return;
			if (i + 2 < slices.size())
				issue(i + 2);
			std::memcpy(slices[i].dst, drain.pinned + slices[i].offset, slices[i].bytes);
		}
	});
	sim.nLogEntriesTransferred = entries;
}

CudaBuffer<PersistentCluster>& Engine::OffloadPclusterState(size_t simulationId) {
	const auto& sim = batch->simulations.at(simulationId);
	const auto range = sim.device.pclusters;
	StoreIntegrationStates(cudaStreams[0]);
	Synchronize();
	const size_t count = sim.finalized ? sim.simulation->box->persistentClusters.size() : range.count;
	batch->pdataCopyBuffer.Expand(count);
	if (sim.finalized) cudaMemcpy(batch->pdataCopyBuffer.Get(), sim.simulation->box->persistentClusters.data(), sizeof(PersistentCluster) * count, cudaMemcpyHostToDevice);
	else cudaMemcpy(batch->pdataCopyBuffer.Get(), batch->pClusterDevice.Get() + range.offset, sizeof(PersistentCluster) * count, cudaMemcpyDeviceToDevice);
	return batch->pdataCopyBuffer;
}

CudaBuffer<float>& Engine::OffloadForcesMagnitudeBuffer(size_t simulationId) {
	const auto& sim = batch->simulations.at(simulationId);
	const auto range = sim.device.particles;
	StoreIntegrationStates(cudaStreams[0]);
	Synchronize();
	const size_t count = sim.finalized ? sim.finalForcesMagnitudeSquared.size() : range.count;
	batch->forcesMagnitudeCopyBuffer.Expand(count);
	if (sim.finalized) cudaMemcpy(batch->forcesMagnitudeCopyBuffer.Get(), sim.finalForcesMagnitudeSquared.data(), sizeof(float) * count, cudaMemcpyHostToDevice);
	else cudaMemcpy(batch->forcesMagnitudeCopyBuffer.Get(), batch->forcesMagnitudeSquareDevice.Get() + range.offset, sizeof(float) * count, cudaMemcpyDeviceToDevice);
	CubWrappers::SqrtInPlace(batch->forcesMagnitudeCopyBuffer.Get(), count);
	return batch->forcesMagnitudeCopyBuffer;
}

namespace {
	template<typename T>
	void SetParticleBuffer(std::optional<CudaBuffer<T>>& buffer, const std::vector<T>& data, BatchRange range, int totalParticles) {
		if (data.empty()) {
			buffer.reset();
			return;
		}
		if (data.size() != range.count) throw std::invalid_argument("Particle buffer size does not match simulation");
		if (!buffer) {
			buffer.emplace();
			buffer->Expand(totalParticles);
		}
		cudaMemcpy(buffer->Get() + range.offset, data.data(), sizeof(T) * data.size(), cudaMemcpyHostToDevice);
	}
}

void Engine::SetFixedParticleMovementBuffer(const std::vector<Float3>& movement, size_t simulationId) {
	Synchronize();
	InvalidateStepGraph();
	auto& sim = batch->simulations.at(simulationId);
	SetParticleBuffer(batch->fixedParticleMovementBuffer, movement, sim.device.particles, batch->nParticles);
}

void Engine::SetFixedParticleRotationBuffer(const std::vector<Rotation>& rotation, size_t simulationId) {
	Synchronize();
	InvalidateStepGraph();
	auto& sim = batch->simulations.at(simulationId);
	SetParticleBuffer(batch->fixedParticleRotationBuffer, rotation, sim.device.particles, batch->nParticles);
}

void Engine::SetForceMask(const std::vector<Float3>& mask, size_t simulationId) {
	Synchronize();
	InvalidateStepGraph();
	auto& sim = batch->simulations.at(simulationId);
	SetParticleBuffer(batch->forceMaskBuffer, mask, sim.device.particles, batch->nParticles);
}

void Engine::SetElasticPositions(const std::vector<Float3>& positions, size_t simulationId) {
	Synchronize();
	InvalidateStepGraph();
	auto& sim = batch->simulations.at(simulationId);
	SetParticleBuffer(batch->elasticPositionsBuffer, positions, sim.device.particles, batch->nParticles);
}

void Engine::BootstrapTrajbufferWithCoords(EngineSimulationData& sim) {
	if (batch->params.data_logging_interval == 0 || !sim.simulation->traj_buffer) return;
	for (int pc = 0; pc < sim.device.pclusters.count; ++pc)
		for (int lane = 0; lane < PersistentCluster::maxParticles; ++lane)
			sim.simulation->traj_buffer->GetDatapoint(pc, lane, 0) = sim.simulation->box->persistentClusters[pc].pqd[lane].position;
}

void Engine::HandleEarlyStoppingInEM(EngineSimulationData& sim) {
	if (!batch->params.em_variant || sim.step == sim.simulation->simParams.n_steps) return;
	const EM::SimState& state = batch->emStatesHost[&sim - batch->simulations.data()];
	sim.runstatus.greatestForce = state.maxForce / KILO;
	sim.simulation->maxForceBuffer.emplace_back(sim.step, sim.runstatus.greatestForce);
	Simulation::EmLogEntry& entry = sim.simulation->emLog.emplace_back();
	entry.step = sim.step;
	entry.maxForce = sim.runstatus.greatestForce;
	entry.dt = state.dt;
	if (!std::isfinite(sim.runstatus.greatestForce))
		throw std::runtime_error("Energy minimization produced non-finite forces at step " + std::to_string(sim.step));
	if (sim.runstatus.greatestForce <= batch->params.em_force_tolerance) sim.runstatus.simulation_finished = true;
}

void Engine::ResetEnergyMinimization() {
	Synchronize();
	cudaMemset(batch->emParticles.Get(), 0, sizeof(EM::ParticleState) * batch->nPclusters * PersistentCluster::maxParticles);
	cudaMemset(batch->emStates.Get(), 0, sizeof(EM::SimState) * batch->simulations.size());
}

// Moves the particles one preconditioned FIRE step, using the forces SuperclusterIntegrateKernel stored in emForces
template <typename BoundaryCondition>
void Engine::UpdateEnergyMinimization(Float3 boxSize) {
	const int nSimulations = static_cast<int>(batch->simulations.size());
	int maxSlots = 0;
	for (const auto& sim : batch->simulations)
		maxSlots = std::max(maxSlots, sim.device.pclusters.count * PersistentCluster::maxParticles);
	if (maxSlots == 0) return;

	EM::PreconditionKernel<BoundaryCondition><<<(batch->nPclusters + 63) / 64, 64, 0, cudaStreams[0]>>>(
		batch->emConfig, batch->pClusterDevice.Get(), batch->pClusterMetaDevice.Get(), batch->emForces.Get(),
		batch->emInverseStiffness.Get(), batch->emWholeMolecule.Get(), batch->emPreconditionedForce.Get(), batch->nPclusters, boxSize);
	constexpr int blockSize = 256;
	const int nBlocks = (maxSlots + blockSize - 1) / blockSize;
	batch->emBlockSums.Expand(size_t(nBlocks) * nSimulations);
	EM::ReduceAndDecideKernel<blockSize><<<dim3(nBlocks, nSimulations), blockSize, 0, cudaStreams[0]>>>(
		batch->emConfig, batch->integrationSimulationDataDevice.Get(), batch->pClusterMetaDevice.Get(), batch->emForces.Get(),
		batch->emPreconditionedForce.Get(), batch->emParticles.Get(), batch->emBlockSums.Get(), batch->emBlocksDone.Get(), batch->emStates.Get());
	EM::UpdateKernel<BoundaryCondition><<<(batch->nSuperclusters + 3) / 4, dim3(16, 4, 1), 0, cudaStreams[0]>>>(
		batch->emConfig, batch->emStates.Get(), batch->superClustersControl->scData, batch->superClustersControl->scMeta,
		batch->pClusterDevice.Get(), batch->emPreconditionedForce.Get(), batch->emParticles.Get(), batch->nSuperclusters, boxSize);
	cudaMemcpyAsync(batch->emStatesHost.data(), batch->emStates.Get(), sizeof(EM::SimState) * nSimulations, cudaMemcpyDeviceToHost, cudaStreams[0]);
}

template <typename BoundaryCondition, bool emvariant, bool logData, bool computeVirial>
void Engine::_deviceMaster(bool measureTemperature) {
	if (emvariant != tasksBuiltForEm)
		MakeSuperClusterTasksGPU(cudaStreams[0]);
	const Float3 boxSize = batch->boxLength;
	const int nScs = batch->nSuperclusters;
	if constexpr (!emvariant) {
		if (batch->hasConstraints && !batch->positionsConstrained && nScs > 0)
			ConstrainPositionsKernel<BoundaryCondition><<<batch->nGridnodes, dim3(SuperCluster::maxParticles, SuperClustersControl::maxClustersPerBlock, 1), 0, cudaStreams[0]>>>(
				*batch->superClustersControl, batch->integrationStates.Get(), batch->pclusterConstraintsDevice.Get(), batch->pclusterParticleSlots.Get(), boxSize);
		batch->positionsConstrained = true;
	}
	// Allocations can't be captured, so the buffers are expanded first. All are noops once allocated
	if constexpr (emvariant)
		batch->scResultsDevice.Expand(batch->nResults, 1.2); // Not allocated at all in MD, where it would be ~1GB for large systems
	if constexpr (computeVirial) {
		batch->pmeVirialPartials.Expand(batch->pmeController ? batch->pmeController->VirialPartialCount() : 0);
		batch->nbVirialPartials.Expand(nScs * batch->nbMaxChunks);
		batch->bondVirialPartials.Expand(batch->nBondgroups);
		if (batch->hasConstraints)
			batch->constraintVirialPartials.Expand(nScs);
	}
	// EM synchronizes within the step. SNFs are not checked for capturability
	const bool captureStep = useStepGraphs && !emvariant && nScs > 0 && batch->params.snf_select.empty();
	if (captureStep && !logData && !computeVirial && replayableStepGraph[measureTemperature]) {
		LaunchStepGraph(stepGraphs[0][0][measureTemperature], false);
		LIMA_UTILS::genericErrorCheckNoSync("Error during batch timestep");
		return;
	}
	if (captureStep)
		BeginStepCapture();
	// The previous step and any rebuild are queued on cudaStreams[0], and the host does not wait for them
	ForkStreamsFromMainStream();
	if constexpr (!emvariant) {
		AdvanceStepCounterKernel<<<1, 1, 0, cudaStreams[1]>>>(stepCounterDevice, batch->step);
		// After the counter, on its stream, overlapping the force kernels. Integration joins it
		if (measureTemperature && nScs > 0)
			deviceThermostat->Queue(batch->integrationStates.Get(), nScs * SuperCluster::maxParticles, batch->pClusterMetaDevice.Get(),
				batch->thermostat->KineticEnergyValues(), stepCounterDevice, cudaStreams[1]);
	}
	// MD: NB, PME and SNF add to forceAcc, which integration left zeroed. Bonded sums have one writer per group
	// and are gathered by integration, so all force kernels can run concurrently.
	ForceAccumulator forceAcc{};
	ForceAccumulator primaryBondForces{};
	ulonglong4* const extraBondForces = batch->extraBondForceResults.Get();
	if constexpr (!emvariant) {
		const size_t n = size_t(nScs) * SuperCluster::maxParticles;
		unsigned long long* const base = batch->forceAccumulatorDevice.Get();
		forceAcc = ForceAccumulator{ base, base + n, base + 2 * n, logData ? base + 3 * n : nullptr };
		primaryBondForces = ForceAccumulator{ base + 4 * n, base + 5 * n, base + 6 * n, logData ? base + 7 * n : nullptr };
	}
	if (ENABLE_ES_LR && batch->params.enable_electrostatics)
		batch->pmeController->CalcCharges(batch->superClustersControl->scData, batch->superClustersControl->scMeta,
			nScs, batch->forceEnergyInterims->pme, forceAcc, computeVirial ? batch->pmeVirialPartials.Get() : nullptr);
	if (nScs > 0) {
		const auto* scData = batch->superClustersControl->scData;
		const Float3 boxSizeInv = boxSize.Inv();
		// Blocksize must be 64, 2 warps per block, see the kernel
		const int nChunks = emvariant ? 1 : batch->nbMaxChunks;
		NbNonlocalKernel<BoundaryCondition, emvariant, logData, computeVirial><<<(nScs * nChunks + 1) / 2, 64, 0, cudaStreams[0]>>>(scData, batch->quarterEntryTasksDevice.Get(),
			batch->quarterEntriesDevice.Get(), forceAcc, batch->scResultsDevice.Get(), batch->quarterEntryResultIndicesDevice.Get(),
			batch->superClustersControl->scMeta, boxSize, boxSizeInv, batch->ewaldKappa, PhysicsUtilsDevice::EwaldRealSpace::Make(batch->ewaldKappa),
			batch->params.cutoff_nm * batch->params.cutoff_nm, nScs, nChunks, batch->nbSimulationChunks.Get(),
			computeVirial ? batch->nbVirialPartials.Get() : nullptr);
	}
	if (!batch->params.snf_select.empty()) SnfHandler<BoundaryCondition, emvariant>(cudaStreams[2], forceAcc);
	if (batch->nBondgroups > 0) {
		const CompactBondGroupsDevice bondGroups{
			batch->bondgroupDescriptors.Get(), batch->bondgroupSinglebonds.Get(), batch->bondgroupPairbonds.Get(),
			batch->bondgroupAnglebonds.Get(), batch->bondgroupDihedralbonds.Get(), batch->bondgroupImproperdihedralbonds.Get(),
			batch->bondgroupCmapbonds.Get(), batch->cmapGridPoints.Get(), batch->cmapGridSpacing
		};
		const Float3 boxSizeInv = boxSize.Inv();
		const auto kernel = batch->cmapGridSpacing > 0 ? BondgroupsKernel<BoundaryCondition, emvariant, logData, computeVirial, true>
			: BondgroupsKernel<BoundaryCondition, emvariant, logData, computeVirial, false>;
		// Warps per bondgroup, see the kernel. The split does not change the forces, so it may depend on the batch
		constexpr int targetWarps = 2048;
		int warpsPerGroup = emvariant ? 1 : std::clamp(targetWarps / batch->nBondgroups, 1, BONDGROUPS_PER_BLOCK);
		if (warpsPerGroup == 3) warpsPerGroup = 2; // Must divide BONDGROUPS_PER_BLOCK
		const int groupsPerBlock = BONDGROUPS_PER_BLOCK / warpsPerGroup;
		kernel<<<(batch->nBondgroups + groupsPerBlock - 1) / groupsPerBlock, THREADS_PER_BONDSGROUPSKERNEL, 0, cudaStreams[4]>>>(
			bondGroups, batch->nBondgroups, warpsPerGroup, batch->boxState, batch->forceEnergyInterims->forceEnergiesBondgroups, primaryBondForces, extraBondForces, batch->bondgroupExtraResultIndices.Get(), batch->superClustersControl->scData,
			batch->bondgroupParticleSlots.Get(), boxSize, boxSizeInv, ENABLE_ES_LR && batch->params.enable_electrostatics ? batch->ewaldKappa : 0.f,
			computeVirial ? batch->bondVirialPartials.Get() : nullptr);
	}
	if (nScs == 0) {
		Synchronize();
	}
	else {
		// Integration waits for the force kernels on the GPU, so the host can queue it without a roundtrip
		JoinStreamsIntoMainStream();

		const LogBuffers log{ batch->dataBuffersDevice->traj_buffer, batch->dataBuffersDevice->potE_buffer, batch->dataBuffersDevice->vel_buffer,
			batch->dataBuffersDevice->forceBuffer, batch->params.data_logging_interval, batch->step };
		const int nBlocks = (nScs + 4 - 1) / 4;
		if constexpr (emvariant) {
			EmCollectForcesKernel<BoundaryCondition, logData><<<nBlocks, dim3(16, 4, 1), 0, cudaStreams[0]>>>(
				*batch->forceEnergyInterims, batch->scResultsDevice.Get(), batch->emForces.Get(), batch->superClustersControl->scData,
				batch->superClustersControl->scMeta, batch->pClusterDevice.Get(), batch->pClusterMetaDevice.Get(),
				batch->integrationSimulationDataDevice.Get(), nScs, boxSize, batch->forcesMagnitudeSquareDevice.Get(), log);
			batch->forceMagnitudesInStates = false;
			UpdateEnergyMinimization<BoundaryCondition>(boxSize);
		}
		else {
			const LiveEditBuffers liveEdit{ batch->fixedParticleMovementBuffer ? batch->fixedParticleMovementBuffer->Get() : nullptr,
				batch->forceMaskBuffer ? batch->forceMaskBuffer->Get() : nullptr,
				batch->fixedParticleRotationBuffer ? batch->fixedParticleRotationBuffer->Get() : nullptr };
			if (batch->hasConstraints) {
				// forcePrev includes the constraint force, so the kernel writes the measured force magnitudes itself
				batch->forceMagnitudesInStates = false;
				SuperclusterIntegrateConstrainedKernel<BoundaryCondition, logData, computeVirial><<<batch->nGridnodes, dim3(SuperCluster::maxParticles, SuperClustersControl::maxClustersPerBlock, 1), 0, cudaStreams[0]>>>(
					forceAcc, primaryBondForces, extraBondForces, batch->slotExtraBondReferences.Get(), *batch->superClustersControl, batch->integrationStates.Get(),
					batch->integrationSimulationDataDevice.Get(), nScs, boxSize, batch->forcesMagnitudeSquareDevice.Get(),
					liveEdit, log, batch->pclusterConstraintsDevice.Get(), batch->pclusterParticleSlots.Get(), batch->pclusterTransfermodule->overflow, stepCounterDevice,
					computeVirial ? batch->constraintVirialPartials.Get() : nullptr, deviceThermostat->impulses.Get());
			}
			else {
				// Without live edits, which may change the force after it is measured, forcePrev is exactly the measured force
				batch->forceMagnitudesInStates = !liveEdit.Any();
				SuperclusterIntegrateKernel<BoundaryCondition, logData><<<nBlocks, dim3(16, 4, 1), 0, cudaStreams[0]>>>(
					forceAcc, primaryBondForces, extraBondForces, batch->slotExtraBondReferences.Get(), batch->superClustersControl->scData, batch->integrationStates.Get(), batch->superClustersControl->scMeta,
					batch->integrationSimulationDataDevice.Get(), nScs, boxSize, batch->forcesMagnitudeSquareDevice.Get(),
					liveEdit, log, batch->pclusterTransfermodule->overflow, stepCounterDevice, deviceThermostat->impulses.Get());
			}
		}
		// MD queues the next step without waiting, the host reads nothing from this one. EM reads emStatesHost in hostMaster
		if constexpr (emvariant)
			cudaStreamSynchronize(cudaStreams[0]);
		if (captureStep)
			LaunchCapturedStep(logData, computeVirial, measureTemperature);
	}
	LIMA_UTILS::genericErrorCheckNoSync("Error during batch timestep");
}

void Engine::deviceMaster() {
	const bool logData = batch->params.data_logging_interval != 0 && batch->step % batch->params.data_logging_interval == 0;
	const bool measureTemperature = DeviceThermostatStep();
	const bool pressureStep = batch->params.apply_barostat
		&& (batch->step + 1) % batch->params.steps_per_pressure_measurement == 0;
	switch (batch->params.bc_select) {
	case NoBC:
		if (batch->params.em_variant) {
			if (logData) {
				_deviceMaster<NoBoundaryCondition, true, true, false>(measureTemperature);
			}
			else {
				_deviceMaster<NoBoundaryCondition, true, false, false>(measureTemperature);
			}
		}
		else {
			if (logData) {
				_deviceMaster<NoBoundaryCondition, false, true, false>(measureTemperature);
			}
			else {
				_deviceMaster<NoBoundaryCondition, false, false, false>(measureTemperature);
			}
		}
		break;
	case PBC:
		if (batch->params.em_variant) {
			if (logData) {
				_deviceMaster<PeriodicBoundaryCondition, true, true, false>(measureTemperature);
			}
			else {
				_deviceMaster<PeriodicBoundaryCondition, true, false, false>(measureTemperature);
			}
		}
		else {
			if (logData) {
				if (pressureStep) _deviceMaster<PeriodicBoundaryCondition, false, true, true>(measureTemperature);
				else _deviceMaster<PeriodicBoundaryCondition, false, true, false>(measureTemperature);
			}
			else {
				if (pressureStep) _deviceMaster<PeriodicBoundaryCondition, false, false, true>(measureTemperature);
				else _deviceMaster<PeriodicBoundaryCondition, false, false, false>(measureTemperature);
			}
		}
		break;
	default:
		throw std::runtime_error("Unsupported boundary condition in LAUNCH_GENERIC_KERNEL");
	}
}


template <typename BoundaryCondition, bool emvariant>
void Engine::SnfHandler(cudaStream_t& stream, const ForceAccumulator& forceAcc) {
	const int nPcs = batch->nPclusters;
	if (nPcs == 0) return;
	const SnfForceOutput output = emvariant ? SnfForceOutput{ batch->forceEnergyInterims->snf } : SnfForceOutput{ nullptr, forceAcc };
	const ParticleSlots particleSlots{ batch->superClustersControl->scData, batch->pclusterParticleSlots.Get() };
	for (const auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		const int count = sim.device.pclusters.count;
		if (batch->params.snf_select.contains(HorizontalChargeField)) {
			PclusterSnfKernel<BoundaryCondition, emvariant><<<(count + 31) / 32, 32, 0, stream>>>(
				batch->pClusterDevice.Get(), batch->pClusterMetaDevice.Get(), sim.device.uniformElectricField,
				output, particleSlots, sim.device.pclusters.offset, count);
		}
		if (batch->params.snf_select.contains(ElasticPosition) && batch->elasticPositionsBuffer) {
			ElasticPositionsForceKernel<<<(count + 31) / 32, 32, 0, stream>>>(
				batch->pClusterMetaDevice.Get(), batch->elasticPositionsBuffer->Get(), output, particleSlots,
				sim.device.pclusters.offset, count, batch->boxLength);
		}
	}
}

// -------------------------- TESTS: Shouldn't be here permanently ---------------------------- //




template <int nBins, int nValuesPerBin>
__global__ void TestSortKernel32(float* keys, int* ids)
{
	//static_assert(nBins * nValuesPerBin == 32);

	// exactly one warp
	LAL::SortBins<nBins, nValuesPerBin>(keys, ids);
}
bool Engine::TestAlgorithms() {

	constexpr int totalValues = 64;


	bool success = true;
	auto runCase = [&](int nBins, int nValuesPerBin)
		{
			std::vector<float> hKeys(totalValues);
			std::vector<int> hIds(totalValues);
			std::vector<float> refKeys(totalValues);
			std::vector<int> refIds(totalValues);

			float* dKeys;
			int* dIds;
			cudaMalloc(&dKeys, totalValues * sizeof(float));
			cudaMalloc(&dIds, totalValues * sizeof(int));

			std::mt19937 rng(12345);
			std::uniform_real_distribution<float> dist(-1e6, 1e6);

			//hKeys = { 4,3,2,1,4,1,2,3,1,2,3,4,4,2,3,1,4,3,2,1,4,1,2,3,1,2,3,4,4,2,3,1 };
			for (int i = 0; i < totalValues; ++i) {
				hKeys[i] = dist(rng);
				hIds[i] = i;
			}


			{
				refKeys = hKeys;
				refIds = hIds;

				std::vector<std::pair<float, int>> keyIdPairs;
				for (int i = 0; i < totalValues; ++i) {
					keyIdPairs.push_back({ refKeys[i], refIds[i] });
				}

				// CPU reference:
				for (int i = 0; i < nBins; i++) {
					std::sort(
						keyIdPairs.begin() + i * nValuesPerBin,
						keyIdPairs.begin() + (i + 1) * nValuesPerBin,
						[](const std::pair<float, int>& a, const std::pair<float, int>& b) {
							return a.first < b.first;
						}
					);
				}
				for (int i = 0; i < totalValues; ++i) {
					refKeys[i] = keyIdPairs[i].first;
					refIds[i] = keyIdPairs[i].second;
				}
			}





			cudaMemcpy(dKeys, hKeys.data(), totalValues * sizeof(float), cudaMemcpyHostToDevice);
			cudaMemcpy(dIds, hIds.data(), totalValues * sizeof(int), cudaMemcpyHostToDevice);

			if (nBins == 1)
				TestSortKernel32<1, 64> <<<1, 32 >>> (dKeys, dIds);
			else if (nBins == 4)
				TestSortKernel32<4, 16> <<<1, 32 >>> (dKeys, dIds);
			else if (nBins == 8)
				TestSortKernel32<8, 8> <<<1, 32 >>> (dKeys, dIds);


			cudaStreamSynchronize(nullptr);

			cudaMemcpy(hKeys.data(), dKeys, totalValues * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(hIds.data(), dIds, totalValues * sizeof(int), cudaMemcpyDeviceToHost);

			for (int i = 0; i < totalValues; ++i) {
				if (hKeys[i] != refKeys[i]) {
					success = false;
				}
			if (hIds[i] != refIds[i]) {
				success = false;
				}
			}
			

			cudaFree(dKeys);
			cudaFree(dIds);
		};

	// 64 values total, tested via 2x32
	runCase(1, 64);   // effectively 1x64
	runCase(4, 16);    // effectively 4x16
	runCase(8, 8);   // effectively 16x4

	return success;
}

namespace EngineLimitTesting {
	__global__ void DeviceFailureKernel() { DeviceIntrinsics::Trap(); }
}

void Engine::TestLimit(EngineLimitProbe probe, int count) {
	using namespace EngineLimitTesting;
	switch (probe) {
	case EngineLimitProbe::ClusterTransfer:
		Require(count >= 1 && count <= PClusterTransfermodule::maxOutgoingClusters + 1, "Invalid transfer fixture size");
		ClusterTransfer(count);
		break;
	case EngineLimitProbe::ClusterOccupancy:
		Require(count >= 1 && count <= PClusterTransfermodule::maxClustersPerBlock + 1, "Invalid occupancy fixture size");
		ClusterOccupancy(count);
		break;
	case EngineLimitProbe::ChargeBlock:
		Require(count >= 1 && count <= 385, "Invalid PME fixture size");
		EngineLimitTesting::ChargeBlock(count);
		break;
	case EngineLimitProbe::DeviceFailure:
		DeviceFailureKernel<<<1, 1>>>();
		break;
	}
}
