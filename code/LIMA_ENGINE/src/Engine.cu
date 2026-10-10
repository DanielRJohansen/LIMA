#include "Engine.cuh"
#include <cstring>
#include <thread>

#include "BoundaryCondition.cuh"
#include "EngineBodies.cuh"
#include "EngineKernels.cuh"
#include "EnergyMinimization.cuh"
#include "LimaPositionSystem.cuh"
#include "PME.cuh"
#include "SimulationData.h"
#include "SupernaturalForces.cuh"
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
	try {
		if (UsesEnergyMinimization()) UploadEnergyMinimizationPreconditioner();
		for (auto& stream : cudaStreams) cudaStreamCreate(&stream);
		cudaStreamCreate(&pmeStream);
		cudaStreamCreateWithFlags(&logCopyStream, cudaStreamNonBlocking); // Not synchronized with the legacy default stream, which would wait for the copy
		for (auto& event : streamJoinEvents) cudaEventCreateWithFlags(&event, cudaEventDisableTiming);
		cudaEventCreateWithFlags(&stepStartEvent, cudaEventDisableTiming);
		for (size_t simulationId = 0; simulationId < renderDataPipes.size(); ++simulationId) {
			if (!renderDataPipes[simulationId]) continue;
			const auto count = batch->simulations[simulationId].device.pclusters.count * PersistentCluster::maxParticles;
			renderDataPipes[simulationId]->Initialize(count);
		}
		batch->dataBuffersDevice = std::make_unique<DatabuffersDeviceController>(batch->nPclusters, batch->params.data_logging_interval);
		batch->superClustersControl = std::make_unique<SuperClustersControl>(batch->nGridnodes, batch->nPclusters);
		batch->pclusterTransfermodule = std::make_unique<PClusterTransfermodule>(PClusterTransfermodule::Create(batch->nGridnodes));
		batch->thermostat = std::make_unique<Thermostat>(batch->nPclusters);
		BootstrapClustering(cudaStreams[0]);
		MakeSuperClusterTasksGPU(cudaStreams[0]);
		InitializePME();
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

const RunStatus& Engine::GetRunStatus(size_t simulationId) const {
	return batch->simulations.at(simulationId).runstatus;
}

void Engine::StopSimulation(size_t simulationId) {
	FinalizeSimulation(batch->simulations.at(simulationId));
}

bool Engine::IsFinished() const {
	return std::none_of(batch->simulations.begin(), batch->simulations.end(), [](const auto& sim) { return sim.device.active; });
}

void Engine::InitializePME() {
	if (ENABLE_ES_LR && batch->params.enable_electrostatics) {
		batch->pmeController = std::make_unique<PME::Controller>(batch->simulations, batch->params.cutoff_nm, pmeStream);
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
	for (size_t id = 0; id < simulations.size(); ++id) {
		const auto& oldSim = batch->simulations[id];
		auto& newSim = rebuilt->simulations[id];
		newSim.runstatus = oldSim.runstatus;
		newSim.step = oldSim.step;
		newSim.nLogEntriesTransferred = oldSim.nLogEntriesTransferred;
		newSim.finalized = oldSim.finalized;
		newSim.finalForcesMagnitudeSquared = oldSim.finalForcesMagnitudeSquared;
		newSim.device.thermostatScalar = oldSim.device.thermostatScalar;
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
	if (UsesEnergyMinimization()) UploadEnergyMinimizationPreconditioner();
	batch->dataBuffersDevice = std::make_unique<DatabuffersDeviceController>(batch->nPclusters, batch->params.data_logging_interval);
	batch->superClustersControl = std::make_unique<SuperClustersControl>(batch->nGridnodes, batch->nPclusters);
	batch->pclusterTransfermodule = std::make_unique<PClusterTransfermodule>(PClusterTransfermodule::Create(batch->nGridnodes));
	batch->thermostat = std::make_unique<Thermostat>(batch->nPclusters);
	if (batch->step == 0) BootstrapClustering(cudaStreams[0]);
	else RunClustering(cudaStreams[0]);
	MakeSuperClusterTasksGPU(cudaStreams[0]);
	InitializePME();
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

void Engine::step() {
	if (IsFinished()) return;
	LIMA_UTILS::genericErrorCheckNoSync("Error before step");
	if (mode == EngineRunMode::Interactive) {
		auto& sim = batch->simulations.front();
		// Live editing changes EM/MD and force selections between steps.
		const bool wasEnergyMinimizing = batch->params.em_variant;
		batch->params = sim.simulation->simParams;
		if (batch->params.em_variant && !wasEnergyMinimizing)
			ResetEnergyMinimization();
		if (sim.device.dt != batch->params.dt) {
			sim.device.dt = batch->params.dt;
			UploadIntegrationSimulationData();
		}
	}
	deviceMaster();
	++batch->step;
	for (auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		++sim.simulation->step;
		sim.step = sim.simulation->getStep();
	}
	PublishRenderData();
	const bool rebuilt = hostMaster();
	if (!rebuilt && !IsFinished() && batch->step % batch->params.stepsPerNlistupdate == 0) {
		batch->superClustersControl->Reset(batch->nGridnodes, cudaStreams[0]);
		RunClustering(cudaStreams[0]);
		MakeSuperClusterTasksGPU(cudaStreams[0]);
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
	// Temperature cadence is independent of logging, otherwise the logging interval changes the thermostat's dynamics
	const bool measureTemperature = batch->step % batch->params.steps_per_temperature_measurement == 0;
	if (measureTemperature) {
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
		if (measureTemperature) {
			auto [temperature, scalar] = batch->thermostat->Temperature(
				sim.simulation->box->boxparams, params, sim.device.pclusters, cudaStreams[0]);
			sim.simulation->temperature_buffer.push_back(temperature);
			sim.runstatus.current_temperature = temperature;
			if (params.apply_thermostat) {
				thermostatScalarsChanged |= sim.device.thermostatScalar != scalar;
				sim.device.thermostatScalar = scalar;
			}
		}
		HandleEarlyStoppingInEM(sim);
		sim.runstatus.current_step = step;
		if (step >= params.n_steps || sim.runstatus.critical_error_occured) sim.runstatus.simulation_finished = true;
		if (mode == EngineRunMode::Simulation && sim.runstatus.simulation_finished) {
			FinalizeSimulation(sim);
			retired = true;
		}
	}
	if (retired && !IsFinished()) {
		RebuildActiveBatch(); // Uploads the integration data itself
		return true;
	}
	if (thermostatScalarsChanged)
		UploadIntegrationSimulationData();
	return false;
}

void Engine::UploadIntegrationSimulationData() {
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
	std::thread unpacker;
	cudaError_t error = cudaSuccess;

	LogDrain() {
		cudaEventCreateWithFlags(&snapshotted, cudaEventDisableTiming);
		for (auto& event : copied) cudaEventCreateWithFlags(&event, cudaEventDisableTiming);
	}
	~LogDrain() {
		if (unpacker.joinable()) unpacker.join();
		cudaFree(device);
		cudaFreeHost(pinned);
		cudaEventDestroy(snapshotted);
		for (auto event : copied) cudaEventDestroy(event);
	}
	void Join() {
		if (unpacker.joinable()) unpacker.join();
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
	drain.unpacker = std::thread([&drain, slices = std::move(slices), stream = logCopyStream] {
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
	auto& sim = batch->simulations.at(simulationId);
	SetParticleBuffer(batch->fixedParticleMovementBuffer, movement, sim.device.particles, batch->nParticles);
}

void Engine::SetFixedParticleRotationBuffer(const std::vector<Rotation>& rotation, size_t simulationId) {
	Synchronize();
	auto& sim = batch->simulations.at(simulationId);
	SetParticleBuffer(batch->fixedParticleRotationBuffer, rotation, sim.device.particles, batch->nParticles);
}

void Engine::SetForceMask(const std::vector<Float3>& mask, size_t simulationId) {
	Synchronize();
	auto& sim = batch->simulations.at(simulationId);
	SetParticleBuffer(batch->forceMaskBuffer, mask, sim.device.particles, batch->nParticles);
}

void Engine::SetElasticPositions(const std::vector<Float3>& positions, size_t simulationId) {
	Synchronize();
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

template <typename BoundaryCondition, bool emvariant, bool logData>
void Engine::_deviceMaster() {
	if (emvariant != tasksBuiltForEm)
		MakeSuperClusterTasksGPU(cudaStreams[0]);
	// The previous step and any rebuild are queued on cudaStreams[0], and the host does not wait for them
	ForkStreamsFromMainStream();
	const Float3 boxSize = NodeIndex(batch->boxSize).toFloat3();
	const int nScs = batch->nSuperclusters;
	// MD: NB, PME and SNF add to forceAcc, which integration left zeroed. Bonded sums have one writer per group
	// and are gathered by integration, so all force kernels can run concurrently.
	ForceAccumulator forceAcc{};
	ForceAccumulator primaryBondForces{};
	ulonglong4* const extraBondForces = batch->extraBondForceResults.Get();
	if constexpr (emvariant) {
		batch->scResultsDevice.Expand(batch->nResults, 1.2); // Noop once allocated. Not allocated at all in MD, where it would be ~1GB for large systems
	}
	else {
		const size_t n = size_t(nScs) * SuperCluster::maxParticles;
		unsigned long long* const base = batch->forceAccumulatorDevice.Get();
		forceAcc = ForceAccumulator{ base, base + n, base + 2 * n, logData ? base + 3 * n : nullptr };
		primaryBondForces = ForceAccumulator{ base + 4 * n, base + 5 * n, base + 6 * n, logData ? base + 7 * n : nullptr };
	}
	if (ENABLE_ES_LR && batch->params.enable_electrostatics)
		batch->pmeController->CalcCharges(batch->superClustersControl->scData, batch->superClustersControl->scMeta,
			nScs, batch->forceEnergyInterims->pme, forceAcc);
	if (nScs > 0) {
		const auto* scData = batch->superClustersControl->scData;
		const Float3 boxSizeInv = boxSize.Inv();
		// Blocksize must be 64, 2 superclusters per block, see the kernel
		NbNonlocalKernel<BoundaryCondition, emvariant, logData><<<(nScs + 1) / 2, 64, 0, cudaStreams[0]>>>(scData, batch->quarterEntryTasksDevice.Get(),
			batch->quarterEntriesDevice.Get(), forceAcc, batch->scResultsDevice.Get(), batch->quarterEntryResultIndicesDevice.Get(),
			batch->superClustersControl->scMeta, boxSize, boxSizeInv, batch->ewaldKappa, batch->params.cutoff_nm * batch->params.cutoff_nm, nScs);
	}
	if (!batch->params.snf_select.empty()) SnfHandler<BoundaryCondition, emvariant>(cudaStreams[2], forceAcc);
	if (batch->nBondgroups > 0) {
		const CompactBondGroupsDevice bondGroups{
			batch->bondgroupDescriptors.Get(), batch->bondgroupSinglebonds.Get(), batch->bondgroupPairbonds.Get(),
			batch->bondgroupAnglebonds.Get(), batch->bondgroupDihedralbonds.Get(), batch->bondgroupImproperdihedralbonds.Get()
		};
		const Float3 boxSizeInv = boxSize.Inv();
		BondgroupsKernel<BoundaryCondition, emvariant, logData><<<(batch->nBondgroups + BONDGROUPS_PER_BLOCK - 1) / BONDGROUPS_PER_BLOCK, THREADS_PER_BONDSGROUPSKERNEL, 0, cudaStreams[4]>>>(
			bondGroups, batch->nBondgroups, batch->boxState, batch->forceEnergyInterims->forceEnergiesBondgroups, primaryBondForces, extraBondForces, batch->bondgroupExtraResultIndices.Get(), batch->superClustersControl->scData,
			batch->bondgroupParticleSlots.Get(), boxSize, boxSizeInv);
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
			// Without live edits, which may change the force after it is measured, forcePrev is exactly the measured force
			batch->forceMagnitudesInStates = !liveEdit.Any();
			SuperclusterIntegrateKernel<BoundaryCondition, logData><<<nBlocks, dim3(16, 4, 1), 0, cudaStreams[0]>>>(
				forceAcc, primaryBondForces, extraBondForces, batch->slotExtraBondReferences.Get(), batch->superClustersControl->scData, batch->integrationStates.Get(), batch->superClustersControl->scMeta,
				batch->integrationSimulationDataDevice.Get(), nScs, boxSize, batch->forcesMagnitudeSquareDevice.Get(),
				liveEdit, log);
		}
		// MD queues the next step without waiting, the host reads nothing from this one. EM reads emStatesHost in hostMaster
		if constexpr (emvariant)
			cudaStreamSynchronize(cudaStreams[0]);
	}
	LIMA_UTILS::genericErrorCheckNoSync("Error during batch timestep");
}

void Engine::deviceMaster() {
	const bool logData = batch->params.data_logging_interval != 0 && batch->step % batch->params.data_logging_interval == 0;
	switch (batch->params.bc_select) {
	case NoBC:
		if (batch->params.em_variant) {
			if (logData) {
				_deviceMaster<NoBoundaryCondition, true, true>();
			}
			else {
				_deviceMaster<NoBoundaryCondition, true, false>();
			}
		}
		else {
			if (logData) {
				_deviceMaster<NoBoundaryCondition, false, true>();
			}
			else {
				_deviceMaster<NoBoundaryCondition, false, false>();
			}
		}
		break;
	case PBC:
		if (batch->params.em_variant) {
			if (logData) {
				_deviceMaster<PeriodicBoundaryCondition, true, true>();
			}
			else {
				_deviceMaster<PeriodicBoundaryCondition, true, false>();
			}
		}
		else {
			if (logData) {
				_deviceMaster<PeriodicBoundaryCondition, false, true>();
			}
			else {
				_deviceMaster<PeriodicBoundaryCondition, false, false>();
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
				sim.device.pclusters.offset, count, NodeIndex(batch->boxSize).toFloat3());
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
				TestSortKernel32<1, 64> << <1, 32 >> > (dKeys, dIds);
			else if (nBins == 4)
				TestSortKernel32<4, 16> << <1, 32 >> > (dKeys, dIds);
			else if (nBins == 8)
				TestSortKernel32<8, 8> << <1, 32 >> > (dKeys, dIds);


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
	__global__ void DeviceFailureKernel() { asm("trap;"); }
}

void Engine::TestLimit(EngineLimitProbe probe, int count) {
	using namespace EngineLimitTesting;
	switch (probe) {
	case EngineLimitProbe::ClusterTransfer:
		Require(count >= 1 && count <= 9, "Invalid transfer fixture size");
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
