#include "Engine.cuh"

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
#include "DebugUtils.h"
#include "EngineHostside.h"
#include "ParticleClusters.cuh"
#include "SuperclusterTaskBuilder.cuh"
#include "BatchData.cuh"
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
		throw;
	}
}

Engine::~Engine() {
	Synchronize();
	StopRenderDataPipes();
	batch.reset(); // Controllers reference streams, so destroy them first.
	for (auto stream : cudaStreams) cudaStreamDestroy(stream);
	cudaStreamDestroy(pmeStream);
	LIMA_UTILS::genericErrorCheckNoSync("Error during Engine destruction");
}

void Engine::Synchronize() {
	cudaStreamSynchronize(pmeStream);
	for (auto stream : cudaStreams) cudaStreamSynchronize(stream);
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
	const bool measureTemperature = DatabuffersDeviceController::IsBufferFull(batch->step, batch->params.data_logging_interval)
		&& batch->step % batch->params.steps_per_temperature_measurement == 0;
	if (measureTemperature)
		batch->thermostat->ComputeKineticEnergy(batch->boxState.pclusterInterimStates, batch->pClusterMetaDevice.Get(),
			batch->nPclusters, cudaStreams[0]);
	for (auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		const auto step = sim.step;
		const auto& params = sim.simulation->simParams;
		if (DatabuffersDeviceController::IsBufferFull(step, params.data_logging_interval)) {
			OffloadLoggingData(sim);
			// Preserve the existing measurement cadence during the storage migration.
			if (measureTemperature) {
				auto [temperature, scalar] = batch->thermostat->Temperature(
					sim.simulation->box->boxparams, params, sim.device.pclusters, cudaStreams[0]);
				sim.simulation->temperature_buffer.push_back(temperature);
				sim.runstatus.current_temperature = temperature;
				if (params.apply_thermostat) {
					sim.device.thermostatScalar = scalar;
				}
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
		RebuildActiveBatch();
		return true;
	}
	return false;
}

void Engine::FinalizeSimulation(EngineSimulationData& sim) {
	if (sim.finalized) return;
	Synchronize();
	OffloadLoggingData(sim);
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

void Engine::OffloadLoggingData(EngineSimulationData& sim) {
	const int interval = batch->params.data_logging_interval;
	if (interval == 0 || sim.step == 0) return;
	const size_t entries = (sim.step - 1) / interval + 1;
	const size_t count = entries - sim.nLogEntriesTransferred;
	if (count == 0) return;
	if (count > DatabuffersDeviceController::nStepsInBuffer) throw std::runtime_error("Logging buffer was not drained");
	Synchronize();
	const size_t width = sim.device.pclusters.count * PersistentCluster::maxParticles;
	// A full ring or its final partial segment needs only four contiguous copies.
	for (size_t entry = sim.nLogEntriesTransferred; entry < entries;) {
		const size_t slot = entry % DatabuffersDeviceController::nStepsInBuffer;
		const size_t chunk = std::min(entries - entry, DatabuffersDeviceController::nStepsInBuffer - slot);
		const size_t src = sim.device.logOffset + slot * width;
		cudaMemcpyAsync(sim.simulation->traj_buffer->getBufferAtIndex(entry), batch->dataBuffersDevice->traj_buffer + src, sizeof(Float3) * width * chunk, cudaMemcpyDeviceToHost, cudaStreams[0]);
		cudaMemcpyAsync(sim.simulation->potE_buffer->getBufferAtIndex(entry), batch->dataBuffersDevice->potE_buffer + src, sizeof(float) * width * chunk, cudaMemcpyDeviceToHost, cudaStreams[0]);
		cudaMemcpyAsync(sim.simulation->vel_buffer->getBufferAtIndex(entry), batch->dataBuffersDevice->vel_buffer + src, sizeof(float) * width * chunk, cudaMemcpyDeviceToHost, cudaStreams[0]);
		cudaMemcpyAsync(sim.simulation->forceBuffer->getBufferAtIndex(entry), batch->dataBuffersDevice->forceBuffer + src, sizeof(Float3) * width * chunk, cudaMemcpyDeviceToHost, cudaStreams[0]);
		entry += chunk;
	}
	cudaStreamSynchronize(cudaStreams[0]);
	sim.nLogEntriesTransferred = entries;
}

CudaBuffer<PersistentCluster>& Engine::OffloadPclusterState(size_t simulationId) {
	const auto& sim = batch->simulations.at(simulationId);
	const auto range = sim.device.pclusters;
	Synchronize();
	const size_t count = sim.finalized ? sim.simulation->box->persistentClusters.size() : range.count;
	batch->pdataCopyBuffer.Expand(count);
	if (sim.finalized) cudaMemcpy(batch->pdataCopyBuffer.Get(), sim.simulation->box->persistentClusters.data(), sizeof(PersistentCluster) * count, cudaMemcpyHostToDevice);
	else cudaMemcpy(batch->pdataCopyBuffer.Get(), batch->pClusterDevice.Get() + range.offset, sizeof(PersistentCluster) * count, cudaMemcpyDeviceToDevice);
	return batch->pdataCopyBuffer;
}

struct SqrtFloat {
	__device__ float operator()(float x) const { return sqrtf(x); }
};

CudaBuffer<float>& Engine::OffloadForcesMagnitudeBuffer(size_t simulationId) {
	const auto& sim = batch->simulations.at(simulationId);
	const auto range = sim.device.particles;
	Synchronize();
	const size_t count = sim.finalized ? sim.finalForcesMagnitudeSquared.size() : range.count;
	batch->forcesMagnitudeCopyBuffer.Expand(count);
	if (sim.finalized) cudaMemcpy(batch->forcesMagnitudeCopyBuffer.Get(), sim.finalForcesMagnitudeSquared.data(), sizeof(float) * count, cudaMemcpyHostToDevice);
	else cudaMemcpy(batch->forcesMagnitudeCopyBuffer.Get(), batch->forcesMagnitudeSquareDevice.Get() + range.offset, sizeof(float) * count, cudaMemcpyDeviceToDevice);
	thrust::device_ptr<float> begin(batch->forcesMagnitudeCopyBuffer.Get());
	thrust::transform(thrust::device, begin, begin + count, begin, SqrtFloat{});
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
	const Float3 boxSize = NodeIndex(batch->boxSize).toFloat3();
	const int nScs = batch->nSuperclusters;
	const int nPcs = batch->nPclusters;
	if (ENABLE_ES_LR && batch->params.enable_electrostatics)
		batch->pmeController->CalcCharges(batch->superClustersControl->scData, batch->superClustersControl->scMeta,
			nScs, batch->forceEnergyInterims->pme);
	NbForceAccumulator nbForceAcc{};
	if constexpr (emvariant) {
		batch->scResultsDevice.Expand(batch->nResults, 1.2); // Noop once allocated. Not allocated at all in MD, where it would be ~1GB for large systems
	}
	else {
		const size_t n = size_t(nScs) * SuperCluster::maxParticles;
		unsigned long long* const base = batch->nbForceAccumulatorDevice.Get();
		nbForceAcc = NbForceAccumulator{ base, base + n, base + 2 * n, base + 3 * n };
		if (nScs > 0)
			cudaMemsetAsync(base, 0, sizeof(unsigned long long) * n * (logData ? 4 : 3), cudaStreams[0]);
	}
	if (nScs > 0) {
		const auto* scData = batch->superClustersControl->scData;
		const auto* scMeta = batch->superClustersControl->scMeta;
		const Float3 boxSizeInv = boxSize.Inv();
		// Blocksize must not exceed 64 threads, see __launch_bounds__ on the kernel
		NbNonlocalKernel<BoundaryCondition, emvariant, logData, true><<<nScs, dim3(16,4,1), 0, cudaStreams[0]>>>(
			scData, batch->scscTasksDevice.Get(), batch->scResultsDevice.Get(), nbForceAcc, batch->idsOfQuerySuperclustersDevice.Get(),
			batch->resultIndicesDevice.Get(), batch->noInteractionMatricesDevice.Get(), scMeta, boxSize, boxSizeInv, batch->ewaldKappa);
	}
	if (!batch->params.snf_select.empty()) SnfHandler<BoundaryCondition, emvariant>(cudaStreams[2]);
	if (batch->nBondgroups > 0) {
		const BondGroupsDevice bondGroups{
			batch->bondgroupDescriptors.Get(), batch->bondgroupParticles.Get(), batch->bondgroupSinglebonds.Get(), batch->bondgroupPairbonds.Get(),
			batch->bondgroupAnglebonds.Get(), batch->bondgroupDihedralbonds.Get(), batch->bondgroupImproperdihedralbonds.Get()
		};
		const Float3 boxSizeInv = boxSize.Inv();
		BondgroupsKernel<BoundaryCondition, emvariant><<<batch->nBondgroups, THREADS_PER_BONDSGROUPSKERNEL, 0, cudaStreams[4]>>>(
			bondGroups, batch->boxState, batch->forceEnergyInterims->forceEnergiesBondgroups, batch->pClusterDevice.Get(), boxSize, boxSizeInv);
		PclusterBondgroupsGather<<<(nPcs + 31) / 32, 32, 0, cudaStreams[4]>>>(
			batch->pClusterMetaDevice.Get(), nPcs, *batch->forceEnergyInterims);
	}
	Synchronize();
	if (nScs > 0) {
		std::vector<IntegrationSimulationData> simulationData;
		simulationData.reserve(batch->simulations.size());
		for (const auto& sim : batch->simulations) {
			simulationData.push_back({
				sim.device.particles, sim.device.pclusters, sim.device.logOffset, sim.device.dt, sim.device.thermostatScalar});
		}
		batch->integrationSimulationDataDevice.Expand(simulationData.size());
		cudaMemcpyAsync(batch->integrationSimulationDataDevice.Get(), simulationData.data(),
			sizeof(IntegrationSimulationData) * simulationData.size(), cudaMemcpyHostToDevice, cudaStreams[0]);

		const int nBlocks = (nScs + 4 - 1) / 4;
		SuperclusterIntegrateKernel<BoundaryCondition, emvariant, logData><<<nBlocks, dim3(16, 4, 1), 0, cudaStreams[0]>>>(
			*batch->forceEnergyInterims, batch->emForces.Get(), batch->params.data_logging_interval, batch->scResultsDevice.Get(), nbForceAcc,
			batch->superClustersControl->scData, batch->superClustersControl->scMeta, batch->pClusterDevice.Get(), batch->pClusterMetaDevice.Get(),
			batch->boxState.pclusterInterimStates, batch->step, batch->integrationSimulationDataDevice.Get(), nScs,
			batch->forcesMagnitudeSquareDevice.Get(), boxSize, batch->fixedParticleMovementBuffer ? batch->fixedParticleMovementBuffer->Get() : nullptr,
			batch->forceMaskBuffer ? batch->forceMaskBuffer->Get() : nullptr,
			batch->fixedParticleRotationBuffer ? batch->fixedParticleRotationBuffer->Get() : nullptr,
			batch->dataBuffersDevice->traj_buffer, batch->dataBuffersDevice->potE_buffer, batch->dataBuffersDevice->vel_buffer,
			batch->dataBuffersDevice->forceBuffer);
		if constexpr (emvariant)
			UpdateEnergyMinimization<BoundaryCondition>(boxSize);
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
void Engine::SnfHandler(cudaStream_t& stream) {
	const int nPcs = batch->nPclusters;
	if (nPcs == 0) return;
	for (const auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		const int count = sim.device.pclusters.count;
		if (batch->params.snf_select.contains(HorizontalChargeField)) {
			PclusterSnfKernel<BoundaryCondition, emvariant><<<(count + 31) / 32, 32, 0, stream>>>(
				batch->pClusterDevice.Get(), batch->pClusterMetaDevice.Get(), sim.device.uniformElectricField,
				batch->forceEnergyInterims->snf, sim.device.pclusters.offset, count);
		}
		if (batch->params.snf_select.contains(ElasticPosition) && batch->elasticPositionsBuffer) {
			ElasticPositionsForceKernel<<<(count + 31) / 32, 32, 0, stream>>>(
				batch->pClusterDevice.Get(), batch->pClusterMetaDevice.Get(), batch->elasticPositionsBuffer->Get(), batch->forceEnergyInterims->snf,
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
