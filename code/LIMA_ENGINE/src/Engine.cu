#include "Engine.cuh"

#include "BoundaryCondition.cuh"
#include "EngineBodies.cuh"
#include "EngineKernels.cuh"
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

#include <random>
#include <numeric>

EngineBatchData::~EngineBatchData() {
	pmeController.reset();
	boxState.FreeMembers();
	cudaFree(adamState);
	if (forceEnergyInterims) forceEnergyInterims->Free();
	if (superClustersControl) superClustersControl->Free();
	if (pclusterTransfermodule) pclusterTransfermodule->Free();
	if (superclusterStagingControl) superclusterStagingControl->Free();
}

Engine::Engine(const std::vector<Simulation*>& simulations, EngineRunMode mode)
	: mode(mode), batch(std::make_unique<EngineBatchData>())
{
	if (mode == EngineRunMode::Interactive && simulations.size() != 1)
		throw std::invalid_argument("Interactive engines require one simulation");
	EngineBatch::Pack(*batch, simulations);
	if (mode == EngineRunMode::Interactive) {
		batch->simulations.front().device.active = true;
		batch->simulations.front().runstatus.simulation_finished = false;
	}
	verifyEngine();
	try {
		for (auto& stream : cudaStreams) cudaStreamCreate(&stream);
		cudaStreamCreate(&pmeStream);
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
		batch.reset();
		for (auto stream : cudaStreams) if (stream) cudaStreamDestroy(stream);
		if (pmeStream) cudaStreamDestroy(pmeStream);
		throw;
	}
}

Engine::~Engine() {
	Synchronize();
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
		newSim.stepAtLastEarlystopCheck = oldSim.stepAtLastEarlystopCheck;
		newSim.nLogEntriesTransferred = oldSim.nLogEntriesTransferred;
		newSim.finalized = oldSim.finalized;
		newSim.finalForcesMagnitudeSquared = oldSim.finalForcesMagnitudeSquared;
		newSim.device.thermostatScalar = oldSim.device.thermostatScalar;
		if (!active[id]) continue;
		cudaMemcpy(rebuilt->adamState + newSim.device.pclusters.offset * PersistentCluster::maxParticles,
			batch->adamState + oldSim.device.pclusters.offset * PersistentCluster::maxParticles,
			sizeof(AdamState) * newSim.device.pclusters.count * PersistentCluster::maxParticles, cudaMemcpyDeviceToDevice);
	}
	CopyParticleBuffer(batch->fixedParticleMovementBuffer, rebuilt->fixedParticleMovementBuffer, *batch, *rebuilt);
	CopyParticleBuffer(batch->fixedParticleRotationBuffer, rebuilt->fixedParticleRotationBuffer, *batch, *rebuilt);
	CopyParticleBuffer(batch->forceMaskBuffer, rebuilt->forceMaskBuffer, *batch, *rebuilt);
	CopyParticleBuffer(batch->elasticPositionsBuffer, rebuilt->elasticPositionsBuffer, *batch, *rebuilt);

	auto oldBatch = std::move(batch);
	batch = std::move(rebuilt);
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

void Engine::step() {
	if (IsFinished()) return;
	LIMA_UTILS::genericErrorCheckNoSync("Error before step");
	if (mode == EngineRunMode::Interactive) {
		auto& sim = batch->simulations.front();
		// Live editing changes EM/MD and force selections between steps.
		batch->params = sim.simulation->simParams;
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
	const bool rebuilt = hostMaster();
	if (!rebuilt && !IsFinished() && batch->step % batch->params.stepsPerNlistupdate == 0) {
		batch->superClustersControl->Reset(batch->nGridnodes, cudaStreams[0]);
		RunClustering(cudaStreams[0]);
		MakeSuperClusterTasksGPU(cudaStreams[0]);
	}
	LIMA_UTILS::genericErrorCheckNoSync("Error after step");
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
	sim.runstatus.stepForMostRecentData = sim.step;
	sim.runstatus.most_recent_positions = sim.simulation->traj_buffer->getBufferAtIndex(entries - 1);
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
	sim.runstatus.most_recent_positions = sim.simulation->traj_buffer->getBufferAtIndex(0);
}

void Engine::HandleEarlyStoppingInEM(EngineSimulationData& sim) {
	if (!batch->params.em_variant || sim.step == sim.simulation->simParams.n_steps) return;
	if (sim.step <= sim.stepAtLastEarlystopCheck + 100) return;
	Synchronize();
	const auto range = sim.device.particles;
	const auto forces = GenericCopyToHost(batch->forcesMagnitudeSquareDevice.Get() + range.offset, range.count);
	sim.runstatus.greatestForce = std::sqrt(Statistics::Max(forces.data(), forces.size())) / KILO;
	sim.simulation->maxForceBuffer.emplace_back(sim.step, sim.runstatus.greatestForce);
	if (sim.runstatus.greatestForce <= batch->params.em_force_tolerance) sim.runstatus.simulation_finished = true;
	sim.stepAtLastEarlystopCheck = sim.step;
}

template <typename BoundaryCondition, bool emvariant, bool logData>
void Engine::_deviceMaster() {
	const Float3 boxSize = NodeIndex(batch->boxSize).toFloat3();
	const int nScs = batch->nSuperclusters;
	const int nPcs = batch->nPclusters;
	if (ENABLE_ES_LR && batch->params.enable_electrostatics)
		batch->pmeController->CalcCharges(batch->superClustersControl->scData, batch->superClustersControl->scMeta,
			nScs, batch->forceEnergyInterims->pme);
	if (nScs > 0) {
		const auto* scData = batch->superClustersControl->scData;
		const auto* scMeta = batch->superClustersControl->scMeta;
		const Float3 boxSizeInv = boxSize.Inv();
		NbNonlocalKernel<BoundaryCondition, emvariant, logData, true><<<nScs, dim3(16,4,1), 0, cudaStreams[0]>>>(
			scData, batch->scscTasksDevice.Get(), batch->scResultsDevice.Get(), batch->idsOfQuerySuperclustersDevice.Get(),
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
		for (const auto& sim : batch->simulations) {
			if (!sim.device.active) continue;
			const int simScs = sim.superclusters.count;
			if (simScs == 0) continue;
			const auto& device = sim.device;
			const int simParticles = device.pclusters.count * PersistentCluster::maxParticles;
			Float3* fixedMovement = batch->fixedParticleMovementBuffer ? batch->fixedParticleMovementBuffer->Get() : nullptr;
			Float3* forceMask = batch->forceMaskBuffer ? batch->forceMaskBuffer->Get() : nullptr;
			const Rotation* fixedRotation = batch->fixedParticleRotationBuffer ? batch->fixedParticleRotationBuffer->Get() : nullptr;
			Float3* trajectoryLog = batch->dataBuffersDevice->traj_buffer + device.logOffset;
			float* potentialEnergyLog = batch->dataBuffersDevice->potE_buffer + device.logOffset;
			float* velocityLog = batch->dataBuffersDevice->vel_buffer + device.logOffset;
			Float3* forceLog = batch->dataBuffersDevice->forceBuffer + device.logOffset;

			const int nBlocks = (simScs + 4 - 1) / 4;
			const dim3 blockDim(16, 4, 1);

			SuperclusterIntegrateKernel<BoundaryCondition, emvariant, logData><<<nBlocks, blockDim, 0, cudaStreams[0]>>>(
				*batch->forceEnergyInterims, batch->adamState, batch->params.data_logging_interval, batch->scResultsDevice.Get(),
				batch->superClustersControl->scData, batch->superClustersControl->scMeta, batch->pClusterDevice.Get(), batch->pClusterMetaDevice.Get(),
				batch->boxState.pclusterInterimStates, batch->step, device.dt, simParticles, device.thermostatScalar,
				device.particles.offset, device.pclusters.offset, sim.superclusters.offset, simScs, batch->forcesMagnitudeSquareDevice.Get(), boxSize,
				fixedMovement, forceMask, fixedRotation, trajectoryLog, potentialEnergyLog, velocityLog, forceLog);
		}
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
						[](const std::pair<int, int>& a, const std::pair<int, int>& b) {
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
