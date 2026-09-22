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
		UploadSimulationData();
		BootstrapClustering(cudaStreams[0]);
		MakeSuperClusterTasksGPU(cudaStreams[0]);
		RefreshActiveWork();
		for (auto& sim : batch->simulations) {
			BootstrapTrajbufferWithCoords(sim);
			if (!sim.device.active) FinalizeSimulation(sim);
		}
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

void Engine::UploadSimulationData() {
	std::vector<SimulationDeviceData> data;
	for (const auto& sim : batch->simulations) data.push_back(sim.device);
	batch->simulationsDevice.SetData(data);
}

void Engine::RefreshActiveWork() {
	Synchronize();
	std::vector<int> pcs, groups, scs;
	for (const auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		for (int i = 0; i < sim.device.pclusters.count; ++i) pcs.push_back(sim.device.pclusters.offset + i);
		for (int i = 0; i < sim.device.bondgroups.count; ++i) groups.push_back(sim.device.bondgroups.offset + i);
	}
	for (const auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		for (int i = 0; i < sim.superclusters.count; ++i) scs.push_back(sim.superclusters.offset + i);
	}
	batch->activePclusterIds.SetData(pcs);
	batch->activeBondgroupIds.SetData(groups);
	batch->activeSuperclusterIds.SetData(scs);
	batch->nActivePclusters = static_cast<int>(pcs.size());
	batch->nActiveBondgroups = static_cast<int>(groups.size());
	batch->nActiveSuperclusters = static_cast<int>(scs.size());
	// PME uses densely packed slots for active simulations; stable simulation IDs
	// map to these slots, independently of the persistent particle ranges.
	if (ENABLE_ES_LR && batch->params.enable_electrostatics) {
		if (!batch->pmeController)
			batch->pmeController = std::make_unique<PME::Controller>(batch->simulations, batch->params.cutoff_nm, pmeStream);
		batch->pmeController->SetActiveSimulations(batch->simulations);
	}
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
			UploadSimulationData();
		}
	}
	deviceMaster();
	++batch->step;
	for (auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		++sim.simulation->step;
		sim.step = sim.simulation->getStep();
	}
	hostMaster();
	if (!IsFinished() && batch->step % batch->params.stepsPerNlistupdate == 0) {
		batch->superClustersControl->Reset(batch->nGridnodes, cudaStreams[0]);
		RunClustering(cudaStreams[0]);
		MakeSuperClusterTasksGPU(cudaStreams[0]);
		RefreshActiveWork();
	}
	LIMA_UTILS::genericErrorCheckNoSync("Error after step");
}

void Engine::hostMaster() {
	bool retired = false;
	bool thermostatChanged = false;
	const bool measureTemperature = DatabuffersDeviceController::IsBufferFull(batch->step, batch->params.data_logging_interval)
		&& batch->step % batch->params.steps_per_temperature_measurement == 0;
	if (measureTemperature)
		batch->thermostat->ComputeKineticEnergy(batch->boxState.pclusterInterimStates, batch->pClusterMetaDevice.Get(),
			batch->activePclusterIds.Get(), batch->nActivePclusters, cudaStreams[0]);
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
					thermostatChanged = true;
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
	if (retired || thermostatChanged) UploadSimulationData();
	if (retired) RefreshActiveWork();
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
	sim.device.active = false;
	sim.runstatus.simulation_finished = true;
	sim.finalized = true;
}

void Engine::terminateSimulation() {
	for (auto& sim : batch->simulations) FinalizeSimulation(sim);
	UploadSimulationData();
	RefreshActiveWork();
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
	const auto range = batch->simulations.at(simulationId).device.pclusters;
	Synchronize();
	batch->pdataCopyBuffer.Expand(range.count);
	cudaMemcpy(batch->pdataCopyBuffer.Get(), batch->pClusterDevice.Get() + range.offset, sizeof(PersistentCluster) * range.count, cudaMemcpyDeviceToDevice);
	return batch->pdataCopyBuffer;
}

struct SqrtFloat {
	__device__ float operator()(float x) const { return sqrtf(x); }
};

CudaBuffer<float>& Engine::OffloadForcesMagnitudeBuffer(size_t simulationId) {
	const auto range = batch->simulations.at(simulationId).device.particles;
	Synchronize();
	batch->forcesMagnitudeCopyBuffer.Expand(range.count);
	cudaMemcpy(batch->forcesMagnitudeCopyBuffer.Get(), batch->forcesMagnitudeSquareDevice.Get() + range.offset, sizeof(float) * range.count, cudaMemcpyDeviceToDevice);
	thrust::device_ptr<float> begin(batch->forcesMagnitudeCopyBuffer.Get());
	thrust::transform(thrust::device, begin, begin + range.count, begin, SqrtFloat{});
	return batch->forcesMagnitudeCopyBuffer;
}

namespace {
	template<typename T>
	void SetParticleBuffer(std::optional<CudaBuffer<T>>& buffer, const std::vector<T>& data, BatchRange range, int totalParticles) {
		if (data.empty()) return;
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
	sim.device.hasFixedMovement = !movement.empty();
	UploadSimulationData();
}

void Engine::SetFixedParticleRotationBuffer(const std::vector<Rotation>& rotation, size_t simulationId) {
	Synchronize();
	auto& sim = batch->simulations.at(simulationId);
	SetParticleBuffer(batch->fixedParticleRotationBuffer, rotation, sim.device.particles, batch->nParticles);
	sim.device.hasFixedRotation = !rotation.empty();
	UploadSimulationData();
}

void Engine::SetForceMask(const std::vector<Float3>& mask, size_t simulationId) {
	Synchronize();
	auto& sim = batch->simulations.at(simulationId);
	SetParticleBuffer(batch->forceMaskBuffer, mask, sim.device.particles, batch->nParticles);
	sim.device.hasForceMask = !mask.empty();
	UploadSimulationData();
}

void Engine::SetElasticPositions(const std::vector<Float3>& positions, size_t simulationId) {
	Synchronize();
	auto& sim = batch->simulations.at(simulationId);
	SetParticleBuffer(batch->elasticPositionsBuffer, positions, sim.device.particles, batch->nParticles);
	sim.device.hasElasticPositions = !positions.empty();
	UploadSimulationData();
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
	const auto* simulations = batch->simulationsDevice.Get();
	const int nScs = batch->nActiveSuperclusters;
	const int nPcs = batch->nActivePclusters;
	if (ENABLE_ES_LR && batch->params.enable_electrostatics)
		batch->pmeController->CalcCharges(batch->superClustersControl->scData, batch->superClustersControl->scMeta,
			batch->activeSuperclusterIds.Get(), nScs, batch->forceEnergyInterims->pme);
	if (nScs > 0) {
		NbNonlocalKernel<BoundaryCondition, emvariant, logData, true><<<nScs, dim3(16,4,1), 0, cudaStreams[0]>>>(
			batch->superClustersControl->scData, batch->scscTasksDevice.Get(), batch->scResultsDevice.Get(), batch->idsOfQuerySuperclustersDevice.Get(),
			batch->resultIndicesDevice.Get(), batch->noInteractionMatricesDevice.Get(), batch->superClustersControl->scMeta,
			batch->activeSuperclusterIds.Get(), simulations, boxSize, boxSize.Inv(), batch->ewaldKappa);
	}
	if (!batch->params.snf_select.empty()) SnfHandler<BoundaryCondition, emvariant>(cudaStreams[2]);
	if (batch->nActiveBondgroups > 0) {
		BondgroupsKernel<BoundaryCondition, emvariant><<<batch->nActiveBondgroups, THREADS_PER_BONDSGROUPSKERNEL, 0, cudaStreams[4]>>>(
			BondGroupsDevice{ batch->bondgroupDescriptors.Get(), batch->bondgroupParticles.Get(), batch->bondgroupSinglebonds.Get(), batch->bondgroupPairbonds.Get(),
				batch->bondgroupAnglebonds.Get(), batch->bondgroupDihedralbonds.Get(), batch->bondgroupImproperdihedralbonds.Get() },
			batch->boxState, batch->forceEnergyInterims->forceEnergiesBondgroups, batch->pClusterDevice.Get(), boxSize, boxSize.Inv(), batch->activeBondgroupIds.Get());
		PclusterBondgroupsGather<<<(nPcs + 31) / 32, 32, 0, cudaStreams[4]>>>(
			batch->pClusterMetaDevice.Get(), nPcs, *batch->forceEnergyInterims, batch->activePclusterIds.Get());
	}
	Synchronize();
	if (nScs > 0) {
		SuperclusterIntegrateKernel<BoundaryCondition, emvariant, logData><<<(nScs + 3) / 4, dim3(16,4,1), 0, cudaStreams[0]>>>(
			*batch->forceEnergyInterims, batch->adamState, batch->params.data_logging_interval, batch->scResultsDevice.Get(),
			batch->superClustersControl->scData, batch->superClustersControl->scMeta, batch->pClusterDevice.Get(), batch->pClusterMetaDevice.Get(),
			batch->boxState.pclusterInterimStates, simulations, batch->step, batch->activeSuperclusterIds.Get(), nScs, batch->forcesMagnitudeSquareDevice.Get(), boxSize,
			batch->fixedParticleMovementBuffer ? batch->fixedParticleMovementBuffer->Get() : nullptr,
			batch->forceMaskBuffer ? batch->forceMaskBuffer->Get() : nullptr,
			batch->fixedParticleRotationBuffer ? batch->fixedParticleRotationBuffer->Get() : nullptr,
			batch->dataBuffersDevice->traj_buffer, batch->dataBuffersDevice->potE_buffer, batch->dataBuffersDevice->vel_buffer, batch->dataBuffersDevice->forceBuffer);
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
	const int nPcs = batch->nActivePclusters;
	if (nPcs == 0) return;
	if (batch->params.snf_select.contains(HorizontalChargeField)) {
		PclusterSnfKernel<BoundaryCondition, emvariant><<<(nPcs + 31) / 32, 32, 0, stream>>>(
			batch->pClusterDevice.Get(), batch->pClusterMetaDevice.Get(), batch->simulationsDevice.Get(), batch->pclusterSimulationIds.Get(),
			batch->forceEnergyInterims->snf, nPcs, batch->activePclusterIds.Get());
	}
	if (batch->params.snf_select.contains(ElasticPosition) && batch->elasticPositionsBuffer) {
		ElasticPositionsForceKernel<<<(nPcs + 31) / 32, 32, 0, stream>>>(
			batch->pClusterDevice.Get(), batch->pClusterMetaDevice.Get(), batch->elasticPositionsBuffer->Get(), batch->forceEnergyInterims->snf,
			nPcs, NodeIndex(batch->boxSize).toFloat3(), batch->activePclusterIds.Get(), batch->simulationsDevice.Get(), batch->pclusterSimulationIds.Get());
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
