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
#include "SuperClusterTaskBuilder.cuh"

#include <random>

EngineSimulationData::EngineSimulationData(Simulation* simulation)
	: simulation(simulation)
	, forceEnergyInterims(std::make_unique<ForceEnergyInterims>(
		simulation->box->bondgroups.particles.size(),
		simulation->box->boxparams.totalParticles,
		simulation->box->persistentClusters.size()))
{
}

EngineSimulationData::~EngineSimulationData() = default;

Engine::Engine(Simulation* _sim, BoundaryConditionSelect bc)
	: simData(std::make_unique<EngineSimulationData>(_sim))
	, bc_select(bc)
{
	simData->ewaldKappa = PhysicsUtils::CalcEwaldkappa(simData->simulation->simParams.cutoff_nm);
	for (cudaStream_t& stream : cudaStreams)
		cudaStreamCreate(&stream);
	cudaStreamCreate(&pmeStream);

    verifyEngine();

	const BoxParams& boxparams = simData->simulation->box->boxparams;

	simData->dataBuffersDevice = std::make_unique<DatabuffersDeviceController>(simData->simulation->box->persistentClusters.size(), simData->simulation->simParams.data_logging_interval);

	simData->superClustersControl = std::make_unique<SuperClustersControl>(boxparams.boxSize, simData->simulation->box->persistentClusters.size());
	simData->pclusterTransfermodule = std::make_unique<PClusterTransfermodule>(PClusterTransfermodule::Create(boxparams.boxSize));
	simData->pClusterDevice.SetData(simData->simulation->box->persistentClusters);
	simData->pClusterMetaDevice.SetData(simData->simulation->box->persistentClustersMetadata);

	simData->forcesMagnitudeSquareDevice.Expand(boxparams.totalParticles);



	simData->boxState = BoxState::Create(*simData->simulation->box);
	cudaMalloc(&simData->adamState, sizeof(AdamState) * simData->simulation->box->persistentClusters.size() * PersistentCluster::maxParticles);
	cudaMemset(simData->adamState, 0, sizeof(AdamState) * simData->simulation->box->persistentClusters.size() * PersistentCluster::maxParticles);
	/*
	// Precomputed LUT initialization is disabled with the LUT declarations. Keep this for potential reuse.
	const float cutoffNM = simData->simulation->simParams.cutoff_nm;
	cudaMemcpyToSymbol(DeviceConstants::bsplineTable, PrecomputeBsplineTable().data(), sizeof(float) * PrecomputeBsplineTable().size(), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(DeviceConstants::erfcForcescalarTable, PrecomputeErfcForcescalarTable(cutoffNM).data(), sizeof(float) * PrecomputeErfcForcescalarTable(cutoffNM).size(), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(DeviceConstants::erfcPotentialscalarTable, PrecomputeErfcPotentialscalarTable(cutoffNM).data(), sizeof(float) * PrecomputeErfcPotentialscalarTable(cutoffNM).size(), 0, cudaMemcpyHostToDevice);
	*/
	BootstrapClustering(cudaStreams[0]);
	MakeSuperClusterTasksGPU(cudaStreams[0]);

	simData->pmeController = std::make_unique<PME::Controller>(*simData->simulation->box, simData->simulation->simParams.cutoff_nm, pmeStream);

	simData->bondgroupDescriptors.SetData(simData->simulation->box->bondgroups.groups);
	simData->bondgroupParticles.SetData(simData->simulation->box->bondgroups.particles);
	simData->bondgroupSinglebonds.SetData(simData->simulation->box->bondgroups.singlebonds);
	simData->bondgroupPairbonds.SetData(simData->simulation->box->bondgroups.pairbonds);
	simData->bondgroupAnglebonds.SetData(simData->simulation->box->bondgroups.anglebonds);
	simData->bondgroupDihedralbonds.SetData(simData->simulation->box->bondgroups.dihedralbonds);
	simData->bondgroupImproperdihedralbonds.SetData(simData->simulation->box->bondgroups.improperdihedralbonds);

	simData->thermostat = std::make_unique<Thermostat>(simData->simulation->box->persistentClusters.size());

	// To create the NLists we need to bootstrap the traj_buffer, since it has no data yet
	bootstrapTrajbufferWithCoords();
}

Engine::~Engine() {
	Synchronize();
	simData->pmeController.reset();
	simData->boxState.FreeMembers();
	if (simData->adamState != nullptr)
		cudaFree(simData->adamState);
	simData->forceEnergyInterims->Free();

	for (cudaStream_t& stream : cudaStreams) {
		cudaStreamDestroy(stream);
	}
	cudaStreamDestroy(pmeStream);

	if (simData->superClustersControl)
		simData->superClustersControl->Free();

	LIMA_UTILS::genericErrorCheckNoSync("Error during Engine destruction");
}

void Engine::Synchronize() {
	cudaStreamSynchronize(pmeStream);
	for (cudaStream_t stream : cudaStreams)
		cudaStreamSynchronize(stream);
}


void Engine::step() {
	LIMA_UTILS::genericErrorCheckNoSync("Error before step!");

	deviceMaster();	// Device first, otherwise offloading data always needs the last datapoint!
	simData->simulation->step++;

	hostMaster();
	

	if (simData->simulation->step % simData->simulation->simParams.stepsPerNlistupdate == 0) {
		simData->superClustersControl->Reset(simData->simulation->box->boxparams.boxSize, cudaStreams[0]);
		RunClustering(cudaStreams[0]);
		MakeSuperClusterTasksGPU(cudaStreams[0]);
	}

	LIMA_UTILS::genericErrorCheckNoSync("Error after step!");
}

void Engine::hostMaster() {						// This is and MUST ALWAYS be called after the deviceMaster, and AFTER incStep()!
	auto t0 = std::chrono::high_resolution_clock::now();
	if (DatabuffersDeviceController::IsBufferFull(simData->simulation->getStep(), simData->simulation->simParams.data_logging_interval)) {
		offloadLoggingData(DatabuffersDeviceController::nStepsInBuffer);
		runstatus.stepForMostRecentData = simData->simulation->getStep();

		if ((simData->simulation->getStep() % simData->simulation->simParams.steps_per_temperature_measurement) == 0 && simData->simulation->getStep() > 0) {
			auto [temperature, newThermostatScalar] = simData->thermostat->Temperature(simData->boxState.pclusterInterimStates, simData->simulation->box->boxparams,
				simData->simulation->simParams, simData->simulation->getStep(), simData->pClusterMetaDevice.Get(), cudaStreams[0]);
			simData->simulation->temperature_buffer.push_back(temperature);
			runstatus.current_temperature = temperature;

			if (simData->simulation->simParams.apply_thermostat)
				simData->thermostatScalar = newThermostatScalar;
		}		
	}
	HandleEarlyStoppingInEM();
	/*if (simData->simulation->getStep() % simData->simulation->simParams.stepsPerNlistupdate == simData->simulation->simParams.stepsPerNlistupdate-1)
		nlistController->UpdateNlist(sim_dev, simData->simulation->box->boxparams, simData->simulation->simParams.bc_select, cudaStreams);*/

	// Handle status
	runstatus.current_step = simData->simulation->getStep();
	if (runstatus.current_step >= simData->simulation->simParams.n_steps || runstatus.critical_error_occured)
		runstatus.simulation_finished = true;


	const auto t1 = std::chrono::high_resolution_clock::now();
	const int cpu_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
}

void Engine::terminateSimulation() {
	const int64_t stepsReadyToTransfer = DatabuffersDeviceController::StepsReadyToTransfer(simData->simulation->getStep(), simData->simulation->simParams.data_logging_interval);
	offloadLoggingData(stepsReadyToTransfer);

	simData->boxState.CopyDataToHost(*simData->simulation->box);

	Synchronize();
	LIMA_UTILS::genericErrorCheckNoSync("Error during TerminateSimulation");
}

//--------------------------------------------------------------------------	CPU workload --------------------------------------------------------------//

void Engine::offloadLoggingData(const int64_t steps_to_transfer) {
	assert(steps_to_transfer <= simData->simulation->getStep());
	if (steps_to_transfer == 0) { return; }

	cudaStreamSynchronize(cudaStreams[0]);

	const int64_t startstep = simData->simulation->getStep() - steps_to_transfer * simData->simulation->simParams.data_logging_interval;
	const int64_t startindex = LIMALOGSYSTEM::getMostRecentDataentryIndex(startstep, simData->simulation->simParams.data_logging_interval);
	const int64_t indices_to_transfer = LIMALOGSYSTEM::getNIndicesBetweenSteps(startstep, simData->simulation->getStep(), simData->simulation->simParams.data_logging_interval);
	//const int particlesUpperbound = simData->simulation->box->boxparams.total_particles_upperbound;
	const int nParticlesUpperbound = simData->simulation->box->persistentClusters.size() * PersistentCluster::maxParticles;
	
	cudaMemcpyAsync(
		simData->simulation->potE_buffer->getBufferAtIndex(startindex),
		simData->dataBuffersDevice->potE_buffer,
		sizeof(float) * nParticlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost, cudaStreams[0]);
	
	cudaMemcpyAsync(
		simData->simulation->vel_buffer->getBufferAtIndex(startindex),
		simData->dataBuffersDevice->vel_buffer,
		sizeof(float) * nParticlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost, cudaStreams[0]);

	cudaMemcpyAsync(
		simData->simulation->forceBuffer->getBufferAtIndex(startindex),
		simData->dataBuffersDevice->forceBuffer,
		sizeof(Float3) * nParticlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost, cudaStreams[0]);

	cudaMemcpyAsync(
		simData->simulation->traj_buffer->getBufferAtIndex(startindex),
		simData->dataBuffersDevice->traj_buffer,
		sizeof(Float3) * nParticlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost, cudaStreams[0]);
	cudaStreamSynchronize(cudaStreams[0]);

	simData->step_at_last_traj_transfer = simData->simulation->getStep();
	runstatus.most_recent_positions = simData->simulation->traj_buffer->getBufferAtIndex(LIMALOGSYSTEM::getMostRecentDataentryIndex(simData->simulation->getStep() - 1, simData->simulation->simParams.data_logging_interval));
}

void Engine::offloadTrainData() {
#ifdef GENERATETRAINDATA
	uint64_t values_per_step = N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simData->simulation->boxparams_host.n_compounds;
	if (values_per_step == 0) {
		return;	// No data to transfer
	}

	uint64_t step_offset = (simData->simulation->getStep() - STEPS_PER_TRAINDATATRANSFER) * values_per_step;	// fix max_compound to the actual count save LOTS of space!. Might need a file in simout that specifies cnt for loading in other programs...
	cudaMemcpy(&simData->simulation->trainingdata[step_offset], simData->dataBuffersDevice->data_GAN, sizeof(Float3) * values_per_step * STEPS_PER_TRAINDATATRANSFER, cudaMemcpyDeviceToHost);
	LIMA_UTILS::genericErrorCheckNoSync("Cuda error during traindata offloading\n");
#endif
}

//
CudaBuffer<PersistentCluster>& Engine::OffloadPclusterState() {
	simData->pdataCopyBuffer.Expand(simData->simulation->box->persistentClusters.size());
	cudaMemcpy(simData->pdataCopyBuffer.Get(), simData->pClusterDevice.Get(), sizeof(PersistentCluster) * simData->simulation->box->persistentClusters.size(), cudaMemcpyDeviceToDevice);
	return simData->pdataCopyBuffer;
}


__device__ struct SqrtFloat {
	__device__ float operator()(float x) const
	{
		return sqrtf(x);
	}
};

CudaBuffer<float>& Engine::OffloadForcesMagnitudeBuffer() {
	// Copy to offloading buffer
	const int nParticles = simData->simulation->box->boxparams.totalParticles;
	simData->forcesMagnitudeCopyBuffer.Expand(nParticles);
	cudaMemcpy(simData->forcesMagnitudeCopyBuffer.Get(), simData->forcesMagnitudeSquareDevice.Get(), sizeof(float) * nParticles, cudaMemcpyDeviceToDevice);

	// Apply sqrt 
	thrust::device_ptr<float> begin(simData->forcesMagnitudeCopyBuffer.Get());
	thrust::transform(
		thrust::device,
		begin,
		begin + nParticles,
		begin,
		SqrtFloat{}
	);

	return simData->forcesMagnitudeCopyBuffer;
}

void Engine::SetFixedParticleMovementBuffer(const std::vector<Float3>& movement) {
	if (movement.empty()) {
		simData->fixedParticleMovementBuffer.reset();
		return;
	}
	assert(movement.size() == simData->simulation->box->boxparams.totalParticles);
	if (!simData->fixedParticleMovementBuffer.has_value())
		simData->fixedParticleMovementBuffer.emplace();
	simData->fixedParticleMovementBuffer->SetData(movement);
}
void Engine::SetFixedParticleRotationBuffer(const std::vector<Rotation>& rotation) {
	if (rotation.empty()) {
		simData->fixedParticleRotationBuffer.reset();
		return;
	}
	assert(rotation.size() == simData->simulation->box->boxparams.totalParticles);
	if (!simData->fixedParticleRotationBuffer.has_value())
		simData->fixedParticleRotationBuffer.emplace();
	simData->fixedParticleRotationBuffer->SetData(rotation);
}

void Engine::SetForceMask(const std::vector<Float3>& mask) {
	if (mask.empty()) {
		simData->forceMaskBuffer.reset();
		return;
	}
	assert(mask.size() == simData->simulation->box->boxparams.totalParticles);
	if (!simData->forceMaskBuffer.has_value())
		simData->forceMaskBuffer.emplace();
	simData->forceMaskBuffer->SetData(mask);
}

void Engine::SetElasticPositions(const std::vector<Float3>& positions) {
	if (positions.empty()) {
		simData->elasticPositionsBuffer.reset();
		return;
	}
	assert(positions.size() == simData->simulation->box->boxparams.totalParticles);
	if (!simData->elasticPositionsBuffer.has_value())
		simData->elasticPositionsBuffer.emplace();
	simData->elasticPositionsBuffer->SetData(positions);
}

void Engine::bootstrapTrajbufferWithCoords() {
	if (simData->simulation->simParams.n_steps == 0) return;

	for (int pcid = 0; pcid < simData->simulation->box->persistentClusters.size(); pcid++) {
		const PersistentCluster& pc = simData->simulation->box->persistentClusters[pcid];
		for (int pid = 0; pid < 4; pid++) {
			simData->simulation->traj_buffer->GetDatapoint(pcid, pid, 0) = pc.pqd[pid].position;
		}
	}

	simData->step_at_last_traj_transfer = 0.f;
	runstatus.most_recent_positions = simData->simulation->traj_buffer->getBufferAtIndex(0);

	LIMA_UTILS::genericErrorCheck(cudaStreams[0], "Error during bootstrapTrajbufferWithCoords");
}

void Engine::HandleEarlyStoppingInEM() {
	if (!simData->simulation->simParams.em_variant || simData->simulation->getStep() == simData->simulation->simParams.n_steps)
		return;
	
	const int minStepsPerCheck = 100;
	if (simData->simulation->getStep() > simData->stepAtLastEarlystopCheck + minStepsPerCheck) {
		auto forceMagSquared = simData->forcesMagnitudeSquareDevice.GetData(); // [(J/mol/nm)^2]
		const float greatestForce = std::sqrt(Statistics::Max(forceMagSquared.data(), forceMagSquared.size()));
		//const float greatestForce = Statistics::MaxLen(simData->simulation->forceBuffer->GetBufferAtStep(simData->simulation->getStep()-1), simData->simulation->forceBuffer->EntriesPerStep());
		runstatus.greatestForce = greatestForce / KILO; // Convert [J/mol/nm] to [kJ/mol/nm]
		simData->simulation->maxForceBuffer.emplace_back(std::pair<int64_t,float>{ simData->simulation->getStep(), runstatus.greatestForce });

		if (runstatus.greatestForce <= simData->simulation->simParams.em_force_tolerance) {
			runstatus.simulation_finished = true;
		}

		simData->stepAtLastEarlystopCheck = simData->simulation->getStep();
	}
	LIMA_UTILS::genericErrorCheck(cudaStreams[0], "HandleEarlyStoppingInEM");
}




//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//
template <typename BoundaryCondition, bool emvariant, bool logData>
void Engine::_deviceMaster() {
	
	const BoxParams& boxparams = simData->simulation->box->boxparams;
	const int step = simData->simulation->getStep();
	const Float3 boxSize = boxparams.BoxSizeFloat();


	// #### Initial round of force computations
    if (ENABLE_ES_LR && simData->simulation->simParams.enable_electrostatics) {
        simData->pmeController->CalcCharges(simData->superClustersControl->scData, simData->superClustersControl->scMeta, simData->nSuperclusters, simData->forceEnergyInterims->pme, step);
    }


	if (simData->nSuperclusters > 0) {
		const bool useNointeractionMatrix = true;
		dim3 blockDim(SuperCluster::maxParticles, 4, 1);
		NbNonlocalKernel<BoundaryCondition, emvariant, logData, useNointeractionMatrix>
			<<<simData->nSuperclusters, blockDim, 0, cudaStreams[0]>>>
			(simData->superClustersControl->scData, simData->scscTasksDevice.Get(), simData->scResultsDevice.Get(), simData->idsOfQuerySuperclustersDevice.Get(), simData->resultIndicesDevice.Get(),
				simData->noInteractionMatricesDevice.Get(), simData->superClustersControl->scMeta, step, boxSize, boxSize.Inv(), simData->ewaldKappa);
		LIMA_UTILS::genericErrorCheckNoSync("Error after NBNonlocalKernel");

		//nbGatherForceenergy.Expand(simData->nSuperclusters * SuperCluster::nParticles, 1.2);
		//NBGather<<<simData->nSuperclusters, dim3(16,4,1), 0, cudaStreams[0]>>>
		//	(simData->superClustersControl->scMeta, simData->scResultsDevice.Get()/*, nbGatherForceenergy.Get()*/);
		//LIMA_UTILS::genericErrorCheckNoSync("Error after NBGather");
	}
	if (!simData->simulation->simParams.snf_select.empty()) {
		SnfHandler<BoundaryCondition, emvariant>(cudaStreams[2]);
		LIMA_UTILS::genericErrorCheckNoSync("Error after SupernaturalForces");
	}

	if (!simData->simulation->box->bondgroups.empty()) {
		BondgroupsKernel<BoundaryCondition, emvariant> << < simData->simulation->box->bondgroups.size(), THREADS_PER_BONDSGROUPSKERNEL, 0, cudaStreams[4]>>>
			(BondGroupsDevice{ simData->bondgroupDescriptors.Get(), simData->bondgroupParticles.Get(), simData->bondgroupSinglebonds.Get(), simData->bondgroupPairbonds.Get(), simData->bondgroupAnglebonds.Get(), simData->bondgroupDihedralbonds.Get(), simData->bondgroupImproperdihedralbonds.Get() }, simData->boxState, simData->forceEnergyInterims->forceEnergiesBondgroups, simData->pClusterDevice.Get(), boxSize, boxSize.Inv());
		LIMA_UTILS::genericErrorCheckNoSync("Error after BondgroupsKernel");

		// Gather bondgroup ordered forces into particle ordered
		const int nPclusters = simData->simulation->box->persistentClusters.size();
		const int nBlocks = (nPclusters + 31) / 32;
		PclusterBondgroupsGather << <nBlocks, 32, 0, cudaStreams[4] >> >
			(simData->pClusterMetaDevice.Get(), nPclusters, *simData->forceEnergyInterims);
		LIMA_UTILS::genericErrorCheckNoSync("Error after PclusterBondgroupsGather");
	}

	// #### Integration and Transfer kernels
	cudaStreamSynchronize(pmeStream);
	for (int i = 0; i < cudaStreams.size(); i++) {
		cudaStreamSynchronize(cudaStreams[i]);
	}



	const bool updateNlistsAfterThisStep = (simData->simulation->getStep()+1) % simData->simulation->simParams.stepsPerNlistupdate == simData->simulation->simParams.stepsPerNlistupdate-1;


	if (simData->nSuperclusters > 0) {
		int totalParticlesUpperbound = simData->simulation->box->persistentClusters.size() * PersistentCluster::maxParticles;
		const int nBlocks = (simData->nSuperclusters + 4 - 1) / 4;
		const dim3 blockDim(16, 4, 1);

		Float3* fixedParticleMovementBufferPtr = simData->fixedParticleMovementBuffer.has_value() ? simData->fixedParticleMovementBuffer->Get() : nullptr;
		Float3* forcesMaskBufferPtr = simData->forceMaskBuffer.has_value() ? simData->forceMaskBuffer->Get() : nullptr;
		Rotation* fixedParticleRotationBufferPtr = simData->fixedParticleRotationBuffer.has_value() ? simData->fixedParticleRotationBuffer->Get() : nullptr;
		SuperclusterIntegrateKernel<BoundaryCondition, emvariant, logData>
			<<<nBlocks, blockDim, 0, cudaStreams[0]>>>
			(*simData->forceEnergyInterims, simData->adamState, simData->simulation->simParams.data_logging_interval, simData->scResultsDevice.Get(), simData->superClustersControl->scData, simData->superClustersControl->scMeta, simData->pClusterDevice.Get(), simData->pClusterMetaDevice.Get(),
				simData->boxState.pclusterInterimStates, step, simData->simulation->simParams.dt, totalParticlesUpperbound, simData->nSuperclusters, simData->forcesMagnitudeSquareDevice.Get(),
				boxSize, simData->thermostatScalar, fixedParticleMovementBufferPtr, forcesMaskBufferPtr, fixedParticleRotationBufferPtr,
				simData->dataBuffersDevice->traj_buffer, simData->dataBuffersDevice->potE_buffer, simData->dataBuffersDevice->vel_buffer, simData->dataBuffersDevice->forceBuffer);
		LIMA_UTILS::genericErrorCheckNoSync("Error after SuperclusterIntegrateKernel");
		cudaStreamSynchronize(cudaStreams[0]);
	}

	//DebugUtils::VerifyIdentical(simData->superClustersControl->scData, simData->nSuperclusters, "Engine_SCData", step);
	//DebugUtils::VerifyIdentical(simData->pClusterDevice, simData->simulation->box->persistentClusters.size(), "Engine_PClusterDevice", step);
}



void Engine::deviceMaster() {

	const bool logData = simData->simulation->simParams.data_logging_interval != 0 && simData->simulation->getStep() % simData->simulation->simParams.data_logging_interval == 0;// TODO maybe log at the final step, not 0th?

	switch (simData->simulation->simParams.bc_select) {
	case NoBC:
		if (simData->simulation->simParams.em_variant) {
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
		if (simData->simulation->simParams.em_variant) {
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





// This function must not have changing template or normal arguments for it's kernels, or it will break cudaGraph
template <typename BoundaryCondition, bool emvariant>
void Engine::SnfHandler(cudaStream_t& stream) {
	if (simData->simulation->simParams.snf_select.contains(HorizontalSqueeze)) {
		//SupernaturalForces::ApplyHorizontalSqueeze << < simData->simulation->box->boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, stream >> > (sim_dev, simData->simulation->getStep());
		//break;
	}
	if (simData->simulation->simParams.snf_select.contains(HorizontalChargeField))
	{
		const int nPclusters = simData->simulation->box->persistentClusters.size();
		const int nCudablocks = (nPclusters + 31) / 32;
		PclusterSnfKernel<BoundaryCondition, emvariant>
			<<<nCudablocks, 32, 0, stream >> >
			(simData->pClusterDevice.Get(), simData->pClusterMetaDevice.Get(), simData->simulation->box->uniformElectricField, simData->forceEnergyInterims->snf, nPclusters);
	}
	
	if (simData->simulation->simParams.snf_select.contains(SupernaturalForcesSelect::ElasticPosition) && simData->elasticPositionsBuffer.has_value()) {
		const int nPclusters = simData->simulation->box->persistentClusters.size();
		const int nCudablocks = (nPclusters + 31) / 32;
		ElasticPositionsForceKernel << <nCudablocks, 32, 0, stream >> >
			(simData->pClusterDevice.Get(), simData->pClusterMetaDevice.Get(), simData->elasticPositionsBuffer->Get(), simData->forceEnergyInterims->snf, nPclusters, simData->simulation->box->boxparams.BoxSizeFloat());
	}
		
		
	//case BoxEdgePotential:
	//	if (simData->simulation->box->boxparams.n_compounds > 0)
	//		SupernaturalForces::BoxEdgeForceCompounds << < simData->simulation->box->boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, stream >> > (sim_dev, simData->simulation->getStep());
	//	if (simData->simulation->box->boxparams.nTinymols > 0)
	//		SupernaturalForces::BoxEdgeForceSolvents<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(simData->simulation->box->boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, stream>>>(sim_dev, simData->simulation->getStep());
	//	break;
	
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
