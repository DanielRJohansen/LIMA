#include "Engine.cuh"

#include "BoundaryCondition.cuh"
#include "EngineBodies.cuh"
#include "EngineKernels.cuh"
#include "LimaPositionSystem.cuh"
#include "PME.cuh"
#include "SimulationDevice.cuh"
#include "SupernaturalForces.cuh"
#include "Thermostat.cuh"
#include "Statistics.h"
#include "Utilities.h"
#include "DebugUtils.h"
#include "EngineHostside.h"
#include "ParticleClusters.cuh"
#include "SuperClusterTaskBuilder.cuh"

#include <random>

Engine::Engine(Simulation* _sim, BoundaryConditionSelect bc)
	: bc_select(bc)
	, forceEnergyInterims(std::make_unique<ForceEnergyInterims>(_sim->box->bondgroups.size(), _sim->box->boxparams.totalParticles, _sim->box->persistentClusters.size()))
{
	simulation = _sim;

    verifyEngine();

	const BoxParams& boxparams = simulation->box->boxparams;

	dataBuffersDevice = std::make_unique<DatabuffersDeviceController>(simulation->box->persistentClusters.size(), simulation->simParams.data_logging_interval);

	superClustersControl = std::make_unique<SuperClustersControl>(boxparams.boxSize, simulation->box->persistentClusters.size());
	pclusterTransfermodule = std::make_unique<PClusterTransfermodule>(PClusterTransfermodule::Create(boxparams.boxSize));
	pClusterDevice.SetData(simulation->box->persistentClusters);
	pClusterMetaDevice.SetData(simulation->box->persistentClustersMetadata);

	forcesMagnitudeSquareDevice.Expand(boxparams.totalParticles);



	// Create the Sim_dev {
	{
		if (sim_dev != nullptr) { throw std::runtime_error("Expected simdev to be null to move sim to device"); }
		SimulationDevice simdevTemp(simulation->simParams, simulation->box.get(), BoxConfig::Create(*simulation->box), BoxState::Create(*simulation->box), *dataBuffersDevice);
		sim_dev = GenericCopyToDevice(&simdevTemp, 1);
	}
	setDeviceConstantMemory();
	boxStateCopy = std::make_unique<BoxState>(); // TODO, just plain copy it now
	boxConfigCopy = std::make_unique<BoxConfig>();
	cudaMemcpy(boxStateCopy.get(), &sim_dev->boxState, sizeof(BoxState), cudaMemcpyDeviceToHost);
	cudaMemcpy(boxConfigCopy.get(), &sim_dev->boxConfig, sizeof(BoxConfig), cudaMemcpyDeviceToHost);	

	BootstrapClustering();
	MakeSuperClusterTasksGPU();


	for (cudaStream_t& stream : cudaStreams) {
		cudaStreamCreate(&stream);
	}
	cudaStreamCreate(&pmeStream);

	pmeController = std::make_unique<PME::Controller>(*simulation->box, simulation->simParams.cutoff_nm, pmeStream);

	bondgroups.SetData(simulation->box->bondgroups);

	thermostat = std::make_unique<Thermostat>(simulation->box->persistentClusters.size());			

	// To create the NLists we need to bootstrap the traj_buffer, since it has no data yet
	bootstrapTrajbufferWithCoords();
}

Engine::~Engine() {
	if (sim_dev != nullptr) {
		sim_dev->FreeMembers();
		cudaFree(sim_dev);
	}
	forceEnergyInterims->Free();

	for (cudaStream_t& stream : cudaStreams) {
		cudaStreamDestroy(stream);
	}

	if (superClustersControl)
		superClustersControl->Free();

	LIMA_UTILS::genericErrorCheck("Error during Engine destruction");
}


void Engine::setDeviceConstantMemory() {
	cudaMemcpyToSymbol(DeviceConstants::forcefield, &simulation->forcefield, sizeof(ForceField_NB), 0, cudaMemcpyHostToDevice);	// So there should not be a & before the device __constant__


	BoxSize boxSize_host;
	boxSize_host.Set(simulation->box->boxparams.boxSize);
	cudaMemcpyToSymbol(DeviceConstants::boxSize, &boxSize_host, sizeof(BoxSize), 0, cudaMemcpyHostToDevice);

	cudaMemcpyToSymbol(DeviceConstants::cutoffNM, &simulation->simParams.cutoff_nm, sizeof(float), 0, cudaMemcpyHostToDevice);
	const float cutoffNmReciprocal = 1.f / simulation->simParams.cutoff_nm;
	cudaMemcpyToSymbol(DeviceConstants::cutoffNmReciprocal, &cutoffNmReciprocal, sizeof(float), 0, cudaMemcpyHostToDevice);
	const float cutoffNmSquaredReciprocal = 1.f / (simulation->simParams.cutoff_nm * simulation->simParams.cutoff_nm );
	cudaMemcpyToSymbol(DeviceConstants::cutoffNmSquaredReciprocal, &cutoffNmSquaredReciprocal, sizeof(float), 0, cudaMemcpyHostToDevice);	
	const float ewaldKappa = PhysicsUtils::CalcEwaldkappa(simulation->simParams.cutoff_nm);
	cudaMemcpyToSymbol(DeviceConstants::ewaldKappa, &ewaldKappa, sizeof(float), 0, cudaMemcpyHostToDevice);
	const float cutoffNmSquared = simulation->simParams.cutoff_nm * simulation->simParams.cutoff_nm;
	cudaMemcpyToSymbol(DeviceConstants::cutoffNMSquared, &cutoffNmSquared, sizeof(float), 0, cudaMemcpyHostToDevice);

	const float initialThermostatScalar = 1.f;
	cudaMemcpyToSymbol(DeviceConstants::thermostatScalar, &initialThermostatScalar, sizeof(float), 0, cudaMemcpyHostToDevice);

	/*assert(simulation->forcefieldTest.size() == ForceField_NB::MAX_TYPES * ForceField_NB::MAX_TYPES);
	cudaMemcpyToSymbol(DeviceConstants::nonbondedinteractionParams, simulation->forcefieldTest.data(), sizeof(NonbondedInteractionParams) * simulation->forcefieldTest.size(), 0, cudaMemcpyHostToDevice);*/

	// Prepare precomputed values on device
	const float cutoffNM = simulation->simParams.cutoff_nm;
	cudaMemcpyToSymbol(DeviceConstants::bsplineTable, PrecomputeBsplineTable().data(), sizeof(float) * PrecomputeBsplineTable().size(), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(DeviceConstants::erfcForcescalarTable, PrecomputeErfcForcescalarTable(cutoffNM).data(), sizeof(float) * PrecomputeErfcForcescalarTable(cutoffNM).size(), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(DeviceConstants::erfcPotentialscalarTable, PrecomputeErfcPotentialscalarTable(cutoffNM).data(), sizeof(float) * PrecomputeErfcPotentialscalarTable(cutoffNM).size(), 0, cudaMemcpyHostToDevice);

	LIMA_UTILS::genericErrorCheck("Error while setting CUDA __constant__ memory\n");
}




void Engine::step() {
	LIMA_UTILS::genericErrorCheckNoSync("Error before step!");

	deviceMaster();	// Device first, otherwise offloading data always needs the last datapoint!
	simulation->step++;

	hostMaster();
	

	if (simulation->step % simulation->simParams.stepsPerNlistupdate == 0) {
		superClustersControl->Reset(simulation->box->boxparams.boxSize);
		RunClustering();
		MakeSuperClusterTasksGPU();
	}

	LIMA_UTILS::genericErrorCheckNoSync("Error after step!");
}

void Engine::hostMaster() {						// This is and MUST ALWAYS be called after the deviceMaster, and AFTER incStep()!
	auto t0 = std::chrono::high_resolution_clock::now();
	if (DatabuffersDeviceController::IsBufferFull(simulation->getStep(), simulation->simParams.data_logging_interval)) {
		offloadLoggingData(DatabuffersDeviceController::nStepsInBuffer);
		runstatus.stepForMostRecentData = simulation->getStep();

		if ((simulation->getStep() % simulation->simParams.steps_per_temperature_measurement) == 0 && simulation->getStep() > 0) {
			auto [temperature, thermostatScalar] = thermostat->Temperature(sim_dev, simulation->box->boxparams, simulation->simParams, simulation->getStep(), pClusterMetaDevice.Get());
			simulation->temperature_buffer.push_back(temperature);
			runstatus.current_temperature = temperature;

			if (simulation->simParams.apply_thermostat)
				cudaMemcpyToSymbol(DeviceConstants::thermostatScalar, &thermostatScalar, sizeof(float), 0, cudaMemcpyHostToDevice);
		}		
	}
	HandleEarlyStoppingInEM();
	/*if (simulation->getStep() % simulation->simParams.stepsPerNlistupdate == simulation->simParams.stepsPerNlistupdate-1)
		nlistController->UpdateNlist(sim_dev, simulation->box->boxparams, simulation->simParams.bc_select, cudaStreams);*/

	// Handle status
	runstatus.current_step = simulation->getStep();
	if (runstatus.current_step >= simulation->simParams.n_steps || runstatus.critical_error_occured)
		runstatus.simulation_finished = true;


	const auto t1 = std::chrono::high_resolution_clock::now();
	const int cpu_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
}

void Engine::terminateSimulation() {
	const int64_t stepsReadyToTransfer = DatabuffersDeviceController::StepsReadyToTransfer(simulation->getStep(), simulation->simParams.data_logging_interval);
	offloadLoggingData(stepsReadyToTransfer);

	sim_dev->boxState.CopyDataToHost(*simulation->box);

	LIMA_UTILS::genericErrorCheck("Error during TerminateSimulation");
}

//--------------------------------------------------------------------------	CPU workload --------------------------------------------------------------//

void Engine::offloadLoggingData(const int64_t steps_to_transfer) {
	assert(steps_to_transfer <= simulation->getStep());
	if (steps_to_transfer == 0) { return; }

	cudaDeviceSynchronize();

	const int64_t startstep = simulation->getStep() - steps_to_transfer * simulation->simParams.data_logging_interval;
	const int64_t startindex = LIMALOGSYSTEM::getMostRecentDataentryIndex(startstep, simulation->simParams.data_logging_interval);
	const int64_t indices_to_transfer = LIMALOGSYSTEM::getNIndicesBetweenSteps(startstep, simulation->getStep(), simulation->simParams.data_logging_interval);
	//const int particlesUpperbound = simulation->box->boxparams.total_particles_upperbound;
	const int nParticlesUpperbound = simulation->box->persistentClusters.size() * PersistentCluster::maxParticles;
	
	cudaMemcpyAsync(
		simulation->potE_buffer->getBufferAtIndex(startindex),
		dataBuffersDevice->potE_buffer,
		sizeof(float) * nParticlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost);
	
	cudaMemcpyAsync(
		simulation->vel_buffer->getBufferAtIndex(startindex),
		dataBuffersDevice->vel_buffer,
		sizeof(float) * nParticlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost);

	cudaMemcpyAsync(
		simulation->forceBuffer->getBufferAtIndex(startindex),
		dataBuffersDevice->forceBuffer,
		sizeof(Float3) * nParticlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost);

	cudaMemcpyAsync(
		simulation->traj_buffer->getBufferAtIndex(startindex),
		dataBuffersDevice->traj_buffer,
		sizeof(Float3) * nParticlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost);

	step_at_last_traj_transfer = simulation->getStep();
	runstatus.most_recent_positions = simulation->traj_buffer->getBufferAtIndex(LIMALOGSYSTEM::getMostRecentDataentryIndex(simulation->getStep() - 1, simulation->simParams.data_logging_interval));
}

void Engine::offloadTrainData() {
#ifdef GENERATETRAINDATA
	uint64_t values_per_step = N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->boxparams_host.n_compounds;
	if (values_per_step == 0) {
		return;	// No data to transfer
	}

	uint64_t step_offset = (simulation->getStep() - STEPS_PER_TRAINDATATRANSFER) * values_per_step;	// fix max_compound to the actual count save LOTS of space!. Might need a file in simout that specifies cnt for loading in other programs...
	cudaMemcpy(&simulation->trainingdata[step_offset], dataBuffersDevice->data_GAN, sizeof(Float3) * values_per_step * STEPS_PER_TRAINDATATRANSFER, cudaMemcpyDeviceToHost);
	LIMA_UTILS::genericErrorCheck("Cuda error during traindata offloading\n");
#endif
}

//
CudaBuffer<PersistentCluster>& Engine::OffloadPclusterState() {
	pdataCopyBuffer.Expand(simulation->box->persistentClusters.size());
	cudaMemcpy(pdataCopyBuffer.Get(), pClusterDevice.Get(), sizeof(PersistentCluster) * simulation->box->persistentClusters.size(), cudaMemcpyDeviceToDevice);
	return pdataCopyBuffer;
}


__device__ struct SqrtFloat {
	__device__ float operator()(float x) const
	{
		return sqrtf(x);
	}
};

CudaBuffer<float>& Engine::OffloadForcesMagnitudeBuffer() {
	// Copy to offloading buffer
	const int nParticles = simulation->box->boxparams.totalParticles;
	forcesMagnitudeCopyBuffer.Expand(nParticles);
	cudaMemcpy(forcesMagnitudeCopyBuffer.Get(), forcesMagnitudeSquareDevice.Get(), sizeof(float) * nParticles, cudaMemcpyDeviceToDevice);

	// Apply sqrt 
	thrust::device_ptr<float> begin(forcesMagnitudeCopyBuffer.Get());
	thrust::transform(
		thrust::device,
		begin,
		begin + nParticles,
		begin,
		SqrtFloat{}
	);

	return forcesMagnitudeCopyBuffer;
}

void Engine::SetFixedParticleMovementBuffer(const std::vector<Float3>& movement) {
	if (movement.empty()) {
		fixedParticleMovementBuffer.reset();
		return;
	}
	assert(movement.size() == simulation->box->boxparams.totalParticles);
	if (!fixedParticleMovementBuffer.has_value())
		fixedParticleMovementBuffer.emplace();
	fixedParticleMovementBuffer->SetData(movement);
}
void Engine::SetFixedParticleRotationBuffer(const std::vector<Rotation>& rotation) {
	if (rotation.empty()) {
		fixedParticleRotationBuffer.reset();
		return;
	}
	assert(rotation.size() == simulation->box->boxparams.totalParticles);
	if (!fixedParticleRotationBuffer.has_value())
		fixedParticleRotationBuffer.emplace();
	fixedParticleRotationBuffer->SetData(rotation);
}

void Engine::SetForceMask(const std::vector<Float3>& mask) {
	if (mask.empty()) {
		forceMaskBuffer.reset();
		return;
	}
	assert(mask.size() == simulation->box->boxparams.totalParticles);
	if (!forceMaskBuffer.has_value())
		forceMaskBuffer.emplace();
	forceMaskBuffer->SetData(mask); 
}

void Engine::SetElasticPositions(const std::vector<Float3>& positions) {
	if (positions.empty()) {
		elasticPositionsBuffer.reset();
		return;
	}
	assert(positions.size() == simulation->box->boxparams.totalParticles);
	if (!elasticPositionsBuffer.has_value())
		elasticPositionsBuffer.emplace();
	elasticPositionsBuffer->SetData(positions);
}

void Engine::bootstrapTrajbufferWithCoords() {
	if (simulation->simParams.n_steps == 0) return;

	for (int pcid = 0; pcid < simulation->box->persistentClusters.size(); pcid++) {
		const PersistentCluster& pc = simulation->box->persistentClusters[pcid];
		for (int pid = 0; pid < 4; pid++) {
			simulation->traj_buffer->GetDatapoint(pcid, pid, 0) = pc.pqd[pid].position;
		}
	}

	step_at_last_traj_transfer = 0.f;
	runstatus.most_recent_positions = simulation->traj_buffer->getBufferAtIndex(0);

	LIMA_UTILS::genericErrorCheck("Error during bootstrapTrajbufferWithCoords");
}

void Engine::HandleEarlyStoppingInEM() {
	if (!simulation->simParams.em_variant || simulation->getStep() == simulation->simParams.n_steps)
		return;
	
	const int minStepsPerCheck = 100;
	if (simulation->getStep() > stepAtLastEarlystopCheck + minStepsPerCheck) {
		auto forceMagSquared = forcesMagnitudeSquareDevice.GetData(); // [(J/mol/nm)^2]
		const float greatestForce = std::sqrt(Statistics::Max(forceMagSquared.data(), forceMagSquared.size()));
		//const float greatestForce = Statistics::MaxLen(simulation->forceBuffer->GetBufferAtStep(simulation->getStep()-1), simulation->forceBuffer->EntriesPerStep());
		runstatus.greatestForce = greatestForce / KILO; // Convert [J/mol/nm] to [kJ/mol/nm]
		simulation->maxForceBuffer.emplace_back(std::pair<int64_t,float>{ simulation->getStep(), runstatus.greatestForce });

		if (runstatus.greatestForce <= simulation->simParams.em_force_tolerance) {
			runstatus.simulation_finished = true;
		}

		stepAtLastEarlystopCheck = simulation->getStep();
	}
	LIMA_UTILS::genericErrorCheck("HandleEarlyStoppingInEM");
}




//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//
template <typename BoundaryCondition, bool emvariant, bool logData>
void Engine::_deviceMaster() {
	
	const BoxParams& boxparams = simulation->box->boxparams;
	const int step = simulation->getStep();
	const Float3 boxSize = boxparams.BoxSizeFloat();


	// #### Initial round of force computations
	//cudaDeviceSynchronize();

    if (ENABLE_ES_LR && simulation->simParams.enable_electrostatics) {
        pmeController->CalcCharges(superClustersControl->scData, superClustersControl->scMeta, nSuperclusters, forceEnergyInterims->pme, step);
    }


	if (nTasks > 0) {
		const bool useNointeractionMatrix = true;
		dim3 blockDim(SuperCluster::maxParticles, 2, 2);
		NbNonlocalKernel<BoundaryCondition, emvariant, logData, useNointeractionMatrix>
			<<<nTasks, blockDim, 0, cudaStreams[0]>>>
			(superClustersControl->scData, scscTasksDevice.Get(), scResultsDevice.Get(), noInteractionMatricesDevice.Get(), superClustersControl->scMeta, step, boxSize, boxSize.Inv());
		LIMA_UTILS::genericErrorCheckNoSync("Error after NBNonlocalKernel");

		//nbGatherForceenergy.Expand(nSuperclusters * SuperCluster::nParticles, 1.2);
		//NBGather<<<nSuperclusters, dim3(16,4,1), 0, cudaStreams[0]>>>
		//	(superClustersControl->scMeta, scResultsDevice.Get()/*, nbGatherForceenergy.Get()*/);
		//LIMA_UTILS::genericErrorCheckNoSync("Error after NBGather");
	}
	if (!simulation->simParams.snf_select.empty()) {
		SnfHandler<BoundaryCondition, emvariant>(cudaStreams[2]);
		LIMA_UTILS::genericErrorCheckNoSync("Error after SupernaturalForces");
	}

	if (!simulation->box->bondgroups.empty()) {
		BondgroupsKernel<BoundaryCondition, emvariant> << < simulation->box->bondgroups.size(), THREADS_PER_BONDSGROUPSKERNEL, 0, cudaStreams[4]>>>
			(bondgroups.Get(), *boxStateCopy, forceEnergyInterims->forceEnergiesBondgroups, pClusterDevice.Get(), boxSize, boxSize.Inv());
		LIMA_UTILS::genericErrorCheckNoSync("Error after BondgroupsKernel");

		// Gather bondgroup ordered forces into particle ordered
		const int nPclusters = simulation->box->persistentClusters.size();
		const int nBlocks = (nPclusters + 31) / 32;
		PclusterBondgroupsGather << <nBlocks, 32, 0, cudaStreams[4] >> >
			(pClusterMetaDevice.Get(), nPclusters, *forceEnergyInterims);
		LIMA_UTILS::genericErrorCheckNoSync("Error after PclusterBondgroupsGather");
	}

	// #### Integration and Transfer kernels
	cudaStreamSynchronize(pmeStream);
	for (int i = 0; i < cudaStreams.size(); i++) {
		cudaStreamSynchronize(cudaStreams[i]);
	}



	const bool updateNlistsAfterThisStep = (simulation->getStep()+1) % simulation->simParams.stepsPerNlistupdate == simulation->simParams.stepsPerNlistupdate-1;


	if (nSuperclusters > 0) {
		int totalParticlesUpperbound = simulation->box->persistentClusters.size() * PersistentCluster::maxParticles;
		const int nBlocks = (nSuperclusters + 4 - 1) / 4;
		const dim3 blockDim(16, 4, 1);

		Float3* fixedParticleMovementBufferPtr = fixedParticleMovementBuffer.has_value() ? fixedParticleMovementBuffer->Get() : nullptr;
		Float3* forcesMaskBufferPtr = forceMaskBuffer.has_value() ? forceMaskBuffer->Get() : nullptr;
		Rotation* fixedParticleRotationBufferPtr = fixedParticleRotationBuffer.has_value() ? fixedParticleRotationBuffer->Get() : nullptr;
		SuperclusterIntegrateKernel<BoundaryCondition, emvariant, logData>
			<<<nBlocks, blockDim, 0, cudaStreams[0]>>>
			(*forceEnergyInterims, sim_dev, simulation->simParams.data_logging_interval, scResultsDevice.Get(), superClustersControl->scData, superClustersControl->scMeta, pClusterDevice.Get(), pClusterMetaDevice.Get(), 
				boxStateCopy->pclusterInterimStates, step, simulation->simParams.dt, totalParticlesUpperbound, nSuperclusters, forcesMagnitudeSquareDevice.Get(), fixedParticleMovementBufferPtr, 
				forcesMaskBufferPtr, fixedParticleRotationBufferPtr);
		LIMA_UTILS::genericErrorCheckNoSync("Error after SuperclusterIntegrateKernel");
		cudaDeviceSynchronize();
	}

	//DebugUtils::VerifyIdentical(superClustersControl->scData, nSuperclusters, "Engine_SCData", step);
	//DebugUtils::VerifyIdentical(pClusterDevice, simulation->box->persistentClusters.size(), "Engine_PClusterDevice", step);
}



void Engine::deviceMaster() {

	const bool logData = simulation->simParams.data_logging_interval != 0 && simulation->getStep() % simulation->simParams.data_logging_interval == 0;// TODO maybe log at the final step, not 0th?

	switch (simulation->simParams.bc_select) {
	case NoBC:
		if (simulation->simParams.em_variant) {
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
		if (simulation->simParams.em_variant) {
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
	if (simulation->simParams.snf_select.contains(HorizontalSqueeze)) {
		//SupernaturalForces::ApplyHorizontalSqueeze << < simulation->box->boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, stream >> > (sim_dev, simulation->getStep());
		//break;
	}
	if (simulation->simParams.snf_select.contains(HorizontalChargeField))
	{
		const int nPclusters = simulation->box->persistentClusters.size();
		const int nCudablocks = (nPclusters + 31) / 32;
		PclusterSnfKernel<BoundaryCondition, emvariant>
			<<<nCudablocks, 32, 0, stream >> >
			(pClusterDevice.Get(), pClusterMetaDevice.Get(), simulation->box->uniformElectricField, forceEnergyInterims->snf, nPclusters);
	}
	
	if (simulation->simParams.snf_select.contains(SupernaturalForcesSelect::ElasticPosition) && elasticPositionsBuffer.has_value()) {
		const int nPclusters = simulation->box->persistentClusters.size();
		const int nCudablocks = (nPclusters + 31) / 32;
		ElasticPositionsForceKernel << <nCudablocks, 32, 0, stream >> >
			(pClusterDevice.Get(), pClusterMetaDevice.Get(), elasticPositionsBuffer->Get(), forceEnergyInterims->snf, nPclusters, simulation->box->boxparams.BoxSizeFloat());
	}
		
		
	//case BoxEdgePotential:
	//	if (simulation->box->boxparams.n_compounds > 0)
	//		SupernaturalForces::BoxEdgeForceCompounds << < simulation->box->boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, stream >> > (sim_dev, simulation->getStep());
	//	if (simulation->box->boxparams.nTinymols > 0)
	//		SupernaturalForces::BoxEdgeForceSolvents<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(simulation->box->boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, stream>>>(sim_dev, simulation->getStep());
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


			cudaDeviceSynchronize();

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
	runCase(1, 64);   // effectively 1×64
	runCase(4, 16);    // effectively 4×16
	runCase(8, 8);   // effectively 16×4

	return success;
}