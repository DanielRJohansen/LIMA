#include "Engine.cuh"

#include "BoundaryCondition.cuh"
#include "EngineBodies.cuh"
#include "EngineKernels.cuh"
#include "LimaPositionSystem.cuh"
#include "Neighborlists.cuh"
#include "PME.cuh"
#include "SimulationDevice.cuh"
#include "SupernaturalForces.cuh"
#include "Thermostat.cuh"

#include "Statistics.h"
#include "Utilities.h"
#include "DebugUtils.h"

#include "EngineHostside.h"

#include <random>
#include <numeric>
#include "ParticleClusters.cuh"

#include <set>
#include <execution>
#include "EngineCore.h"

Engine::Engine(std::unique_ptr<Simulation> _sim, BoundaryConditionSelect bc, std::unique_ptr<LimaLogger> logger)
	: bc_select(bc)
	, m_logger(std::move(logger))
	, forceEnergyInterims(std::make_unique<ForceEnergyInterims>(_sim->box_host->boxparams.n_compounds, _sim->box_host->boxparams.nTinymols, BoxGrid::BlocksTotal(_sim->box_host->boxparams.boxSize), _sim->box_host->bondgroups.size(), _sim->box_host->boxparams.total_particles))
{
	simulation = std::move(_sim);

    verifyEngine();

	const BoxParams boxparams = simulation->box_host->boxparams;

	dataBuffersDevice = std::make_unique<DatabuffersDeviceController>(boxparams.total_particles_upperbound, 
		boxparams.n_compounds, simulation->simparams_host.data_logging_interval);

	superClustersControl = std::make_unique<SuperClustersControl>(SuperClustersControl::Create(boxparams.boxSize, simulation->box_host->persistentClusters.size()));
	pclusterTransfermodule = std::make_unique<PClusterTransfermodule>(PClusterTransfermodule::Create(boxparams.boxSize));
	//cudaMalloc(&pClusterDevice, sizeof(PersistentCluster) * simulation->box_host->persistentClusters.size());// We will never have more sc than pc
	pClusterDevice = GenericCopyToDevice(simulation->box_host->persistentClusters);
	pClusterMetaDevice = GenericCopyToDevice(simulation->box_host->persistentClustersMetadata);
	particleToCompoundOrSolventMappingDevice = GenericCopyToDevice(simulation->box_host->particleToCompoundOrSolventMapping);





	// Create the Sim_dev {
	{
		if (sim_dev != nullptr) { throw std::runtime_error("Expected simdev to be null to move sim to device"); }
		sim_dev = new SimulationDevice(simulation->simparams_host, simulation->box_host.get(), BoxConfig::Create(*simulation->box_host), BoxState::Create(*simulation->box_host), *dataBuffersDevice);
		sim_dev = genericMoveToDevice(sim_dev, 1);
	}
	setDeviceConstantMemory();
	boxStateCopy = std::make_unique<BoxState>(); // TODO, just plain copy it now
	boxConfigCopy = std::make_unique<BoxConfig>();
	cudaMemcpy(boxStateCopy.get(), &sim_dev->boxState, sizeof(BoxState), cudaMemcpyDeviceToHost);
	cudaMemcpy(boxConfigCopy.get(), &sim_dev->boxConfig, sizeof(BoxConfig), cudaMemcpyDeviceToHost);	
	nParticlesInCompoundsBufferPtr = sim_dev->nParticlesInCompoundsBuffer;

	BootstrapClustering();
	cudaMalloc(&scscTasksDevice, 1); // Bootstrap these, will be reallocated in the next function. TODO: CHange that system so we dont need realloc
	cudaMalloc(&noInteractionMatricesDevice, 1);
	cudaMalloc(&scResultsDevice, 1);
	MakeSuperClusterTasksCPU();

    std::vector<ForceField_NB::ParticleParameters> compoundParticleParams(boxparams.n_compounds * MAX_COMPOUND_PARTICLES, ForceField_NB::ParticleParameters{0,0});
    for (int cid = 0; cid < simulation->box_host->compounds.size(); cid++) {
        const Compound& compound = simulation->box_host->compounds[cid];
        for (int pid = 0; pid < compound.n_particles; pid++) {
            compoundParticleParams[cid*MAX_COMPOUND_PARTICLES + pid] = simulation->forcefield.particle_parameters[compound.atom_types[pid]];
        }

    }

	for (cudaStream_t& stream : cudaStreams) {
		cudaStreamCreate(&stream);
	}
	cudaStreamCreate(&pmeStream);

	pmeController = std::make_unique<PME::Controller>(*simulation->box_host, simulation->simparams_host.cutoff_nm, pmeStream);

	bondgroups = GenericCopyToDevice(simulation->box_host->bondgroups);

	compoundQuickData = CompoundQuickData::CreateBuffer(*simulation);

	thermostat = std::make_unique<Thermostat>(boxparams.n_compounds, boxparams.nTinymolParticles, boxparams.total_particles_upperbound);

	nlistController = std::make_unique<NeighborList::Controller>(boxparams);
	
	tinymolTransferModule = std::make_unique<TinymolTransferModule>(TinymolTransferModule::Create(BoxGrid::BlocksTotal(boxparams.boxSize)));
		




	// To create the NLists we need to bootstrap the traj_buffer, since it has no data yet
	bootstrapTrajbufferWithCoords();

	BootstrapSolventblockDistributeFromDensity();

	nlistController->UpdateNlist(sim_dev, boxparams, simulation->simparams_host.bc_select, cudaStreams);
	m_logger->finishSection("Engine Ready");
}

Engine::~Engine() {
	if (sim_dev != nullptr) {
		sim_dev->FreeMembers();
		cudaFree(sim_dev);
	}
	forceEnergyInterims->Free();

	cudaFree(bondgroups);

	for (cudaStream_t& stream : cudaStreams) {
		cudaStreamDestroy(stream);
	}

	if (superClustersControl)
		superClustersControl->Free();

	LIMA_UTILS::genericErrorCheck("Error during Engine destruction");
	//assert(simulation == nullptr);
}


void Engine::setDeviceConstantMemory() {
	//const int forcefield_bytes = sizeof(ForceField_NB);
	cudaMemcpyToSymbol(DeviceConstants::forcefield, &simulation->forcefield, sizeof(ForceField_NB), 0, cudaMemcpyHostToDevice);	// So there should not be a & before the device __constant__
	cudaMemcpyToSymbol(DeviceConstants::tinymolForcefield, &simulation->forcefieldTinymol, sizeof(ForcefieldTinymol), 0, cudaMemcpyHostToDevice);

	{
		auto t0 = simulation->forcefieldTinymol.types[0];
		auto t1 = simulation->forcefieldTinymol.types[1];
		NonbondedInteractionParams precomputedParams[3]{
			{LJ::CalcSigma(t0.sigmaHalf, t0.sigmaHalf), LJ::CalcEpsilon(t0.epsilonSqrt, t0.epsilonSqrt), t0.charge*t0.charge},
			{LJ::CalcSigma(t0.sigmaHalf, t1.sigmaHalf), LJ::CalcEpsilon(t0.epsilonSqrt, t1.epsilonSqrt), t0.charge*t1.charge},
			{LJ::CalcSigma(t1.sigmaHalf, t1.sigmaHalf), LJ::CalcEpsilon(t1.epsilonSqrt, t1.epsilonSqrt), t1.charge*t1.charge}
		};
		cudaMemcpyToSymbol(DeviceConstants::tinymolPrecomputedParams, precomputedParams, sizeof(NonbondedInteractionParams) * 3, 0, cudaMemcpyHostToDevice);
	}

	BoxSize boxSize_host;
	boxSize_host.Set(simulation->box_host->boxparams.boxSize);
	cudaMemcpyToSymbol(DeviceConstants::boxSize, &boxSize_host, sizeof(BoxSize), 0, cudaMemcpyHostToDevice);

	cudaMemcpyToSymbol(DeviceConstants::cutoffNM, &simulation->simparams_host.cutoff_nm, sizeof(float), 0, cudaMemcpyHostToDevice);
	const float cutoffNmReciprocal = 1.f / simulation->simparams_host.cutoff_nm;
	cudaMemcpyToSymbol(DeviceConstants::cutoffNmReciprocal, &cutoffNmReciprocal, sizeof(float), 0, cudaMemcpyHostToDevice);
	const float cutoffNmSquaredReciprocal = 1.f / (simulation->simparams_host.cutoff_nm * simulation->simparams_host.cutoff_nm );
	cudaMemcpyToSymbol(DeviceConstants::cutoffNmSquaredReciprocal, &cutoffNmSquaredReciprocal, sizeof(float), 0, cudaMemcpyHostToDevice);	
	const float ewaldKappa = PhysicsUtils::CalcEwaldkappa(simulation->simparams_host.cutoff_nm);
	cudaMemcpyToSymbol(DeviceConstants::ewaldKappa, &ewaldKappa, sizeof(float), 0, cudaMemcpyHostToDevice);
	const float cutoffNmSquared = simulation->simparams_host.cutoff_nm * simulation->simparams_host.cutoff_nm;
	cudaMemcpyToSymbol(DeviceConstants::cutoffNMSquared, &cutoffNmSquared, sizeof(float), 0, cudaMemcpyHostToDevice);

	const float initialThermostatScalar = 1.f;
	cudaMemcpyToSymbol(DeviceConstants::thermostatScalar, &initialThermostatScalar, sizeof(float), 0, cudaMemcpyHostToDevice);

	/*assert(simulation->forcefieldTest.size() == ForceField_NB::MAX_TYPES * ForceField_NB::MAX_TYPES);
	cudaMemcpyToSymbol(DeviceConstants::nonbondedinteractionParams, simulation->forcefieldTest.data(), sizeof(NonbondedInteractionParams) * simulation->forcefieldTest.size(), 0, cudaMemcpyHostToDevice);*/

	// Prepare precomputed values on device
	const float cutoffNM = simulation->simparams_host.cutoff_nm;
	cudaMemcpyToSymbol(DeviceConstants::bsplineTable, PrecomputeBsplineTable().data(), sizeof(float) * PrecomputeBsplineTable().size(), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(DeviceConstants::erfcForcescalarTable, PrecomputeErfcForcescalarTable(cutoffNM).data(), sizeof(float) * PrecomputeErfcForcescalarTable(cutoffNM).size(), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(DeviceConstants::erfcPotentialscalarTable, PrecomputeErfcPotentialscalarTable(cutoffNM).data(), sizeof(float) * PrecomputeErfcPotentialscalarTable(cutoffNM).size(), 0, cudaMemcpyHostToDevice);

	LIMA_UTILS::genericErrorCheck("Error while setting CUDA __constant__ memory\n");
}




void Engine::step() {
	LIMA_UTILS::genericErrorCheckNoSync("Error before step!");

	deviceMaster();	// Device first, otherwise offloading data always needs the last datapoint!
	assert(simulation);
	assert(sim_dev);
	simulation->step++;

	hostMaster();
	

	if (true) {
		superClustersControl->Reset();
		RunClustering();
		MakeSuperClusterTasksCPU();
	}

	LIMA_UTILS::genericErrorCheckNoSync("Error after step!");
}

void Engine::hostMaster() {						// This is and MUST ALWAYS be called after the deviceMaster, and AFTER incStep()!
	auto t0 = std::chrono::high_resolution_clock::now();
	if (DatabuffersDeviceController::IsBufferFull(simulation->getStep(), simulation->simparams_host.data_logging_interval)) {
		offloadLoggingData(DatabuffersDeviceController::nStepsInBuffer);
		runstatus.stepForMostRecentData = simulation->getStep();

		if ((simulation->getStep() % simulation->simparams_host.steps_per_temperature_measurement) == 0 && simulation->getStep() > 0) {
			auto [temperature, thermostatScalar] = thermostat->Temperature(sim_dev, simulation->box_host->boxparams, simulation->simparams_host, simulation->getStep());
			simulation->temperature_buffer.push_back(temperature);
			runstatus.current_temperature = temperature;

			if (simulation->simparams_host.apply_thermostat)
				cudaMemcpyToSymbol(DeviceConstants::thermostatScalar, &thermostatScalar, sizeof(float), 0, cudaMemcpyHostToDevice);
		}
		
		HandleEarlyStoppingInEM();
	}
	if (simulation->getStep() % simulation->simparams_host.stepsPerNlistupdate == simulation->simparams_host.stepsPerNlistupdate-1)
		nlistController->UpdateNlist(sim_dev, simulation->box_host->boxparams, simulation->simparams_host.bc_select, cudaStreams);

	// Handle status
	runstatus.current_step = simulation->getStep();
	if (runstatus.current_step >= simulation->simparams_host.n_steps || runstatus.critical_error_occured)
		runstatus.simulation_finished = true;


	const auto t1 = std::chrono::high_resolution_clock::now();
	const int cpu_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
}

void Engine::terminateSimulation() {
	const int64_t stepsReadyToTransfer = DatabuffersDeviceController::StepsReadyToTransfer(simulation->getStep(), simulation->simparams_host.data_logging_interval);
	offloadLoggingData(stepsReadyToTransfer);

	sim_dev->boxState.CopyDataToHost(*simulation->box_host);

	const float greatestForce = Statistics::MaxLen(simulation->forceBuffer->GetBufferAtStep(simulation->getStep() - 1), simulation->forceBuffer->EntriesPerStep());
	runstatus.greatestForce = greatestForce / KILO; // Convert [J/mol/nm] to [kJ/mol/nm]
	simulation->maxForceBuffer.emplace_back(std::pair<int64_t, float>{ simulation->getStep(), runstatus.greatestForce });

	LIMA_UTILS::genericErrorCheck("Error during TerminateSimulation");
}

//--------------------------------------------------------------------------	CPU workload --------------------------------------------------------------//

void Engine::offloadLoggingData(const int64_t steps_to_transfer) {
	assert(steps_to_transfer <= simulation->getStep());
	if (steps_to_transfer == 0) { return; }

	cudaDeviceSynchronize();

	const int64_t startstep = simulation->getStep() - steps_to_transfer * simulation->simparams_host.data_logging_interval;
	const int64_t startindex = LIMALOGSYSTEM::getMostRecentDataentryIndex(startstep, simulation->simparams_host.data_logging_interval);
	const int64_t indices_to_transfer = LIMALOGSYSTEM::getNIndicesBetweenSteps(startstep, simulation->getStep(), simulation->simparams_host.data_logging_interval);
	const int particlesUpperbound = simulation->box_host->boxparams.total_particles_upperbound;
	
	cudaMemcpyAsync(
		simulation->potE_buffer->getBufferAtIndex(startindex),
		dataBuffersDevice->potE_buffer,
		sizeof(float) * particlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost);
	
	cudaMemcpyAsync(
		simulation->vel_buffer->getBufferAtIndex(startindex),
		dataBuffersDevice->vel_buffer,
		sizeof(float) * particlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost);

	cudaMemcpyAsync(
		simulation->forceBuffer->getBufferAtIndex(startindex),
		dataBuffersDevice->forceBuffer,
		sizeof(Float3) * particlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost);

	cudaMemcpyAsync(
		simulation->traj_buffer->getBufferAtIndex(startindex),
		dataBuffersDevice->traj_buffer,
		sizeof(Float3) * particlesUpperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost);

	step_at_last_traj_transfer = simulation->getStep();
	runstatus.most_recent_positions = simulation->traj_buffer->getBufferAtIndex(LIMALOGSYSTEM::getMostRecentDataentryIndex(simulation->getStep() - 1, simulation->simparams_host.data_logging_interval));
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


void Engine::bootstrapTrajbufferWithCoords() {
	if (simulation->simparams_host.n_steps == 0) return;

	LIMA_UTILS::genericErrorCheck("Error during bootstrapTrajbufferWithCoords");

	// We need to bootstrap step-0 which is used for traj-buffer
	for (int compound_id = 0; compound_id < simulation->box_host->boxparams.n_compounds; compound_id++) {
		for (int particle_id = 0; particle_id < MAX_COMPOUND_PARTICLES; particle_id++) {
			const Float3 particle_abspos = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(simulation->box_host->compoundCoordsBuffer[compound_id].origo, simulation->box_host->compoundCoordsBuffer[compound_id].rel_positions[particle_id]);
			simulation->traj_buffer->getCompoundparticleDatapointAtIndex(compound_id, particle_id, 0) = particle_abspos;
		}
	}

	for (int blockId = 0; blockId < BoxGrid::BlocksTotal(simulation->box_host->boxparams.boxSize); blockId++) {
		const SolventBlock& solventBlock = *SolventBlocksCircularQueue::getBlockPtr(simulation->box_host->solventblockgrid_circularqueue.data(), BoxGrid::NodesPerDim(simulation->box_host->boxparams.boxSize), blockId, 0);
		const NodeIndex origo = BoxGrid::Get3dIndexWithNNodes(blockId, BoxGrid::NodesPerDim(simulation->box_host->boxparams.boxSize));

		for (int pid = 0; pid < solventBlock.nParticles; pid++) {
			const int solventId = solventBlock.ids[pid];
			const Float3 pos = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(origo, solventBlock.rel_pos[pid].ToRelpos());
			simulation->traj_buffer->getSolventparticleDatapointAtIndex(solventId, 0) = pos;
		}		
	}
	step_at_last_traj_transfer = 0.f;
	runstatus.most_recent_positions = simulation->traj_buffer->getBufferAtIndex(0);

	LIMA_UTILS::genericErrorCheck("Error during bootstrapTrajbufferWithCoords");
}

void Engine::BootstrapSolventblockDistributeFromDensity() {
	Int3 boxSize = simulation->box_host->boxparams.boxSize;

	// Bootstrap compressed positions
	int nGridblocks = BoxGrid::NodesPerDim(boxSize.y) * BoxGrid::NodesPerDim(boxSize.z);
	SolventPositionsBufferCompress << <nGridblocks, 32, 0, cudaStreams[1] >> >
		(*boxStateCopy, *boxConfigCopy, simulation->box_host->boxparams);

	SolventBlockAdjacencySequenceUpdate << <BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxSize)), 64, 0, cudaStreams[1] >> >
		(*boxStateCopy, *boxConfigCopy, simulation->box_host->boxparams);
}



void Engine::HandleEarlyStoppingInEM() {
	if (!simulation->simparams_host.em_variant || simulation->getStep() == simulation->simparams_host.n_steps)
		return;
	
	const int minStepsPerCheck = 100;
	if (simulation->getStep() > stepAtLastEarlystopCheck + minStepsPerCheck) {
		const float greatestForce = Statistics::MaxLen(simulation->forceBuffer->GetBufferAtStep(simulation->getStep()-1), simulation->forceBuffer->EntriesPerStep());
		runstatus.greatestForce = greatestForce / KILO; // Convert [J/mol/nm] to [kJ/mol/nm]
		simulation->maxForceBuffer.emplace_back(std::pair<int64_t,float>{ simulation->getStep(), runstatus.greatestForce });

		if (runstatus.greatestForce <= simulation->simparams_host.em_force_tolerance) {
			runstatus.simulation_finished = true;
		}

		stepAtLastEarlystopCheck = simulation->getStep();
	}
}




//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//
template <typename BoundaryCondition, bool emvariant, bool computePotE>
void Engine::_deviceMaster() {
	
	const BoxParams& boxparams = simulation->box_host->boxparams;
	const int step = simulation->getStep();

	// #### Initial round of force computations
	cudaDeviceSynchronize();

    if (ENABLE_ES_LR && simulation->simparams_host.enable_electrostatics) {
        pmeController->CalcCharges(*boxConfigCopy, *boxStateCopy, boxparams.n_compounds, forceEnergyInterims->forceEnergiesPME, forceEnergyInterims->solvents.pmeInteraction, pmeStream);
        LIMA_UTILS::genericErrorCheckNoSync("Error after HandleElectrostatics");
    }

	bool newAlg = true;

	//if (newAlg)
	if (nTasks > 0) {
		const bool useNointeractionMatrix = true;
		dim3 blockDim(16, 1, 1); // TEMP
		NbNonlocalKernel<BoundaryCondition, emvariant, computePotE, useNointeractionMatrix>
			<<<nTasks, blockDim, 0, cudaStreams[0]>>>
			(superClustersControl->scData, scscTasksDevice, scResultsDevice, noInteractionMatricesDevice, superClustersControl->scMeta, step);
		LIMA_UTILS::genericErrorCheckNoSync("Error after NBNonlocalKernel");

		/*std::vector<SCResult> results = GenericCopyToHost(scResultsDevice, nResults);
		DebugUtils::VerifyIdentical(results, "SCresults" + std::to_string(simulation->getStep()));*/


		//std::vector<PersistentClusterMeta> pcMetaTemp = GenericCopyToHost(pClusterMetaDevice, simulation->box_host->persistentClusters.size());

		SuperclusterForceenergyReduce<<<nSuperclusters, 16, 0, cudaStreams[0] >> >
			(superClustersControl->scMeta, pClusterMetaDevice, scResultsDevice, forceEnergyInterims->nbNonlocal);
		LIMA_UTILS::genericErrorCheckNoSync("Error after SuperclusterForceenergyReduce");		

		/*std::vector<ForceEnergy> feNonlocal = GenericCopyToHost(forceEnergyInterims->nbNonlocal, boxparams.total_particles);
		DebugUtils::VerifyIdentical(feNonlocal, "FeNonlocal" + std::to_string(simulation->getStep()));*/
	}
	/*else*/ 
	//{
	//	cudaDeviceSynchronize();
	//	if (boxparams.n_compounds > 0) {
	//	compoundFarneighborShortrangeInteractionsKernel<BoundaryCondition, emvariant, computePotE> 
	//		<<<boxparams.n_compounds, MAX_COMPOUND_PARTICLES>>>
 //           (simulation->simparams_host.enable_electrostatics,
 //               forceEnergyInterims->forceEnergyFarneighborShortrange, compoundQuickData, nlistController->GetBuffers().compoundsNNeighborNonbondedCompounds, 
	//			nlistController->GetBuffers().compoundsNeighborNonbondedCompounds, nParticlesInCompoundsBufferPtr, 
	//			sim_dev,
	//			//nullptr,
 //               step
	//			);
	//	cudaDeviceSynchronize();
	//	int a = 0;
	//	LIMA_UTILS::genericErrorCheckNoSync("Error after compoundFarneighborShortrangeInteractionsKernel");

	//	compoundImmediateneighborAndSelfShortrangeInteractionsKernel<BoundaryCondition, emvariant, computePotE> 
	//		<<<boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, cudaStreams[1] >>> 
	//		(sim_dev, step, forceEnergyInterims->forceEnergyImmediateneighborShortrange, nlistController->GetBuffers());
	//	LIMA_UTILS::genericErrorCheckNoSync("Error after compoundImmediateneighborAndSelfShortrangeInteractionsKernel");
	//	}
	//}

	

	if (boxparams.nTinymols > 0) {
		const int nSolventblocks = BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize));
		// Should only use max_compound_particles threads here. and let 1 thread handle multiple solvents
		TinymolCompoundinteractionsKernel<BoundaryCondition, emvariant>
			<<<nSolventblocks, SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, cudaStreams[2]>>>
			(*boxStateCopy, *boxConfigCopy, nlistController->GetBuffers(), step, forceEnergyInterims->solvents.compoundsInteractions);
		LIMA_UTILS::genericErrorCheckNoSync("Error after TinymolCompoundinteractionsKernel");
	 

		constexpr auto occRanges = SolventBlockOccupancy::ranges;
		solventForceKernel<BoundaryCondition, emvariant, computePotE, occRanges[0].batchsize, occRanges[0].min, occRanges[0].max>
			<<<nSolventblocks, occRanges[0].max, 0, cudaStreams[3] >> >
			(*boxStateCopy, forceEnergyInterims->solvents.solventsInteractions);
		LIMA_UTILS::genericErrorCheckNoSync("Error after solventForceKernel");		

		solventForceKernel<BoundaryCondition, emvariant, computePotE, occRanges[1].batchsize, occRanges[1].min, occRanges[1].max>
			<<<nSolventblocks, occRanges[1].max, 0, cudaStreams[3] >> >
			(*boxStateCopy, forceEnergyInterims->solvents.solventsInteractions);
		LIMA_UTILS::genericErrorCheckNoSync("Error after solventForceKernel");

		solventForceKernel<BoundaryCondition, emvariant, computePotE, occRanges[2].batchsize, occRanges[2].min, occRanges[2].max>
			<<<nSolventblocks, occRanges[2].max, 0, cudaStreams[3] >> >
			(*boxStateCopy, forceEnergyInterims->solvents.solventsInteractions);
		LIMA_UTILS::genericErrorCheckNoSync("Error after solventForceKernel");

		TinymolBondgroupsKernel<emvariant>
			<< <nSolventblocks, dim3(SolventBlock::maxBondgroups, 1, 1), 0, cudaStreams[2] >> >
			(sim_dev, step, forceEnergyInterims->solvents.bondgroupsInteractions);
		LIMA_UTILS::genericErrorCheckNoSync("Error after TinymolBondgroupsKernel");
	}
	


	if (simulation->simparams_host.snf_select != None) {
		SnfHandler<BoundaryCondition, emvariant>(cudaStreams[2]);
		LIMA_UTILS::genericErrorCheckNoSync("Error after SupernaturalForces");
	}

	if (!simulation->box_host->bondgroups.empty()) {
		BondgroupsKernel<BoundaryCondition, emvariant> << < simulation->box_host->bondgroups.size(), THREADS_PER_BONDSGROUPSKERNEL, 0, cudaStreams[4]>>> 
			(bondgroups, *boxStateCopy, forceEnergyInterims->forceEnergiesBondgroups, pClusterDevice);
		LIMA_UTILS::genericErrorCheckNoSync("Error after BondgroupsKernel");

		// Gather bondgroup ordered forces into particle ordered
		{
			const int nPclusters = simulation->box_host->persistentClusters.size();
			const int nBlocks = (nPclusters + 31) / 32;
			PclusterBondgroupsGather << <nBlocks, 32, 0, cudaStreams[0] >> >
				(pClusterMetaDevice, nPclusters, *forceEnergyInterims);
			LIMA_UTILS::genericErrorCheckNoSync("Error after PclusterBondgroupsGather");
		}
	}

	// #### Integration and Transfer kernels
	cudaStreamSynchronize(pmeStream);
	for (int i = 0; i < cudaStreams.size(); i++) {
		cudaStreamSynchronize(cudaStreams[i]);
	}

	cudaDeviceSynchronize();
	DistributePlcusterForceenergyToCompoundsAndSolvents << <(boxparams.total_particles + 31) / 32, 32, 0, cudaStreams[0] >> >
		(*forceEnergyInterims, particleToCompoundOrSolventMappingDevice, boxparams.total_particles);
	LIMA_UTILS::genericErrorCheckNoSync("Error after DistributePlcusterForceenergyToCompoundsAndSolvents");
	cudaDeviceSynchronize();
	// TODO: Do i need sync before and after this? Yes, right?? Which is why i dont want this kernel at all, the logic should be inside the integration kernel


	//auto bondsPresort = GenericCopyToHost(forceEnergyInterims->forceEnergiesBondgroups, boxparams.total_particles);
	//auto bondForces = GenericCopyToHost(forceEnergyInterims->bonded, boxparams.total_particles);
	//auto sumForces = GenericCopyToHost(forceEnergyInterims->fromSuperclusters, boxparams.total_particles);
	//if (step >= 100)
	//	int a = 0;

	const bool updateNlistsAfterThisStep = (simulation->getStep()+1) % simulation->simparams_host.stepsPerNlistupdate == simulation->simparams_host.stepsPerNlistupdate-1;
	if (boxparams.n_compounds > 0) {
		CompoundIntegrationKernel<BoundaryCondition, emvariant> 
			<<<boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, cudaStreams[0] >> >
			(sim_dev, step, *forceEnergyInterims, compoundQuickData, updateNlistsAfterThisStep);
		LIMA_UTILS::genericErrorCheckNoSync("Error after CompoundIntegrationKernel");
	}

	if (boxparams.nTinymols > 0) {	
		TinymolIntegrateAndLogKernel<BoundaryCondition, emvariant>
			<< <BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, cudaStreams[1] >> >
			(sim_dev, step, *forceEnergyInterims);
		LIMA_UTILS::genericErrorCheckNoSync("Error after TinymolIntegrateAndLogKernel");

		if (SolventBlocksCircularQueue::isTransferStep(step)) {
			SolventPretransferKernel<BoundaryCondition> 
				<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, cudaStreams[1]>>> 
				(sim_dev, step, *tinymolTransferModule);
			LIMA_UTILS::genericErrorCheckNoSync("Error after SolventPretransferKernel");

			SolventTransferKernel<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, cudaStreams[1] >>> (sim_dev, step, *tinymolTransferModule);
			LIMA_UTILS::genericErrorCheckNoSync("Error after SolventTransferKernel");

			int nGridblocks = BoxGrid::NodesPerDim(boxparams.boxSize.y) * BoxGrid::NodesPerDim(boxparams.boxSize.z);
			SolventPositionsBufferCompress << <nGridblocks, 32, 0, cudaStreams[1] >> >
				(*boxStateCopy, *boxConfigCopy, boxparams);
			LIMA_UTILS::genericErrorCheckNoSync("Error after SolventPositionsBufferCompress");

			SolventBlockAdjacencySequenceUpdate << <BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), 32, 0, cudaStreams[1] >> >
				(*boxStateCopy, *boxConfigCopy, boxparams);
			LIMA_UTILS::genericErrorCheckNoSync("Error after SolventBlockAdjacencySequenceUpdate");
		}
	}

	UpdatePdataPositions<<<nSuperclusters, 16>>>
		(*boxStateCopy, particleToCompoundOrSolventMappingDevice, superClustersControl->scMeta, pClusterMetaDevice, superClustersControl->scData, pClusterDevice);
	LIMA_UTILS::genericErrorCheckNoSync("Error after SolventBlockAdjacencySequenceUpdate");
}



void Engine::deviceMaster() {

	const bool logData = simulation->getStep() % simulation->simparams_host.data_logging_interval == 0;// TODO maybe log at the final step, not 0th?

	switch (simulation->simparams_host.bc_select) {
	case NoBC:
		if (simulation->simparams_host.em_variant) {
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
		if (simulation->simparams_host.em_variant) {
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
	switch (simulation->simparams_host.snf_select) {
	case None:
		break;
	case HorizontalSqueeze:
		SupernaturalForces::ApplyHorizontalSqueeze << < simulation->box_host->boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, stream >> > (sim_dev, simulation->getStep());
		break;
	case HorizontalChargeField:
		CompoundSnfKernel<BoundaryCondition, emvariant>
			<< <simulation->box_host->boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, stream >>>
			(sim_dev, simulation->box_host->uniformElectricField, forceEnergyInterims->forceEnergySNF);
		break;
	case BoxEdgePotential:
		if (simulation->box_host->boxparams.n_compounds > 0)
			SupernaturalForces::BoxEdgeForceCompounds << < simulation->box_host->boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, stream >> > (sim_dev, simulation->getStep());
		if (simulation->box_host->boxparams.nTinymols > 0)
			SupernaturalForces::BoxEdgeForceSolvents<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(simulation->box_host->boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, stream>>>(sim_dev, simulation->getStep());
		break;
	}
}























bool ScAreBonded(const SuperClusterMeta& sc0, const SuperClusterMeta& sc1, const std::vector<std::set<int>>& bondedToPclusters) {
	for (auto pclusterId0 : sc0.pclusterIds) {
		if (pclusterId0 == -1) break;
		for (auto pclusterId1 : sc1.pclusterIds) {
			if (pclusterId1 == -1) break;

			if (bondedToPclusters[pclusterId0].contains(pclusterId1)) {
				return true;
			}
		}
	}
	return false;
}

std::array<int, 16> GetParticleIdsOfSuperCluster(const std::vector<PersistentClusterMeta>& pClusterMeta, const SuperClusterMeta& scMeta) {
	std::array<int, 16> particleIds{};
	for (auto& e : particleIds) { e = -1; }
	int cnt = 0;
	for (auto pcId : scMeta.pclusterIds) {
		if (pcId == -1) break;
		for (int particleId : pClusterMeta[pcId].particleIdsGlobal) {
			particleIds[cnt++] = particleId;
		}
	}
	return particleIds;
}

struct ReservedTask {
	int queryScId;
	//bool areBonded;
	int resultIndexRelativeSelf = -1;
	int resultIndexRelativeQuery = -1;
	int nointeractionMatrixIndexRelative = -1;
};

template <typename T>
std::vector<size_t> ExlusivePrefixsum(const std::vector<T>& counts) {
	static_assert(std::is_integral<T>::value, "ExlusivePrefixsum only supports integral types");
	std::vector<size_t> prefixsum(counts.size());
	//std::exclusive_scan(std::execution::par, counts.begin(), counts.end(), prefixsum.begin(), 0);
	std::exclusive_scan(counts.begin(), counts.end(), prefixsum.begin(), size_t{ 0 });
    return prefixsum;
}
template <typename T>
std::vector<size_t> ExlusivePrefixsum(const std::vector<std::vector<T>>& sizes) {
	std::vector<size_t> prefixsum(sizes.size());
	std::transform_exclusive_scan(std::execution::par, sizes.begin(), sizes.end(), prefixsum.begin(), size_t{ 0 }, std::plus<>{},
		[](const std::vector<T>& v) { return v.size(); }
	);
	return prefixsum;
}

std::vector<BoolMatrix16x16> BuildNointeractionMatrices(const std::vector<SuperClusterMeta>& superClusterMetas, const std::vector<int>& nBondedmatricesReserved, const std::vector<size_t>& nBondedMatricesPrefixsum, 
	const std::vector<PersistentClusterMeta>& pClustersMeta, const std::vector<std::vector<ReservedTask>>& workPerSc, const std::vector<std::set<int>>& particleBondedToParticle) {
	const size_t numNointeractionMatricesTotal = nBondedMatricesPrefixsum.back() + nBondedmatricesReserved.back();
	std::vector<BoolMatrix16x16> nointeractionMatrices(numNointeractionMatricesTotal);

	// Build all the nointeractionMatrices
	for (int scId = 0; scId < superClusterMetas.size(); ++scId) {
		std::array<int, 16> particleIdsSelf = GetParticleIdsOfSuperCluster(pClustersMeta, superClusterMetas[scId]);
		for (int i = 0; i < workPerSc[scId].size(); i++) {
			if (workPerSc[scId][i].nointeractionMatrixIndexRelative == -1)
				continue;

			const int queryScId = workPerSc[scId][i].queryScId;
			BoolMatrix16x16 nointeractionMatrix{};

			std::array<int, 16> particleIdsQuery = GetParticleIdsOfSuperCluster(pClustersMeta, superClusterMetas[queryScId]);
			const bool isSelfInteractionTask = scId == queryScId;

			for (int col = 0; col < 16; ++col) {
				for (int row = 0; row < 16; ++row) {
					if (particleIdsSelf[row] == -1)
						continue;
					int pidSelf = particleIdsSelf[row];
					int pidQuery = particleIdsQuery[col];
					if (pidSelf == pidQuery && pidSelf == 0)
						int a = 0;

					bool noInteraction = particleBondedToParticle[particleIdsSelf[row]].contains(particleIdsQuery[col]);
					if (isSelfInteractionTask && row == col) {
						noInteraction = true;
					}
			//		noInteraction = true;
					nointeractionMatrix.Set(row, col, noInteraction);
				}
			}
			

			const int matrixIndex = workPerSc[scId][i].nointeractionMatrixIndexRelative + nBondedMatricesPrefixsum[scId];
			nointeractionMatrices[matrixIndex] = nointeractionMatrix;
		}
	}

	return nointeractionMatrices;
}

//float MinDistanceBetweenPclustersInSupercluster(const SuperCluster& sc0, const SuperCluster& sc1, const Float3& boxSize) {
//	float minDist = FLT_MAX;
//	for (int pcid0 = 0; pcid0 < 4; pcid0++) {
//		if (!sc0.pData->Valid())
//			continue;
//		for (int pcid1 = 0; pcid1 < 4; pcid1++) {
//			if (!sc1.pData->Valid())
//				continue;
//			const float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(sc0.pData[pcid0].position, sc1.pData[pcid1].position, boxSize, BoundaryConditionSelect::PBC);
//			if (dist < minDist) {
//				minDist = dist;
//			}
//		}
//	}
//	return minDist;
//}

std::vector<std::array<float4, 4>> ComputeMeanposAndRadiiForEachPclusterInEachSupercluster(const std::vector<SuperCluster>& superclusters) {
	std::vector<std::array<float4, 4>> out(superclusters.size());

	// Debugging
	float maxRadius = 0;
	float maxIntraScDistance = 0;

	for (int scId = 0; scId < superclusters.size(); scId++) {
		for (int pcid = 0; pcid < 4; pcid++) {
			Float3 sum{};
			int cnt = 0;
			for (int pid = 0; pid < 4; pid++){
				const PData& pData = superclusters[scId].pData[pcid * 4 +pid];
				if (pData.Valid()) {
					sum += pData.position;
					cnt++;
				}
			}

			const Float3 meanPos = sum * (1.0f / static_cast<float>(cnt));
			float radius = 0;
			for (int pid = 0; pid < cnt; pid++) {
				const PData& pData = superclusters[scId].pData[pcid * 4 + pid];
				radius = std::max(radius, (pData.position - meanPos).len());
			}
			out[scId][pcid] = float4{ meanPos.x, meanPos.y, meanPos.z, radius };	

			// Debug
			maxRadius = std::max(maxRadius, radius);
			if (pcid != 0)
				maxIntraScDistance = std::max(maxIntraScDistance, (meanPos - Float3{ out[scId][pcid - 1] }).len());
			if (radius > .8f || maxIntraScDistance > 1.2f)
				int a = 0;
			//
		}
	}

	return out;
}

bool DoesSuperclustersInteract(const std::vector<std::array<float4, 4>>& superclusterPositionSpheres, int scId0, int scId1, float cutoffDistance, Float3 boxSize) {
	for (int pcid0 = 0; pcid0 < 4; pcid0++) {
		for (int pcid1 = 0; pcid1 < 4; pcid1++) {
			float4 p0 = superclusterPositionSpheres[scId0][pcid0];
			float4 p1 = superclusterPositionSpheres[scId1][pcid1];

			float distance = LIMAPOSITIONSYSTEM::calcHyperDistNM(Float3{ p0 }, Float3{ p1 }, boxSize, BoundaryConditionSelect::PBC);
			float radiusSum = p0.w + p1.w;

			if (distance + radiusSum <= cutoffDistance) {	// optim use LenSq
				return true;
			}
		}
	}
	//return true;
	return false;
}

bool Engine::MakeSuperClusterTasksCPU() {
	if (nSuperclusters == 0)
		return true;

	const Box& box = *simulation->box_host;
	Float3 boxSizeF = simulation->box_host->boxparams.BoxSizeFloat();

	//simulation->box_host->persistentClusters;
	const std::vector<PersistentClusterMeta>& pClustersMeta = simulation->box_host->persistentClustersMetadata;
	//const int nSuperClusters = GenericCopyToHost(superClustersControl->nSuperclustersAtomic);
	const std::vector<SuperCluster> superClusters = GenericCopyToHost(superClustersControl->scData, nSuperclusters);
	std::vector<SuperClusterMeta> superClusterMetas = GenericCopyToHost(superClustersControl->scMeta, nSuperclusters);
	
	//// Debug
	//std::vector<int> pclustersMissing(pClustersMeta.size(), 1);
	//std::vector<int> particlesMissing(box.boxparams.total_particles, 1);
	//for (const auto& scm : superClusterMetas) {
	//	for (int pid : scm.particlesIds)
	//		if (pid != -1)
	//			particlesMissing[pid] = 0;
	//	for (int pcid : scm.pclusterIds)
	//		if (pcid != -1)
	//			pclustersMissing[pcid] = 0;
	//}
	//const int nMissingPclusters = std::accumulate(pclustersMissing.begin(), pclustersMissing.end(), 0);
	//const int nMissingParticles = std::accumulate(particlesMissing.begin(), particlesMissing.end(), 0);
	//if (nMissingPclusters > 0 || nMissingParticles > 0) {
	//	int a = 0;
	//}


	const std::vector<std::array<float4, 4>> superclusterPositionSpheres = ComputeMeanposAndRadiiForEachPclusterInEachSupercluster(superClusters);


	std::vector<std::vector<ReservedTask>> workPerSc(superClusters.size());
	std::vector<int> nResultsReserved(superClusters.size(), 0);
	std::vector<int> nBondedmatricesReserved(superClusters.size(), 0);

	for (int scId = 0; scId < superClusterMetas.size(); ++scId) {
		for (int queryScId = scId; queryScId < superClusterMetas.size(); ++queryScId) {

			if (DoesSuperclustersInteract(superclusterPositionSpheres, scId, queryScId, simulation->simparams_host.cutoff_nm, boxSizeF)){
				const bool useNointeractionMatrix = scId == queryScId || ScAreBonded(superClusterMetas[scId], superClusterMetas[queryScId], box.pclusterBondedToPcluster);

				workPerSc[scId].emplace_back(ReservedTask{
					queryScId,
					nResultsReserved[scId],
					scId != queryScId ? nResultsReserved[queryScId] : nResultsReserved[queryScId],
					useNointeractionMatrix ? nBondedmatricesReserved[scId] : -1
					});

				nResultsReserved[scId]++;
				if (scId != queryScId)
					nResultsReserved[queryScId]++;				
				if (useNointeractionMatrix) {
					nBondedmatricesReserved[scId]++;
				}
			}
		}
	}

	// Make prefixsums
	const std::vector<size_t> nResultsPrefixsum = ExlusivePrefixsum(nResultsReserved);
	const std::vector<size_t> nBondedMatricesPrefixsum = ExlusivePrefixsum(nBondedmatricesReserved);
	const std::vector<size_t> nTasksPrefixsum = ExlusivePrefixsum(workPerSc);

	const size_t numTasksTotal = nTasksPrefixsum.back() + workPerSc.back().size();
	std::vector<ScScTask> tasks(numTasksTotal);	

	// Build all the tasks and update the scMeta
	for (int scId = 0; scId < superClusterMetas.size(); ++scId) {
		for (int i = 0; i < workPerSc[scId].size(); i++) {
			//const bool bondedTask = workPerSc[scId][i].areBonded;
			ScScTask task;
			task.nointeractionMatrixIndex = workPerSc[scId][i].nointeractionMatrixIndexRelative != -1 ? workPerSc[scId][i].nointeractionMatrixIndexRelative + nBondedMatricesPrefixsum[scId] : -1;
			task.scIds[0] = scId;
			task.scIds[1] = workPerSc[scId][i].queryScId;
			task.resultIndices[0] = workPerSc[scId][i].resultIndexRelativeSelf + nResultsPrefixsum[scId];
			task.resultIndices[1] = scId != workPerSc[scId][i].queryScId ? (workPerSc[scId][i].resultIndexRelativeQuery + nResultsPrefixsum[workPerSc[scId][i].queryScId]) : -1;
			//task.resultIndices[1] = (workPerSc[scId][i].resultIndexRelativeQuery + nResultsPrefixsum[workPerSc[scId][i].queryScId]);
			tasks[nTasksPrefixsum[scId] + i] = task;
		}

		superClusterMetas[scId].resultsStartIndex = nResultsPrefixsum[scId];
		superClusterMetas[scId].nResults = nResultsReserved[scId];
	}

	//const size_t numNointeractionMatricesTotal = nBondedMatricesPrefixsum.back() + nBondedmatricesReserved.back();
	//nointeractionMatrices.resize(numNointeractionMatricesTotal);



	// Build all the nointeractionMatrices
	const std::vector<BoolMatrix16x16> nointeractionMatrices = BuildNointeractionMatrices(superClusterMetas, nBondedmatricesReserved, nBondedMatricesPrefixsum, pClustersMeta, workPerSc, box.particleBondedToParticle);
	//for (const auto& mat : nointeractionMatrices) {
	//	mat.Print();
	//}
	//{
	//	std::vector<std::set<int>> expectedLjInteractions(16);
	//	for (int row = 0; row < 16; row++) {
	//		for (int col = 0; col < 16; col++) {
	//			if (col == 8 && row == 8)
	//				int aa = 0;
	//			auto _row = nointeractionMatrices[0].GetRow(row);
	//			if (!nointeractionMatrices[0].Get(_row, col)) {
	//				int pid0 = superClusterMetas[0].particlesIds[row];
	//				int pid1 = superClusterMetas[0].particlesIds[col];
	//				expectedLjInteractions[pid0].insert(pid1);
	//			}
	//		}
	//	}
	//	for (int pid = 0; pid < 16; pid++) {
	//		for (auto& interactPid : expectedLjInteractions[pid]) {
	//			printf("%d ", interactPid);
	//		}
	//		printf("\n");
	//	}
	//}


	//DebugUtils::VerifyIdentical(tasks, "ScScTasks" + std::to_string(simulation->getStep()));
	//DebugUtils::VerifyIdentical


	// Push back to device
	cudaMemcpy(superClustersControl->scMeta, superClusterMetas.data(), superClusterMetas.size() * sizeof(SuperClusterMeta), cudaMemcpyHostToDevice);
	cudaFree(scscTasksDevice);
	cudaFree(noInteractionMatricesDevice);	
	scscTasksDevice = GenericCopyToDevice(tasks);
	noInteractionMatricesDevice = GenericCopyToDevice(nointeractionMatrices);

	nResults = nResultsPrefixsum.back() + nResultsReserved.back();
	cudaFree(scResultsDevice);
	cudaMalloc(&scResultsDevice, sizeof(SCResult)* nResults);
	cudaMemset(scResultsDevice, 0, sizeof(SCResult)* nResults);

	nTasks = numTasksTotal;
	if (simulation->getStep() == 787) {
		int a = 0;
	}
	return true;
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