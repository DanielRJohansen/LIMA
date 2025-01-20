#include "Engine.cuh"

#include "EngineBodies.cuh"
#include "Neighborlists.cuh"
#include "BoundaryCondition.cuh"
#include "SimulationDevice.cuh"
#include "LimaPositionSystem.cuh"
#include "EngineKernels.cuh"
#include "Thermostat.cuh"
#include "SupernaturalForces.cuh"
#include "PME.cuh"

#include "Statistics.h"
#include "Utilities.h"

#include "EngineHostside.h"


Engine::Engine(std::unique_ptr<Simulation> _sim, BoundaryConditionSelect bc, std::unique_ptr<LimaLogger> logger)
	: bc_select(bc), m_logger(std::move(logger)), compoundForceEnergyInterims(_sim->box_host->boxparams.n_compounds)
{
	simulation = std::move(_sim);

	verifyEngine();

	dataBuffersDevice = std::make_unique<DatabuffersDeviceController>(simulation->box_host->boxparams.total_particles_upperbound, 
		simulation->box_host->boxparams.n_compounds, simulation->simparams_host.data_logging_interval);


	// Create the Sim_dev {
	{
		if (sim_dev != nullptr) { throw std::runtime_error("Expected simdev to be null to move sim to device"); }
		sim_dev = new SimulationDevice(simulation->simparams_host, simulation->box_host.get(), BoxConfig::Create(*simulation->box_host), BoxState::Create(*simulation->box_host), *dataBuffersDevice);
		sim_dev = genericMoveToDevice(sim_dev, 1);
	}
	setDeviceConstantMemory();
	boxStateCopy = std::make_unique<BoxState>(nullptr, nullptr, nullptr, nullptr, nullptr);
	boxConfigCopy = std::make_unique<BoxConfig>(nullptr, nullptr, nullptr, nullptr, nullptr);
	cudaMemcpy(boxStateCopy.get(), sim_dev->boxState, sizeof(BoxState), cudaMemcpyDeviceToHost);
	cudaMemcpy(boxConfigCopy.get(), &sim_dev->boxConfig, sizeof(BoxConfig), cudaMemcpyDeviceToHost);
	neighborlistsPtr = sim_dev->compound_neighborlists;
	compoundgridPtr = sim_dev->compound_grid;



    std::vector<ForceField_NB::ParticleParameters> compoundParticleParams(simulation->box_host->boxparams.n_compounds * MAX_COMPOUND_PARTICLES, ForceField_NB::ParticleParameters{0,0});
    for (int cid = 0; cid < simulation->box_host->compounds.size(); cid++) {
        const Compound& compound = simulation->box_host->compounds[cid];
        for (int pid = 0; pid < compound.n_particles; pid++) {
            compoundParticleParams[cid*MAX_COMPOUND_PARTICLES + pid] = simulation->forcefield.particle_parameters[compound.atom_types[pid]];
        }

    }
    compoundLjParameters = GenericCopyToDevice(compoundParticleParams);

	cudaStreamCreate(&masterStream);
	cudaStreamCreate(&pmeStream);

	pmeController = std::make_unique<PME::Controller>(simulation->box_host->boxparams.boxSize, *simulation->box_host, simulation->simparams_host.cutoff_nm, pmeStream);
	cudaMalloc(&forceEnergiesPME, sizeof(ForceEnergy) * simulation->box_host->boxparams.n_compounds * MAX_COMPOUND_PARTICLES); // TODO: make cudaFree ...
	cudaMemset(forceEnergiesPME, 0, sizeof(ForceEnergy) * simulation->box_host->boxparams.n_compounds * MAX_COMPOUND_PARTICLES);

	bondgroups = GenericCopyToDevice(simulation->box_host->bondgroups);
	cudaMalloc(&forceEnergiesBondgroups, sizeof(ForceEnergy) * simulation->box_host->bondgroups.size() * BondGroup::maxParticles);






	auto boxparams = simulation->box_host->boxparams;
	thermostat = std::make_unique<Thermostat>(boxparams.n_compounds, boxparams.n_solvents, boxparams.total_particles_upperbound);

	// To create the NLists we need to bootstrap the traj_buffer, since it has no data yet
	bootstrapTrajbufferWithCoords();

	NeighborLists::updateNlists(sim_dev, simulation->getStep(), simulation->simparams_host.bc_select, simulation->box_host->boxparams, timings.nlist);

	// Create kernels graphs. We need for for the steps where we log data, and 1 for all other
	int a = 0;
	cudaMemcpyToSymbol(step, &a, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(transferOutThisStep, &a, sizeof(bool), 0, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	BuildKernelgraph(kernelsGraph_log, true);
	auto err1 = cudaGraphInstantiate(&kernelsGraphExec_log, kernelsGraph_log, nullptr, nullptr, 0);
	LIMA_UTILS::genericErrorCheck("2");
	BuildKernelgraph(kernelsGraph_nolog, false);
	auto err2 = cudaGraphInstantiate(&kernelsGraphExec_nolog, kernelsGraph_nolog, nullptr, nullptr, 0);
	LIMA_UTILS::genericErrorCheck("4");


	m_logger->finishSection("Engine Ready");
}

Engine::~Engine() {
	if (sim_dev != nullptr) {
		sim_dev->FreeMembers();
		cudaFree(sim_dev);
	}
	compoundForceEnergyInterims.Free();

	cudaFree(compoundLjParameters);
	cudaFree(forceEnergiesPME);
	cudaFree(forceEnergiesBondgroups);
	cudaFree(bondgroups);



	LIMA_UTILS::genericErrorCheck("Error during Engine destruction");
	assert(simulation == nullptr);
}


void Engine::setDeviceConstantMemory() {
	//const int forcefield_bytes = sizeof(ForceField_NB);
	cudaMemcpyToSymbol(DeviceConstants::forcefield, &simulation->forcefield, sizeof(ForceField_NB), 0, cudaMemcpyHostToDevice);	// So there should not be a & before the device __constant__
	cudaMemcpyToSymbol(DeviceConstants::tinymolForcefield, &simulation->forcefieldTinymol, sizeof(ForcefieldTinymol), 0, cudaMemcpyHostToDevice);

	BoxSize boxSize_host;
	boxSize_host.Set(simulation->box_host->boxparams.boxSize);
	cudaMemcpyToSymbol(DeviceConstants::boxSize, &boxSize_host, sizeof(BoxSize), 0, cudaMemcpyHostToDevice);
	//SetConstantMem(simulation->boxparams_host.boxSize);
	//BoxSize bs;
	//cudaMemcpyFromSymbol(&bs, DeviceConstants::boxSize, sizeof(BoxSize));

//	cudaDeviceSynchronize();


	//BoxSize bs1;
	//cudaMemcpyFromSymbol(&bs1, DeviceConstants::boxSize, sizeof(BoxSize));

	cudaMemcpyToSymbol(DeviceConstants::cutoffNM, &simulation->simparams_host.cutoff_nm, sizeof(float), 0, cudaMemcpyHostToDevice);
	const float cutoffNmReciprocal = 1.f / simulation->simparams_host.cutoff_nm;
	cudaMemcpyToSymbol(DeviceConstants::cutoffNmReciprocal, &cutoffNmReciprocal, sizeof(float), 0, cudaMemcpyHostToDevice);
	const float cutoffNmSquaredReciprocal = 1.f / (simulation->simparams_host.cutoff_nm * simulation->simparams_host.cutoff_nm );
	cudaMemcpyToSymbol(DeviceConstants::cutoffNmSquaredReciprocal, &cutoffNmSquaredReciprocal, sizeof(float), 0, cudaMemcpyHostToDevice);	
	const float ewaldKappa = PhysicsUtils::CalcEwaldkappa(simulation->simparams_host.cutoff_nm);
	cudaMemcpyToSymbol(DeviceConstants::ewaldKappa, &ewaldKappa, sizeof(float), 0, cudaMemcpyHostToDevice);

	const float initialThermostatScalar = 1.f;
	cudaMemcpyToSymbol(DeviceConstants::thermostatScalar, &initialThermostatScalar, sizeof(float), 0, cudaMemcpyHostToDevice);

	assert(simulation->forcefieldTest.size() == ForceField_NB::MAX_TYPES * ForceField_NB::MAX_TYPES);
	cudaMemcpyToSymbol(DeviceConstants::nonbondedinteractionParams, simulation->forcefieldTest.data(), sizeof(NonbondedInteractionParams) * simulation->forcefieldTest.size(), 0, cudaMemcpyHostToDevice);

	// Prepare precomputed values on device
	const float cutoffNM = simulation->simparams_host.cutoff_nm;
	cudaMemcpyToSymbol(DeviceConstants::bsplineTable, PrecomputeBsplineTable().data(), sizeof(float) * PrecomputeBsplineTable().size(), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(DeviceConstants::erfcForcescalarTable, PrecomputeErfcForcescalarTable(cutoffNM).data(), sizeof(float) * PrecomputeErfcForcescalarTable(cutoffNM).size(), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(DeviceConstants::erfcPotentialscalarTable, PrecomputeErfcPotentialscalarTable(cutoffNM).data(), sizeof(float) * PrecomputeErfcPotentialscalarTable(cutoffNM).size(), 0, cudaMemcpyHostToDevice);

	LIMA_UTILS::genericErrorCheck("Error while setting CUDA __constant__ memory\n");
}




void Engine::Step() {
	LIMA_UTILS::genericErrorCheckNoSync("Error before step!");

	deviceMaster();	// Device first, otherwise offloading data always needs the last datapoint!
	assert(simulation);
	assert(sim_dev);
	simulation->step++;

	hostMaster();

	LIMA_UTILS::genericErrorCheckNoSync("Error after step!");
}

void Engine::hostMaster() {						// This is and MUST ALWAYS be called after the deviceMaster, and AFTER incStep()!
	auto t0 = std::chrono::high_resolution_clock::now();
	if (DatabuffersDeviceController::IsBufferFull(simulation->getStep(), simulation->simparams_host.data_logging_interval)) {
		offloadLoggingData(DatabuffersDeviceController::nStepsInBuffer);
		runstatus.stepForMostRecentData = simulation->getStep();

		if ((simulation->getStep() % simulation->simparams_host.steps_per_temperature_measurement) == 0 && simulation->getStep() > 0) {
			auto [temperature, thermostatScalar] = thermostat->Temperature(sim_dev, simulation->box_host->boxparams, simulation->simparams_host);
			simulation->temperature_buffer.push_back(temperature);
			runstatus.current_temperature = temperature;

			if (simulation->simparams_host.apply_thermostat)
				cudaMemcpyToSymbol(DeviceConstants::thermostatScalar, &thermostatScalar, sizeof(float), 0, cudaMemcpyHostToDevice);
		}
		
		HandleEarlyStoppingInEM();

		NeighborLists::updateNlists(sim_dev, simulation->getStep(), simulation->simparams_host.bc_select, simulation->box_host->boxparams, timings.nlist);
	}

	// Handle status
	runstatus.current_step = simulation->getStep();
	if (runstatus.current_step >= simulation->simparams_host.n_steps || runstatus.critical_error_occured)
		runstatus.simulation_finished = true;


	const auto t1 = std::chrono::high_resolution_clock::now();
	const int cpu_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	timings.cpu_master += cpu_duration;
}

void Engine::terminateSimulation() {
	const int64_t stepsReadyToTransfer = DatabuffersDeviceController::StepsReadyToTransfer(simulation->getStep(), simulation->simparams_host.data_logging_interval);
	offloadLoggingData(stepsReadyToTransfer);

	sim_dev->boxState->CopyDataToHost(*simulation->box_host);

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

	step_at_last_traj_transfer = 0.f;
	runstatus.most_recent_positions = simulation->traj_buffer->getBufferAtIndex(0);

	LIMA_UTILS::genericErrorCheck("Error during bootstrapTrajbufferWithCoords");
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


void Engine::deviceMaster() {
	const BoxParams& boxparams = simulation->box_host->boxparams;
	//const int step = simulation->getStep();
	// #### Initial round of force computations
	cudaDeviceSynchronize();

	const bool logData = simulation->getStep() % simulation->simparams_host.data_logging_interval == 0;// TODO maybe log at the final step, not 0th?

	// Overwrite the "step" argument for all kernels, somehow????
	const int _step = simulation->getStep();
	const bool isTransferStep = SolventBlocksCircularQueue::isTransferStep(step);

	cudaMemcpyToSymbol(step, &_step, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(transferOutThisStep, &isTransferStep, sizeof(bool), 0, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	cudaError_t err;

	if (logData)
		err = cudaGraphLaunch(kernelsGraphExec_log, masterStream);
	else
		err = cudaGraphLaunch(kernelsGraphExec_nolog, masterStream);

	if (ENABLE_ES_LR && simulation->simparams_host.enable_electrostatics) {
		pmeController->CalcCharges(*boxConfigCopy, *boxStateCopy, boxparams.n_compounds, forceEnergiesPME, pmeStream);
	}

	cudaStreamSynchronize(masterStream);
	cudaStreamSynchronize(pmeStream);


	if (boxparams.n_compounds > 0) {
		CompoundIntegrationKernel<PeriodicBoundaryCondition, false>
			<< <boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, masterStream >> >
			(sim_dev, compoundForceEnergyInterims, forceEnergiesBondgroups, forceEnergiesPME);
	}

	if (boxparams.n_solvents > 0) {
		TinymolIntegrationLoggingAndTransferout<PeriodicBoundaryCondition, false>
			<< <BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, pmeStream >> >
			(sim_dev);
		solventTransferKernel<PeriodicBoundaryCondition>
			<< <BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlockTransfermodule::max_queue_size, 0, pmeStream >> >
			(sim_dev);
	}

	cudaStreamSynchronize(masterStream);
	cudaStreamSynchronize(pmeStream);


	/*cudaDeviceSynchronize();
	LIMA_UTILS::genericErrorCheckNoSync("Error after step!");*/

}

template <typename BoundaryCondition, bool emvariant, bool computePotE>
void Engine::_BuildKernelgraph(cudaGraph_t& graph) {

	const BoxParams boxparams = simulation->box_host->boxparams;

	cudaGraphCreate(&graph, 0);

	cudaGraphNode_t nodeFarneighbor = nullptr;
	cudaGraphNode_t nodeBonds = nullptr;
	cudaGraphNode_t nodeTinymolCompoundinteractions = nullptr;
	cudaGraphNode_t nodeTinymolForceKernel = nullptr;
	cudaGraphNode_t nodeSnf = nullptr;

	cudaGraphNode_t nodeBarrierSection = nullptr;

	cudaGraphNode_t nodeCompoundIntegration = nullptr;
	cudaGraphNode_t nodeTinymolIntegration = nullptr;
	cudaGraphNode_t nodeSolventTransfer = nullptr;

	cudaGraphNode_t pmeNode1;

	std::vector<cudaGraphNode_t> forceDependencies; // Add all nodes to this vector, so they are sync before moving to integration


	if (boxparams.n_compounds > 0) {
		{
			void* kernelArgs[] = {
				(void*)&(*boxStateCopy),
				(void*)&(*boxConfigCopy),
				(void*)&neighborlistsPtr,
				(void*)&simulation->simparams_host.enable_electrostatics,
				(void*)&compoundForceEnergyInterims.forceEnergyFarneighborShortrange,
				(void*)&compoundLjParameters
			};
			cudaKernelNodeParams kernelnodeParams{
				(void*)compoundFarneighborShortrangeInteractionsKernel<BoundaryCondition, emvariant, computePotE>,
				dim3(boxparams.n_compounds),
				dim3(MAX_COMPOUND_PARTICLES),
				0,
				kernelArgs,
				nullptr
			};

			LIMA_UTILS::genericErrorCheck(cudaGraphAddKernelNode(&nodeFarneighbor, graph, nullptr, 0, &kernelnodeParams));
			forceDependencies.emplace_back(nodeFarneighbor);
		}

		{
			cudaGraphNode_t node;
			void* kernelArgs[] = { &sim_dev, &compoundForceEnergyInterims.forceEnergyImmediateneighborShortrange };
			cudaKernelNodeParams kernelnodeParams = {
				(void*)compoundImmediateneighborAndSelfShortrangeInteractionsKernel<BoundaryCondition, emvariant, computePotE>,
				boxparams.n_compounds,
				MAX_COMPOUND_PARTICLES,
				0,
				kernelArgs,
				nullptr 
			};
						
			LIMA_UTILS::genericErrorCheck(cudaGraphAddKernelNode(&node, graph, nullptr, 0, &kernelnodeParams));
			forceDependencies.emplace_back(node);
		}
	}


	if (!simulation->box_host->bondgroups.empty()) {
		void* kernelArgs[] = { &bondgroups, &(*boxStateCopy), &forceEnergiesBondgroups };
		cudaKernelNodeParams kernelnodeParams = {
			(void*)BondgroupsKernel<BoundaryCondition, emvariant>,
			simulation->box_host->bondgroups.size(),
			THREADS_PER_BONDSGROUPSKERNEL,
			0,
			kernelArgs,
			nullptr
		};
		LIMA_UTILS::genericErrorCheck(cudaGraphAddKernelNode(&nodeBonds, graph, nullptr, 0, &kernelnodeParams));
		forceDependencies.emplace_back(nodeBonds);
	}

	if (boxparams.n_solvents > 0) {
		{
			void* kernelArgs[] = {
				(void*)&(*boxStateCopy),
				(void*)&(*boxConfigCopy),
				&compoundgridPtr
			};
			cudaKernelNodeParams kernelnodeParams = {
				(void*)TinymolCompoundinteractionsKernel<BoundaryCondition, emvariant>,
				BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)),
				SolventBlock::MAX_SOLVENTS_IN_BLOCK,
				0,
				kernelArgs,
				nullptr
			};
			LIMA_UTILS::genericErrorCheck(cudaGraphAddKernelNode(&nodeTinymolCompoundinteractions, graph, nullptr, 0, &kernelnodeParams));
			forceDependencies.emplace_back(nodeTinymolCompoundinteractions);
		}
		{
			void* kernelArgs[] = { (void*)&(*boxStateCopy), (void*)&(*boxConfigCopy) };
			cudaKernelNodeParams kernelnodeParams = {
				(void*)solventForceKernel<BoundaryCondition, emvariant>,
				BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)),
				SolventBlock::MAX_SOLVENTS_IN_BLOCK,
				0,
				kernelArgs,
				nullptr
			};
			LIMA_UTILS::genericErrorCheck(cudaGraphAddKernelNode(&nodeTinymolForceKernel, graph, nullptr, 0, &kernelnodeParams));
			forceDependencies.emplace_back(nodeTinymolForceKernel);
		}
	}
	//if (ENABLE_ES_LR && simulation->simparams_host.enable_electrostatics) {
	//	pmeController->AddToGraph(*boxConfigCopy, *boxStateCopy, boxparams.n_compounds, forceEnergiesPME, graph, forceDependencies, pmeNode1);
	//}
	switch (simulation->simparams_host.snf_select)
	{

	case None:
		break;
	case HorizontalSqueeze:	{
		//SupernaturalForces::ApplyHorizontalSqueeze << < simulation->box_host->boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, stream >> > (sim_dev, simulation->getStep());
		throw std::runtime_error("Not implemented");

		//void* kernelArgs[] = { &sim_dev, &compoundForceEnergyInterims.forceEnergyBonds };
		//cudaKernelNodeParams kernelnodeParams = {
		//	(void*)SupernaturalForces::ApplyHorizontalSqueeze
		//	boxparams.n_compounds,
		//	MAX_COMPOUND_PARTICLES,
		//	0,
		//	kernelArgs,
		//	nullptr
		//};
		//LIMA_UTILS::genericErrorCheck(cudaGraphAddKernelNode(&nodeSnf, graph, nullptr, 0, &kernelnodeParams));
		//forceDependencies.emplace_back(nodeSnf);
		break;
	}		
	case HorizontalChargeField: {
		void* kernelArgs[] = { &sim_dev, &simulation->box_host->uniformElectricField, &compoundForceEnergyInterims.forceEnergyBonds };
		cudaKernelNodeParams kernelnodeParams = {
			(void*)CompoundSnfKernel<BoundaryCondition, emvariant>,
			boxparams.n_compounds,
			MAX_COMPOUND_PARTICLES,
			0,
			kernelArgs,
			nullptr
		};
		LIMA_UTILS::genericErrorCheck(cudaGraphAddKernelNode(&nodeSnf, graph, nullptr, 0, &kernelnodeParams));
		forceDependencies.emplace_back(nodeSnf);
		break;

		//CompoundSnfKernel<BoundaryCondition, emvariant>
		//	<< <simulation->box_host->boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, stream >> >
		//	(sim_dev, simulation->box_host->uniformElectricField, compoundForceEnergyInterims.forceEnergyBonds);
	}

	case BoxEdgePotential: {
		throw std::runtime_error("Not implemented");
		/*if (simulation->box_host->boxparams.n_compounds > 0)
			SupernaturalForces::BoxEdgeForceCompounds << < simulation->box_host->boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, stream >> > (sim_dev, simulation->getStep());
		if (simulation->box_host->boxparams.n_solvents > 0)
			SupernaturalForces::BoxEdgeForceSolvents << <BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(simulation->box_host->boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, stream >> > (sim_dev, simulation->getStep());
		break;*/
		break;
	}
		
	default:
		break;
	}













	// Sync all prev nodes
	{
		auto err = cudaGraphAddEmptyNode(&nodeBarrierSection, graph, forceDependencies.data(), forceDependencies.size());
		LIMA_UTILS::genericErrorCheck("e");
	}

	//if (boxparams.n_compounds > 0) {
	//	void* compoundIntegrationArgs[] = {
	//		&sim_dev,
	//		&compoundForceEnergyInterims,
	//		&forceEnergiesBondgroups,
	//		&forceEnergiesPME
	//	};
	//	cudaKernelNodeParams compoundIntegrationParams = {
	//		(void*)CompoundIntegrationKernel<BoundaryCondition, emvariant>,
	//		dim3(boxparams.n_compounds),
	//		dim3(MAX_COMPOUND_PARTICLES),
	//		0,
	//		compoundIntegrationArgs,
	//		nullptr
	//	};
	//	cudaGraphAddKernelNode(&nodeCompoundIntegration, graph, &nodeBarrierSection, 1, &compoundIntegrationParams);
	//	LIMA_UTILS::genericErrorCheck("e");
	//}

	//if (boxparams.n_solvents > 0) {
	//	{
	//		void* tinymolIntegrationArgs[] = { &sim_dev };
	//		cudaKernelNodeParams tinymolIntegrationParams = {
	//			(void*)TinymolIntegrationLoggingAndTransferout<BoundaryCondition, emvariant>,
	//			dim3(BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize))),
	//			dim3(SolventBlock::MAX_SOLVENTS_IN_BLOCK),
	//			0,
	//			tinymolIntegrationArgs,
	//			nullptr
	//		};
	//		cudaGraphAddKernelNode(&nodeTinymolIntegration, graph, &nodeBarrierSection, 1, &tinymolIntegrationParams);
	//		LIMA_UTILS::genericErrorCheck("e");
	//	}
	//	{
	//		void* solventTransferArgs[] = { &sim_dev };
	//		cudaKernelNodeParams solventTransferParams = {
	//			(void*)solventTransferKernel<BoundaryCondition>,
	//			dim3(BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize))),
	//			dim3(SolventBlockTransfermodule::max_queue_size),
	//			0,
	//			solventTransferArgs,
	//			nullptr
	//		};
	//		cudaGraphAddKernelNode(&nodeSolventTransfer, graph, &nodeTinymolIntegration, 1, &solventTransferParams);
	//		LIMA_UTILS::genericErrorCheck("e");
	//	}
	//}
}
template void Engine::_BuildKernelgraph<PeriodicBoundaryCondition, true, true>(cudaGraph_t&);
template void Engine::_BuildKernelgraph<PeriodicBoundaryCondition, true, false>(cudaGraph_t&);
template void Engine::_BuildKernelgraph<PeriodicBoundaryCondition, false, true>(cudaGraph_t&);
template void Engine::_BuildKernelgraph<PeriodicBoundaryCondition, false, false>(cudaGraph_t&);
template void Engine::_BuildKernelgraph<NoBoundaryCondition, true, true>(cudaGraph_t&);
template void Engine::_BuildKernelgraph<NoBoundaryCondition, true, false>(cudaGraph_t&);
template void Engine::_BuildKernelgraph<NoBoundaryCondition, false, true>(cudaGraph_t&);
template void Engine::_BuildKernelgraph<NoBoundaryCondition, false, false>(cudaGraph_t&);



void Engine::BuildKernelgraph(cudaGraph_t& graph, bool logData) {

	switch (simulation->simparams_host.bc_select) {
	case NoBC:
		if (simulation->simparams_host.em_variant) {
			if (logData) {
				_BuildKernelgraph<NoBoundaryCondition, true, true>(graph);
			}
			else {
				_BuildKernelgraph<NoBoundaryCondition, true, false>(graph);
			}
		}
		else {
			if (logData) {
				_BuildKernelgraph<NoBoundaryCondition, false, true>(graph);
			}
			else {
				_BuildKernelgraph<NoBoundaryCondition, false, false>(graph);
			}
		}
		break;
	case PBC:
		if (simulation->simparams_host.em_variant) {
			if (logData) {
				_BuildKernelgraph<PeriodicBoundaryCondition, true, true>(graph);
			}
			else {
				_BuildKernelgraph<PeriodicBoundaryCondition, true, false>(graph);
			}
		}
		else {
			if (logData) {
				_BuildKernelgraph<PeriodicBoundaryCondition, false, true>(graph);
			}
			else {
				_BuildKernelgraph<PeriodicBoundaryCondition, false, false>(graph);
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
			(sim_dev, simulation->box_host->uniformElectricField, compoundForceEnergyInterims.forceEnergyBonds);
		break;
	case BoxEdgePotential:
		if (simulation->box_host->boxparams.n_compounds > 0)
			SupernaturalForces::BoxEdgeForceCompounds << < simulation->box_host->boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, stream >> > (sim_dev, simulation->getStep());
		if (simulation->box_host->boxparams.n_solvents > 0)
			SupernaturalForces::BoxEdgeForceSolvents<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(simulation->box_host->boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, stream>>>(sim_dev, simulation->getStep());
		break;
	}
}





//cudaKernelNodeParams kFarParams{};
//// For demonstration, we show typical fields:
//kFarParams.func = (void*)compoundFarneighborShortrangeInteractionsKernel<BoundaryCondition, emvariant, computePotE>;
//kFarParams.gridDim = dim3(boxparams.n_compounds);
//kFarParams.blockDim = dim3(MAX_COMPOUND_PARTICLES);
//kFarParams.sharedMemBytes = 0;
//{
//	// Put real arguments in a vector<void*>
//	// Adjust these to match your kernel signature
//	static bool cpe = true; // or false, per your usage
//	void* args[] = {
//		&boxStateCopy,
//		&boxConfigCopy,
//		&neighborlistsPtr,
//		&compoundForceEnergyInterims.forceEnergyFarneighborShortrange,
//		&compoundLjParameters
//	};
//	kFarParams.kernelParams = args;
//}
//kFarParams.extra = nullptr;

//// No dependencies  can be added with a null dependency array
//auto err = cudaGraphAddKernelNode(&nodeFarneighbor, graph, nullptr, 0, &kFarParams);







// ------------------------------------------------------------------
// 3) Add a barrier node so that Section 2 will start only after
//    both kernels in Section 1 have completed.
// ------------------------------------------------------------------
   //{
   //	cudaGraphNode_t dependencies[2] = { nodeFarneighbor, nodeImmediateneighbor };
   //	cudaGraphAddEmptyNode(&nodeBarrierSection1, graph, dependencies, 2);
   //}





   // ------------------------------------------------------------------
  // 4) Section 2
  //    - CompoundIntegrationKernel (if n_compounds > 0)
  //    - TinymolIntegrationLoggingAndTransferout + solventTransferKernel (in sequence)
  // ------------------------------------------------------------------
//if (boxparams.n_compounds > 0)
//{
	//cudaKernelNodeParams compoundIntegrationParams = {};
	//compoundIntegrationParams.func = (void*)CompoundIntegrationKernel<BoundaryCondition, emvariant>;
	//compoundIntegrationParams.gridDim = dim3(boxparams.n_compounds);
	//compoundIntegrationParams.blockDim = dim3(MAX_COMPOUND_PARTICLES);
	//compoundIntegrationParams.sharedMemBytes = 0;
	//void* compoundIntegrationArgs[] = {
	//	&sim_dev,
	//	&compoundForceEnergyInterims,
	//	&forceEnergiesBondgroups,
	//	&forceEnergiesPME
	//};
	//compoundIntegrationParams.kernelParams = compoundIntegrationArgs;
	//compoundIntegrationParams.extra = nullptr;

	//cudaGraphAddKernelNode(&nodeCompoundIntegration, graph, &nodeBarrierSection1, 1, &compoundIntegrationParams);
//}
//
//// ------------------------------------------------------------------
//// 4b) If n_solvents > 0, run 2 kernels in sequence
//// ------------------------------------------------------------------
//if (boxparams.n_solvents > 0)
//{
//	//cudaKernelNodeParams tinymolIntegrationParams = {};
//	//tinymolIntegrationParams.func = (void*)TinymolIntegrationLoggingAndTransferout<BoundaryCondition, emvariant>;
//	//tinymolIntegrationParams.gridDim = dim3(BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)));
//	//tinymolIntegrationParams.blockDim = dim3(SolventBlock::MAX_SOLVENTS_IN_BLOCK);
//	//tinymolIntegrationParams.sharedMemBytes = 0;
//	//void* tinymolIntegrationArgs[] = { &sim_dev };
//	//tinymolIntegrationParams.kernelParams = tinymolIntegrationArgs;
//	//tinymolIntegrationParams.extra = nullptr;
//
//	//cudaGraphAddKernelNode(&nodeTinymolIntegration, graph, &nodeBarrierSection1, 1, &tinymolIntegrationParams);
//
//	//cudaKernelNodeParams solventTransferParams = {};
//	//solventTransferParams.func = (void*)solventTransferKernel<BoundaryCondition>;
//	//solventTransferParams.gridDim = dim3(BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)));
//	//solventTransferParams.blockDim = dim3(SolventBlockTransfermodule::max_queue_size);
//	//solventTransferParams.sharedMemBytes = 0;
//	//void* solventTransferArgs[] = { &sim_dev };
//	//solventTransferParams.kernelParams = solventTransferArgs;
//	//solventTransferParams.extra = nullptr;
//
//	//cudaGraphAddKernelNode(&nodeSolventTransfer, graph, &nodeTinymolIntegration, 1, &solventTransferParams);
//}

// #### Force kernels

//if (boxparams.n_compounds > 0) {
//	compoundFarneighborShortrangeInteractionsKernel<BoundaryCondition, emvariant, computePotE> 
//		<<<boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, cudaStreams[0] >> >
//           (*boxStateCopy, *boxConfigCopy, neighborlistsPtr, simulation->simparams_host.enable_electrostatics, 
   //			compoundForceEnergyInterims.forceEnergyFarneighborShortrange, compoundLjParameters);

   //	/*compoundImmediateneighborAndSelfShortrangeInteractionsKernel<BoundaryCondition, emvariant, computePotE> 
   //		<<<boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, cudaStreams[1] >>> 
   //		(sim_dev, compoundForceEnergyInterims.forceEnergyImmediateneighborShortrange);*/
   //}

   //if (boxparams.n_solvents > 0) {
   //	// Should only use max_compound_particles threads here. and let 1 thread handle multiple solvents
   //	TinymolCompoundinteractionsKernel<BoundaryCondition, emvariant>
   //		<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, cudaStreams[2]>>>
   //		(*boxStateCopy, *boxConfigCopy, compoundgridPtr);

   //	// TODO: Too many threads, we rarely get close to filling the block
   //	solventForceKernel<BoundaryCondition, emvariant> 
   //		<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, cudaStreams[3]>>>
   //		(*boxStateCopy, *boxConfigCopy);
   //}

   //if (ENABLE_ES_LR && simulation->simparams_host.enable_electrostatics) {
   //	pmeController->CalcCharges(*boxConfigCopy, *boxStateCopy, boxparams.n_compounds, forceEnergiesPME, pmeStream);
   //}

   /*if (simulation->simparams_host.snf_select != None) {
	   SnfHandler<BoundaryCondition, emvariant>(cudaStreams[2]);
   }*/

   //if (!simulation->box_host->bondgroups.empty()) {
   //	BondgroupsKernel<BoundaryCondition, emvariant> << < simulation->box_host->bondgroups.size(), THREADS_PER_BONDSGROUPSKERNEL, 0, cudaStreams[4]>>> 
   //		(bondgroups, *boxStateCopy, forceEnergiesBondgroups);
   //}

   //// #### Sync all streams
   //for (int i = 0; i < 5; i++)
   //	cudaEventRecord(streamSync[i], cudaStreams[i]);
   //cudaEventRecord(pmeStreamSync, pmeStream);

   //for (int i = 0; i < 5; i++) {
   //	if (i != 0)
   //		cudaStreamWaitEvent(cudaStreams[0], streamSync[i], 0);
   //	if (i != 1)
   //		cudaStreamWaitEvent(cudaStreams[1], streamSync[i], 0);
   //}
   //cudaStreamWaitEvent(cudaStreams[0], pmeStreamSync, 0);// solvents dont need to wait for pme, as they are only neutral currently
   //cudaStreamWaitEvent(cudaStreams[1], pmeStreamSync, 0);// solvents dont need to wait for pme, as they are only neutral currently

   //// #### Integration and Transfer kernels

   //if (boxparams.n_compounds > 0) {
   //	CompoundIntegrationKernel<BoundaryCondition, emvariant> 
   //		<<<boxparams.n_compounds, MAX_COMPOUND_PARTICLES, 0, cudaStreams[0] >> >
   //		(sim_dev, compoundForceEnergyInterims, forceEnergiesBondgroups, forceEnergiesPME);
   //}

   //if (boxparams.n_solvents > 0) {		
   //	TinymolIntegrationLoggingAndTransferout<BoundaryCondition, emvariant>
   //		<< <BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, cudaStreams[1] >> >
   //		(sim_dev);
   //	solventTransferKernel<BoundaryCondition>
   //		<< <BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlockTransfermodule::max_queue_size, 0, cudaStreams[1] >> >
   //		(sim_dev);
   //}	


   //// End graph capture
   //auto c = cudaStreamEndCapture(masterStream, &graph);

   //for (int i = 0; i < 5; i++)
   //	cudaEventDestroy(streamSync[i]);
   //cudaEventDestroy(pmeStreamSync);

   //for (cudaStream_t& stream : cudaStreams) {
   //	cudaStreamDestroy(stream);
   //}
   //cudaStreamDestroy(pmeStream);