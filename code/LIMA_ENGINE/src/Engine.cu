#include "Engine.cuh"

#include "EngineBodies.cuh"
#include "Neighborlists.cuh"
#include "BoundaryCondition.cuh"
#include "SimulationDevice.cuh"
#include "LimaPositionSystem.cuh"
#include "ChargeOcttree.cuh"
#include "EngineKernels.cuh"
#include "Thermostat.cuh"
#include "SupernaturalForces.cuh"
//#include "PME.cuh"
#include "PMEfast.cuh"

#include "Statistics.h"
#include "Utilities.h"

#include <unordered_set>


// Optim: make sure any ptr sent to a kernel is 128 byte aligned

//const int compound_size = sizeof(CompoundCompact);
//const int nlsit_size = sizeof(NeighborList);
//const int sssize = (sizeof(Float3) + sizeof(float)) * THREADS_PER_COMPOUNDBLOCK;
//const int Ckernel_shared_mem = sizeof(CompoundCompact) + sizeof(NeighborList) +
//	(2* sizeof(Float3)) * THREADS_PER_COMPOUNDBLOCK + sizeof(Coord) + sizeof(Float3) + clj_utilitybuffer_bytes;
//static_assert(Ckernel_shared_mem < 45000, "Not enough shared memory for CompoundKernel");

//const int sbsize = sizeof(SolventBlock);
//const int Skernel_shared_mem = (sizeof(Float3) + 1) * SolventBlock::MAX_SOLVENTS_IN_BLOCK + sizeof(SolventBlock)
//	+ sizeof(SolventTransferqueue<SolventBlockTransfermodule::max_queue_size>) * 6
//	+ 4 + 4 * 3 * 2;

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
    //cudaMalloc(&compoundLjParameters, sizeof(ForceField_NB::ParticleParameters) * MAX_COMPOUND_PARTICLES * simulation->box_host->boxparams.n_compounds);

	for (cudaStream_t& stream : cudaStreams) {
		cudaStreamCreate(&stream);
	}


	pmeController = std::make_unique<PME::Controller>(simulation->box_host->boxparams.boxSize, *simulation->box_host, simulation->simparams_host.cutoff_nm);
	cudaMalloc(&forceEnergiesPME, sizeof(ForceEnergy) * simulation->box_host->boxparams.n_compounds * MAX_COMPOUND_PARTICLES); // TODO: make cudaFree ...
	cudaMemset(forceEnergiesPME, 0, sizeof(ForceEnergy) * simulation->box_host->boxparams.n_compounds * MAX_COMPOUND_PARTICLES);
	//std::unordered_set<std::string> unique_compounds;
	//for (int i = 0; i < simulation->box_host->boxparams.n_compounds; i++) {
	//	char types[64];
	//	memcpy(types, simulation->box_host->compounds[i].atom_types, 64);
	//	unique_compounds.insert(std::string(types));
	//}
	//int a = 0;

	bondgroups = GenericCopyToDevice(simulation->box_host->bondgroups);
	cudaMalloc(&forceEnergiesBondgroups, sizeof(ForceEnergy) * simulation->box_host->bondgroups.size() * BondGroup::maxParticles);






	auto boxparams = simulation->box_host->boxparams;
	thermostat = std::make_unique<Thermostat>(boxparams.n_compounds, boxparams.n_solvents, boxparams.total_particles_upperbound);

	// To create the NLists we need to bootstrap the traj_buffer, since it has no data yet
	bootstrapTrajbufferWithCoords();

	NeighborLists::updateNlists(sim_dev, simulation->getStep(), simulation->simparams_host.bc_select, simulation->box_host->boxparams, timings.nlist);
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

	for (cudaStream_t& stream : cudaStreams) {
		cudaStreamDestroy(stream);
	}

	LIMA_UTILS::genericErrorCheck("Error during Engine destruction");
	assert(simulation == nullptr);
}


void Engine::setDeviceConstantMemory() {
	//const int forcefield_bytes = sizeof(ForceField_NB);
	cudaMemcpyToSymbol(forcefield_device, &simulation->forcefield, sizeof(ForceField_NB), 0, cudaMemcpyHostToDevice);	// So there should not be a & before the device __constant__
	cudaMemcpyToSymbol(tinymolForcefield_device, &simulation->forcefieldTinymol, sizeof(ForcefieldTinymol), 0, cudaMemcpyHostToDevice);

	BoxSize boxSize_host;
	boxSize_host.Set(simulation->box_host->boxparams.boxSize);
	cudaMemcpyToSymbol(boxSize_device, &boxSize_host, sizeof(BoxSize), 0, cudaMemcpyHostToDevice);
	//SetConstantMem(simulation->boxparams_host.boxSize);
	//BoxSize bs;
	//cudaMemcpyFromSymbol(&bs, boxSize_device, sizeof(BoxSize));

//	cudaDeviceSynchronize();


	//BoxSize bs1;
	//cudaMemcpyFromSymbol(&bs1, boxSize_device, sizeof(BoxSize));

	cudaMemcpyToSymbol(cutoffNm_device, &simulation->simparams_host.cutoff_nm, sizeof(float), 0, cudaMemcpyHostToDevice);
	const float cutoffNmSquaredReciprocal = 1.f / (simulation->simparams_host.cutoff_nm * simulation->simparams_host.cutoff_nm );
	cudaMemcpyToSymbol(cutoffNmSquaredReciprocal_device, &cutoffNmSquaredReciprocal, sizeof(float), 0, cudaMemcpyHostToDevice);	
	const float ewaldKappa = PhysicsUtils::CalcEwaldkappa(simulation->simparams_host.cutoff_nm);
	cudaMemcpyToSymbol(ewaldkappa_device, &ewaldKappa, sizeof(float), 0, cudaMemcpyHostToDevice);

	const float initialThermostatScalar = 1.f;
	cudaMemcpyToSymbol(thermostatScalar_device, &initialThermostatScalar, sizeof(float), 0, cudaMemcpyHostToDevice);

	assert(simulation->forcefieldTest.size() == ForceField_NB::MAX_TYPES * ForceField_NB::MAX_TYPES);
	cudaMemcpyToSymbol(nonbondedInteractionParams_device, simulation->forcefieldTest.data(), sizeof(NonbondedInteractionParams) * simulation->forcefieldTest.size(), 0, cudaMemcpyHostToDevice);

	LIMA_UTILS::genericErrorCheck("Error while setting CUDA __constant__ memory\n");
}


std::unique_ptr<Simulation> Engine::takeBackSim() {
	assert(sim_dev);
	sim_dev->boxState->CopyDataToHost(*simulation->box_host);
	return std::move(simulation);
}

void Engine::verifyEngine() {
	LIMA_UTILS::genericErrorCheck("Error before engine initialization.\n");

	const int nBlocks = simulation->box_host->boxparams.boxSize;
	assert(nBlocks* nBlocks* nBlocks < INT32_MAX && "Neighborlist cannot handle such large gridnode_ids");

	if constexpr (ENABLE_ES_LR) {
		if (simulation->simparams_host.enable_electrostatics && simulation->simparams_host.bc_select != PBC) {
			throw std::invalid_argument("Electrostatics only supported with PBC at the current time");
		}
	}
}

void Engine::step() {
	LIMA_UTILS::genericErrorCheck("Error before step!");

	deviceMaster();	// Device first, otherwise offloading data always needs the last datapoint!
	assert(simulation);
	assert(sim_dev);
	simulation->step++;

	hostMaster();

	LIMA_UTILS::genericErrorCheck("Error after step!");
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
				cudaMemcpyToSymbol(thermostatScalar_device, &thermostatScalar, sizeof(float), 0, cudaMemcpyHostToDevice);
		}
		
		HandleEarlyStoppingInEM();

		NeighborLists::updateNlists(sim_dev, simulation->getStep(), simulation->simparams_host.bc_select, simulation->box_host->boxparams, timings.nlist);
	}

	// Handle status
	runstatus.current_step = simulation->getStep();
	runstatus.critical_error_occured = sim_dev->signals->critical_error_encountered;	// TODO: Can i get this from simparams_host? UNSAFE
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

	cudaDeviceSynchronize();
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
template <typename BoundaryCondition, bool emvariant, bool computePotE>
void Engine::_deviceMaster() {
	cudaDeviceSynchronize();
	
	const BoxParams& boxparams = simulation->box_host->boxparams;
	const bool logData = simulation->getStep() % simulation->simparams_host.data_logging_interval == 0;// TODO maybe log at the final step, not 0th?


	// #### Pre force kernels
	if (ENABLE_ES_LR && simulation->simparams_host.enable_electrostatics && boxparams.n_compounds > 0) {
		//Electrostatics::DistributeCompoundchargesToGridKernel<<<boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, 0, cudaStreams[0]>>>(sim_dev);
		LIMA_UTILS::genericErrorCheckNoSync("Error after DistributeCompoundchargesToGridKernel");
	}

	// #### Initial round of force computations
	cudaDeviceSynchronize();
	if (boxparams.n_compounds > 0) {
		compoundFarneighborShortrangeInteractionsKernel<BoundaryCondition, emvariant, computePotE> 
			<<<boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, 0, cudaStreams[0]>>>
            (simulation->getStep(), *boxStateCopy, *boxConfigCopy, neighborlistsPtr, simulation->simparams_host.enable_electrostatics, compoundForceEnergyInterims.forceEnergyFarneighborShortrange, compoundLjParameters);
		LIMA_UTILS::genericErrorCheckNoSync("Error after compoundFarneighborShortrangeInteractionsKernel");

		compoundImmediateneighborAndSelfShortrangeInteractionsKernel<BoundaryCondition, emvariant, computePotE> 
			<<<boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, 0, cudaStreams[1] >>> 
			(sim_dev, simulation->getStep(), compoundForceEnergyInterims.forceEnergyImmediateneighborShortrange);
		LIMA_UTILS::genericErrorCheckNoSync("Error after compoundImmediateneighborAndSelfShortrangeInteractionsKernel");
	}

	if (boxparams.n_solvents > 0) {
		// Should only use max_compound_particles threads here. and let 1 thread handle multiple solvents
		TinymolCompoundinteractionsKernel<BoundaryCondition, emvariant>
			<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, cudaStreams[2]>>>
			(*boxStateCopy, *boxConfigCopy, compoundgridPtr, simulation->getStep());
		LIMA_UTILS::genericErrorCheckNoSync("Error after TinymolCompoundinteractionsKernel");

		// TODO: Too many threads, we rarely get close to filling the block
		solventForceKernel<BoundaryCondition, emvariant> 
			<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, cudaStreams[3]>>>
			(*boxStateCopy, *boxConfigCopy, compoundgridPtr, simulation->getStep());
		LIMA_UTILS::genericErrorCheckNoSync("Error after solventForceKernel");
	}
	
	if (ENABLE_ES_LR && simulation->simparams_host.enable_electrostatics) {
		// Must occur after DistributeCompoundchargesToGridKernel
		//timings.electrostatics += Electrostatics::HandleElectrostatics(sim_dev, boxparams, cudaStreams[2]);
		pmeController->CalcCharges(*boxConfigCopy, *boxStateCopy, boxparams.n_compounds, forceEnergiesPME);
		LIMA_UTILS::genericErrorCheckNoSync("Error after HandleElectrostatics");
	}

	if (simulation->simparams_host.snf_select != None) {
		SnfHandler<BoundaryCondition, emvariant>(cudaStreams[2]);
		LIMA_UTILS::genericErrorCheckNoSync("Error after SupernaturalForces");
	}

	if (!simulation->box_host->bondgroups.empty()) {
		BondgroupsKernel<BoundaryCondition, emvariant> << < simulation->box_host->bondgroups.size(), THREADS_PER_BONDSGROUPSKERNEL, 0, cudaStreams[4]>>> 
			(bondgroups, *boxStateCopy, forceEnergiesBondgroups);
		LIMA_UTILS::genericErrorCheckNoSync("Error after BondgroupsKernel");
	}

	// #### Integration and Transfer kernels
	cudaDeviceSynchronize();

	if (boxparams.n_compounds > 0) {
		CompoundIntegrationKernel<BoundaryCondition, emvariant> 
			<<<boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, 0, cudaStreams[0] >> >
			(sim_dev, simulation->getStep(), compoundForceEnergyInterims, forceEnergiesBondgroups, forceEnergiesPME);
		LIMA_UTILS::genericErrorCheckNoSync("Error after CompoundIntegrationKernel");
	}

	if (boxparams.n_solvents > 0) {
		const bool isTransferStep = SolventBlocksCircularQueue::isTransferStep(simulation->getStep());
		if (isTransferStep)
			TinymolIntegrationLoggingAndTransferout<BoundaryCondition, emvariant, true> 
				<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, cudaStreams[1] >>>
					(sim_dev, simulation->getStep());
		else
			TinymolIntegrationLoggingAndTransferout<BoundaryCondition, emvariant, false>
				<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, cudaStreams[1] >>>
					(sim_dev, simulation->getStep());
		LIMA_UTILS::genericErrorCheckNoSync("Error after TinymolIntegrationLoggingAndTransferout");

		if (isTransferStep) {
			cudaDeviceSynchronize();
			solventTransferKernel<BoundaryCondition> 
				<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlockTransfermodule::max_queue_size, 0, cudaStreams[1]>>> 
				(sim_dev, simulation->getStep());
			LIMA_UTILS::genericErrorCheckNoSync("Error after solventTransferKernel");
		}
	}	

	cudaDeviceSynchronize();
}
template void Engine::_deviceMaster<PeriodicBoundaryCondition, true, true>();
template void Engine::_deviceMaster<PeriodicBoundaryCondition, true, false>();
template void Engine::_deviceMaster<PeriodicBoundaryCondition, false, true>();
template void Engine::_deviceMaster<PeriodicBoundaryCondition, false, false>();
template void Engine::_deviceMaster<NoBoundaryCondition, true, true>();
template void Engine::_deviceMaster<NoBoundaryCondition, true, false>();
template void Engine::_deviceMaster<NoBoundaryCondition, false, true>();
template void Engine::_deviceMaster<NoBoundaryCondition, false, false>();



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






template <typename BoundaryCondition, bool emvariant>
void Engine::SnfHandler(cudaStream_t& stream) {
	switch (simulation->simparams_host.snf_select) {
	case None:
		break;
	case HorizontalSqueeze:
		SupernaturalForces::ApplyHorizontalSqueeze << < simulation->box_host->boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, 0, stream >> > (sim_dev, simulation->getStep());
		break;
	case HorizontalChargeField:
		CompoundSnfKernel<BoundaryCondition, emvariant>
			<< <simulation->box_host->boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, 0, stream >>>
			(sim_dev, simulation->getStep(), simulation->box_host->uniformElectricField, compoundForceEnergyInterims.forceEnergyBonds);
		break;
	case BoxEdgePotential:
		if (simulation->box_host->boxparams.n_compounds > 0)
			SupernaturalForces::BoxEdgeForceCompounds << < simulation->box_host->boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, 0, stream >> > (sim_dev, simulation->getStep());
		if (simulation->box_host->boxparams.n_solvents > 0)
			SupernaturalForces::BoxEdgeForceSolvents<<<BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(simulation->box_host->boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, 0, stream>>>(sim_dev, simulation->getStep());
		break;
	}
}








CompoundForceEnergyInterims::CompoundForceEnergyInterims(int nCompounds) {
	const size_t byteSize = sizeof(ForceEnergy) * nCompounds * MAX_COMPOUND_PARTICLES;
	cudaMalloc(&forceEnergyFarneighborShortrange, byteSize);
	cudaMalloc(&forceEnergyImmediateneighborShortrange, byteSize);
	cudaMalloc(&forceEnergyBonds, byteSize);

	cudaMemset(forceEnergyFarneighborShortrange, 0, byteSize);
	cudaMemset(forceEnergyImmediateneighborShortrange, 0, byteSize);
	cudaMemset(forceEnergyBonds, 0, byteSize);
}

void CompoundForceEnergyInterims::Free() {
	cudaFree(forceEnergyFarneighborShortrange);
	cudaFree(forceEnergyImmediateneighborShortrange);
	cudaFree(forceEnergyBonds);

	LIMA_UTILS::genericErrorCheck("Error during CompoundForceEnergyInterims destruction");
}
