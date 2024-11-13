#include "Engine.cuh"
#include "Utilities.h"
#include "Neighborlists.cuh"
#include "Statistics.h"

#include "BoundaryCondition.cuh"
#include "EngineBodies.cuh"
#include "SimulationDevice.cuh"
#include "LimaPositionSystem.cuh"

#include "ChargeOcttree.cuh"

#include "EngineKernels.cuh"

#include "Thermostat.cuh"

#include "SupernaturalForces.cuh"


#include <unordered_set>

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

Engine::Engine(std::unique_ptr<Simulation> sim, BoundaryConditionSelect bc, std::unique_ptr<LimaLogger> logger)
	: bc_select(bc), m_logger(std::move(logger))
{
	simulation = std::move(sim);

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





	//std::unordered_set<std::string> unique_compounds;
	//for (int i = 0; i < simulation->box_host->boxparams.n_compounds; i++) {
	//	char types[64];
	//	memcpy(types, simulation->box_host->compounds[i].atom_types, 64);
	//	unique_compounds.insert(std::string(types));
	//}
	//int a = 0;






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
	const float cutoffLmSquaredReciprocal = 1.f / (simulation->simparams_host.cutoff_nm * NANO_TO_LIMA * simulation->simparams_host.cutoff_nm * NANO_TO_LIMA);
	cudaMemcpyToSymbol(cutoffLmSquaredReciprocal_device, &cutoffLmSquaredReciprocal, sizeof(float), 0, cudaMemcpyHostToDevice);


	const float initialThermostatScalar = 1.f;
	cudaMemcpyToSymbol(thermostatScalar_device, &initialThermostatScalar, sizeof(float), 0, cudaMemcpyHostToDevice);

	LIMA_UTILS::genericErrorCheck("Error while setting Global Constants\n");
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
		runstatus.greatestForce = greatestForce / LIMA * NANO / KILO; // Convert to [kJ/mol/nm]
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
	const auto t0 = std::chrono::high_resolution_clock::now();
	cudaDeviceSynchronize();

	const BoxParams& boxparams = simulation->box_host->boxparams;

	const bool logData = simulation->getStep() % simulation->simparams_host.data_logging_interval == 0;// TODO maybe log at the final step, not 0th?

	if (boxparams.n_compounds > 0) {
		compoundFarneighborShortrangeInteractionsKernel<BoundaryCondition, emvariant, computePotE> << <boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK >> > (simulation->getStep(), *boxStateCopy, *boxConfigCopy, neighborlistsPtr, simulation->simparams_host.enable_electrostatics);
		//LAUNCH_GENERIC_KERNEL_3(compoundFarneighborShortrangeInteractionsKernel, boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, bc_select, simulation->simparams_host.em_variant, logData, sim_dev, simulation->getStep());
	}

	LIMA_UTILS::genericErrorCheck("Error after compoundForceKernel");

	if (boxparams.n_compounds > 0) {
		//LAUNCH_GENERIC_KERNEL_3(compoundLJKernel, boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, bc_select, simulation->simparams_host.em_variant, logData, sim_dev, simulation->getStep());		
		//LAUNCH_GENERIC_KERNEL_3(compoundImmediateneighborAndSelfShortrangeInteractionsKernel, boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, bc_select, simulation->simparams_host.em_variant, logData, sim_dev, simulation->getStep());
		compoundImmediateneighborAndSelfShortrangeInteractionsKernel<BoundaryCondition, emvariant, computePotE> << <boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK >> > (sim_dev, simulation->getStep());
	}

	LIMA_UTILS::genericErrorCheck("Error after compoundForceKernel");


	const auto t0a = std::chrono::high_resolution_clock::now();
	cudaDeviceSynchronize();

	if constexpr (ENABLE_ES_LR) {
		if (simulation->simparams_host.enable_electrostatics) {
			timings.electrostatics += Electrostatics::HandleElectrostatics(sim_dev, boxparams);
		}
	}
	const auto t0b = std::chrono::high_resolution_clock::now();

	if (boxparams.n_bridges > 0) {
		//LAUNCH_GENERIC_KERNEL(compoundBridgeKernel, boxparams.n_bridges, MAX_PARTICLES_IN_BRIDGE, bc_select, sim_dev, simulation->getStep());
		compoundBridgeKernel<BoundaryCondition> << <boxparams.n_bridges, MAX_PARTICLES_IN_BRIDGE >> > (sim_dev, simulation->getStep());
	}

	if (simulation->simparams_host.snf_select != None) {
		SupernaturalForces::SnfHandler(simulation.get(), sim_dev, simulation->getStep());
	}

	cudaDeviceSynchronize();
	if (boxparams.n_compounds > 0) {
		//LAUNCH_GENERIC_KERNEL_2(compoundBondsAndIntegrationKernel, boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, bc_select, simulation->simparams_host.em_variant, sim_dev, simulation->getStep());
		compoundBondsAndIntegrationKernel<BoundaryCondition, emvariant> << <boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK >> > (sim_dev, simulation->getStep(), simulation->box_host->uniformElectricField);
	}
	LIMA_UTILS::genericErrorCheck("Error after compoundForceKernel");
	const auto t1 = std::chrono::high_resolution_clock::now();


#ifdef ENABLE_SOLVENTS
	if (boxparams.n_solvents > 0) {
		//LAUNCH_GENERIC_KERNEL_2(solventForceKernel, BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK, bc_select, simulation->simparams_host.em_variant, sim_dev, simulation->getStep());
		solventForceKernel<BoundaryCondition, emvariant> << <BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlock::MAX_SOLVENTS_IN_BLOCK >> > (sim_dev, simulation->getStep());


		cudaDeviceSynchronize();
		LIMA_UTILS::genericErrorCheck("Error after solventForceKernel");
		if (SolventBlocksCircularQueue::isTransferStep(simulation->getStep())) {
			//LAUNCH_GENERIC_KERNEL(solventTransferKernel, BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlockTransfermodule::max_queue_size, bc_select, sim_dev, simulation->getStep());
			solventTransferKernel<BoundaryCondition> << <BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), SolventBlockTransfermodule::max_queue_size >> > (sim_dev, simulation->getStep());
		}
	}
	cudaDeviceSynchronize();
	LIMA_UTILS::genericErrorCheck("Error after solventTransferKernel");
#endif

	if (boxparams.n_compounds > 0) {
		//LAUNCH_GENERIC_KERNEL_2(CompoundIntegrationKernel, boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK, bc_select, simulation->simparams_host.em_variant, sim_dev, simulation->getStep());
		CompoundIntegrationKernel<BoundaryCondition, emvariant> << <boxparams.n_compounds, THREADS_PER_COMPOUNDBLOCK >> > (sim_dev, simulation->getStep());
	}
	LIMA_UTILS::genericErrorCheck("Error after CompoundIntegrationKernel");

	const auto t2 = std::chrono::high_resolution_clock::now();

	const int compounds_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0b + t0a - t0).count();
	const int solvents_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

	timings.compound_kernels += compounds_duration;
	timings.solvent_kernels += solvents_duration;
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
