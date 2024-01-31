

#include "Engine.cuh"
#include "Utilities.h"
#include "EngineUtils.cuh"
#include "Neighborlists.cuh"

#include "BoundaryCondition.cuh"
#include "EngineBodies.cuh"
#include "SimulationDevice.cuh"

#include <algorithm>






Engine::Engine(std::unique_ptr<Simulation> sim, BoundaryConditionSelect bc, std::unique_ptr<LimaLogger> logger)
	: bc_select(bc), m_logger(std::move(logger))
{
	simulation = std::move(sim);

	verifyEngine();

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

	// Create the Sim_dev
	if (sim_dev != nullptr) { throw std::runtime_error("Expected simdev to be null to move sim to device"); }
	sim_dev = new SimulationDevice(simulation->simparams_host, std::move(simulation->box_host));
	sim_dev = genericMoveToDevice(sim_dev, 1);




	//this->forcefield_host = forcefield_host;
	setDeviceConstantMemory();

	LIMA_UTILS::genericErrorCheck("Error during bootstrapTrajbufferWithCoords");

	// To create the NLists we need to bootstrap the traj_buffer, since it has no data yet
	bootstrapTrajbufferWithCoords();

	NeighborLists::updateNlists(sim_dev, simulation->simparams_host.bc_select, simulation->boxparams_host, timings.nlist);
	m_logger->finishSection("Engine Ready");
}

Engine::~Engine() {
	if (sim_dev != nullptr) {
		sim_dev->deleteMembers();
		cudaFree(sim_dev);
	}

	assert(simulation == nullptr);
}

std::unique_ptr<Simulation> Engine::takeBackSim() {
	assert(sim_dev);
	simulation->box_host = SimUtils::copyToHost(sim_dev->box);
	return std::move(simulation);
}

void Engine::verifyEngine() {
	LIMA_UTILS::genericErrorCheck("Error before engine initialization.\n");

	if (simulation->boxparams_host.dims.x != BOX_LEN_NM) {
		//throw std::runtime_error(std::format("This simulations box_size of {} did not match the size the engine is compiled with {}", simulation->boxparams_host.dims.x, BOX_LEN_NM));
		throw std::runtime_error("This simulations box_size of "+ std::to_string(simulation->boxparams_host.dims.x)
		+ "did not match the size the engine is compiled with" + std::to_string(BOX_LEN_NM));
	}
}

void Engine::step() {
	LIMA_UTILS::genericErrorCheck("Error before step!");

	deviceMaster();	// Device first, otherwise offloading data always needs the last datapoint!
	//simulation->incStep();
	assert(simulation);
	assert(sim_dev);
	simulation->simsignals_host.step++;
	sim_dev->signals->step++;	// UNSAFE

	hostMaster();

	LIMA_UTILS::genericErrorCheck("Error after step!");
}

void Engine::hostMaster() {						// This is and MUST ALWAYS be called after the deviceMaster, and AFTER incStep()!
	auto t0 = std::chrono::high_resolution_clock::now();
	if ((simulation->getStep() % STEPS_PER_LOGTRANSFER) == 0) {
		offloadLoggingData();
		offloadTrajectory();


		if ((simulation->getStep() % STEPS_PER_THERMOSTAT) == 0 && ENABLE_BOXTEMP) {
			handleBoxtemp();
		}
		if (simulation->simparams_host.em_variant) {
			sim_dev->signals->thermostat_scalar = 0.9995f;	// The current implementation of em is only this, overwriting the velocity scaler to always be 0.5
		}

		//nlist_manager->handleNLISTS(simulation.get(), ALLOW_ASYNC_NLISTUPDATE, false, &timings.nlist);
		NeighborLists::updateNlists(sim_dev, simulation->simparams_host.bc_select, simulation->boxparams_host, timings.nlist);
	}
	if ((simulation->getStep() % STEPS_PER_TRAINDATATRANSFER) == 0) {
		offloadTrainData();
	}

	// Handle status
	runstatus.current_step = simulation->getStep();
	runstatus.critical_error_occured = sim_dev->signals->critical_error_encountered;	// TODO: Can i get this from simparams_host? UNSAFE
	// most recent positions are handled automaticall by transfer_traj
	runstatus.simulation_finished = runstatus.current_step >= simulation->simparams_host.n_steps || runstatus.critical_error_occured;

	//if ((simulation->getStep() % STEPS_PER_THERMOSTAT) == 1) {	// So this runs 1 step AFTER handleBoxtemp
	//	simulation->box->thermostat_scalar = 1.f;
	//}

	const auto t1 = std::chrono::high_resolution_clock::now();
	const int cpu_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	timings.cpu_master += cpu_duration;
}

void Engine::terminateSimulation() {
	const auto steps_since_transfer = simulation->getStep() % STEPS_PER_LOGTRANSFER;
	if ((steps_since_transfer) > LOG_EVERY_N_STEPS) {
		offloadLoggingData(steps_since_transfer);
		offloadTrajectory(steps_since_transfer);
	}
}
#include <assert.h>

//--------------------------------------------------------------------------	CPU workload --------------------------------------------------------------//

void Engine::offloadLoggingData(const int steps_to_transfer) {
	assert(steps_to_transfer <= simulation->getStep());

	const int64_t startstep = simulation->getStep() - steps_to_transfer;
	const int64_t startindex = LIMALOGSYSTEM::getMostRecentDataentryIndex(startstep);
	const int64_t indices_to_transfer = LIMALOGSYSTEM::getNIndicesBetweenSteps(startstep, simulation->getStep());

	cudaMemcpy(
		simulation->potE_buffer->getBufferAtIndex(startindex),
		//&simulation->potE_buffer[step_relative * simulation->boxparams_host.total_particles_upperbound],
		sim_dev->databuffers->potE_buffer, 
		sizeof(float) * simulation->boxparams_host.total_particles_upperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost);

	cudaMemcpy(
		simulation->vel_buffer->getBufferAtIndex(startindex),
		sim_dev->databuffers->vel_buffer,
		sizeof(float) * simulation->boxparams_host.total_particles_upperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost);

#ifdef GENERATETRAINDATA
	cudaMemcpy(	// THIS IS PROLLY WRONG NOW
		&simulation->loggingdata[startindex * 10],
		sim_dev->databuffers->outdata, 
		sizeof(float) * 10 * indices_to_transfer,
		cudaMemcpyDeviceToHost);
#endif
	cudaDeviceSynchronize();
}

void Engine::offloadTrajectory(const int steps_to_transfer) {
#ifndef DONTGENDATA

	const int64_t startstep = simulation->getStep() - steps_to_transfer;
	const int64_t startindex = LIMALOGSYSTEM::getMostRecentDataentryIndex(startstep);
	const int64_t indices_to_transfer = LIMALOGSYSTEM::getNIndicesBetweenSteps(startstep, simulation->getStep());

	cudaMemcpy(
		//&simulation->traj_buffer[step_relative * simulation->total_particles_upperbound],
		simulation->traj_buffer->getBufferAtIndex(startindex),
		sim_dev->databuffers->traj_buffer,
		sizeof(Float3) * simulation->boxparams_host.total_particles_upperbound * indices_to_transfer,
		cudaMemcpyDeviceToHost
	);

	cudaDeviceSynchronize();
#endif
	step_at_last_traj_transfer = simulation->getStep();
	runstatus.most_recent_positions = simulation->traj_buffer->getBufferAtIndex(LIMALOGSYSTEM::getMostRecentDataentryIndex(simulation->getStep()-1));
}


void Engine::offloadTrainData() {
#ifdef GENERATETRAINDATA
	uint64_t values_per_step = N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->boxparams_host.n_compounds;
	if (values_per_step == 0) {
		return;	// No data to transfer
	}

	uint64_t step_offset = (simulation->getStep() - STEPS_PER_TRAINDATATRANSFER) * values_per_step;	// fix max_compound to the actual count save LOTS of space!. Might need a file in simout that specifies cnt for loading in other programs...
	cudaMemcpy(&simulation->trainingdata[step_offset], sim_dev->databuffers->data_GAN, sizeof(Float3) * values_per_step * STEPS_PER_TRAINDATATRANSFER, cudaMemcpyDeviceToHost);
	LIMA_UTILS::genericErrorCheck("Cuda error during traindata offloading\n");
#endif
}


void Engine::bootstrapTrajbufferWithCoords() {
	if (simulation->simparams_host.n_steps == 0) return;
	LIMA_UTILS::genericErrorCheck("Error during bootstrapTrajbufferWithCoords");

	std::vector<CompoundCoords> compoundcoords_array(simulation->boxparams_host.n_compounds);
	auto error = cudaMemcpy(compoundcoords_array.data(), sim_dev->box->coordarray_circular_queue, sizeof(CompoundCoords) * simulation->boxparams_host.n_compounds, cudaMemcpyDeviceToHost);
	LIMA_UTILS::genericErrorCheck(error);
	

	// We need to bootstrap step-0 which is used for traj-buffer
	for (int compound_id = 0; compound_id < simulation->boxparams_host.n_compounds; compound_id++) {
		for (int particle_id = 0; particle_id < MAX_COMPOUND_PARTICLES; particle_id++) {

			const Float3 particle_abspos = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compoundcoords_array[compound_id].origo, compoundcoords_array[compound_id].rel_positions[particle_id]);
			simulation->traj_buffer->getCompoundparticleDatapointAtIndex(compound_id, particle_id, 0) = particle_abspos;
		}
	}

	LIMA_UTILS::genericErrorCheck("Error during bootstrapTrajbufferWithCoords");
}




//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//

void Engine::deviceMaster() {
	const auto t0 = std::chrono::high_resolution_clock::now();
	cudaDeviceSynchronize();


	if (simulation->boxparams_host.n_compounds > 0) {
		LAUNCH_GENERIC_KERNEL(compoundLJKernel, simulation->boxparams_host.n_compounds, THREADS_PER_COMPOUNDBLOCK, bc_select, sim_dev);
		//compoundLJKernel<BoundaryCondition> << < simulation->boxparams_host.n_compounds, THREADS_PER_COMPOUNDBLOCK >> > (sim_dev);
	}

	cudaDeviceSynchronize();

	if (simulation->boxparams_host.n_bridges > 0) {
		LAUNCH_GENERIC_KERNEL(compoundBridgeKernel, simulation->boxparams_host.n_bridges, MAX_PARTICLES_IN_BRIDGE, bc_select, sim_dev);
		//compoundBridgeKernel<BoundaryCondition> <<< simulation->boxparams_host.n_bridges, MAX_PARTICLES_IN_BRIDGE >> > (sim_dev);	// Must come before compoundKernel()
	}

	cudaDeviceSynchronize();
	if (simulation->boxparams_host.n_compounds > 0) {
		LAUNCH_GENERIC_KERNEL(compoundBondsAndIntegrationKernel, simulation->boxparams_host.n_compounds, THREADS_PER_COMPOUNDBLOCK, bc_select, sim_dev);
		//compoundBondsAndIntegrationKernel<BoundaryCondition> << <simulation->boxparams_host.n_compounds, THREADS_PER_COMPOUNDBLOCK >> > (sim_dev);
	}
	LIMA_UTILS::genericErrorCheck("Error after compoundForceKernel");
	const auto t1 = std::chrono::high_resolution_clock::now();


#ifdef ENABLE_SOLVENTS
	if (simulation->boxparams_host.n_solvents > 0) {
		LAUNCH_GENERIC_KERNEL(solventForceKernel, SolventBlocksCircularQueue::blocks_per_grid, SolventBlock::MAX_SOLVENTS_IN_BLOCK, bc_select, sim_dev);
		//solventForceKernel<BoundaryCondition> << < SolventBlocksCircularQueue::blocks_per_grid, SolventBlock::MAX_SOLVENTS_IN_BLOCK>> > (sim_dev);


		cudaDeviceSynchronize();
		LIMA_UTILS::genericErrorCheck("Error after solventForceKernel");
		if (SolventBlocksCircularQueue::isTransferStep(simulation->getStep())) {
			LAUNCH_GENERIC_KERNEL(solventTransferKernel, SolventBlocksCircularQueue::blocks_per_grid, SolventBlockTransfermodule::max_queue_size, bc_select, sim_dev);
			//solventTransferKernel<BoundaryCondition> << < SolventBlocksCircularQueue::blocks_per_grid, SolventBlockTransfermodule::max_queue_size >> > (sim_dev);
		}
	}
	cudaDeviceSynchronize();
	LIMA_UTILS::genericErrorCheck("Error after solventTransferKernel");
#endif
	const auto t2 = std::chrono::high_resolution_clock::now();

	const int compounds_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	const int solvents_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

	timings.compound_kernels += compounds_duration;
	timings.solvent_kernels += solvents_duration;
}
