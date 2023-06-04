#include "LIMA_ENGINE/include/Engine.cuh"
#include <algorithm>

Engine::Engine() {}
Engine::Engine(Simulation* simulation, ForceField_NB forcefield_host) {
	EngineUtils::genericErrorCheck("Error before engine initialization.\n");
	this->simulation = simulation;

	const int Ckernel_shared_mem = sizeof(Compound) + sizeof(CompoundState) + sizeof(CompoundCoords) + sizeof(NeighborList) + sizeof(BondedParticlesLUT) + sizeof(Float3) * THREADS_PER_COMPOUNDBLOCK + sizeof(Coord) * 2;
	static_assert(Ckernel_shared_mem < 45000, "Not enough shared memory for CompoundKernel");
	//int Skernel_shared_mem = sizeof(Float3) * MAX_COMPOUND_PARTICLES + sizeof(uint8_t) * MAX_COMPOUND_PARTICLES + sizeof(Solvent) * THREADS_PER_SOLVENTBLOCK;
	printf("Compoundkernel shared mem. size: %d B\n", Ckernel_shared_mem);
	//printf("Solventkernel shared mem. size: %d B\n", Skernel_shared_mem);


	this->forcefield_host = forcefield_host;
	setDeviceConstantMemory();

	EngineUtils::genericErrorCheck("Error during bootstrapTrajbufferWithCoords");

	// To create the NLists we need to bootstrap the traj_buffer, since it has no data yet
	nlist_manager = std::make_unique<NListManager>(simulation);
				//nlist_manager->bootstrapCompoundgrid(simulation);
	bootstrapTrajbufferWithCoords();
	nlist_manager->handleNLISTS(simulation, true, true, &timings.z);


	printf("Engine ready\n\n\n");
}





template<bool em_variant> void Engine::step() {
	EngineUtils::genericErrorCheck("Error before step!");

	deviceMaster<em_variant>();	// Device first, otherwise offloading data always needs the last datapoint!
	simulation->incStep();
	hostMaster();

	EngineUtils::genericErrorCheck("Error after step!");
}

void Engine::hostMaster() {						// This is and MUST ALWAYS be called after the deviceMaster, and AFTER incStep()!
	auto t0 = std::chrono::high_resolution_clock::now();
	if ((simulation->getStep() % STEPS_PER_LOGTRANSFER) == 0) {
		offloadLoggingData();
		offloadTrajectory();


		if ((simulation->getStep() % STEPS_PER_THERMOSTAT) == 0 && ENABLE_BOXTEMP) {
			handleBoxtemp();
		}
		nlist_manager->handleNLISTS(simulation, ALLOW_ASYNC_NLISTUPDATE, false, &timings.z);
	}
	if ((simulation->getStep() % STEPS_PER_TRAINDATATRANSFER) == 0) {
		offloadTrainData();
	}


	//if ((simulation->getStep() % STEPS_PER_THERMOSTAT) == 1) {	// So this runs 1 step AFTER handleBoxtemp
	//	simulation->box->thermostat_scalar = 1.f;
	//}

	auto t1 = std::chrono::high_resolution_clock::now();

	int cpu_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	timings = timings + Int3(0,0,cpu_duration);
}

void Engine::terminateSimulation() {
	const int steps_since_transfer = simulation->getStep() % STEPS_PER_LOGTRANSFER;
	if ((steps_since_transfer) > 0) {
		offloadLoggingData(steps_since_transfer);
		offloadTrajectory(steps_since_transfer);
	}
}
#include <assert.h>

//--------------------------------------------------------------------------	CPU workload --------------------------------------------------------------//

void Engine::offloadLoggingData(const int steps_to_transfer) {
	assert(steps_to_transfer <= simulation->getStep());
	uint64_t step_relative = (simulation->getStep() - steps_to_transfer) ;	// Tongue in cheek here, i think this is correct...

	cudaMemcpy(
		&simulation->potE_buffer[step_relative * simulation->total_particles_upperbound], 
		simulation->box->potE_buffer, 
		sizeof(float) * simulation->total_particles_upperbound * steps_to_transfer, 
		cudaMemcpyDeviceToHost);

	cudaMemcpy(
		&simulation->logging_data[step_relative * 10], 
		simulation->box->outdata, 
		sizeof(float) * 10 * steps_to_transfer, 
		cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();
}

void Engine::offloadTrajectory(const int steps_to_transfer) {
	uint64_t step_relative = (simulation->getStep() - steps_to_transfer);	// Tongue in cheek here, i think this is correct...

	cudaMemcpy(
		&simulation->traj_buffer[step_relative * simulation->total_particles_upperbound],
		simulation->box->traj_buffer,
		sizeof(Float3) * simulation->total_particles_upperbound * steps_to_transfer,
		cudaMemcpyDeviceToHost
	);

	cudaDeviceSynchronize();
	step_at_last_traj_transfer = simulation->getStep();
}


void Engine::offloadTrainData() {
	uint64_t values_per_step = N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->n_compounds;
	uint64_t step_offset = (simulation->getStep() - STEPS_PER_TRAINDATATRANSFER) * values_per_step;	// fix max_compound to the actual count save LOTS of space!. Might need a file in simout that specifies cnt for loading in other programs...
	cudaMemcpy(&simulation->traindata_buffer[step_offset], simulation->box->data_GAN, sizeof(Float3) * values_per_step * STEPS_PER_TRAINDATATRANSFER, cudaMemcpyDeviceToHost);
	EngineUtils::genericErrorCheck("Cuda error during traindata offloading\n");
}



void Engine::bootstrapTrajbufferWithCoords() {
	EngineUtils::genericErrorCheck("Error during bootstrapTrajbufferWithCoords");


	CompoundCoords* compoundcoords_array = new CompoundCoords[simulation->n_compounds];
	auto error = cudaMemcpy(compoundcoords_array, simulation->box->coordarray_circular_queue, sizeof(CompoundCoords) * simulation->n_compounds, cudaMemcpyDeviceToHost);
	EngineUtils::genericErrorCheck(error);
	

	// We need to bootstrap step-0 which is used for traj-buffer
	for (int compound_id = 0; compound_id < simulation->n_compounds; compound_id++) {
		for (int particle_id = 0; particle_id < MAX_COMPOUND_PARTICLES; particle_id++) {
			const size_t buffer_index = EngineUtils::getAlltimeIndexOfParticle(0, simulation->total_particles_upperbound, compound_id, particle_id);

			//simulation->traj_buffer[index] = compoundcoords_array[compound_id].getAbsolutePositionLM(particle_id);
			simulation->traj_buffer[buffer_index] = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compoundcoords_array[compound_id].origo, compoundcoords_array[compound_id].rel_positions[particle_id]);
		}
	}

	EngineUtils::genericErrorCheck("Error during bootstrapTrajbufferWithCoords");

	delete[] compoundcoords_array;
}




//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//
template <bool em_variant>
void Engine::deviceMaster() {
	auto t0 = std::chrono::high_resolution_clock::now();
	cudaDeviceSynchronize();



	if (simulation->box->bridge_bundle->n_bridges > 0) {																		// TODO: Illegal access to device mem!!
		compoundBridgeKernel<<< simulation->box->bridge_bundle->n_bridges, MAX_PARTICLES_IN_BRIDGE >> > (simulation->box, simulation->simparams_device);	// Must come before compoundKernel()		// DANGER
	}

	cudaDeviceSynchronize();
	if (simulation->n_compounds > 0) {
		compoundKernel<em_variant><< < simulation->n_compounds, THREADS_PER_COMPOUNDBLOCK >> > (simulation->box, simulation->simparams_device);
	}

	cudaDeviceSynchronize();	// Prolly not necessary
	EngineUtils::genericErrorCheck("Error after compoundForceKernel");
#ifdef ENABLE_SOLVENTS
	if (simulation->n_solvents > 0) { 
		solventForceKernel<em_variant> << < SolventBlockGrid::blocks_total, MAX_SOLVENTS_IN_BLOCK>> > (simulation->box, simulation->simparams_device);
//		return;

		cudaDeviceSynchronize();
		EngineUtils::genericErrorCheck("Error after solventForceKernel");
		if (SolventBlockHelpers::isTransferStep(simulation->getStep())) {
			solventTransferKernel << < SolventBlockGrid::blocks_total, SolventBlockTransfermodule::max_queue_size >> > (simulation->box, simulation->simparams_device);
		}
	}
	cudaDeviceSynchronize();
	EngineUtils::genericErrorCheck("Error after solventTransferKernel");
#endif
	return;

	auto t1 = std::chrono::high_resolution_clock::now();

	EngineUtils::genericErrorCheck("Error during step\n");		// Temp, we want to do host stuff while waiting for async GPU operations...	// SLOW


	cudaDeviceSynchronize();
	EngineUtils::genericErrorCheck("Error during step or state_transfer\n");		// Temp, we want to do host stuff while waiting for async GPU operations...	// SLOW
	auto t2 = std::chrono::high_resolution_clock::now();

	int force_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	int copy_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	timings = timings + Int3(force_duration, copy_duration, 0);
}