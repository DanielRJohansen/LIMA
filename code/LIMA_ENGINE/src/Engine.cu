#include "Engine.cuh"

#include <algorithm>

Engine::Engine() {}
Engine::Engine(Simulation* simulation, ForceField_NB forcefield_host) {
	EngineUtils::genericErrorCheck("Error before engine initialization.\n");
	this->simulation = simulation;


	int Ckernel_shared_mem = sizeof(Compound) + sizeof(CompoundState) + sizeof(NeighborList) + sizeof(Float3) * NEIGHBORLIST_MAX_SOLVENTS + sizeof(uint8_t) * NEIGHBORLIST_MAX_SOLVENTS;	
	//int Skernel_shared_mem = sizeof(Float3) * MAX_COMPOUND_PARTICLES + sizeof(uint8_t) * MAX_COMPOUND_PARTICLES + sizeof(Solvent) * THREADS_PER_SOLVENTBLOCK;
	printf("Compoundkernel shared mem. size: %d B\n", Ckernel_shared_mem);
	//printf("Solventkernel shared mem. size: %d B\n", Skernel_shared_mem);


	this->forcefield_host = forcefield_host;
	setDeviceConstantMemory();


	// To create the NLists we need to bootstrap the traj_buffer, since it has no data yet
	nlist_manager = new NListManager(simulation);
	nlist_manager->bootstrapCompoundgrid(simulation);
	bootstrapTrajbufferWithCoords();
	handleNLISTS(simulation, true, true);


	printf("Engine ready\n\n\n");
}





void Engine::runOnce() {
	deviceMaster();	// Device first, otherwise offloading data always needs the last datapoint!
	hostMaster();
}

void Engine::deviceMaster() {
	EngineUtils::genericErrorCheck("Error before step!");
	step();
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
		handleNLISTS(simulation, true);
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


//--------------------------------------------------------------------------	CPU workload --------------------------------------------------------------//

// maybe rename forceupdate, to wait for update and push?
void Engine::handleNLISTS(Simulation* simulation, bool async, const bool force_update) {
	if (neighborlistUpdateRequired() && !updatenlists_mutexlock) {
		updatenlists_mutexlock = 1;

		nlist_manager->updateNeighborLists(simulation, &updatenlists_mutexlock, force_update, async, &timings.z, &critical_error);
	}

	if (nlist_manager->updated_neighborlists_ready) {
		nlist_manager->pushNlistsToDevice(simulation);
	}
}


void Engine::offloadLoggingData(const int steps_to_transfer) {
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

bool Engine::neighborlistUpdateRequired() const {
	const auto step = simulation->getStep();
	if (nlist_manager->stepsSinceUpdate(step) >= STEPS_PER_NLIST_UPDATE
		|| step == 0) {
		return true;
	}
	return false;
}


void Engine::bootstrapTrajbufferWithCoords() {
	CompoundCoords* compoundcoords_array = new CompoundCoords[simulation->n_compounds];
	cudaMemcpy(compoundcoords_array, simulation->box->coordarray_circular_queue, sizeof(CompoundCoords) * simulation->n_compounds, cudaMemcpyDeviceToHost);

	//SolventCoord* solventcoord_array = new SolventCoord[simulation->n_solvents];
	//cudaMemcpy(solventcoord_array, simulation->box->solventcoordarray_circular_queue, sizeof(SolventCoord) * simulation->n_solvents, cudaMemcpyDeviceToHost);

	// We need to bootstrap step-0 which is used for traj-buffer
	for (int compound_id = 0; compound_id < simulation->n_compounds; compound_id++) {
		for (int particle_id = 0; particle_id < MAX_COMPOUND_PARTICLES; particle_id++) {
			const int index = EngineUtils::getAlltimeIndexOfParticle(0, simulation->total_particles_upperbound, compound_id, particle_id);
			simulation->traj_buffer[index] = compoundcoords_array[compound_id].getAbsolutePositionLM(particle_id);
		}
	}

	for (int solvent_id = 0; solvent_id < simulation->n_solvents; solvent_id++) {
		const int index = EngineUtils::getAlltimeIndexOfParticle(0, simulation->total_particles_upperbound, simulation->n_compounds, solvent_id);
		//simulation->traj_buffer[index] = solventcoord_array[solvent_id].getAbsolutePositionLM();
	}

	delete[] compoundcoords_array;
	//delete[] solventcoord_array;
}




//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//
void Engine::step() {
	auto t0 = std::chrono::high_resolution_clock::now();
	cudaDeviceSynchronize();

	if (simulation->box->bridge_bundle->n_bridges > 0) {																		// TODO: Illegal access to device mem!!
		compoundBridgeKernel << < simulation->box->bridge_bundle->n_bridges, MAX_PARTICLES_IN_BRIDGE >> > (simulation->box);	// Must come before compoundKernel()		// DANGER
	}

	cudaDeviceSynchronize();
	if (simulation->n_compounds > 0) {
		compoundKernel << < simulation->n_compounds, THREADS_PER_COMPOUNDBLOCK >> > (simulation->box);
	}
	cudaDeviceSynchronize();	// Prolly not necessary

#ifdef ENABLE_SOLVENTS
	if (simulation->n_solvents > 0) { 
		solventForceKernel << < SolventBlockGrid::blocks_total, MAX_SOLVENTS_IN_BLOCK>> > (simulation->box);
	}
	cudaDeviceSynchronize();
	EngineUtils::genericErrorCheck("Error after solventForceKernel");
	if (SolventBlockHelpers::isTransferStep(simulation->getStep())) {
		solventTransferKernel <<< SolventBlockGrid::blocks_total, SolventBlockTransfermodule::max_queue_size>>> (simulation->box);
	}
	cudaDeviceSynchronize();
	EngineUtils::genericErrorCheck("Error after solventTransferKernel");
#endif

	auto t1 = std::chrono::high_resolution_clock::now();

	EngineUtils::genericErrorCheck("Error during step\n");		// Temp, we want to do host stuff while waiting for async GPU operations...	// SLOW


	cudaDeviceSynchronize();
	EngineUtils::genericErrorCheck("Error during step or state_transfer\n");		// Temp, we want to do host stuff while waiting for async GPU operations...	// SLOW
	auto t2 = std::chrono::high_resolution_clock::now();


	simulation->incStep();


	int force_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	int copy_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	timings = timings + Int3(force_duration, copy_duration, 0);
}