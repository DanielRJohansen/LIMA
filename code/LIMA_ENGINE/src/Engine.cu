#include "Engine.cuh"

#include <algorithm>

Engine::Engine() {}
Engine::Engine(Simulation* simulation, ForceField_NB forcefield_host) {
	EngineUtils::genericErrorCheck("Error before engine initialization.\n");
	this->simulation = simulation;


	int Ckernel_shared_mem = sizeof(Compound) + sizeof(CompoundState) + sizeof(NeighborList) + sizeof(Float3) * NEIGHBORLIST_MAX_SOLVENTS + sizeof(uint8_t) * NEIGHBORLIST_MAX_SOLVENTS;	
	int Skernel_shared_mem = sizeof(Float3) * MAX_COMPOUND_PARTICLES + sizeof(uint8_t) * MAX_COMPOUND_PARTICLES + sizeof(Solvent) * THREADS_PER_SOLVENTBLOCK;
	printf("Compoundkernel shared mem. size: %d B\n", Ckernel_shared_mem);
	printf("Solventkernel shared mem. size: %d B\n", Skernel_shared_mem);


	this->forcefield_host = forcefield_host;
	setDeviceConstantMemory();


	nlist_manager = new NListManager(simulation);
	handleNLISTS(simulation, false, true);				// Fix neighborlists before running


	printf("Engine ready\n\n\n");
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
		//offloadPositionData();

		if ((simulation->getStep() % STEPS_PER_THERMOSTAT) == 0 && ENABLE_BOXTEMP) {
			handleBoxtemp();
		}
	}
	if ((simulation->getStep() % STEPS_PER_TRAINDATATRANSFER) == 0) {
		offloadTrainData();
	}

	handleNLISTS(simulation, false);
	

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
	}
}


//--------------------------------------------------------------------------	CPU workload --------------------------------------------------------------//


void Engine::handleNLISTS(Simulation* simulation, bool async, bool force_update) {

	if (((nlist_manager->stepsSinceUpdate(simulation->getStep() ) >= STEPS_PER_NLIST_UPDATE) || force_update) && !updatenlists_mutexlock) {
		updatenlists_mutexlock = 1;
		
		nlist_manager->offloadPositionDataNLIST(simulation);
		// Lots of waiting time here...
		cudaDeviceSynchronize();

		nlist_manager->updateNeighborLists(simulation, &updatenlists_mutexlock, force_update, async, &timings.z);
	}

	if (nlist_manager->updated_neighborlists_ready) {
		nlist_manager->pushNlistsToDevice(simulation);
	}
}


void Engine::offloadLoggingData(const int steps_to_transfer) {
	uint64_t step_offset = (simulation->getStep() - STEPS_PER_LOGTRANSFER) ;	// Tongue in cheek here, i think this is correct...

	cudaMemcpy(&simulation->potE_buffer[step_offset * simulation->total_particles_upperbound], simulation->box->potE_buffer, sizeof(float) * simulation->total_particles_upperbound * steps_to_transfer, cudaMemcpyDeviceToHost);
	
	cudaMemcpy(&simulation->traj_buffer[step_offset * simulation->total_particles_upperbound], simulation->box->traj_buffer, sizeof(Float3) * simulation->total_particles_upperbound * steps_to_transfer, cudaMemcpyDeviceToHost);

	cudaMemcpy(&simulation->logging_data[step_offset * 10], simulation->box->outdata, sizeof(float) * 10 * steps_to_transfer, cudaMemcpyDeviceToHost);
}

void Engine::offloadPositionData() {
	//uint64_t step_offset = (simulation->getStep() - STEPS_PER_LOGTRANSFER) * simulation->total_particles_upperbound;	// Tongue in cheek here, i think this is correct...
	//cudaMemcpy(&simulation->traj_buffer[step_offset], simulation->box->traj_buffer, sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER, cudaMemcpyDeviceToHost);
}

void Engine::offloadTrainData() {
	uint64_t values_per_step = N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->n_compounds;
	uint64_t step_offset = (simulation->getStep() - STEPS_PER_TRAINDATATRANSFER) * values_per_step;	// fix max_compound to the actual count save LOTS of space!. Might need a file in simout that specifies cnt for loading in other programs...
	cudaMemcpy(&simulation->traindata_buffer[step_offset], simulation->box->data_GAN, sizeof(Float3) * values_per_step * STEPS_PER_TRAINDATATRANSFER, cudaMemcpyDeviceToHost);
	EngineUtils::genericErrorCheck("Cuda error during traindata offloading\n");
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
		compoundKernel << < simulation->box->n_compounds, THREADS_PER_COMPOUNDBLOCK >> > (simulation->box);
	}
	cudaDeviceSynchronize();

#ifdef ENABLE_SOLVENTS
	if (simulation->n_solvents > 0) { 
		solventForceKernel << < simulation->blocks_per_solventkernel, THREADS_PER_SOLVENTBLOCK >> > (simulation->box); 
	}
#endif
	cudaDeviceSynchronize();

	auto t1 = std::chrono::high_resolution_clock::now();

	EngineUtils::genericErrorCheck("Error during step\n");		// Temp, we want to do host stuff while waiting for async GPU operations...	// SLOW


	
	
	Solvent* temp_s = simulation->box->solvents;
	simulation->box->solvents = simulation->box->solvents_next;
	simulation->box->solvents_next = temp_s;


	
	
	
	
	
	//cudaMemcpy(simulation->box->solvents, simulation->box->solvents_next, sizeof(Solvent) * MAX_SOLVENTS, cudaMemcpyDeviceToDevice);
	cudaDeviceSynchronize();
	EngineUtils::genericErrorCheck("Error during step or state_transfer\n");		// Temp, we want to do host stuff while waiting for async GPU operations...	// SLOW
	auto t2 = std::chrono::high_resolution_clock::now();


	simulation->incStep();


	int force_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	int copy_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	timings = timings + Int3(force_duration, copy_duration, 0);
}