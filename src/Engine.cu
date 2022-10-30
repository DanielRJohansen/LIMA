#include "Engine.cuh"



Engine::Engine() {}
Engine::Engine(Simulation* simulation, ForceField forcefield_host) {
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

void Engine::hostMaster() {						// This is and MUST ALWAYS be called after the deviceMaster, and AFTER intStep()!
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

	handleNLISTS(simulation);
	

	if ((simulation->getStep() % STEPS_PER_THERMOSTAT) == 1) {	// So this runs 1 step AFTER handleBoxtemp
		simulation->box->thermostat_scalar = 1.f;
	}

	auto t1 = std::chrono::high_resolution_clock::now();

	int cpu_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	timings = timings + Int3(0,0,cpu_duration);
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


void Engine::offloadLoggingData() {
	uint64_t step_offset = (simulation->getStep() - STEPS_PER_LOGTRANSFER) ;	// Tongue in cheek here, i think this is correct...

	cudaMemcpy(&simulation->potE_buffer[step_offset * simulation->total_particles_upperbound], simulation->box->potE_buffer, sizeof(float) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER, cudaMemcpyDeviceToHost);
	
	cudaMemcpy(&simulation->traj_buffer[step_offset * simulation->total_particles_upperbound], simulation->box->traj_buffer, sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER, cudaMemcpyDeviceToHost);

	cudaMemcpy(&simulation->logging_data[step_offset * 10], simulation->box->outdata, sizeof(float) * 10 * STEPS_PER_LOGTRANSFER, cudaMemcpyDeviceToHost);
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



																																			// THIS fn requires mallocmanaged!!   // HARD DISABLED HERE
void Engine::handleBoxtemp() {
	const float target_temp = 310.f;				// [k]
	Float3 temp_package = EngineUtils::getBoxTemperature(simulation, forcefield_host);
	float temp = temp_package.x;
	float biggest_contribution = temp_package.y;

	simulation->temperature_buffer[simulation->n_temp_values++] = temp;

	// So we avoid dividing by 0
	float temp_safe = temp == 0.f ? 1 : temp;
	float temp_scalar = target_temp / temp;
		// I just added this to not change any temperatures too rapidly
	temp_scalar = min(temp_scalar, 1.01f);
	temp_scalar = max(temp_scalar, 0.99f);

	if (PRINT_TEMP || temp > 500.f || temp < 100.f) { printf("\n %llu Temperature: %.1f Biggest contrib: %.0f avg kinE %.0f\n", (simulation->getStep() - 1) / STEPS_PER_THERMOSTAT, temp, biggest_contribution, temp_package.z); }
		
	if (temp > target_temp/4.f && temp < target_temp*4.f || true) {
		if (APPLY_THERMOSTAT && simulation->getStep() > 10) {
			
			simulation->box->thermostat_scalar = temp_scalar;

			if (temp_scalar != temp_scalar ){//} || abs(temp_scalar) == "inf") {
				printf("Scalar: %f\n", simulation->box->thermostat_scalar);
				exit(0);
			}			
		}		
	}
	else {
		printf("Critical temperature encountered (%0.02f [k])\n", temp);
		simulation->box->critical_error_encountered = true;
	}
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


	
	CompoundState* temp = simulation->box->compound_state_array;
	simulation->box->compound_state_array = simulation->box->compound_state_array_next;
	simulation->box->compound_state_array_next = temp;
	
	
	Solvent* temp_s = simulation->box->solvents;
	simulation->box->solvents = simulation->box->solvents_next;
	simulation->box->solvents_next = temp_s;
	
	
	//cudaMemcpy(simulation->box->compound_state_array, simulation->box->compound_state_array_next, sizeof(CompoundState) * MAX_COMPOUNDS, cudaMemcpyDeviceToDevice);	// Update all positions, after all forces have been calculated
	
	
	
	//cudaMemcpy(simulation->box->solvents, simulation->box->solvents_next, sizeof(Solvent) * MAX_SOLVENTS, cudaMemcpyDeviceToDevice);
	cudaDeviceSynchronize();
	EngineUtils::genericErrorCheck("Error during step or state_transfer\n");		// Temp, we want to do host stuff while waiting for async GPU operations...	// SLOW
	auto t2 = std::chrono::high_resolution_clock::now();


	simulation->incStep();


	int force_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	int copy_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	timings = timings + Int3(force_duration, copy_duration, 0);
}