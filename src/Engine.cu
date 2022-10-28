#include "Engine.cuh"


__constant__ ForceField forcefield_device;
//__constant__ ForceField forcefield_nb_device;


Engine::Engine() {}
Engine::Engine(Simulation* simulation, ForceField forcefield_host) {
	LIMAENG::genericErrorCheck("Error before engine initialization.\n");
	this->simulation = simulation;
	nlist_data_collection = new NListDataCollection(simulation);




	int Ckernel_shared_mem = sizeof(Compound) + sizeof(CompoundState) + sizeof(NeighborList) + sizeof(Float3) * NEIGHBORLIST_MAX_SOLVENTS + sizeof(uint8_t) * NEIGHBORLIST_MAX_SOLVENTS;	
	int Skernel_shared_mem = sizeof(Float3) * MAX_COMPOUND_PARTICLES + sizeof(uint8_t) * MAX_COMPOUND_PARTICLES + sizeof(Solvent) * THREADS_PER_SOLVENTBLOCK;
	printf("Compoundkernel shared mem. size: %d B\n", Ckernel_shared_mem);
	printf("Solventkernel shared mem. size: %d B\n", Skernel_shared_mem);

	LIMAENG::genericErrorCheck("Error before moving forcefield to device\n");



	//ForceField forcefield_host = FFM.getForcefield();
	//ForceField forcefield_host = FFM.getNBForcefield();
	this->forcefield_host = forcefield_host;
	cudaMemcpyToSymbol(forcefield_device, &forcefield_host, sizeof(ForceField), 0, cudaMemcpyHostToDevice);	// So there should not be a & before the device __constant__

	printf("Forcefield size: %d bytes\n", sizeof(ForceField));

	//ForceField forcefield_nb_host = FFM.getNBForcefield();
	//cudaMemcpyToSymbol(forcefield_nb_device, &forcefield_nb_host, sizeof(ForceField), 0, cudaMemcpyHostToDevice);
	//cudaDeviceSynchronize();
	LIMAENG::genericErrorCheck("Error while moving forcefield to device\n");

	for (int i = 0; i < nlist_data_collection->n_compounds; i++) {
		nlist_data_collection->compound_neighborlists[i].associated_id = i;
	}
	for (int i = 0; i < nlist_data_collection->n_solvents; i++) {
		nlist_data_collection->solvent_neighborlists[i].associated_id = i;
	}
	handleNLISTS(simulation, false, true);				// Fix neighborlists before running

	printf("Engine ready\n\n\n");
}





void Engine::deviceMaster() {
	LIMAENG::genericErrorCheck("Error before step!");
	step();
	LIMAENG::genericErrorCheck("Error after step!");
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



	if (((simulation->getStep() - prev_nlist_update_step >= STEPS_PER_NLIST_UPDATE) || force_update) && !updatenlists_mutexlock) {
		updatenlists_mutexlock = 1;
		offloadPositionDataNLIST(simulation);
		// Lots of waiting time here...
		cudaDeviceSynchronize();


		if (async) {
			std::thread nlist_worker(Engine::updateNeighborLists, simulation, nlist_data_collection, &updated_neighborlists_ready, &timings.z, &updatenlists_mutexlock);
			nlist_worker.detach();
		}
		else {
			updateNeighborLists(simulation, nlist_data_collection, &updated_neighborlists_ready, &timings.z, &updatenlists_mutexlock);
		}

		prev_nlist_update_step = simulation->getStep();
	}

	if (updated_neighborlists_ready) {
		pushNlistsToDevice();
	}

	if (force_update) {
		Int3 n_data(nlist_data_collection->compound_neighborlists[0].n_compound_neighbors, nlist_data_collection->compound_neighborlists[0].n_solvent_neighbors, 0);
		//Int3 after(nlist_data_collection->solvent_neighborlists[193].n_compound_neighbors, nlist_data_collection->solvent_neighborlists[193].n_solvent_neighbors, 0);

		printf("\nEntity neighbors: %d %d\n", n_data.x, n_data.y);
	}
}

void Engine::offloadPositionDataNLIST(Simulation* simulation) {
	//cudaMemcpyAsync(compoundstatearray_host, simulation->box->compound_state_array, sizeof(CompoundState) * simulation->box->n_compounds, cudaMemcpyDeviceToHost);
	cudaMemcpy(nlist_data_collection->compoundstates, simulation->box->compound_state_array, sizeof(CompoundState) * simulation->n_compounds, cudaMemcpyDeviceToHost);
	if (simulation->n_solvents > 0)
		cudaMemcpy(nlist_data_collection->solvents, simulation->box->solvents, sizeof(Solvent) * simulation->n_solvents, cudaMemcpyDeviceToHost);
}

void Engine::pushNlistsToDevice() {
	cudaMemcpy(simulation->box->compound_neighborlists, nlist_data_collection->compound_neighborlists, sizeof(NeighborList) * simulation->n_compounds, cudaMemcpyHostToDevice);
	cudaMemcpy(simulation->box->solvent_neighborlists, nlist_data_collection->solvent_neighborlists, sizeof(NeighborList) * simulation->n_solvents, cudaMemcpyHostToDevice);
	updated_neighborlists_ready = 0;
}

/*
void Engine::offloadPositionData(Simulation* simulation) {
	//cudaMemcpyAsync(compoundstatearray_host, simulation->box->compound_state_array, sizeof(CompoundState) * simulation->box->n_compounds, cudaMemcpyDeviceToHost);

	const int step_offset = (simulation->getStep() - STEPS_PER_LOGTRANSFER) * simulation->total_particles_upperbound;	// Tongue in cheek here, i think this is correct...
	cudaMemcpy(&nlist_data_collection->compoundstates[step_offset], simulation->box->compound_state_array, sizeof(CompoundState) * simulation->n_compounds, cudaMemcpyDeviceToHost);
	cudaMemcpy(&nlist_data_collection->solvents[step_offset], simulation->box->solvents, sizeof(Solvent) * simulation->n_solvents, cudaMemcpyDeviceToHost);
}
*/



bool Engine::neighborWithinCutoff(Float3* pos_a, Float3* pos_b, float cutoff_offset) {		// This is used for compounds with a confining_particle_sphere from key_particle BEFORE CUTOFF begins
	Float3 pos_b_temp = *pos_b;
	LIMAENG::applyHyperpos(pos_a, &pos_b_temp);
	float dist = (*pos_a - pos_b_temp).len();
	return (dist < (CUTOFF + cutoff_offset));
}



void Engine::cullDistantNeighbors(Simulation* simulation, NListDataCollection* nlist_data_collection) {	
	for (int id_self = 0; id_self < nlist_data_collection->n_compounds; id_self++) {
		NeighborList* nlist_self = &nlist_data_collection->compound_neighborlists[id_self];
		float cutoff_add_self = simulation->compounds_host[id_self].confining_particle_sphere;



		for (int j = 0; j < nlist_self->n_compound_neighbors; j++) {		// Cull compound-compound
			int id_neighbor = nlist_self->neighborcompound_ids[j];
			NeighborList* nlist_neighbor = &nlist_data_collection->compound_neighborlists[id_neighbor];
			float cutoff_add_neighbor = simulation->compounds_host[id_neighbor].confining_particle_sphere;

			if (id_self < id_neighbor) {
				if (!neighborWithinCutoff(&nlist_data_collection->compound_key_positions[id_self], &nlist_data_collection->compound_key_positions[id_neighbor], cutoff_add_self + cutoff_add_neighbor + CUTOFF)) {
					nlist_self->removeId(id_neighbor, NeighborList::NEIGHBOR_TYPE::COMPOUND);
					nlist_neighbor->removeId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);				
					j--;	// Decrement, as the removeId puts the last element at the current and now vacant spot.
				}
			}
		}


		for (int j = 0; j < nlist_self->n_solvent_neighbors; j++) {			// Cull compound-solvent
			int id_neighbor = nlist_self->neighborsolvent_ids[j];
			NeighborList* nlist_neighbor = &nlist_data_collection->solvent_neighborlists[id_neighbor];

			//printf("Dist: %f\n", (nlist_data_collection->compound_key_positions[id_self] - nlist_data_collection->solvent_positions[id_neighbor]).len());
			if (!neighborWithinCutoff(&nlist_data_collection->compound_key_positions[id_self], &nlist_data_collection->solvent_positions[id_neighbor], cutoff_add_self + CUTOFF) && false ) {
				nlist_self->removeId(id_neighbor, NeighborList::NEIGHBOR_TYPE::SOLVENT);
				//	printf("J: %d\n", j);
				nlist_neighbor->removeId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);
				j--;	// Decrement, as the removeId puts the last element at the current and now va-cant spot.
			}
		}
	}


	for (int id_self = 0; id_self < nlist_data_collection->n_solvents; id_self++) {																// Cull solvent-solvent
		NeighborList* nlist_self = &nlist_data_collection->solvent_neighborlists[id_self];

		int cnt = 0;

		for (int j = 0; j < nlist_self->n_solvent_neighbors; j++) {			/// NOT FINISHED HERE
			int id_neighbor = nlist_self->neighborsolvent_ids[j];
			NeighborList* nlist_neighbor = &nlist_data_collection->solvent_neighborlists[id_neighbor];

			if (!neighborWithinCutoff(&nlist_data_collection->solvent_positions[id_self], &nlist_data_collection->solvent_positions[id_neighbor], CUTOFF)) {
				cnt++;
				if (!nlist_self->removeId(id_neighbor, NeighborList::NEIGHBOR_TYPE::SOLVENT))
					printf("J1: %d id_self %d id_neighbor %d    cnt %d\n", j, id_self, id_neighbor, cnt);
				if (!nlist_neighbor->removeId(id_self, NeighborList::NEIGHBOR_TYPE::SOLVENT)) {
					printf("J2: %d of %d.   id_self %d id_neighbor %d count: %d\n", j, nlist_self->n_solvent_neighbors, id_self, id_neighbor, cnt);
					for (int i = 0; i < nlist_self->n_solvent_neighbors; i++) {
						printf("neighbor %d\n", nlist_self->neighborsolvent_ids[i]);
					}
					printf("\n\n\n");
					exit(1);
				}
					

				j--;	// Decrement, as the removeId puts the last element at the current and now vacant spot.
			}			
		}
	}
}


void Engine::updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection, volatile bool* finished, int* timing, bool* mutex_lock) {	// This is a thread worker-function, so it can't own the object, thus i pass a ref to the engine object..
	auto t0 = std::chrono::high_resolution_clock::now();
	Int3 before(nlist_data_collection->compound_neighborlists[0].n_compound_neighbors, nlist_data_collection->compound_neighborlists[0].n_solvent_neighbors, 0);
	//Int3 before(nlist_data_collection->solvent_neighborlists[193].n_compound_neighbors, nlist_data_collection->solvent_neighborlists[193].n_solvent_neighbors, 0);
	//nlist_data_collection->preparePositionData();


	nlist_data_collection->preparePositionData(simulation->compounds_host);		// Makes key positions addressable in arrays: compound_key_positions and solvent_positions
	// First do culling of neighbors that has left CUTOFF
	cullDistantNeighbors(simulation, nlist_data_collection);


	// First add compound->solvent, compound->compound
	for (uint16_t id_self = 0; id_self < simulation->n_compounds; id_self++) {

		NeighborList* nlist_self = &nlist_data_collection->compound_neighborlists[id_self];
		HashTable hashtable_compoundneighbors(nlist_self->neighborcompound_ids, nlist_self->n_compound_neighbors, NEIGHBORLIST_MAX_COMPOUNDS * 2);
		HashTable hashtable_solventneighbors(nlist_self->neighborsolvent_ids, nlist_self->n_solvent_neighbors, NEIGHBORLIST_MAX_SOLVENTS * 2);
		float cutoff_add_self = simulation->compounds_host[id_self].confining_particle_sphere;

		//printf("\nSize to hashtable: %d\n", nlist_self->n_solvent_neighbors);
		//printf("cutoff %f\n", cutoff_offset);

		// Go through all solvents in box!
		for (uint16_t id_candidate = 0; id_candidate < simulation->n_solvents; id_candidate++) {
			NeighborList* nlist_candidate = &nlist_data_collection->solvent_neighborlists[id_candidate];
			//continue;
			if (neighborWithinCutoff(&nlist_data_collection->compound_key_positions[id_self], &nlist_data_collection->solvent_positions[id_candidate], cutoff_add_self + CUTOFF)) {
				if (hashtable_solventneighbors.insert(id_candidate)) {
					nlist_self->addId(id_candidate, NeighborList::NEIGHBOR_TYPE::SOLVENT);
					nlist_candidate->addId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);
				}
			}

		}

		// Go through all compounds in box, with higher ID than self!
		for (uint16_t id_candidate = id_self + 1; id_candidate < simulation->n_compounds; id_candidate++) {	// For finding new nearby compounds, it is faster and simpler to just check all compounds, since there are so few
			NeighborList* nlist_candidate = &nlist_data_collection->compound_neighborlists[id_candidate];
			float cutoff_add_candidate = simulation->compounds_host[id_self].confining_particle_sphere;
			//continue;
			//printf("Distance to neighbor compound: %f\n", (nlist_data_collection->compound_key_positions[id_self] - nlist_data_collection->compound_key_positions[id_candidate]).len());
			//if (id_self == 0)
				//printf("Adding compound %d to %d\n", id_candidate, id_self);

			if (neighborWithinCutoff(&nlist_data_collection->compound_key_positions[id_self], &nlist_data_collection->compound_key_positions[id_candidate], cutoff_add_self + cutoff_add_candidate + CUTOFF)) {
				if (hashtable_compoundneighbors.insert(id_candidate)) {
					nlist_self->addId(id_candidate, NeighborList::NEIGHBOR_TYPE::COMPOUND);
					nlist_candidate->addId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);
				}
			}
		}
	}

	// Finally add all solvent->solvent candidates
	for (int id_self = 0; id_self < simulation->n_solvents; id_self++) {
		NeighborList* nlist_self = &nlist_data_collection->solvent_neighborlists[id_self];
		HashTable hashtable_solventneighbors(nlist_self->neighborsolvent_ids, (int)nlist_self->n_solvent_neighbors, NEIGHBORLIST_MAX_SOLVENTS * 2);



		for (int id_candidate = id_self + 1; id_candidate < simulation->n_solvents; id_candidate++) {
			NeighborList* nlist_candidate = &nlist_data_collection->solvent_neighborlists[id_candidate];
			if (neighborWithinCutoff(&nlist_data_collection->solvent_positions[id_self], &nlist_data_collection->solvent_positions[id_candidate], CUTOFF)) {
				if (hashtable_solventneighbors.insert(id_candidate)) {
					nlist_self->addId(id_candidate, NeighborList::NEIGHBOR_TYPE::SOLVENT);
					nlist_candidate->addId(id_self, NeighborList::NEIGHBOR_TYPE::SOLVENT);
					//printf("\n");
					if (id_self == 193) {
						//printf("Adding id %d with dist %f\n", id_candidate, (nlist_data_collection->solvent_positions[id_self] - nlist_data_collection->solvent_positions[id_candidate]).len());
					}
				}
			}
			
		}
	}

	Int3 after(nlist_data_collection->compound_neighborlists[0].n_compound_neighbors, nlist_data_collection->compound_neighborlists[0].n_solvent_neighbors, 0);
	//Int3 after(nlist_data_collection->solvent_neighborlists[193].n_compound_neighbors, nlist_data_collection->solvent_neighborlists[193].n_solvent_neighbors, 0);

	//printf("\nEntity went from %d %d neighbors to %d %d\n", before.x, before.y, after.x, after.y);
	
	auto t1 = std::chrono::high_resolution_clock::now();
	*timing = (int) std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

	*finished = 1;		// Thread terminates here!
	*mutex_lock = 0;	// Unlock
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
	LIMAENG::genericErrorCheck("Cuda error during traindata offloading\n");
}



Float3 Engine::getBoxTemperature() {
	
	const int step = simulation->getStep() - 1;

	long double temp_sum = 0;				// [k]

	long double kinE_sum = 0;
	float biggest_contribution = 0;


	const int step_offset_a = step * simulation->total_particles_upperbound;
	const int step_offset_b = (step - 2) * simulation->total_particles_upperbound;
	const int solvent_offset = MAX_COMPOUND_PARTICLES * simulation->n_compounds;


	for (int c = 0; c < simulation->n_compounds; c++) {
		int compound_offset = c * MAX_COMPOUND_PARTICLES;
		for (int i = 0; i < simulation->compounds_host[c].n_particles; i++) {	// i gotta move this somewhere else....

			Float3 posa = simulation->traj_buffer[i + compound_offset + step_offset_a];
			Float3 posb = simulation->traj_buffer[i + compound_offset + step_offset_b];
			float kinE = LIMAENG::calcKineticEnergy(&posa, &posb, forcefield_host.particle_parameters[simulation->compounds_host[c].atom_types[i]].mass, simulation->dt);			// Doesnt work, use forcefield_host!!
			//float kinE = LIMAENG::calcKineticEnergy(&posa, &posb, forcefield_device.particle_parameters[simulation->box->compounds[c].atom_types[i]].mass, simulation->dt);			// Doesnt work, use forcefield_host!!
			//printf("kinE %f\n", kinE);
			//printf("mass %f\n", forcefield_host.particle_parameters[simulation->compounds_host[c].atom_types[i]].mass);
			biggest_contribution = max(biggest_contribution, kinE);
				
			kinE_sum += kinE;
			//particles_total++;
		}
	}
	//printf("\nKin e %Lf\n", kinE_sum);
	//kinE_sum = 0;
	for (int i = 0; i < simulation->n_solvents; i++) {
		Float3 posa = simulation->traj_buffer[i + solvent_offset + step_offset_a];
		Float3 posb = simulation->traj_buffer[i + solvent_offset + step_offset_b];
		float kinE = LIMAENG::calcKineticEnergy(&posa, &posb, SOLVENT_MASS, simulation->dt);
		biggest_contribution = max(biggest_contribution, kinE);
		kinE_sum += kinE;
	}
	//double avg_kinE = kinE_sum / (long double)particles_total;
	double avg_kinE = kinE_sum / (long double)simulation->total_particles;
	float temperature = avg_kinE * 2.f / (3.f * 8.3145);
	//printf("\nTemp: %f\n", temperature);
	return Float3(temperature, biggest_contribution, avg_kinE);
}

																																			// THIS fn requires mallocmanaged!!   // HARD DISABLED HERE
void Engine::handleBoxtemp() {
	const float target_temp = 310.f;				// [k]
	Float3 temp_package = getBoxTemperature();
	float temp = temp_package.x;
	float biggest_contribution = temp_package.y;

	temp = temp == 0.f ? 1 : temp;																			// So we avoid dividing by 0
	simulation->temperature_buffer[simulation->n_temp_values++] = temp;

	float temp_scalar = target_temp / temp;
		// I just added this to not change any temperatures too rapidly
	temp_scalar = min(temp_scalar, 1.01f);
	temp_scalar = max(temp_scalar, 0.99f);

	if (PRINT_TEMP || temp > 500.f || temp < 100.f) { printf("\n %d Temperature: %.1f Biggest contrib: %.0f avg kinE %.0f\n", (simulation->getStep() - 1) / STEPS_PER_THERMOSTAT, temp, biggest_contribution, temp_package.z); }
		
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
	LIMAENG::genericErrorCheck("Error during step or state_transfer\n");		// Temp, we want to do host stuff while waiting for async GPU operations...	// SLOW
	auto t2 = std::chrono::high_resolution_clock::now();


	simulation->incStep();


	int force_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	int copy_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	timings = timings + Int3(force_duration, copy_duration, 0);
}





// ------------------------------------------------------------------------------------------- DEVICE FUNCTIONS -------------------------------------------------------------------------------------------//



__device__ void applyPBC(Float3* current_position) {	// Only changes position if position is outside of box;
	for (int dim = 0; dim < 3; dim++) {
		*current_position->placeAt(dim) += BOX_LEN * (current_position->at(dim) < 0.f);
		*current_position->placeAt(dim) -= BOX_LEN * (current_position->at(dim) > BOX_LEN);
	}
}



/*
constexpr double sigma = 0.3923f;										// nm, point at which d/dE or force is 0.
constexpr double epsilon = 0.5986 * 1'000.f;							// J/mol
*/

__device__ inline float calcSigma(uint8_t atomtype1, uint8_t atomtype2) {
	return (forcefield_device.particle_parameters[atomtype1].sigma + forcefield_device.particle_parameters[atomtype2].sigma) * 0.5;
}
__device__ inline float calcEpsilon(uint8_t atomtype1, uint8_t atomtype2) {
	return sqrtf(forcefield_device.particle_parameters[atomtype1].epsilon * forcefield_device.particle_parameters[atomtype2].epsilon);
}

__device__ Float3 calcLJForce(Float3* pos0, Float3* pos1, float* data_ptr, float* potE, float sigma, float epsilon, int type1=-1, int type2=-1) {
	// Calculates LJ force on p0	(attractive to p1. Negative values = repulsion )//
	// input positions in cartesian coordinates [nm]
	// sigma [nm]
	// epsilon [J/mol]->[(kg*nm^2)/(ns^2*mol)]
	//	Returns force in J/mol*M		?????????????!?!?//


	//return Float3(sigma);	// Test for forcing LJ calc, but always small force, even with bonded particles!


	// Directly from book
	float dist_sq = (*pos1 - *pos0).lenSquared();
	float s = sigma * sigma / dist_sq;								// [nm^2]/[nm^2] -> unitless
	s = s * s * s;
	float force_scalar = 24.f * epsilon * s / dist_sq * (1.f - 2.f * s);	// Attractive. Negative, when repulsive		[(kg*nm^2)/(nm^2*ns^2*mol)] ->----------------------	[(kg)/(ns^2*mol)]	
	//float force_scalar = 1.f;

	//if (force_scalar != force_scalar) {
	//	printf("Error is here\n\n");
	//}



	//if (blockIdx.x == 12 && threadIdx.x == 27 && ddd.len() > 1000000)
		//ddd.print('f');
#ifndef ENABLE_BLJV
	return (*pos1 - *pos0) * sigma * epsilon;				// This line forces the CUDA compile to include all calculations, but can never produce a large force.
#endif

#ifdef LIMA_VERBOSE
	Float3 ddd = (*pos1 - *pos0) * force_scalar;
	if (ddd.x != ddd.x) {
		printf("\nError is here\n");
		printf("Scalar %f dist_sq %f\n", force_scalar, dist_sq);
		pos0->print('0');
		pos0->print('1');
	}
	//if (threadIdx.x == 32 && blockIdx.x == 28 && force < -600000.f && force > -700000) {
	//if (abs(force) > 20e+8 && type1 == 200 && type2 == 1) {
	if (abs(force_scalar) > 20e+9 || *potE != *potE ) {
		//printf("\nDist %f   D2 %f    force %.1f [MN]     sigma: %.3f \tepsilon %.0f  t1 %d t2 %d\n", (float)sqrt(dist_sq), (float) dist_sq, (float)force_scalar *1e-6, sigma, epsilon, type1, type2);
		//pos0->print('1');
		//pos1->print('2');
		//printf("Block %d Thread %d\n", blockIdx.x, threadIdx.x);
		//printf("Thread %d Block %d self %f %f %f other %f %f %f\n", threadIdx.x, blockIdx.x, pos0->x, pos0->y, pos0->z, pos1->x, pos1->y, pos1->z);
		//printf("Force %f %f %f\n", force_unit_vector.x * force, force_unit_vector.y * force, force_unit_vector.z * force);
	}

	{
		//Float3 force_vec = force_unit_vector * force;
		/*if (force_vec.x != force_vec.x) {
			//printf("Force: %f\n", force);
			//force_unit_vector.print('u');
			printf("Thread %d Block %d self %f %f %f other %f %f %f\n", threadIdx.x, blockIdx.x, pos0->x, pos0->y, pos0->z, pos1->x, pos1->y, pos1->z);

		}*/
	}
#endif
	return (*pos1 - *pos0) * force_scalar;										// GN/mol [(kg*nm)/(ns^2*mol)]
	//return Float3(0.);
}

//constexpr double kb = 17.5 * 1e+6;		//	J/(mol*nm^2)	/ kg/(ns^2 * mol)
__device__ void calcPairbondForces(Float3* pos_a, Float3* pos_b, PairBond* bondtype, Float3* results, float* potE) {
	// Calculates bond force on both particles					//
	// Calculates forces as J/mol*M								//

	Float3 difference = *pos_a - *pos_b;						//	[nm]
	double error = difference.len() - bondtype->b0;			//	[nm]

	*potE += 0.5 * bondtype->kb * (error * error);					// [J/mol]

	difference = difference.norm();								// dif_unit_vec, but shares register with dif
	double force_scalar = -bondtype->kb * error;				//	[J/(mol*nm^2)]*nm =>	kg*nm^2*ns^-2/(mol*nm^2)*nm = kg*nm/(mol*ns^2)				 [N/mol] directionless, yes?

	results[0] = difference * force_scalar;				// [GN]
	results[1] = difference * force_scalar * -1;		// [GN]

#ifdef LIMA_VERBOSE
	if (force_scalar > 2e+7) {
		printf("thread %d  error %f ref %f force %f\n", threadIdx.x, error, bondtype->b0, force_scalar);
		//pos_a->print('a');
		//pos_b->print('b');
	}
#endif
		
}


//constexpr double ktheta = 65 * 1e+3;	// J/mol
__device__ void calcAnglebondForces(Float3* pos_left, Float3* pos_middle, Float3* pos_right, AngleBond* angletype, Float3* results, float* potE) {
	Float3 v1 = *pos_left - *pos_middle;
	Float3 v2 = *pos_right - *pos_middle;
	Float3 normal1 = v1.cross(v2);
	Float3 normal2 = v2.cross(v1);

	Float3 inward_force_direction1 = (v1.cross(normal2)).norm();
	Float3 inward_force_direction2 = (v2.cross(normal1)).norm();

	double angle = Float3::getAngle(v1, v2);
	double error = angle - angletype->theta_0;// *reference_angle;

	//error = (error / (2.f * PI)) * 360.f;

	//*potE += 0.5 * ktheta * error * error * 0.5;
	//double force_scalar = ktheta * (error);
	*potE += 0.5 * angletype->k_theta * error * error * 0.5;
	double force_scalar = angletype->k_theta * (error);

	results[0] = inward_force_direction1 * force_scalar;
	results[2] = inward_force_direction2 * force_scalar;
	results[1] = (results[0] + results[2]) * -1;


	//printf("\nangle %f error %f force %f t0 %f kt %f\n", angle, error, force_scalar, angletype->theta_0, angletype->k_theta);

}
__device__ void calcDihedralbondForces(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	Float3 normal1 = (*pos_left-*pos_lm).cross((*pos_rm-*pos_lm)).norm();		// Should this not be normalized????
	Float3 normal2 = (*pos_lm-*pos_rm).cross((*pos_right-*pos_rm)).norm();			// Is inward or outward? Refactor!!!
	// Both vectors point "left". If pos_left+n1 is closer to pos_right than pos_left-n1, then n1 is pointing inwards, and n2 outwards. Also vice-versa.
	float torsion = Float3::getAngle(normal1, normal2);




	if (((*pos_left + normal1) - *pos_right).len() < ((*pos_left - normal1) - *pos_right).len())	// Wait, can't i just check if torsion > pi | torsion < 0`???????????????????????
		normal2 *= -1;
	else
		normal1 *= -1;
	// Now both normals are pointing inwards

	//printf("Torsion %f\n", normal1.len());
	//float error = (torsion - dihedral->phi_0);	// *360.f / 2.f * PI;	// Check whether k_phi is in radians or degrees

	//error = (error / (2.f * PI)) * 360.f;

	
	//*potE += 0.5 * dihedral->k_phi * error * error * 0.5;
	//double force_scalar = dihedral->k_phi * (error);

	//double force_scalar = sinf(torsion - dihedral->phi_0);		// Should be -sinf? (cos)' = -sin??!?!

	float force_scalar = dihedral->k_phi * sinf(dihedral->n * torsion - dihedral->phi_0);
	*potE = dihedral->k_phi * (1 - cosf(dihedral->n * torsion - dihedral->phi_0));

	//printf("Torsion %f ref %f force_scalar %f\n", torsion, dihedral->phi_0, force_scalar);
	//force_scalar *= dihedral->k_phi;
	if (abs(force_scalar) > 183300) {
		pos_left->print('L');
		pos_lm->print('l');
		pos_rm->print('r');
		pos_right->print('R');
		//printf("torsion %f      ref %f     error %f     force: %f\n", torsion, dihedral->phi_0, error, force_scalar);
		printf("torsion %f      ref %f     error %f     force: %f\n", torsion, dihedral->phi_0, 0.f, force_scalar);
	}
		



	//results[0] = normal1 * force_scalar * -1.f;
	//results[3] = normal2 * force_scalar;
	results[0] = normal1 * force_scalar;
	results[3] = normal2 * force_scalar;
	// Not sure about the final two forces, for now we'll jsut split the sum of opposite forces between them.
	results[1] = (results[0] + results[3]) * -1.f * 0.5;
	results[2] = (results[0] + results[3]) * -1.f * 0.5;


}




//__device__ Float3 computeIntermolecularLJForces(Float3* self_pos, int n_particles, Float3* positions, float* data_ptr, float* potE_sum, uint8_t atomtype_self, uint8_t* atomtypes_others, LJ_Ignores* lj_ignore_list, uint32_t* global_ids) {	// Assumes all positions are 
__device__ Float3 computerIntermolecularLJForces(Float3* self_pos, uint8_t atomtype_self, LJ_Ignores* lj_ignore_list, float* potE_sum, uint32_t global_id_self, float* data_ptr,
	Compound* neighbor_compound, Float3* neighbor_positions, int neighborcompound_id) {
	
	
	
	Float3 force(0.f);
	//Float3 a;
	for (int neighborparticle_id = 0; neighborparticle_id < neighbor_compound->n_particles; neighborparticle_id++) {

#ifdef ENABLE_BLJV
		//if (!lj_ignore_list->ignore((uint8_t)neighborparticle_id, (uint8_t)neighborcompound_id))	// Check if particle_self is bonded to particle_i in neighbor compound
			//continue;
#endif

		int neighborparticle_atomtype = neighbor_compound->atom_types[neighborparticle_id];	//// TEMPORARY, this is waaaay to many global mem accesses

		
		//continue;	// Tester
		force += calcLJForce(self_pos, &neighbor_positions[neighborparticle_id], data_ptr, potE_sum,
			calcSigma(atomtype_self, neighborparticle_atomtype), calcEpsilon(atomtype_self, neighborparticle_atomtype),
			//(forcefield_device.particle_parameters[atomtype_self].sigma + forcefield_device.particle_parameters[neighborparticle_atomtype].sigma) * 0.5f,
			//(forcefield_device.particle_parameters[atomtype_self].epsilon + forcefield_device.particle_parameters[neighborparticle_atomtype].epsilon) * 0.5,
			//atomtype_self, atomtypes_others[i]
			global_id_self, neighbor_compound->particle_global_ids[neighborparticle_id]
		);

	}	
	return force;// *24.f * 1e-9;
}

__device__ Float3 computeIntramolecularLJForces(Compound* compound, CompoundState* compound_state, float* potE_sum, float* data_ptr, BondedParticlesLUT& bonded_particles_lut) {
	Float3 force(0.f);
	for (int i = 0; i < compound->n_particles; i++) {
#ifndef INWD
		if (i != threadIdx.x && threadIdx.x < compound->n_particles) {
		//if (i != threadIdx.x && !compound->lj_ignore_list[threadIdx.x].ignore((uint8_t)i, (uint8_t)blockIdx.x) && threadIdx.x < compound->n_particles) {
			//printf("Thread %d computing LJ to index %d\n", threadIdx.x, i);
#else
		if (i != threadIdx.x && !*bonded_particles_lut.get(threadIdx.x, i)) {
#endif
			//if (i != threadIdx.x  && threadIdx.x < compound.n_particles) {																											// DANGER



			force += calcLJForce(&compound_state->positions[threadIdx.x], &compound_state->positions[i], data_ptr, potE_sum,
				//(forcefield_device.particle_parameters[compound->atom_types[threadIdx.x]].sigma + forcefield_device.particle_parameters[compound->atom_types[i]].sigma) * 0.5f,		// Don't know how to handle atoms being this close!!!
				//sqrtf(forcefield_device.particle_parameters[compound->atom_types[threadIdx.x]].epsilon * forcefield_device.particle_parameters[compound->atom_types[i]].epsilon),
				// 
				calcSigma(compound->atom_types[threadIdx.x], compound->atom_types[i]), calcEpsilon(compound->atom_types[threadIdx.x], compound->atom_types[i]),
				// 
				//compound.atom_types[threadIdx.x], compound.atom_types[i]
				//compound.particle_global_ids[threadIdx.x], blockIdx.x
				//threadIdx.x, i
				compound->particle_global_ids[threadIdx.x], compound->particle_global_ids[i]
			);// *24.f * 1e-9;
		}
	}
	return force;
}





__device__ Float3 computeSolventToSolventLJForces(Float3* self_pos, NeighborList* nlist, Solvent* solvents, float* data_ptr, float* potE_sum) {	// Specific to solvent kernel
	Float3 force(0.f);
	for (int i = 0; i < nlist->n_solvent_neighbors; i++) {
		Solvent neighbor = solvents[nlist->neighborsolvent_ids[i]];			
		LIMAENG::applyHyperpos(&neighbor.pos, self_pos);
		force += calcLJForce(self_pos, &neighbor.pos, data_ptr, potE_sum, 
			forcefield_device.particle_parameters[ATOMTYPE_SOL].sigma, 
			forcefield_device.particle_parameters[ATOMTYPE_SOL].epsilon);
		/*if (abs(force.len()) > 1000000) {
			(*self_pos - neighbor.pos).print('d');
			printf("F %f %f %f n %d  solvent_id %d neighbor_id %d\n", force.x, force.y, force.z, nlist->n_solvent_neighbors, threadIdx.x + blockIdx.x*blockDim.x, nlist->neighborsolvent_ids[i]);
		}*/
	}
	return force;// *24.f * 1e-9;
}
__device__ Float3 computeSolventToCompoundLJForces(Float3* self_pos, int n_particles, Float3* positions, float* data_ptr, float* potE_sum, uint8_t atomtype_self) {	// Specific to solvent kernel
	Float3 force(0, 0, 0);
	for (int i = 0; i < n_particles; i++) {
		force += calcLJForce(self_pos, &positions[i], data_ptr, potE_sum,
			calcSigma(atomtype_self, ATOMTYPE_SOL), calcEpsilon(atomtype_self, ATOMTYPE_SOL)
			//(forcefield_device.particle_parameters[atomtype_self].sigma + forcefield_device.particle_parameters[ATOMTYPE_SOL].sigma) * 0.5f,
			//sqrtf(forcefield_device.particle_parameters[atomtype_self].epsilon * forcefield_device.particle_parameters[ATOMTYPE_SOL].epsilon)
		);
	}
	return force;// *24.f * 1e-9;
}
__device__ Float3 computeCompoundToSolventLJForces(Float3* self_pos, int n_particles, Float3* positions, float* data_ptr, float* potE_sum, uint8_t atomtype_self, uint8_t* atomtypes_others) {	// Assumes all positions are 
	Float3 force(0.f);
	for (int i = 0; i < n_particles; i++) {
		force += calcLJForce(self_pos, &positions[i], data_ptr, potE_sum,
			(forcefield_device.particle_parameters[atomtype_self].sigma + forcefield_device.particle_parameters[atomtypes_others[i]].sigma) * 0.5f,
			(forcefield_device.particle_parameters[atomtype_self].epsilon + forcefield_device.particle_parameters[atomtypes_others[i]].epsilon) * 0.5,
			atomtype_self, atomtypes_others[i]
		);
	}
	return force;// *24.f * 1e-9;
}

/*
__device__ Float3 computeCompoundToSolventLJForces(Float3* self_pos, int n_particles, Float3* positions, double* data_ptr, double* potE_sum) {	// Specific to solvent kernel
	Float3 force(0, 0, 0);
	for (int i = 0; i < n_particles; i++) {
		Float3 hyperpos = positions[i];			// copy, DONT ref it as all threads will cann applyHyperpos
		LIMAENG::applyHyperpos(self_pos, &hyperpos);
		force += calcLJForce(self_pos, &hyperpos, data_ptr, potE_sum);	
	}
	return force;
}
*/
template <typename T>	// Can either be Compound or CompoundBridgeCompact
__device__ Float3 computePairbondForces(T* entity, Float3* positions, Float3* utility_buffer, float* potE) {	// only works if n threads >= n bonds
//__device__ Float3 computePairbondForces(CompoundBridgeCompact* entity, Float3* positions, Float3* utility_buffer, float* potE) {	// only works if n threads >= n bonds
	utility_buffer[threadIdx.x] = Float3(0.f);
	for (int bond_offset = 0; (bond_offset * blockDim.x) < entity->n_singlebonds; bond_offset++) {
		PairBond* pb = nullptr;
		Float3 forces[2];
		int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < entity->n_singlebonds) {
			pb = &entity->singlebonds[bond_index];

			calcPairbondForces(
				&positions[pb->atom_indexes[0]],
				&positions[pb->atom_indexes[1]],
				//&pb->reference_dist,
				pb,
				forces, potE
			);

			//forces[0] = Float3(0.f); forces[1] = Float3(0.f);
#ifndef INWD
			if (pb->invertLJ) {
				float temp = 0.f;
				//float sigma = (forcefield_device.particle_parameters[entity->atom_types[pb->atom_indexes[0]]].sigma + forcefield_device.particle_parameters[entity->atom_types[pb->atom_indexes[1]]].sigma) * 0.5;
				
				float epsilon = sqrtf(forcefield_device.particle_parameters[entity->atom_types[pb->atom_indexes[0]]].epsilon * forcefield_device.particle_parameters[entity->atom_types[pb->atom_indexes[1]]].epsilon);
				//Float3 anti_lj_force = calcLJForce(&positions[pb->atom_indexes[0]], &positions[pb->atom_indexes[1]], &temp, potE, sigma, epsilon, -1,-1);
				Float3 anti_lj_force = calcLJForce(&positions[pb->atom_indexes[0]], &positions[pb->atom_indexes[1]], &temp, potE, 
					calcSigma(entity->atom_types[pb->atom_indexes[0]], entity->atom_types[pb->atom_indexes[1]]), 
					calcEpsilon(entity->atom_types[pb->atom_indexes[0]], entity->atom_types[pb->atom_indexes[1]]),
					-1, -1);
				forces[0] -= anti_lj_force; forces[1] += anti_lj_force;
				if (blockIdx.x == 0 && (pb->atom_indexes[0] == 1 || pb->atom_indexes[1] == 1)) {
					//printf("ANTI LJ %d %d\n", pb->atom_indexes[0], pb->atom_indexes[1]);
					//anti_lj_force.print('A');
				}

				//if (entity->particle_refs[pb->atom_indexes[0]].global_id == 649 || entity->particle_refs[pb->atom_indexes[1]].global_id == 649) {
					//anti_lj_force.print('B');
					//printf("B dist %.08f\n", (positions[pb->atom_indexes[0]] - positions[pb->atom_indexes[1]]).len());
				//}
					
				
			}
#endif
		}


		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && pb != nullptr) {
				utility_buffer[pb->atom_indexes[0]] += forces[0];
				utility_buffer[pb->atom_indexes[1]] += forces[1];
			}
			__syncthreads();
		}
	}
	return utility_buffer[threadIdx.x];
}

template <typename T>	// Can either be Compound or CompoundBridgeCompact
__device__ Float3 computeAnglebondForces(T* entity, Float3* positions, Float3* utility_buffer, float* potE) {
	utility_buffer[threadIdx.x] = Float3(0, 0, 0);
	for (int bond_offset = 0; (bond_offset * blockDim.x) < entity->n_anglebonds; bond_offset++) {
		AngleBond* ab = nullptr;
		Float3 forces[3];
		int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < entity->n_anglebonds) {
			ab = &entity->anglebonds[bond_index];


			calcAnglebondForces(
				&positions[ab->atom_indexes[0]],
				&positions[ab->atom_indexes[1]],
				&positions[ab->atom_indexes[2]],
				//&ab->reference_angle,
				ab,
				forces, potE
			);

			//forces[0] = Float3(0.f); forces[1] = Float3(0.f); forces[2] = Float3(0.f);
#ifndef INWD
			if (ab->invertLJ) {
				float temp = 0.f;
				//float sigma = (forcefield_device.particle_parameters[entity->atom_types[ab->atom_indexes[0]]].sigma + forcefield_device.particle_parameters[entity->atom_types[ab->atom_indexes[2]]].sigma) * 0.5;
				//float epsilon = sqrtf(forcefield_device.particle_parameters[entity->atom_types[ab->atom_indexes[0]]].epsilon * forcefield_device.particle_parameters[entity->atom_types[ab->atom_indexes[2]]].epsilon);
				//Float3 anti_lj_force = calcLJForce(&positions[ab->atom_indexes[0]], &positions[ab->atom_indexes[2]], &temp, potE, sigma, epsilon, -1, -1);
				Float3 anti_lj_force = calcLJForce(&positions[ab->atom_indexes[0]], &positions[ab->atom_indexes[2]], &temp, potE, 
					calcSigma(entity->atom_types[ab->atom_indexes[0]], entity->atom_types[ab->atom_indexes[2]]), 
					calcEpsilon(entity->atom_types[ab->atom_indexes[0]], entity->atom_types[ab->atom_indexes[2]]),
					-1, -1);
				forces[0] -= anti_lj_force; forces[2] += anti_lj_force;
				if (blockIdx.x == 0 && (ab->atom_indexes[0] == 1 || ab->atom_indexes[2] == 1)) {
					//printf("ANTI LJ %d %d %d\n", ab->atom_indexes[0], ab->atom_indexes[1], ab->atom_indexes[2]);
					//anti_lj_force.print('A');
				}
			}
#endif
		}


		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && ab != nullptr) {
				utility_buffer[ab->atom_indexes[0]] += forces[0];
				utility_buffer[ab->atom_indexes[1]] += forces[1];
				utility_buffer[ab->atom_indexes[2]] += forces[2];
			}
			__syncthreads();
		}
	}

	return utility_buffer[threadIdx.x];
}

template <typename T>	// Can either be Compound or CompoundBridgeCompact
__device__ Float3 computeDihedralForces(T* entity, Float3* positions, Float3* utility_buffer, float* potE) {
	utility_buffer[threadIdx.x] = Float3(0, 0, 0);
	for (int bond_offset = 0; (bond_offset * blockDim.x) < entity->n_dihedrals; bond_offset++) {
		DihedralBond* db = nullptr;
		Float3 forces[4];
		int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < entity->n_dihedrals) {
			db = &entity->dihedrals[bond_index];
			//printf("Firing %d of %d\n", bond_index, entity);
			calcDihedralbondForces(
				&positions[db->atom_indexes[0]],
				&positions[db->atom_indexes[1]],
				&positions[db->atom_indexes[2]],
				&positions[db->atom_indexes[3]],
				db,
				forces, 
				potE
			);

#ifndef INWD
			//forces[0] = Float3(0.f); forces[1] = Float3(0.f); forces[2] = Float3(0.f); forces[3] = Float3(0.f);
			if (db->invertLJ) {
				float temp = 0.f;
				//float sigma = (forcefield_device.particle_parameters[entity->atom_types[db->atom_indexes[0]]].sigma + forcefield_device.particle_parameters[entity->atom_types[db->atom_indexes[3]]].sigma) * 0.5;
				//float epsilon = sqrtf(forcefield_device.particle_parameters[entity->atom_types[db->atom_indexes[0]]].epsilon * forcefield_device.particle_parameters[entity->atom_types[db->atom_indexes[3]]].epsilon);
				//Float3 anti_lj_force = calcLJForce(&positions[db->atom_indexes[0]], &positions[db->atom_indexes[3]], &temp, potE, sigma, epsilon, db->atom_indexes[0], db->atom_indexes[3]);
				Float3 anti_lj_force = calcLJForce(&positions[db->atom_indexes[0]], &positions[db->atom_indexes[3]], &temp, potE, 
					calcSigma(entity->atom_types[db->atom_indexes[0]], entity->atom_types[db->atom_indexes[3]]), 
					calcEpsilon(entity->atom_types[db->atom_indexes[0]], entity->atom_types[db->atom_indexes[3]]),
					db->atom_indexes[0], db->atom_indexes[3]);

				//Float3 anti_lj_force = calcLJForce(&positions[db->atom_indexes[0]], &positions[db->atom_indexes[3]], &temp, potE, sigma, epsilon, );
				forces[0] -= anti_lj_force; forces[3] += anti_lj_force;

			}
#endif

		}

		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && db != nullptr) {
				utility_buffer[db->atom_indexes[0]] += forces[0];
				utility_buffer[db->atom_indexes[1]] += forces[1];
				utility_buffer[db->atom_indexes[2]] += forces[2];
				utility_buffer[db->atom_indexes[3]] += forces[3];
			}
			__syncthreads();
		}
	}
	//if (utility_buffer[threadIdx.x].x != utility_buffer[threadIdx.x].x)
		//utility_buffer[threadIdx.x].print('d');
	return utility_buffer[threadIdx.x];
}

__device__ void integratePosition(Float3* pos, Float3* pos_tsub1, Float3* force, const double mass, const double dt, float* thermostat_scalar, int p_index, bool issolvent =false) {
	// Force is in ??Newton, [kg * nm /(mol*ns^2)] //

	float prev_vel = (*pos - *pos_tsub1).len();
		
	Float3 temp = *pos;
	LIMAENG::applyHyperpos(pos, pos_tsub1);
	//printf("dt %.08f\n", dt);
	//*pos = *pos * 2 - *pos_tsub1 + *force * (1.f / mass) * dt * dt;	// no *0.5?	// [nm] - [nm] + [kg*m*s^-2]/[kg] * [ns]^2
	//double magnitude_equalizer = 1e+4;		// Due to force being very large and dt being very small
	//double masstime_scaling = ((dt*magnitude_equalizer)*(dt * magnitude_equalizer)) / mass;

	//*pos = *pos * 2 - *pos_tsub1 + *force * masstime_scaling;		// [nm] - [nm] + [kg/mol*m*/s^2]/[kg/mol] * [s]^2 * (1e-9)^2	=> [nm]-[nm]+[]
	//pos = *pos * 2 - *pos_tsub1 + (*force * (1./(magnitude_equalizer * magnitude_equalizer))) * masstime_scaling;		// [nm] - [nm] + [kg/mol*m * /s^2]/[kg/mol] * [s]^2 * (1e-9)^2	=> [nm]-[nm]+[]
	*pos = *pos * 2 - *pos_tsub1 + *force * (dt / mass) * dt;		// [nm] - [nm] + [kg/mol*m*/s ^ 2] / [kg / mol] * [s] ^ 2 * (1e-9) ^ 2 = > [nm] - [nm] + []
	*pos_tsub1 = temp;

	
	return;
	Float3 delta_pos = *pos - *pos_tsub1;
	*pos = *pos_tsub1 + delta_pos * *thermostat_scalar;

	if (delta_pos.len() > 0.05)
		printf("\nSol: %d b %d t %d.      Distance/step: %f prev: %f    Force: %f\n", issolvent, blockIdx.x, threadIdx.x, delta_pos.len(), prev_vel, force->len());
#ifdef LIMA_VERBOSE
	if ((*pos-*pos_tsub1).len() > 0.1) {
		printf("\nP_index %d Thread %d blockId %d\tForce %f mass  %f \Dist %f\n", p_index, threadIdx.x, blockIdx.x, force->len(), mass, (*pos - *pos_tsub1).len());
		//printf("\nP_index %d Thread %d blockId %d\tForce %f mass  %f \tFrom %f %f %f\tTo %f %f %f\n", p_index, threadIdx.x, blockIdx.x, force->len(), mass, pos_tsub1->x, pos_tsub1->y, pos_tsub1->z, pos->x, pos->y, pos->z);		
	}
#endif
}

__device__ void integratePositionRampUp(Float3* pos, Float3* pos_tsub1, Float3* force, const double mass, const double dt, float rampup_scalar) {
	// Force is in ??Newton, [kg * nm /(mol*ns^2)] //

	//float prev_vel = (*pos - *pos_tsub1).len();

	Float3 temp = *pos;
	LIMAENG::applyHyperpos(pos, pos_tsub1);
	//printf("dt %.08f\n", dt);
	//*pos = *pos * 2 - *pos_tsub1 + *force * (1.f / mass) * dt * dt;	// no *0.5?	// [nm] - [nm] + [kg*m*s^-2]/[kg] * [ns]^2
	//double magnitude_equalizer = 1e+4;		// Due to force being very large and dt being very small
	//double masstime_scaling = ((dt * magnitude_equalizer) * (dt * magnitude_equalizer)) / mass;

	//Float3 vel_scaled = (*pos - *pos_tsub1) * (1. / rampup_scalar);
	//Float3 acc_scaled = (*force * (1. / (magnitude_equalizer * magnitude_equalizer))) * masstime_scaling * (1./rampup_scalar);

	//*pos = *pos + vel_scaled + acc_scaled;
	//*pos = *pos * 2 - *pos_tsub1 + *force * masstime_scaling;		// [nm] - [nm] + [kg/mol*m*/s^2]/[kg/mol] * [s]^2 * (1e-9)^2	=> [nm]-[nm]+[]
	//*pos = *pos * 2 - *pos_tsub1 + (*force * (1. / (magnitude_equalizer * magnitude_equalizer))) * masstime_scaling;		// [nm] - [nm] + [kg/mol*m*/s^2]/[kg/mol] * [s]^2 * (1e-9)^2	=> [nm]-[nm]+[]


	*pos = *pos * 2. - *pos_tsub1 + *force * (dt / mass) * dt * 0.f;		// [nm] - [nm] + [kg/mol*m*/s ^ 2] / [kg / mol] * [s] ^ 2 * (1e-9) ^ 2 = > [nm] - [nm] + []
	*pos_tsub1 = temp;

	return;

	float vel_scalar = log2f(force->len()/1000.);
	vel_scalar = max(vel_scalar, 1.f);
//	rampup_scalar = min(rampup_scalar, 2.f);
	Float3 delta_pos = *pos - *pos_tsub1;
	*pos = *pos_tsub1 + delta_pos * (1. / vel_scalar);
}




// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//




#define compound_index blockIdx.x
__global__ void compoundKernel(Box* box) {
	__shared__ Compound compound;
	__shared__ CompoundState compound_state;
	__shared__ NeighborList neighborlist;

	__shared__ BondedParticlesLUT bonded_particles_lut;

#ifdef ENABLE_SOLVENTS
	//__shared__ Float3 utility_buffer[NEIGHBORLIST_MAX_SOLVENTS];							// waaaaay too biggg
	//__shared__ uint8_t utility_buffer_small[NEIGHBORLIST_MAX_SOLVENTS];
	__shared__ Float3 utility_buffer[THREADS_PER_COMPOUNDBLOCK];
#else
	//__shared__ Float3 utility_buffer[MAX_COMPOUND_PARTICLES];
	//__shared__ uint8_t utility_buffer_small[MAX_COMPOUND_PARTICLES];
#endif



	if (threadIdx.x == 0) {
		compound.loadMeta(&box->compounds[blockIdx.x]);
		compound_state.setMeta(compound.n_particles);
		neighborlist.loadMeta(&box->compound_neighborlists[blockIdx.x]);

	}
	__syncthreads();
	compound.loadData(&box->compounds[blockIdx.x]);
	compound_state.loadData(&box->compound_state_array[blockIdx.x]);
	neighborlist.loadData(&box->compound_neighborlists[blockIdx.x]);
	__syncthreads();






	float potE_sum = 0;
	float data_ptr[4];
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999;
	data_ptr[3] = box->step + 1;
	

	Float3 force = compound.forces[threadIdx.x];
	
	force = Float3(0.f);
	Float3 force_LJ_sol(0.f);
	


	// ------------------------------------------------------------ Intramolecular Operations ------------------------------------------------------------ //
	{
		//BondedParticlesLUT* lut = box->bonded_particles_lut_manager->get(compound_index, compound_index);
		bonded_particles_lut.load(*box->bonded_particles_lut_manager->get(compound_index, compound_index));
		LIMAENG::applyHyperpos(&compound_state.positions[0], &compound_state.positions[threadIdx.x]);
		__syncthreads();
	/*	if (threadIdx.x == 0 && compound_index == 0) {

			bool match = true;
			for (int i = 0; i < MAX_COMPOUND_PARTICLES; i++) {
				for (int ii = 0; ii < MAX_COMPOUND_PARTICLES; ii++) {
					printf("%d ", *bonded_particles_lut.get(i, ii));
					if (*bonded_particles_lut.get(i, ii) != compound.bondedparticles_lookup[i][ii])
						printf("FAIL!");
				}
				printf("\n");
			}
			printf("\n\n\n\n");
		}*/
		__syncthreads();
		force += computePairbondForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeAnglebondForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeDihedralForces(&compound, compound_state.positions, utility_buffer, &potE_sum);


		force += computeIntramolecularLJForces(&compound, &compound_state, &potE_sum, data_ptr, bonded_particles_lut);
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //
	//force = Float3(0.f);

	// --------------------------------------------------------------- Intermolecular forces --------------------------------------------------------------- //
	
	//for (int neighborcompound_id = 0; neighborcompound_id < box->n_compounds; neighborcompound_id++) {
	for (int i = 0; i < neighborlist.n_compound_neighbors; i++) {
		int neighborcompound_id = neighborlist.neighborcompound_ids[i];
		if (neighborcompound_id == blockIdx.x) {
			printf("i: %d  Block %d neighbors nc %d of %d\n",i, blockIdx.x, neighborcompound_id, neighborlist.n_compound_neighbors);
		}
		//if (neighborcompound_id == blockIdx.x)	// Only needed untill we have proper neighbor lists
//			continue;
		int neighborcompound_particles = box->compound_state_array[neighborcompound_id].n_particles;

		if (threadIdx.x < neighborcompound_particles) {
			utility_buffer[threadIdx.x] = box->compound_state_array[neighborcompound_id].positions[threadIdx.x];
			//utility_buffer_small[threadIdx.x] = box->compounds[neighborcompound_id].atom_types[threadIdx.x];				/// HEEEEY WHY ARE WE NOT USING THIS?!??!?!?!?!
			LIMAENG::applyHyperpos(&compound_state.positions[0], &utility_buffer[threadIdx.x]);
		}
		__syncthreads();				// CRITICAL
		//continue;

		continue;									// DANGER
		if (threadIdx.x < compound.n_particles) {
			force += computerIntermolecularLJForces(&compound_state.positions[threadIdx.x], compound.atom_types[threadIdx.x], &compound.lj_ignore_list[threadIdx.x], &potE_sum, compound.particle_global_ids[threadIdx.x], data_ptr,
			//anti_inter += computerIntermolecularLJForces(&compound_state.positions[threadIdx.x], compound.atom_types[threadIdx.x], &compound.lj_ignore_list[threadIdx.x], &potE_sum, compound.particle_global_ids[threadIdx.x], data_ptr,
				&box->compounds[neighborcompound_id], utility_buffer, neighborcompound_id);
		}
		__syncthreads();				// CRITICAL			
	}
	//force += anti_inter;
	//force = Float3(0.f);
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //




	// --------------------------------------------------------------- Solvation forces --------------------------------------------------------------- //
#ifdef ENABLE_SOLVENTS
	/*for (int i = threadIdx.x; i < neighborlist.n_solvent_neighbors; i += blockDim.x) {
		utility_buffer[i] = box->solvents[neighborlist.neighborsolvent_ids[i]].pos;
		LIMAENG::applyHyperpos(&compound_state.positions[0], &utility_buffer[i]);
	}
	__syncthreads();
	if (threadIdx.x < compound.n_particles) {
		force_LJ_sol = computeSolventToCompoundLJForces(&compound_state.positions[threadIdx.x], neighborlist.n_solvent_neighbors, utility_buffer, data_ptr, &potE_sum, compound.atom_types[threadIdx.x]);
		force += force_LJ_sol;
	}
	*/
	for (int offset = 0; offset * blockDim.x < neighborlist.n_solvent_neighbors; offset += blockDim.x) {
		int solvent_nlist_index = offset + threadIdx.x; // index in neighborlist

		if (solvent_nlist_index < neighborlist.n_solvent_neighbors) {
			utility_buffer[threadIdx.x] = box->solvents[neighborlist.neighborsolvent_ids[solvent_nlist_index]].pos;
			LIMAENG::applyHyperpos(&compound_state.positions[0], &utility_buffer[threadIdx.x]);
		}		
		__syncthreads();

		if (threadIdx.x < compound.n_particles) {
			force_LJ_sol += computeSolventToCompoundLJForces(&compound_state.positions[threadIdx.x], blockDim.x, utility_buffer, data_ptr, &potE_sum, compound.atom_types[threadIdx.x]);			
		}
	}
	force += force_LJ_sol;
#endif
	// ------------------------------------------------------------------------------------------------------------------------------------------------ //


	//box->potE_buffer[threadIdx.x + blockIdx.x * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES + (box->step % STEPS_PER_TRAINDATATRANSFER) * gridDim.x * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES] = force.len();
	box->traj_buffer[threadIdx.x + blockIdx.x * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES + (box->step % STEPS_PER_TRAINDATATRANSFER) * gridDim.x * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES] = Float3(0.);// compound_state.positions[threadIdx.x];	// N DATAGAN???!?!


	// ------------------------------------------------------------ Integration ------------------------------------------------------------ //
	if (threadIdx.x < compound.n_particles) {
		//Float3 force_(force.x, force.y, force.z);
		//integratePosition(&compound_state.positions[threadIdx.x], &compound.prev_positions[threadIdx.x], &force_, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt, &box->thermostat_scalar, threadIdx.x, false);		
		if (box->step >= RAMPUP_STEPS || 1) {
			integratePosition(&compound_state.positions[threadIdx.x], &compound.prev_positions[threadIdx.x], &force, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt, &box->thermostat_scalar, threadIdx.x, false);
		}
		else {
			integratePositionRampUp(&compound_state.positions[threadIdx.x], &compound.prev_positions[threadIdx.x], &force, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt, RAMPUP_STEPS-box->step);
		}
		
		box->compounds[blockIdx.x].prev_positions[threadIdx.x] = compound.prev_positions[threadIdx.x];
	}
	__syncthreads();
	// ------------------------------------------------------------------------------------------------------------------------------------- //





	// ------------------------------------ PERIODIC BOUNDARY CONDITION ------------------------------------------------------------------------------------------------- // 
	if (threadIdx.x == 0) {
		applyPBC(&compound_state.positions[threadIdx.x]);
	}
	LIMAENG::applyHyperpos(&compound_state.positions[0], &compound_state.positions[threadIdx.x]);	// So all particles follows p0
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------ //
	





	

	
	// ------------------------------------ DATA LOG ------------------------------- //	
	{
		if (threadIdx.x < compound.n_particles || 1) {																							// TODO: Remove || 1
			int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound;
			int compound_offset = blockIdx.x * MAX_COMPOUND_PARTICLES;
			//box->potE_buffer[threadIdx.x + compound_offset + step_offset] = potE_sum;			// TODO: Should be += since bridge sets value first

		}		
		__syncthreads();


		
		if (blockIdx.x == 0 && threadIdx.x == 0) {			
			Float3 pos_prev_temp = compound.prev_positions[threadIdx.x];			
			LIMAENG::applyHyperpos(&compound_state.positions[threadIdx.x], &pos_prev_temp);
			
			int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * 10;

			box->outdata[0 + step_offset] = (compound_state.positions[threadIdx.x] - pos_prev_temp).len() / box->dt;
			box->outdata[1 + step_offset] = LIMAENG::calcKineticEnergy(&compound_state.positions[threadIdx.x], &pos_prev_temp, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt);
			box->outdata[2 + step_offset] = potE_sum;																											// This does NOT take bridge potE into account!!!!!!!
			box->outdata[3 + step_offset] = force.len();

			//box->outdata[5 + box->step * 10] = data_ptr[2];// closest particle
			//box->outdata[6 + box->step * 10] = data_ptr[1];// force.len();
		}


		int n_compounds_total = gridDim.x;
		int step_offset = (box->step % STEPS_PER_TRAINDATATRANSFER) * n_compounds_total * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES;
		int compound_offset = blockIdx.x * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES;
		int particle_offset = threadIdx.x * N_DATAGAN_VALUES;
		box->data_GAN[0 + particle_offset + compound_offset + step_offset] = compound_state.positions[threadIdx.x];
		box->data_GAN[1 + particle_offset + compound_offset + step_offset] = force_LJ_sol;
		box->data_GAN[2 + particle_offset + compound_offset + step_offset] = force;

		if (threadIdx.x >= compound.n_particles)
			box->data_GAN[0 + particle_offset + compound_offset + step_offset] = Float3(-1.f);
//		box->data_GAN[1 + threadIdx.x * 6 + step_offset] = force_bond + force_angle;
//		box->data_GAN[2 + threadIdx.x * 6 + step_offset] = force_LJ_com;
//		box->data_GAN[3 + threadIdx.x * 6 + step_offset] = force_LJ_sol;
		//box->outdata[0 + threadIdx.x * 10]
	}
	// ----------------------------------------------------------------------------- //

	if (force.len() > 2e+10) {
		printf("\n\nCritical force %.0f           block %d thread %d\n\n\n", force.len(), blockIdx.x, threadIdx.x);
		box->critical_error_encountered = true;
	}
		
	
	box->compound_state_array_next[blockIdx.x].loadData(&compound_state);
}
#undef compound_id










__global__ void solventForceKernel(Box* box) {
#define solvent_index (blockIdx.x * blockDim.x + threadIdx.x)
#define thread_active (solvent_index < box->n_solvents)
	//__shared__ Float3 solvent_positions[THREADS_PER_SOLVENTBLOCK];
	__shared__ Float3 utility_buffer[MAX_COMPOUND_PARTICLES];
	__shared__ uint8_t utility_buffer_small[MAX_COMPOUND_PARTICLES];
	__shared__ Solvent solvents[THREADS_PER_SOLVENTBLOCK];

	float potE_sum = 0;
	float data_ptr[4];	// Pot, force, closest particle, ?
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999.f;
	Float3 force(0,0,0);


	//Solvent* solvent = &solvents[threadIdx.x];
	Solvent solvent;
	//Float3 solvent_pos;
	if (thread_active) {
		solvent = box->solvents[solvent_index];
	}	

	// --------------------------------------------------------------- Molecule Interactions --------------------------------------------------------------- //
	for (int i = 0; i < box->n_compounds; i++) {
		//continue;		 // DANGER
		int n_compound_particles = box->compound_state_array[i].n_particles;
		// First all threads help loading the molecule
		if (threadIdx.x < n_compound_particles) {
			utility_buffer[threadIdx.x] = box->compound_state_array[i].positions[threadIdx.x];
			utility_buffer_small[threadIdx.x] = box->compounds[i].atom_types[threadIdx.x];
		}
		__syncthreads();
			

		if (thread_active) {
			LIMAENG::applyHyperpos(&utility_buffer[0], &solvent.pos);									// Move own particle in relation to compound-key-position
			force += computeCompoundToSolventLJForces(&solvent.pos, n_compound_particles, utility_buffer, data_ptr, &potE_sum, ATOMTYPE_SOL, utility_buffer_small);
		}
		__syncthreads();


	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //



	// --------------------------------------------------------------- Solvent Interactions --------------------------------------------------------------- //
	if (thread_active) {
		// DANGER
		force += computeSolventToSolventLJForces(&solvent.pos, &box->solvent_neighborlists[solvent_index], box->solvents, data_ptr, &potE_sum);
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //


	if (solvent_index < box->n_solvents)
		box->traj_buffer[box->n_compounds * MAX_COMPOUND_PARTICLES + solvent_index + (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound] = Float3(0.);
	//__syncthreads();
	// ------------------------------------ DATA LOG ------------------------------- //
	if (thread_active) {
		int compounds_offset = box->n_compounds * MAX_COMPOUND_PARTICLES;
		int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound;

		//box->potE_buffer[compounds_offset + solvent_index + step_offset] = (double) potE_sum;	//  data_ptr[0];
		//box->potE_buffer[compounds_offset + solvent_index + step_offset] = 0;
		//printf("pot: %f\n", box->potE_buffer[compounds_offset + solvent_index + (box->step) * box->total_particles]);
		box->potE_buffer[compounds_offset + solvent_index + step_offset] = force.len();
		box->traj_buffer[compounds_offset + solvent_index + step_offset] = solvent.pos;
		/*if ((solvent.pos.x > 7 || solvent.pos.y > 7 || solvent.pos.z > 7) && box->step == 0)
			solvent.pos.print('y');*/
		box->traj_buffer[compounds_offset + solvent_index + step_offset] = solvent.pos;
	}

	// ----------------------------------------------------------------------------- //

	if (thread_active) {
//		printf("%f\n", (solvent.pos - box->solvents[solvent_index].pos).len());

		int p_index = MAX_COMPOUND_PARTICLES + solvent_index;
		if (box->step >= RAMPUP_STEPS || 1) {
			integratePosition(&solvent.pos, &solvent.pos_tsub1, &force, forcefield_device.particle_parameters[ATOMTYPE_SOL].mass, box->dt, &box->thermostat_scalar, p_index, true);
		}
		else {
			integratePositionRampUp(&solvent.pos, &solvent.pos_tsub1, &force, forcefield_device.particle_parameters[ATOMTYPE_SOL].mass, box->dt, RAMPUP_STEPS - box->step);

		}
		float len = (solvent.pos - solvent.pos_tsub1).len();


		applyPBC(&solvent.pos);	


		
	}



	if (thread_active) {
		if (solvent.pos.x != solvent.pos.x) {
			solvent.pos.print('s');
			box->critical_error_encountered = true;
		}
		box->solvents_next[solvent_index] = solvent;
	}




#undef solvent_index
#undef thread_active
}






__global__ void compoundBridgeKernel(Box* box) {
//#define compound_id blockIdx.x
#define particle_id_bridge threadIdx.x
//#define particle_active particle_id_bridge < bridge.n_particles
	__shared__ CompoundBridgeCompact bridge;
	__shared__ Float3 positions[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Float3 utility_buffer[MAX_PARTICLES_IN_BRIDGE];							// waaaaay too biggg
	__shared__ uint8_t utility_buffer_small[MAX_PARTICLES_IN_BRIDGE];

	if (threadIdx.x == 0) {
		bridge.loadMeta(&box->bridge_bundle->compound_bridges[blockIdx.x]);
	}
	__syncthreads();
	bridge.loadData(&box->bridge_bundle->compound_bridges[blockIdx.x]);
	positions[particle_id_bridge] = Float3(0.f);
	if (particle_id_bridge < bridge.n_particles) {
		ParticleRefCompact* p_ref = &bridge.particle_refs[particle_id_bridge];
		positions[particle_id_bridge] = box->compound_state_array[p_ref->compound_id].positions[p_ref->local_id_compound];
	}
	__syncthreads();



	float potE_sum = 0;

	Float3 force(0.f);

	// ------------------------------------------------------------ Intercompund Operations ------------------------------------------------------------ //
	{											// So for the very first step, these ´should all be 0, but they are not??										TODO: Look into this at some point!!!! Also, can cant really apply hyperpos here without breaking stuff, sindce
																																								// The bridge spans such a big volume! :/
		LIMAENG::applyHyperpos(&positions[0], &positions[particle_id_bridge]);
		__syncthreads();	
		force += computePairbondForces(&bridge, positions, utility_buffer, &potE_sum);
		force += computeAnglebondForces(&bridge, positions, utility_buffer, &potE_sum);
		force += computeDihedralForces(&bridge, positions, utility_buffer, &potE_sum);
	}

	__syncthreads();
	// --------------------------------------------------------------------------------------------------------------------------------------------------- //

	if (particle_id_bridge < bridge.n_particles) {
		ParticleRefCompact* p_ref = &bridge.particle_refs[particle_id_bridge];
		box->compounds[p_ref->compound_id].forces[p_ref->local_id_compound] = force;



		int compound_offset = p_ref->compound_id * MAX_COMPOUND_PARTICLES;
		int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound;

		box->potE_buffer[p_ref->local_id_compound + compound_offset + step_offset] = potE_sum;
		//force.print(char((int)threadIdx.x+48));
	}


}



















/*	VELOCITY VERLET STORMER
__device__ void integratePosition(CompactParticle* particle, Float3* particle_pos, Float3* particle_force, double* dt) {
	*particle_pos = *particle_pos + (particle->vel + *particle_force * (0.5 / particle->mass) * *dt) * *dt;
}
__device__ void integrateVelocity(CompactParticle* particle, Float3* particle_force, double* dt) {
	particle->vel = particle->vel + (*particle_force + particle->force_prev) * (0.5 / particle->mass) * *dt;
}
*/