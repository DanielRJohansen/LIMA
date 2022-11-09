#include "Neighborlists.cuh"


NListManager::NListManager(Simulation* simulation) {
	nlist_data_collection = new NListDataCollection(simulation);


	for (int i = 0; i < nlist_data_collection->n_compounds; i++) {
		nlist_data_collection->compound_neighborlists[i].associated_id = i;
	}
	for (int i = 0; i < nlist_data_collection->n_solvents; i++) {
		nlist_data_collection->solvent_neighborlists[i].associated_id = i;
	}
}


void NListManager::updateNeighborLists(Simulation* simulation, bool* updatenlists_mutexlock, bool force_update, bool async, int* timings) {
	if (async && !force_update) {
		std::thread nlist_worker(NListUtils::updateNeighborLists, simulation, nlist_data_collection, &updated_neighborlists_ready, timings, updatenlists_mutexlock);
		nlist_worker.detach();
	}
	else {
		NListUtils::updateNeighborLists(simulation, nlist_data_collection, &updated_neighborlists_ready, timings, updatenlists_mutexlock);
	}

	prev_update_step = simulation->getStep();

	if (force_update) {
		Int3 n_data(nlist_data_collection->compound_neighborlists[0].n_compound_neighbors, nlist_data_collection->compound_neighborlists[0].n_solvent_neighbors, 0);
		printf("\nEntity neighbors: %d %d\n", n_data.x, n_data.y);
	}
}

void NListManager::offloadPositionDataNLIST(Simulation* simulation) {
	//cudaMemcpyAsync(compoundstatearray_host, simulation->box->compound_state_array, sizeof(CompoundState) * simulation->box->n_compounds, cudaMemcpyDeviceToHost);
	cudaMemcpy(nlist_data_collection->compoundstates, simulation->box->compound_state_array, sizeof(CompoundState) * simulation->n_compounds, cudaMemcpyDeviceToHost);
	if (simulation->n_solvents > 0)
		cudaMemcpy(nlist_data_collection->solvents, simulation->box->solvents, sizeof(Solvent) * simulation->n_solvents, cudaMemcpyDeviceToHost);
}

void NListManager::pushNlistsToDevice(Simulation* simulation) {
	cudaMemcpy(simulation->box->compound_neighborlists, nlist_data_collection->compound_neighborlists, sizeof(NeighborList) * simulation->n_compounds, cudaMemcpyHostToDevice);
	cudaMemcpy(simulation->box->solvent_neighborlists, nlist_data_collection->solvent_neighborlists, sizeof(NeighborList) * simulation->n_solvents, cudaMemcpyHostToDevice);
	updated_neighborlists_ready = 0;
}



namespace NListUtils {

	bool neighborWithinCutoff(Float3* pos_a, Float3* pos_b, float cutoff_offset) {		// This is used for compounds with a confining_particle_sphere from key_particle BEFORE CUTOFF begins
		Float3 pos_b_temp = *pos_b;
		EngineUtils::applyHyperpos(pos_a, &pos_b_temp);
		float dist = (*pos_a - pos_b_temp).len();
		return (dist < (CUTOFF + cutoff_offset));
	}


	void cullDistantNeighbors(Simulation* simulation, NListDataCollection* nlist_data_collection) {
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
				if (!neighborWithinCutoff(&nlist_data_collection->compound_key_positions[id_self], &nlist_data_collection->solvent_positions[id_neighbor], cutoff_add_self + CUTOFF) && false) {
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







	void updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection, volatile bool* finished, int* timing, bool* mutex_lock) {	// This is a thread worker-function, so it can't own the object, thus i pass a ref to the engine object..
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
		*timing = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();




		// SIGNALING MAIN THREAD //
		*finished = 1;		// Thread terminates here!
		*mutex_lock = 0;	// Unlock
	}
}
