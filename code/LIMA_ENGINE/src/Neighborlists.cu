#include "Neighborlists.cuh"




bool neighborWithinCutoff(const Float3* pos_a, const Float3* pos_b, const float cutoff_lm) {		// This is used for compounds with a confining_particle_sphere from key_particle BEFORE CUTOFF begins
	const float dist = EngineUtils::calcHyperDist(pos_a, pos_b);
	return dist < cutoff_lm;
}





void inline addNeighborIfEligible(HashTable& currentNeighbors,
	NeighborList& nlist_self, NeighborList& nlist_other,
	const Float3& pos_self, const Float3& pos_other,
	const int& id_self, const int& id_other,
	NeighborList::NEIGHBOR_TYPE neighbortype_self, NeighborList::NEIGHBOR_TYPE neighbortype_other,
	const float cutoff_extension)
{
	if (neighborWithinCutoff(&pos_self, &pos_other, CUTOFF_LM + cutoff_extension)) {
		if (currentNeighbors.insert(id_other)) {
			nlist_self.addId(id_other, neighbortype_other);
			nlist_other.addId(id_self, neighbortype_self);
		}
	}
}

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
//	if (simulation->n_solvents > 0)
//		cudaMemcpy(nlist_data_collection->solvents, simulation->box->solvents, sizeof(Solvent) * simulation->n_solvents, cudaMemcpyDeviceToHost);
}

void NListManager::pushNlistsToDevice(Simulation* simulation) {
	cudaMemcpy(simulation->box->compound_neighborlists, nlist_data_collection->compound_neighborlists, sizeof(NeighborList) * simulation->n_compounds, cudaMemcpyHostToDevice);
	cudaMemcpy(simulation->box->solvent_neighborlists, nlist_data_collection->solvent_neighborlists, sizeof(NeighborList) * simulation->n_solvents, cudaMemcpyHostToDevice);
	updated_neighborlists_ready = 0;
}








namespace NListUtils {





	void cullDistantNeighbors(Simulation* simulation, NListDataCollection* nlist_data_collection) {
		for (int id_self = 0; id_self < nlist_data_collection->n_compounds; id_self++) {
			NeighborList* nlist_self = &nlist_data_collection->compound_neighborlists[id_self];
			float cutoff_add_self = simulation->compounds_host[id_self].confining_particle_sphere;



			for (int j = 0; j < nlist_self->n_compound_neighbors; j++) {		// Cull compound-compound
				int id_neighbor = nlist_self->neighborcompound_ids[j];
				NeighborList* nlist_neighbor = &nlist_data_collection->compound_neighborlists[id_neighbor];
				float cutoff_add_neighbor = simulation->compounds_host[id_neighbor].confining_particle_sphere;

				if (id_self < id_neighbor) {
					if (!neighborWithinCutoff(&nlist_data_collection->compound_key_positions[id_self], &nlist_data_collection->compound_key_positions[id_neighbor], cutoff_add_self + cutoff_add_neighbor + CUTOFF_LM)) {
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
				if (!neighborWithinCutoff(&nlist_data_collection->compound_key_positions[id_self], &nlist_data_collection->solvent_positions[id_neighbor], cutoff_add_self + CUTOFF_LM) && false) {
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

				if (!neighborWithinCutoff(&nlist_data_collection->solvent_positions[id_self], &nlist_data_collection->solvent_positions[id_neighbor], CUTOFF_LM)) {
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
		//Int3 before(nlist_data_collection->compound_neighborlists[0].n_compound_neighbors, nlist_data_collection->compound_neighborlists[0].n_solvent_neighbors, 0);

		// Make key positions addressable in arrays: compound_key_positions and solvent_positions
		//nlist_data_collection->preparePositionData(simulation->compounds_host);		
		nlist_data_collection->preparePositionData(*simulation);

		// First do culling of neighbors that has left CUTOFF
		NListUtils::cullDistantNeighbors(simulation, nlist_data_collection);


		// Now add compound->solvent, compound->compound
		for (uint16_t id_self = 0; id_self < simulation->n_compounds; id_self++) {

			NeighborList* nlist_self = &nlist_data_collection->compound_neighborlists[id_self];
			HashTable hashtable_compoundneighbors(nlist_self->neighborcompound_ids, nlist_self->n_compound_neighbors, NEIGHBORLIST_MAX_COMPOUNDS * 2);
			HashTable hashtable_solventneighbors(nlist_self->neighborsolvent_ids, nlist_self->n_solvent_neighbors, NEIGHBORLIST_MAX_SOLVENTS * 2);
			const float cutoff_add_self = simulation->compounds_host[id_self].confining_particle_sphere;
			const Float3& pos_self = nlist_data_collection->compound_key_positions[id_self];


			// Go through all solvents in box!
			for (uint16_t id_candidate = 0; id_candidate < simulation->n_solvents; id_candidate++) {
				NeighborList* nlist_candidate = &nlist_data_collection->solvent_neighborlists[id_candidate];
				const Float3& pos_other = nlist_data_collection->solvent_positions[id_candidate];

				addNeighborIfEligible(hashtable_solventneighbors, *nlist_self, *nlist_candidate,
					pos_self, pos_other,
					id_self, id_candidate,
					NeighborList::NEIGHBOR_TYPE::COMPOUND, NeighborList::NEIGHBOR_TYPE::SOLVENT,
					cutoff_add_self
				);
			}

			// Go through all compounds in box, with higher ID than self!
			for (uint16_t id_other = id_self + 1; id_other < simulation->n_compounds; id_other++) {	// For finding new nearby compounds, it is faster and simpler to just check all compounds, since there are so few
				NeighborList* nlist_candidate = &nlist_data_collection->compound_neighborlists[id_other];
				const Float3& pos_other = nlist_data_collection->compound_key_positions[id_other];
				const float cutoff_add_candidate = simulation->compounds_host[id_self].confining_particle_sphere;	// THIS IS BORKEN SINCE LIMAMETRES

				addNeighborIfEligible(hashtable_compoundneighbors, *nlist_self, *nlist_candidate,
					pos_self, pos_other,
					id_self, id_other,
					NeighborList::NEIGHBOR_TYPE::COMPOUND, NeighborList::NEIGHBOR_TYPE::COMPOUND,
					cutoff_add_self + cutoff_add_candidate
				);
			}
		}

		// Finally add all solvent->solvent candidates
		/*for (int id_self = 0; id_self < simulation->n_solvents; id_self++) {
			NeighborList* nlist_self = &nlist_data_collection->solvent_neighborlists[id_self];
			HashTable hashtable_solventneighbors(nlist_self->neighborsolvent_ids, (int)nlist_self->n_solvent_neighbors, NEIGHBORLIST_MAX_SOLVENTS * 2);
			const Float3& pos_self = nlist_data_collection->solvent_positions[id_self];


			for (int id_other = id_self + 1; id_other < simulation->n_solvents; id_other++) {
				NeighborList* nlist_candidate = &nlist_data_collection->solvent_neighborlists[id_other];
				const Float3 pos_other = nlist_data_collection->solvent_positions[id_other];

				addNeighborIfEligible(hashtable_solventneighbors, *nlist_self, *nlist_candidate,
					pos_self, pos_other,
					id_self, id_other,
					NeighborList::NEIGHBOR_TYPE::SOLVENT, NeighborList::NEIGHBOR_TYPE::SOLVENT,
					0.f
				);
			}
		}*/

		// New method:
		SolventBlockCollection solventblock_collection(nlist_data_collection->solvent_positions, simulation->n_solvents);
		const auto& neighborCandidatesAll = solventblock_collection.getNeighborSolventForAllSolvents(simulation->n_solvents);

		//for (int id_self = 0; id_self < simulation->n_solvents; id_self++) {
		//	NeighborList* nlist_self = &nlist_data_collection->solvent_neighborlists[id_self];
		//	HashTable hashtable_solventneighbors(nlist_self->neighborsolvent_ids, (int)nlist_self->n_solvent_neighbors, NEIGHBORLIST_MAX_SOLVENTS * 2);
		//	const Float3& pos_self = nlist_data_collection->solvent_positions[id_self];
		//	const auto& neighborCandidates = neighborCandidatesAll[id_self];

		//	for (int i = 0; i < neighborCandidates.n_candidates; i++) {
		//		const auto id_other = neighborCandidates.candidates[i];
		//		NeighborList* nlist_candidate = &nlist_data_collection->solvent_neighborlists[id_other];
		//		const Float3 pos_other = nlist_data_collection->solvent_positions[id_other];
		//		addNeighborIfEligible(hashtable_solventneighbors, *nlist_self, *nlist_candidate,
		//			pos_self, pos_other,
		//			id_self, id_other,
		//			NeighborList::NEIGHBOR_TYPE::SOLVENT, NeighborList::NEIGHBOR_TYPE::SOLVENT,
		//			0.f
		//		);
		//	}

		//}




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

SolventBlockCollection::SolventBlockCollection(const Float3* positions, int n_solvents) {
	for (auto i = 0; i < n_solvents; i++) {
		addSolventId(i, positions[i]);
	}
}

void SolventBlockCollection::addSolventId(uint32_t id, const Float3& pos) {
	const Int3 block_index = getSolventblockIndex(pos);

	m_blocks[block_index.x][block_index.y][block_index.z].insert(id);
}

std::vector<CandidateList> SolventBlockCollection::getNeighborSolventForAllSolvents(uint32_t n_solvents)
{
	std::vector<CandidateList> neighborCandidates(n_solvents);

	const auto& all_indices = getAllIndices();

	for (const auto& index : all_indices) {
		// First add all index combinations inside the block
		SolventBlock& block_self = getBlock(index);
		addAllInsideBlock(neighborCandidates, block_self);

		const auto& query_blockindices = getAdjacentIndicesThatAreGreater(index);
		for (const auto& query_index : query_blockindices) {
			SolventBlock& block_query = getBlock(query_index);

			// Now add all index combinations between blocks
			addAllBetweenBlocks(neighborCandidates, block_self, block_query);
		}
	}
	return neighborCandidates;
}







Int3 SolventBlockCollection::getSolventblockIndex(const Float3& pos) {
	return Int3(
		static_cast<int>(pos.x / block_len),
		static_cast<int>(pos.y / block_len),
		static_cast<int>(pos.z / block_len)
	);
}

constexpr std::array<Int3, SolventBlockCollection::blocks_per_dim> SolventBlockCollection::getAllIndices()
{
	std::array<Int3, blocks_per_dim> indices{};
	int index1d = 0;
	for (int z = 0; z < blocks_per_dim; z++) {
		for (int y = 0; y < blocks_per_dim; y++) {
			for (int x = 0; x < blocks_per_dim; x++) {
				indices[index1d++] = { Int3{x, y, z} };
			}
		}
	}
	return indices;
}

constexpr std::array<Int3, 2*2*2> SolventBlockCollection::getAdjacentIndicesThatAreGreater(Int3 index)
{
	std::array<Int3, 2*2*2> indices{};
	int index1d = 0;
	for (int z = index.z; z <= index.z + 1; z++) {
		for (int y = index.y; y <= index.y + 1; y++) {
			for (int x = index.x; x <= index.x + 1; x++) {
				indices[index1d++] = { Int3{x, y, z} };
			}
		}
	}
	return indices;
}

SolventBlockCollection::SolventBlock& SolventBlockCollection::getBlock(Int3 index)
{
	return m_blocks[index.x][index.y][index.z];
}

void SolventBlockCollection::addAllInsideBlock(std::vector<CandidateList>& neighborCandidates, const SolventBlock& block) {
	for (auto i = 0; i < block.n_elements; i++) {
		for (int ii = i; ii < block.n_elements; ii++) {
			neighborCandidates[i].candidates[neighborCandidates[i].n_candidates++] = ii;
			neighborCandidates[ii].candidates[neighborCandidates[ii].n_candidates++] = i;
		}
	}
}

void SolventBlockCollection::addAllBetweenBlocks(std::vector<CandidateList>& neighborCandidates, const SolventBlock& blocka, const SolventBlock& blockb) {
	for (auto ia = 0; ia < blocka.n_elements; ia++) {
		for (int ib = 0; ib < blockb.n_elements; ib++) {
			neighborCandidates[ia].candidates[neighborCandidates[ia].n_candidates++] = ib;
			neighborCandidates[ib].candidates[neighborCandidates[ib].n_candidates++] = ia;
		}
	}
}