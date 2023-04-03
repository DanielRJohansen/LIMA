#include "Neighborlists.cuh"
#include <algorithm>
#include <execution>
#include <algorithm>
#include <unordered_set>

bool neighborWithinCutoff(const Float3* pos_a, const Float3* pos_b, const float cutoff_lm) {		// This is used for compounds with a confining_particle_sphere from key_particle BEFORE CUTOFF begins
	const float dist = EngineUtils::calcHyperDistNM(pos_a, pos_b);
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


NListDataCollection::NListDataCollection(Simulation* simulation) {
	n_compounds = simulation->n_compounds;
	compoundstates = new CompoundState[n_compounds];
	compound_neighborlists = new NeighborList[MAX_COMPOUNDS];
	cudaMemcpy(compound_neighborlists, simulation->box->compound_neighborlists, sizeof(NeighborList) * n_compounds, cudaMemcpyDeviceToHost);

	compoundgrid = new CompoundGrid{};
	EngineUtils::genericErrorCheck("Error creating NListDataCollection");
}

// This doesn't currently work
void NListDataCollection::preparePositionData(const Simulation& simulation, const uint32_t step_at_update) {

	// Data for the current step has not yet been generated so we need to use the previous step.
	// For the very first step, engine has cheated and already written the traj from the initial setup.	
	const auto step = step_at_update == 0 ? 0 : step_at_update - 1;	

	for (int compound_id = 0; compound_id < n_compounds; compound_id++) {
		const size_t index = EngineUtils::getAlltimeIndexOfParticle(step, simulation.total_particles_upperbound, compound_id, 0);

		compound_key_positions[compound_id] = simulation.traj_buffer[index];
		compound_origos[compound_id] = LIMAPOSITIONSYSTEM::absolutePositionToNodeIndex(simulation.traj_buffer[index]);
		if (compound_origos[compound_id].x >= BOXGRID_N_NODES || compound_origos[compound_id].y >= BOXGRID_N_NODES || compound_origos[compound_id].z >= BOXGRID_N_NODES) {
			int a = 0;
		}
	}
}








namespace NListUtils {
	void cullDistantNeighbors(Simulation* simulation, NListDataCollection* nlist) {
		for (int id_self = 0; id_self < nlist->n_compounds; id_self++) {
			NeighborList* nlist_self = &nlist->compound_neighborlists[id_self];
			float cutoff_add_self = simulation->compounds_host[id_self].confining_particle_sphere;


			// Cull compound-compound
			for (int j = 0; j < nlist_self->n_compound_neighbors; j++) {
				int id_neighbor = nlist_self->neighborcompound_ids[j];
				NeighborList* nlist_neighbor = &nlist->compound_neighborlists[id_neighbor];
				float cutoff_add_neighbor = simulation->compounds_host[id_neighbor].confining_particle_sphere;

				if (id_self < id_neighbor) {
					if (!neighborWithinCutoff(&nlist->compound_key_positions[id_self], &nlist->compound_key_positions[id_neighbor], cutoff_add_self + cutoff_add_neighbor + CUTOFF_LM)) {
						nlist_self->removeId(id_neighbor, NeighborList::NEIGHBOR_TYPE::COMPOUND);
						nlist_neighbor->removeId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);
						j--;	// Decrement, as the removeId puts the last element at the current and now vacant spot.
					}
				}
			}

			nlist_self->n_gridnodes = 0;	// These are completely redone each iteration
		}
	}


	void matchCompoundNeighbors(Simulation* simulation, NListDataCollection* nlist_data_collection) {
		for (uint16_t id_self = 0; id_self < simulation->n_compounds; id_self++) {

			NeighborList* nlist_self = &nlist_data_collection->compound_neighborlists[id_self];
			HashTable hashtable_compoundneighbors(nlist_self->neighborcompound_ids, nlist_self->n_compound_neighbors, NEIGHBORLIST_MAX_COMPOUNDS * 2);
			HashTable hashtable_solventneighbors(nlist_self->neighborsolvent_ids, nlist_self->n_solvent_neighbors, NEIGHBORLIST_MAX_SOLVENTS * 2);
			const float cutoff_add_self = simulation->compounds_host[id_self].confining_particle_sphere;
			const Float3& pos_self = nlist_data_collection->compound_key_positions[id_self];	// abs pos [nm]

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
	}








	void distributeCompoundsInGrid(Simulation* simulation, NListDataCollection& nlist_data_collection) {
		for (int compound_index = 0; compound_index < simulation->n_compounds; compound_index++) {
			const NodeIndex& compound_nodeindex = nlist_data_collection.compound_origos[compound_index];
			CompoundGridNode& node = *nlist_data_collection.compoundgrid->getBlockPtr(compound_nodeindex);
			//const Coord& compound_origo = compoundgrid_host->getOrigosPtr()[compound_index];
			//CompoundGridNode& node = *compoundgrid_host->getBlockPtr(compound_origo);
			node.addAssociatedCompound(compound_index);
		}
	}



	//void NListManager::bootstrapCompoundgrid(Simulation* simulation) {
	//	nlist_data_collection->compoundgrid = new CompoundGrid{};
	//
	//	CompoundCoords* compoundcoords_array = new CompoundCoords[simulation->n_compounds];
	//	cudaMemcpy(compoundcoords_array, simulation->box->coordarray_circular_queue, sizeof(CompoundCoords) * simulation->n_compounds, cudaMemcpyDeviceToHost);
	//
	//	// We need to bootstrap origo's before we can use the normal functionality to find neighbors
	//	for (int compound_id = 0; compound_id < simulation->n_compounds; compound_id++) {
	//		nlist_data_collection->compoundgrid->getOrigosPtr()[compound_id] = compoundcoords_array[compound_id].origo;
	//	}
	//	delete[] compoundcoords_array;
	//	
	//	EngineUtils::genericErrorCheck("Here");
	//
	//	NListUtils::distributeCompoundsInGrid(simulation, nlist_data_collection->compoundgrid);
	//	NListUtils::assignNearbyCompoundsToGridnodes(simulation, nlist_data_collection);
	//	NListUtils::transferCompoundgridToDevice(simulation, nlist_data_collection->compoundgrid);
	//}



	bool isNearby(const Simulation& simulation, const NodeIndex& nodeindex_self, const int querycompound_id, NListDataCollection& nlist_data_collection) {
		//const Float3 querycompound_origo = (compoundgrid_host->getOrigosPtr()[querycompound_id]).toFloat3();
		//const Float3 querycompound_origo = nlist_data_collection.compound_origos[querycompound_id].toFloat3();
		//const float dist = (querycompound_origo - node_origo.toFloat3()).len();
		const NodeIndex nodeindex_querycompound = nlist_data_collection.compound_origos[querycompound_id];
		const float dist_nm = LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(nodeindex_querycompound - nodeindex_self).len();


		const float querycompound_radius = simulation.compounds_host[querycompound_id].confining_particle_sphere;	// Is this nm or lm?=?!??!!

		return dist_nm < (CUTOFF_NM + querycompound_radius);
	}

	void assignNearbyCompoundsToGridnodes(Simulation* simulation, NListDataCollection* nlist_data_collection) {
		//nlist_data_collection->compound_origos[0].print('C');
		for (int z = 0; z < BOXGRID_N_NODES; z++) {
			for (int y = 0; y < BOXGRID_N_NODES; y++) {
				for (int x = 0; x < BOXGRID_N_NODES; x++) {
					const NodeIndex node_origo{ x, y, z };	// Doubles as the 3D index of the block!
					CompoundGridNode* node_self = nlist_data_collection->compoundgrid->getBlockPtr(node_origo);
					int nodeself_id = CompoundGrid::get1dIndex(node_origo);


					const int query_range = 1;
					for (int x = -query_range; x <= query_range; x++) {
						for (int y = -query_range; y <= query_range; y++) {
							for (int z = -query_range; z <= query_range; z++) {
								NodeIndex query_origo = node_origo + NodeIndex{ x,y,z };
								LIMAPOSITIONSYSTEM::applyPBC(query_origo);
								const CompoundGridNode* node_query = nlist_data_collection->compoundgrid->getBlockPtr(query_origo);

								for (int i = 0; i < node_query->n_associated_compounds; i++) {
									const int querycompound_id = node_query->associated_ids[i];

									if (isNearby(*simulation, node_origo, querycompound_id, *nlist_data_collection) || true) {
										node_self->addNearbyCompound(querycompound_id);	// Add compound so solvents can see it
										nlist_data_collection->compound_neighborlists[querycompound_id].addGridnode(nodeself_id);	// Add grid so compound can see solvents

										if (node_origo == NodeIndex{ 3,3,3 }) {
											//printf("\nAdding  to %d\n", querycompound_id); 
											//nlist_data_collection->compound_origos[querycompound_id].print('C');
											//node_origo.print('N');
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}


	void updateCompoundGrid(Simulation* simulation, NListDataCollection* nlist) {
		*nlist->compoundgrid = CompoundGrid{};	// Reset the grid
		//printf("\n clear \n");
		// Load the origo's of each compound
		//cudaMemcpy(
		//	nlist->compoundgrid->getOrigosPtr(),
		//	simulation->box->compound_grid->getOrigosPtr(),
		//	sizeof(Coord) * MAX_COMPOUNDS,
		//	cudaMemcpyDeviceToHost
		//);

		distributeCompoundsInGrid(simulation, *nlist);
		assignNearbyCompoundsToGridnodes(simulation, nlist);
		//NListUtils::transferCompoundgridToDevice(simulation, compoundgrid_host);
	}


	void transferCompoundgridToDevice(Simulation* simulation, CompoundGrid* compoundgrid_host) {
		cudaMemcpy(simulation->box->compound_grid, compoundgrid_host, sizeof(CompoundGrid), cudaMemcpyHostToDevice);
		EngineUtils::genericErrorCheck("Error after transferring CompoundGrid to device");
	}


	// Important: do NOT call getStep during this funciton, as it runs async!!!!
	// This is a thread worker-function, so it can't own the object, thus i pass a ref to the engine object..
	void updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection, volatile bool* finished, int* timing, std::mutex& mutex, const uint32_t step_at_update) {
		auto t0 = std::chrono::high_resolution_clock::now();
		mutex.lock();

		// Make key positions addressable in arrays: compound_key_positions and solvent_positions
		nlist_data_collection->preparePositionData(*simulation, step_at_update);

		// First do culling of neighbors that has left CUTOFF
		NListUtils::cullDistantNeighbors(simulation, nlist_data_collection);
		EngineUtils::genericErrorCheck("2");

		// Add all compound->compound neighbors
		matchCompoundNeighbors(simulation, nlist_data_collection);

		updateCompoundGrid(simulation, nlist_data_collection);
		EngineUtils::genericErrorCheck("4");



		auto t1 = std::chrono::high_resolution_clock::now();
		*timing = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

		// SIGNALING MAIN THREAD //
		*finished = 1;		// Thread terminates here!
		mutex.unlock();	// Unlock
	}
}












NListManager::NListManager(Simulation* simulation) {
	nlist_data_collection = new NListDataCollection(simulation);

	for (int i = 0; i < nlist_data_collection->n_compounds; i++) {
		nlist_data_collection->compound_neighborlists[i].associated_id = i;
	}
}

// Main sim thread enters this block, so make sure it can leave VERY quickly
void NListManager::handleNLISTS(Simulation* simulation, const bool async, const bool force_update, int* timing) {
	const auto step = simulation->getStep();

	// Check whether we are getting too far behind
	if (stepsSinceUpdate(step) > STEPS_PER_NLIST_UPDATE * 3) {
		printf("We are now 3 nlist updates behind!");
		exit(1);
	}

	// If module is busy, return
	if (!m_mutex.try_lock()) { return; }

	// If new data is ready, push it
	if (updated_neighborlists_ready) {
		pushNlistsToDevice(simulation);
	}
	m_mutex.unlock();


	// If we dont need to update nlist, return
	if (!(stepsSinceUpdate(step) >= STEPS_PER_NLIST_UPDATE || step == 0)) { 
		return; 
	}

	if (async && !force_update) {
		std::thread nlist_worker(NListUtils::updateNeighborLists, simulation, nlist_data_collection, &updated_neighborlists_ready, timing, std::ref(m_mutex), step);
		nlist_worker.detach();
	}
	else {
		NListUtils::updateNeighborLists(simulation, nlist_data_collection, &updated_neighborlists_ready, timing, m_mutex, step);
	}
	prev_update_step = step;

	// If we are not async we can update immediately. If force_update, we need to wait. In either case lock
	if (!async || force_update) {
		const std::chrono::microseconds sleep_duration{ 500 };
		while (!updated_neighborlists_ready) { std::this_thread::sleep_for(sleep_duration); }
		m_mutex.lock();
		pushNlistsToDevice(simulation);
		m_mutex.unlock();
	}
}

//
//void NListManager::updateNeighborLists(Simulation* simulation, bool* updatenlists_mutexlock, bool force_update, bool async, int* timings, bool* critical_error) {
//	const uint32_t step_at_update = simulation->getStep();
//
//	if (async && !force_update) {
//		std::thread nlist_worker(NListUtils::updateNeighborLists, simulation, nlist_data_collection, &updated_neighborlists_ready, timings, updatenlists_mutexlock, step_at_update);
//		nlist_worker.detach();
//	}
//	else {
//		NListUtils::updateNeighborLists(simulation, nlist_data_collection, &updated_neighborlists_ready, timings, updatenlists_mutexlock, step_at_update);
//	}
//
//	prev_update_step = step_at_update;
//
//	if (force_update) {
//		Int3 n_data(nlist_data_collection->compound_neighborlists[0].n_compound_neighbors, nlist_data_collection->compound_neighborlists[0].n_solvent_neighbors, 0);
//		printf("\nEntity neighbors: %d %d\n", n_data.x, n_data.y);
//	}
//}

void NListManager::pushNlistsToDevice(Simulation* simulation) {
	cudaMemcpy(simulation->box->compound_neighborlists, nlist_data_collection->compound_neighborlists, sizeof(NeighborList) * simulation->n_compounds, cudaMemcpyHostToDevice);
	EngineUtils::genericErrorCheck("Error after transferring compound neighborlists to device");

	//cudaMemcpy(simulation->box->solvent_neighborlists, nlist_data_collection->solvent_neighborlists, sizeof(NeighborList) * simulation->n_solvents, cudaMemcpyHostToDevice);

	cudaMemcpy(simulation->box->compound_grid, nlist_data_collection->compoundgrid, sizeof(CompoundGrid), cudaMemcpyHostToDevice);
	EngineUtils::genericErrorCheck("Error after transferring CompoundGrid to device");

	updated_neighborlists_ready = 0;
}