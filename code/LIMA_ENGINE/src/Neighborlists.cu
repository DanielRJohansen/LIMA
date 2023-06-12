#include "LIMA_ENGINE/include/Neighborlists.cuh"
#include <algorithm>
#include <execution>
#include <algorithm>
#include <unordered_set>
#include "LIMA_BASE/include/Utilities.h"

// ------------------------------------------------------------------------------------------- PRIVATE HELPERS -------------------------------------------------------------------------------------------//


bool neighborWithinCutoff(const Float3* pos_a, const Float3* pos_b, const float cutoff_nm) {		// This is used for compounds with a confining_particle_sphere from key_particle BEFORE CUTOFF begins
	const float dist = EngineUtils::calcHyperDistNM(pos_a, pos_b);
	return dist < cutoff_nm;
}

void inline addNeighborIfEligible(HashTable& currentNeighbors,
	NeighborList& nlist_self, NeighborList& nlist_other,
	const Float3& pos_self, const Float3& pos_other,
	const int& id_self, const int& id_other,
	NeighborList::NEIGHBOR_TYPE neighbortype_self, NeighborList::NEIGHBOR_TYPE neighbortype_other,
	const float cutoff_extension)
{
	if (neighborWithinCutoff(&pos_self, &pos_other, CUTOFF_NM + cutoff_extension)) {
		if (currentNeighbors.insert(id_other)) {
			nlist_self.addId(id_other, neighbortype_other);
			nlist_other.addId(id_self, neighbortype_self);
		}
	}
}


NListDataCollection::NListDataCollection(Simulation* simulation) {
	compound_neighborlists.resize(MAX_COMPOUNDS);
	cudaMemcpy(compound_neighborlists.data(), simulation->sim_dev->box->compound_neighborlists, sizeof(NeighborList) * simulation->n_compounds, cudaMemcpyDeviceToHost);

	compoundgrid = std::make_unique<CompoundGrid>();
	LIMA_UTILS::genericErrorCheck("Error creating NListDataCollection");
}

void NListDataCollection::preparePositionData(const Simulation& simulation, const uint32_t step_at_update) {

	// Data for the current step has not yet been generated so we need to use the previous step.
	// For the very first step, engine has cheated and already written the traj from the initial setup.	
	const auto step = step_at_update == 0 ? 0 : step_at_update - 1;	

	for (int compound_id = 0; compound_id < simulation.n_compounds; compound_id++) {
		const size_t index = EngineUtils::getAlltimeIndexOfParticle(step, simulation.total_particles_upperbound, compound_id, 0);

		compound_key_positions[compound_id] = simulation.traj_buffer[index]; // Temp?
		// const LimaPosition position = LIMAPOSITIONSYSTEM::createLimaPosition(simulation.traj_buffer[index]);
		//compound_origos[compound_id] = LIMAPOSITIONSYSTEM::absolutePositionToNodeIndex(position);
		//if (compound_origos[compound_id].x >= BOXGRID_N_NODES || compound_origos[compound_id].y >= BOXGRID_N_NODES || compound_origos[compound_id].z >= BOXGRID_N_NODES) {
		//	int a = 0;
		//}
	}
}



namespace NListUtils {
	void cullDistantNeighbors(Simulation* simulation, NListDataCollection* nlist) {
		for (int id_self = 0; id_self < simulation->n_compounds; id_self++) {
			NeighborList* nlist_self = &nlist->compound_neighborlists[id_self];
			float cutoff_add_self = simulation->compounds_host[id_self].radius;


			// Cull compound-compound
			for (int j = 0; j < nlist_self->n_compound_neighbors; j++) {
				int id_neighbor = nlist_self->neighborcompound_ids[j];
				NeighborList* nlist_neighbor = &nlist->compound_neighborlists[id_neighbor];
				float cutoff_add_neighbor = simulation->compounds_host[id_neighbor].radius;

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
			const float cutoff_add_self = simulation->compounds_host[id_self].radius;
			const Float3& pos_self = nlist_data_collection->compound_key_positions[id_self];	// abs pos [nm]

			// Go through all compounds in box, with higher ID than self!
			for (uint16_t id_other = id_self + 1; id_other < simulation->n_compounds; id_other++) {	// For finding new nearby compounds, it is faster and simpler to just check all compounds, since there are so few
				NeighborList* nlist_candidate = &nlist_data_collection->compound_neighborlists[id_other];
				const Float3& pos_other = nlist_data_collection->compound_key_positions[id_other];
				const float cutoff_add_candidate = simulation->compounds_host[id_other].radius;	// THIS IS BORKEN SINCE LIMAMETRES

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
			const Float3& compound_position = nlist_data_collection.compound_key_positions[compound_index];
			const NodeIndex& nearest_nodeindex = LIMAPOSITIONSYSTEM::absolutePositionToNodeIndex(compound_position);

			CompoundGridNode& node = *nlist_data_collection.compoundgrid->getBlockPtr(nearest_nodeindex);

			node.addAssociatedCompound(compound_index);
		}
	}


	bool isNearby(const Simulation& simulation, const NodeIndex& nodeindex_self, const int querycompound_id, NListDataCollection& nlist_data) {
		const Float3& querycompound_pos = nlist_data.compound_key_positions[querycompound_id];

		const Float3 currentnode_pos = LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(nodeindex_self);

		const float dist = EngineUtils::calcHyperDistNM(&querycompound_pos, &currentnode_pos);

		const float querycompound_radius = simulation.compounds_host[querycompound_id].radius;	// Is this nm or lm?=?!??!!

		return dist < (CUTOFF_NM + querycompound_radius);
	}

	void assignNearbyCompoundsToGridnodes(Simulation* simulation, NListDataCollection* nlist_data_collection) {
		for (int z = 0; z < BOXGRID_N_NODES; z++) {
			for (int y = 0; y < BOXGRID_N_NODES; y++) {
				for (int x = 0; x < BOXGRID_N_NODES; x++) {
					const NodeIndex node_origo{ x, y, z };	// Doubles as the 3D index of the block!
					CompoundGridNode* node_self = nlist_data_collection->compoundgrid->getBlockPtr(node_origo);
					int nodeself_id = CompoundGrid::get1dIndex(node_origo);


					const int query_range = 2;
					for (int x = -query_range; x <= query_range; x++) {
						for (int y = -query_range; y <= query_range; y++) {
							for (int z = -query_range; z <= query_range; z++) {
								NodeIndex query_origo = node_origo + NodeIndex{ x,y,z };
								LIMAPOSITIONSYSTEM::applyPBC(query_origo);
								const CompoundGridNode* node_query = nlist_data_collection->compoundgrid->getBlockPtr(query_origo);

								for (int i = 0; i < node_query->n_associated_compounds; i++) {
									const int querycompound_id = node_query->associated_ids[i];

									if (isNearby(*simulation, node_origo, querycompound_id, *nlist_data_collection)) {
										node_self->addNearbyCompound(querycompound_id);	// Add compound so solvents can see it
										nlist_data_collection->compound_neighborlists[querycompound_id].addGridnode(nodeself_id);	// Add grid so compound can see solvents
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
		nlist->compoundgrid = std::make_unique<CompoundGrid>();	// Reset the grid. Maybe there is a way to do this faster?

		distributeCompoundsInGrid(simulation, *nlist);
		assignNearbyCompoundsToGridnodes(simulation, nlist);
	}


	// Important: do NOT call getStep during this funciton, as it runs async!!!!
	// This is a thread worker-function, so it can't own the object, thus i pass a ref to the engine object..
	void updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection, volatile bool* finished, int* timing, std::mutex& mutex, const uint32_t step_at_update) {
		auto t0 = std::chrono::high_resolution_clock::now();
		mutex.lock();

		LIMA_UTILS::genericErrorCheck("24");

		// Make key positions addressable in arrays: compound_key_positions and solvent_positions
		nlist_data_collection->preparePositionData(*simulation, step_at_update);
		LIMA_UTILS::genericErrorCheck("24");
		// First do culling of neighbors that has left CUTOFF
		NListUtils::cullDistantNeighbors(simulation, nlist_data_collection);
		LIMA_UTILS::genericErrorCheck("2");

		// Add all compound->compound neighbors
		matchCompoundNeighbors(simulation, nlist_data_collection);

		updateCompoundGrid(simulation, nlist_data_collection);
		LIMA_UTILS::genericErrorCheck("4");



		auto t1 = std::chrono::high_resolution_clock::now();
		*timing = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

		// SIGNALING MAIN THREAD //
		*finished = 1;		// Thread terminates here!
		mutex.unlock();	// Unlock
	}
}








// ------------------------------------------------------------------------------------------- PUBLIC INTERFACE -------------------------------------------------------------------------------------------//

NListManager::NListManager(Simulation* simulation) {
	nlist_data_collection = std::make_unique<NListDataCollection>(simulation);

	for (int i = 0; i < simulation->n_compounds; i++) {
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
		std::thread nlist_worker(NListUtils::updateNeighborLists, simulation, nlist_data_collection.get(), &updated_neighborlists_ready, timing, std::ref(m_mutex), step);
		nlist_worker.detach();
	}
	else {
		NListUtils::updateNeighborLists(simulation, nlist_data_collection.get(), &updated_neighborlists_ready, timing, m_mutex, step);
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


void NListManager::pushNlistsToDevice(Simulation* simulation) {
	cudaMemcpy(simulation->sim_dev->box->compound_neighborlists, nlist_data_collection->compound_neighborlists.data(), sizeof(NeighborList) * simulation->n_compounds, cudaMemcpyHostToDevice);
	LIMA_UTILS::genericErrorCheck("Error after transferring compound neighborlists to device");

	//cudaMemcpy(simulation->box->solvent_neighborlists, nlist_data_collection->solvent_neighborlists, sizeof(NeighborList) * simulation->n_solvents, cudaMemcpyHostToDevice);

	cudaMemcpy(simulation->sim_dev->box->compound_grid, nlist_data_collection->compoundgrid.get(), sizeof(CompoundGrid), cudaMemcpyHostToDevice);
	LIMA_UTILS::genericErrorCheck("Error after transferring CompoundGrid to device");

	updated_neighborlists_ready = 0;
}