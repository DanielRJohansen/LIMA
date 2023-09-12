
#include <algorithm>
//#include <execution>
#include <algorithm>
#include <unordered_set>
#include "Neighborlists.cuh"
#include "EngineUtils.cuh"
#include "Utilities.h"

// ------------------------------------------------------------------------------------------- PRIVATE HELPERS -------------------------------------------------------------------------------------------//


bool neighborWithinCutoff(const Float3* pos_a, const Float3* pos_b, const float cutoff_nm) {		// This is used for compounds with a confining_particle_sphere from key_particle BEFORE CUTOFF begins
	const float dist = EngineUtils::calcHyperDistNM(pos_a, pos_b);
	return dist < cutoff_nm;
}


NListDataCollection::NListDataCollection(Simulation* simulation) {
	compound_neighborlists.resize(MAX_COMPOUNDS);
	auto cuda_status = cudaMemcpy(compound_neighborlists.data(), simulation->sim_dev->box->compound_neighborlists, sizeof(NeighborList) * simulation->boxparams_host.n_compounds, cudaMemcpyDeviceToHost);
	LIMA_UTILS::genericErrorCheck(cuda_status);

	compoundgrid = std::make_unique<CompoundGrid>();
}

void NListDataCollection::preparePositionData(const Simulation& simulation, const uint32_t step_at_update) {
	// Data for the current step has not yet been generated so we need to use the previous step.
	// For the very first step, engine has cheated and already written the traj from the initial setup.	
	const auto step = step_at_update == 0 ? 0 : step_at_update - 1;	
	const auto entryindex = LIMALOGSYSTEM::getMostRecentDataentryIndex(step);

	for (int compound_id = 0; compound_id < simulation.boxparams_host.n_compounds; compound_id++) {
		const int key_index = simulation.compounds_host[compound_id].key_particle_index;
		compound_key_positions[compound_id] = simulation.traj_buffer->getCompoundparticleDatapointAtIndex(compound_id, key_index, entryindex);
	}
}



namespace NListUtils {

	void matchCompoundNeighbors(Simulation* simulation, NListDataCollection* nlist_data_collection) 
	{
		for (int i = 0; i < simulation->boxparams_host.n_compounds; i++) {
			nlist_data_collection->compound_neighborlists[i].n_compound_neighbors = 0;
		}

		for (uint16_t id_self = 0; id_self < simulation->boxparams_host.n_compounds; id_self++) {

			NeighborList* nlist_self = &nlist_data_collection->compound_neighborlists[id_self];
			const float cutoff_add_self = simulation->compounds_host[id_self].radius;
			const Float3& pos_self = nlist_data_collection->compound_key_positions[id_self];	// abs pos [nm]


			// Go through all compounds in box, with higher ID than self!
			for (uint16_t id_other = id_self + 1; id_other < simulation->boxparams_host.n_compounds; id_other++) {	// For finding new nearby compounds, it is faster and simpler to just check all compounds, since there are so few
				NeighborList* nlist_candidate = &nlist_data_collection->compound_neighborlists[id_other];
				const Float3& pos_other = nlist_data_collection->compound_key_positions[id_other];
				const float cutoff_add_candidate = simulation->compounds_host[id_other].radius;
				
				if (neighborWithinCutoff(&pos_self, &pos_other, CUTOFF_NM + cutoff_add_self + cutoff_add_candidate)) {
					nlist_self->addCompound(id_other);
					nlist_candidate->addCompound(id_self);					
				}
			}
		}
	}

	void assignNearbyCompoundsToGridnodes(Simulation* simulation, NListDataCollection* nlist_data_collection) {
#ifdef ENABLE_SOLVENTS
		for (int compound_id = 0; compound_id < simulation->boxparams_host.n_compounds; compound_id++) 
		{
			nlist_data_collection->compound_neighborlists[compound_id].n_gridnodes = 0;

			const Float3& compound_pos = nlist_data_collection->compound_key_positions[compound_id];
			const NodeIndex& compound_nodeindex = LIMAPOSITIONSYSTEM::absolutePositionToNodeIndex(compound_pos);

			const float compound_radius = simulation->compounds_host[compound_id].radius;
			const float max_dist_nm = CUTOFF_NM + compound_radius;

			for (int x = -GRIDNODE_QUERY_RANGE; x <= GRIDNODE_QUERY_RANGE; x++) {
				for (int y = -GRIDNODE_QUERY_RANGE; y <= GRIDNODE_QUERY_RANGE; y++) {
					for (int z = -GRIDNODE_QUERY_RANGE; z <= GRIDNODE_QUERY_RANGE; z++) {

						NodeIndex query_origo = compound_nodeindex + NodeIndex{ x,y,z };
						LIMAPOSITIONSYSTEM::applyPBC(query_origo);

						CompoundGridNode* querynode = nlist_data_collection->compoundgrid->getBlockPtr(query_origo);
						const int querynode_id = CompoundGrid::get1dIndex(query_origo);

						const Float3 querynode_pos = LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(query_origo);
						const float dist = EngineUtils::calcHyperDistNM(&compound_pos, &querynode_pos);

						if (dist < max_dist_nm) {
							querynode->addNearbyCompound(compound_id);	// Add compound so solvents can see it
							nlist_data_collection->compound_neighborlists[compound_id].addGridnode(querynode_id);	// Add grid so compound can see solvents
						}						
					}
				}
			}
		}
#endif
	}

	// Important: do NOT call getStep during this funciton, as it runs async!!!!
	// This is a thread worker-function, so it can't own the object, thus i pass a ref to the engine object..
	void updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection, volatile bool* finished, int* timing, std::mutex& mutex, const uint32_t step_at_update) {
		try {
			auto t0 = std::chrono::high_resolution_clock::now();
			mutex.lock();

			// Make key positions addressable in arrays: compound_key_positions and solvent_positions
			nlist_data_collection->preparePositionData(*simulation, step_at_update);

			// First do culling of neighbors that has left CUTOFF
			//NListUtils::cullDistantNeighbors(simulation, nlist_data_collection);

			// Add all compound->compound neighbors
			matchCompoundNeighbors(simulation, nlist_data_collection);

			// updateCompoundGrid
			nlist_data_collection->compoundgrid = std::make_unique<CompoundGrid>();	// Reset the grid. Maybe there is a way to do this faster?
			assignNearbyCompoundsToGridnodes(simulation, nlist_data_collection);

			auto t1 = std::chrono::high_resolution_clock::now();
			*timing = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

			// SIGNALING MAIN THREAD //
			*finished = 1;		// Thread terminates here!
			mutex.unlock();	// Unlock
		}
		catch (const std::runtime_error& ex) {
			std::cerr << "\nCaught exception: " << ex.what() << std::endl;	// TODO: Remove before final release
		}
	}
}








// ------------------------------------------------------------------------------------------- PUBLIC INTERFACE -------------------------------------------------------------------------------------------//

NListManager::NListManager(Simulation* simulation) {
	nlist_data_collection = std::make_unique<NListDataCollection>(simulation);

	for (int i = 0; i < simulation->boxparams_host.n_compounds; i++) {
		nlist_data_collection->compound_neighborlists[i].associated_id = i;
	}
}

// Main sim thread enters this block, so make sure it can leave VERY quickly
void NListManager::handleNLISTS(Simulation* simulation, const bool async, const bool force_update, int* timing) {
	const auto step = simulation->getStep();

	// Check whether we are getting too far behind
	if (stepsSinceUpdate(step) / STEPS_PER_NLIST_UPDATE > 4) {
#ifndef LIMAPUSH
		printf("\nWe are now %d nlist updates behind!\n", stepsSinceUpdate(step) / STEPS_PER_NLIST_UPDATE);
#endif
		//exit(1);
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
	cudaMemcpy(simulation->sim_dev->box->compound_neighborlists, nlist_data_collection->compound_neighborlists.data(), sizeof(NeighborList) * simulation->boxparams_host.n_compounds, cudaMemcpyHostToDevice);
	LIMA_UTILS::genericErrorCheck("Error after transferring compound neighborlists to device");

	//cudaMemcpy(simulation->box->solvent_neighborlists, nlist_data_collection->solvent_neighborlists, sizeof(NeighborList) * simulation->n_solvents, cudaMemcpyHostToDevice);

	cudaMemcpy(simulation->sim_dev->box->compound_grid, nlist_data_collection->compoundgrid.get(), sizeof(CompoundGrid), cudaMemcpyHostToDevice);
	LIMA_UTILS::genericErrorCheck("Error after transferring CompoundGrid to device");

	updated_neighborlists_ready = 0;
}
