#include "Neighborlists.cuh"

/*
* const int key_index = simulation.compounds_host[compound_id].key_particle_index;
		compound_key_positions[compound_id] = simulation.traj_buffer->getCompoundparticleDatapointAtIndex(compound_id, key_index, entryindex);
*/



// Assumes the compound is active
__device__ Float3 getCompoundAbspos(SimulationDevice& sim_dev, int compound_id)
{
	const int compound_key_index = sim_dev.box->compounds[compound_id].key_particle_index;

	const CompoundCoords& compound_coords = *CoordArrayQueueHelpers::getCoordarrayRef(sim_dev.box->coordarray_circular_queue, sim_dev.params->step, compound_id);

	const NodeIndex compound_origo = compound_coords.origo;

	const Float3 abspos = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compound_origo, compound_coords.rel_positions[compound_key_index]);

	return abspos;
}

const int threads_in_compoundnlist_kernel = 256;
__global__ void updateCompoundNlistsKernel(SimulationDevice* sim_dev) {

	const int n_compounds = sim_dev->box->boxparams.n_compounds;
	const int compound_id = blockIdx.x * blockDim.x + threadIdx.x;
	const bool compound_active = compound_id < n_compounds;

	NeighborList* nlist = compound_active
		? nlist = &sim_dev->box->compound_neighborlists[compound_id]
		: nullptr;

	const Float3 abspos_self = compound_active
		? getCompoundAbspos(*sim_dev, compound_id)
		: Float3{};

	const float cutoff_add_self = compound_active
		? sim_dev->box->compounds[compound_id].radius
		: 0.f;



	// Clear old nlist
	if (compound_active) {
		nlist->n_compound_neighbors = 0;
		//nlist->n_gridnodes = 0;// Do this in the other kernel
	}

	__shared__ Float3 abs_positions_buffer[threads_in_compoundnlist_kernel];

	// Loop over all compounds
	for (int offset = 0; offset < n_compounds; offset += blockDim.x) {
		
		// All threads help load a batch of compound_positions
		__syncthreads();
		if (threadIdx.x + offset < n_compounds)
			abs_positions_buffer[threadIdx.x] == getCompoundAbspos(*sim_dev, offset + compound_id);
		__syncthreads();

		// All active-compound threads now loop through the batch
		if (compound_active) {
			for (int i = 0; i < threads_in_compoundnlist_kernel; i++) {
				const int query_compound_id = offset + i;

				if (query_compound_id == n_compounds) { break; }

				const Float3& querycompound_pos = abs_positions_buffer[i];
				const float cutoff_add_query = sim_dev->box->compounds[query_compound_id].radius;	// Todo: move to shared mem, no reason for all threads to load the same value

				//if (neighborWithinCutoff(&pos_self, &pos_other, CUTOFF_NM + cutoff_add_self + cutoff_add_candidate)) {


			}
		}
	}




}




void NListManager::handleNlistGPU(Simulation* simulation, bool async, bool force_update, int* timing)
{
	auto t0 = std::chrono::high_resolution_clock::now();

	//NeighborList* nlists_dev;
	//cudaMalloc(&nlists_dev, sizeof(NeighborList) * simulation->boxparams_host.n_compounds);





	//// Make key positions addressable in arrays: compound_key_positions and solvent_positions
	//nlist_data_collection->preparePositionData(*simulation, step_at_update);

	//// Add all compound->compound neighbors
	//matchCompoundNeighbors(simulation, nlist_data_collection);

	//// updateCompoundGrid
	//nlist_data_collection->compoundgrid = std::make_unique<CompoundGrid>();	// Reset the grid. Maybe there is a way to do this faster?
	//assignNearbyCompoundsToGridnodes(simulation, nlist_data_collection);


	//pushNlists();

	//auto t1 = std::chrono::high_resolution_clock::now();
	//*timing = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

	//// SIGNALING MAIN THREAD //
	//*finished = 1;		// Thread terminates here!
}