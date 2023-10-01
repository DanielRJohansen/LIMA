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

// Returns false if an error occured
__device__ bool addAllNearbyCompounds(const SimulationDevice& sim_dev, NeighborList& nlist, const Float3* const abs_positions_buffer, const float* const cutoff_add_buffer,
	int offset, int n_compounds, float cutoff_add_self, const Float3& abspos_self, int compound_id) 
{
	for (int i = 0; i < threads_in_compoundnlist_kernel; i++) {
		const int query_compound_id = offset + i;

		if (query_compound_id == n_compounds) { break; }

		if (query_compound_id == compound_id) { continue; }	// dont add self to self

		const Float3& querycompound_pos = abs_positions_buffer[i];
		const float cutoff_add_query = cutoff_add_buffer[i];

		const float dist = EngineUtils::calcHyperDistNM(&abspos_self, &querycompound_pos);

		if (NListUtils::neighborWithinCutoff(&abspos_self, &querycompound_pos, CUTOFF_NM + cutoff_add_self + cutoff_add_query)) {
			if (!nlist.addCompound(static_cast<uint16_t>(query_compound_id)))
				return false;
		}
	}
	return true;
}


// This kernel creates a new nlist and pushes that to the kernel. Any other kernels that may
// interact with the neighborlist should come AFTER this kernel. Also, they should only run if this
// has run, and thus it is not allowed to comment out this kernel call.
__global__ void updateCompoundNlistsKernel(SimulationDevice* sim_dev) {

	const int n_compounds = sim_dev->box->boxparams.n_compounds;
	const int compound_id = blockIdx.x * blockDim.x + threadIdx.x;
	const bool compound_active = compound_id < n_compounds;
	
	NeighborList nlist;

	const Float3 abspos_self = compound_active
		? getCompoundAbspos(*sim_dev, compound_id)
		: Float3{};

	const float cutoff_add_self = compound_active
		? sim_dev->box->compounds[compound_id].radius
		: 0.f;



	//// Clear old nlist
	//if (compound_active) {
	//	nlist->n_compound_neighbors = 0;
	//	//nlist->n_gridnodes = 0;// Do this in the other kernel
	//}

	__shared__ Float3 abs_positions_buffer[threads_in_compoundnlist_kernel];
	__shared__ float cutoff_add_buffer[threads_in_compoundnlist_kernel];

	// Loop over all compounds and add all nearbys
	for (int offset = 0; offset < n_compounds; offset += blockDim.x) {
		// All threads help load a batch of compound_positions
		const int query_compound_id = threadIdx.x + offset;
		__syncthreads();
		if (query_compound_id < n_compounds) {
			abs_positions_buffer[threadIdx.x] = getCompoundAbspos(*sim_dev, query_compound_id);
			cutoff_add_buffer[threadIdx.x] = sim_dev->box->compounds[compound_id].radius;
		}
		__syncthreads();

		// All active-compound threads now loop through the batch
		if (compound_active) {
			const bool success = addAllNearbyCompounds(*sim_dev, nlist, abs_positions_buffer, cutoff_add_buffer, offset, n_compounds, cutoff_add_self, abspos_self, compound_id);

			if (!success) {
				sim_dev->params->critical_error_encountered = true;
			}
		}
	}

#ifdef ENABLE_SOLVENTS
	// Loop over the nearby gridnodes, and add them if they're within range
	if (compound_active) 
	{
		const CompoundCoords& compound_coords = *CoordArrayQueueHelpers::getCoordarrayRef(sim_dev->box->coordarray_circular_queue, sim_dev->params->step, compound_id);
		const NodeIndex compound_origo = compound_coords.origo;

		for (int x = -GRIDNODE_QUERY_RANGE; x <= GRIDNODE_QUERY_RANGE; x++) {
			for (int y = -GRIDNODE_QUERY_RANGE; y <= GRIDNODE_QUERY_RANGE; y++) {
				for (int z = -GRIDNODE_QUERY_RANGE; z <= GRIDNODE_QUERY_RANGE; z++) {

					NodeIndex query_origo = compound_origo + NodeIndex{ x,y,z };
					LIMAPOSITIONSYSTEM::applyPBC(query_origo);

					const int querynode_id = CompoundGrid::get1dIndex(query_origo);

					const Float3 querynode_pos = LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(query_origo);
					const float dist = EngineUtils::calcHyperDistNM(&abspos_self, &querynode_pos);

					if (dist < CUTOFF_NM + cutoff_add_self) {
						const bool success = nlist.addGridnode(querynode_id);
						// It is NOT handled here, that the gridnode must also have information about this compound;

						if (!success) {
							sim_dev->params->critical_error_encountered = true;
						}
					}
				}
			}
		}
	}
#endif

	// Push the new nlist
	if (compound_active) {
		sim_dev->box->compound_neighborlists[compound_id] = nlist;
	}
}

const int nthreads_in_blockgridkernel = 128;
__global__ void updateBlockgridKernel(SimulationDevice* sim_dev) 
{
	const int n_blocks = CompoundGrid::blocks_total;
	const int block_id = blockIdx.x * blockDim.x + threadIdx.x;
	const bool block_active = block_id < n_blocks;
	const int n_compounds = sim_dev->box->boxparams.n_compounds;

	CompoundGridNode gridnode;

	const NodeIndex block_origo = block_active
		? CompoundGrid::get3dIndex(block_id)
		: NodeIndex{};

	const Float3 block_abspos = LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(block_origo);

	__shared__ Float3 abs_positions_buffer[nthreads_in_blockgridkernel];
	__shared__ float cutoff_add_buffer[nthreads_in_blockgridkernel];

	// Loop over all compounds in batches
	for (int offset = 0; offset < n_compounds; offset += blockDim.x) {
		const int compound_id = offset + threadIdx.x;
		__syncthreads();
		if (compound_id < n_compounds) {
			abs_positions_buffer[threadIdx.x] = getCompoundAbspos(*sim_dev, compound_id);
			cutoff_add_buffer[threadIdx.x] = sim_dev->box->compounds[compound_id].radius;
		}
		__syncthreads();

		if (block_active) {
			for (int i = 0; i < nthreads_in_blockgridkernel; i++) {
				const int querycompound_id = i + offset;

				if (querycompound_id >= n_compounds) { break; }

				const float dist = EngineUtils::calcHyperDistNM(&block_abspos, &abs_positions_buffer[i]);

				if (dist < CUTOFF_NM + cutoff_add_buffer[i]) {
					gridnode.addNearbyCompound(querycompound_id);
					// It is NOT handled here, that the compound must also have information about this gridnode;
				}
			}
		}
	}
	
	if (block_active) {
		CompoundGridNode* gridnode_global = sim_dev->box->compound_grid->getBlockPtr(block_id);
		gridnode_global->loadData(gridnode);
	}
}



void NListManager::handleNlistGPU(Simulation* simulation, int& timing)
{
	const auto t0 = std::chrono::high_resolution_clock::now();

	{
		const int n_blocks = simulation->boxparams_host.n_compounds / threads_in_compoundnlist_kernel + 1;
		updateCompoundNlistsKernel << < n_blocks, threads_in_compoundnlist_kernel >> > (simulation->sim_dev);	// Must come before any other kernel()
	}

	cudaDeviceSynchronize();	// The above kernel overwrites the nlists, while the below fills ut the nlists present, so the above must be completed before progressing

#ifdef ENABLE_SOLVENTS
	{
		const int n_blocks = CompoundGrid::blocks_total / nthreads_in_blockgridkernel + 1;
		updateBlockgridKernel<<<n_blocks, nthreads_in_blockgridkernel>>>(simulation->sim_dev);
	}
#endif

	cudaDeviceSynchronize();
	const auto t1 = std::chrono::high_resolution_clock::now();
	timing += static_cast<int>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());
}