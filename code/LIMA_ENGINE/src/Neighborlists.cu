#include "Neighborlists.cuh"
#include "BoundaryCondition.cuh"
#include "SimulationDevice.cuh"
#include "EngineUtils.cuh"

#include <chrono>


// Assumes the compound is active
template <typename BoundaryCondition>
__device__ void getCompoundAbspositions(SimulationDevice& sim_dev, int compound_id, Float3* result)
{
	const CompoundCoords& compound_coords = *CoordArrayQueueHelpers::getCoordarrayRef(sim_dev.box->coordarray_circular_queue, sim_dev.signals->step, compound_id);
	const NodeIndex compound_origo = compound_coords.origo;

	for (int i = 0; i < CompoundInteractionBoundary::k; i++) {
		const int particle_index = sim_dev.box->compounds[compound_id].interaction_boundary.key_particle_indices[i];
		const Float3 abspos = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compound_origo, compound_coords.rel_positions[particle_index]);
		result[i] = abspos;
	}
}
template __device__ void getCompoundAbspositions<PeriodicBoundaryCondition>(SimulationDevice& sim_dev, int compound_id, Float3* result);

template <typename BoundaryCondition>
__device__ bool canCompoundsInteract(const CompoundInteractionBoundary& left, const CompoundInteractionBoundary& right, const Float3* const posleft, const Float3* const posright) 
{
	for (int ileft = 0; ileft < CompoundInteractionBoundary::k; ileft++) {
		for (int iright = 0; iright < CompoundInteractionBoundary::k; iright++) {

			const float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM<BoundaryCondition>(&posleft[ileft], &posright[iright]);
			const float max_dist = CUTOFF_NM + left.radii[ileft] + right.radii[iright];
			if (dist < max_dist)
				return true;
		}
	}

	return false;
}

template <typename BoundaryCondition>
__device__ bool canCompoundInteractWithPoint(const CompoundInteractionBoundary& boundary, const Float3* const posleft, const Float3& point)
{
	for (int ileft = 0; ileft < CompoundInteractionBoundary::k; ileft++) {

		const float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM<BoundaryCondition>(&posleft[ileft], &point);
		const float max_dist = CUTOFF_NM + boundary.radii[ileft];

		if (dist < max_dist)
			return true;

	}

	return false;
}

// Returns false if an error occured
template <typename BoundaryCondition>
__device__ bool addAllNearbyCompounds(const SimulationDevice& sim_dev, NeighborList& nlist, const Float3* const key_positions_others /*[n_compounds, k]*/,
	const Float3* const key_positions_self, int offset, int n_compounds, int compound_id, const CompoundInteractionBoundary& boundary_self, 
	const CompoundInteractionBoundary* const boundaries_others,
	int n_bonded_compounds, const int* const bonded_compound_ids)
{
	// Now add all compounds nearby we are NOT bonded to. (They were added before this)
	for (int i = 0; i < blockDim.x; i++) {
		const int query_compound_id = offset + i;

		if (query_compound_id == n_compounds) { break; }

		if (query_compound_id == compound_id) { continue; }	// dont add self to self

		// Dont add bonded compounds to list again
		bool is_bonded_to_query = false;
		for (int i = 0; i < n_bonded_compounds; i++) {
			if (query_compound_id == bonded_compound_ids[i]) { 
				is_bonded_to_query = true;
				break;
			}
		}
		if (is_bonded_to_query) { continue; }

	/*	auto lut = sim_dev.box->bonded_particles_lut_manager->get(compound_id, query_compound_id);
		if (lut != nullptr || query_compound_id == nlist.neighborcompound_ids[0]) {
			printf("WTF is going on here %d %d\n", compound_id, query_compound_id);
		}*/

		const Float3* const positionsbegin_other = &key_positions_others[i * CompoundInteractionBoundary::k];
		if (canCompoundsInteract<BoundaryCondition>(boundary_self, boundaries_others[i], key_positions_self, positionsbegin_other))
		{
			if (!nlist.addCompound(static_cast<uint16_t>(query_compound_id)))
				return false;
		}
	}
	return true;
}
template __device__ bool addAllNearbyCompounds<PeriodicBoundaryCondition>(const SimulationDevice&, NeighborList&, const Float3* const, const Float3* const, int, int, int, const CompoundInteractionBoundary&,
	const CompoundInteractionBoundary* const, int, const int* const);



// This kernel creates a new nlist and pushes that to the kernel. Any other kernels that may
// interact with the neighborlist should come AFTER this kernel. Also, they should only run if this
// has run, and thus it is not allowed to comment out this kernel call.
const int threads_in_compoundnlist_kernel = 256;
template <typename BoundaryCondition>
__global__ void updateCompoundNlistsKernel(SimulationDevice* sim_dev) {

	const int n_compounds = sim_dev->box->boxparams.n_compounds;
	const int compound_id = blockIdx.x * blockDim.x + threadIdx.x;
	const bool compound_active = compound_id < n_compounds;
	
	NeighborList nlist;

	Float3 key_positions_self[CompoundInteractionBoundary::k];
	if (compound_active)
		getCompoundAbspositions<BoundaryCondition>(*sim_dev, compound_id, key_positions_self);

	const CompoundInteractionBoundary boundary_self = compound_active
		? sim_dev->box->compounds[compound_id].interaction_boundary
		: CompoundInteractionBoundary{};

	int bonded_compound_ids[Compound::max_bonded_compounds];
	const int n_bonded_compounds = compound_active
		? sim_dev->box->compounds[compound_id].n_bonded_compounds
		: 0;
	for (int i = 0; i < n_bonded_compounds; i++) {
		bonded_compound_ids[i] = sim_dev->box->compounds[compound_id].bonded_compound_ids[i];
	}
	// First add all the bonded compounds to the list
	{
		// First add the compounds that we are bonded to
		for (int i = 0; i < n_bonded_compounds; i++) {
			if (!nlist.addCompound(static_cast<uint16_t>(bonded_compound_ids[i]))) { sim_dev->signals->critical_error_encountered = true; }
		}
	}


	__shared__ Float3 key_positions_buffer[threads_in_compoundnlist_kernel * CompoundInteractionBoundary::k];
	__shared__ CompoundInteractionBoundary boundaries[threads_in_compoundnlist_kernel];

	// Loop over all compounds and add all nearbys
	for (int offset = 0; offset < n_compounds; offset += blockDim.x) {
		// All threads help load a batch of compound_positions
		const int query_compound_id = threadIdx.x + offset;
		__syncthreads();
		if (query_compound_id < n_compounds) {
			Float3* const positionsbegin = &key_positions_buffer[threadIdx.x * CompoundInteractionBoundary::k];
			getCompoundAbspositions<BoundaryCondition>(*sim_dev, query_compound_id, positionsbegin);
			boundaries[threadIdx.x] = sim_dev->box->compounds[query_compound_id].interaction_boundary;
		}
		__syncthreads();

		// All active-compound threads now loop through the batch
		if (compound_active) {
			const bool success = addAllNearbyCompounds<BoundaryCondition>(*sim_dev, nlist, key_positions_buffer, key_positions_self, offset, n_compounds,
				compound_id, boundary_self, boundaries, n_bonded_compounds, bonded_compound_ids);
			if (!success) {
				sim_dev->signals->critical_error_encountered = true;
			}
		}
	}
#ifdef ENABLE_SOLVENTS
	// Loop over the nearby gridnodes, and add them if they're within range
	if (compound_active) 
	{
		const CompoundCoords& compound_coords = *CoordArrayQueueHelpers::getCoordarrayRef(sim_dev->box->coordarray_circular_queue, sim_dev->signals->step, compound_id);
		const NodeIndex compound_origo = compound_coords.origo;

		for (int x = -GRIDNODE_QUERY_RANGE; x <= GRIDNODE_QUERY_RANGE; x++) {
			for (int y = -GRIDNODE_QUERY_RANGE; y <= GRIDNODE_QUERY_RANGE; y++) {
				for (int z = -GRIDNODE_QUERY_RANGE; z <= GRIDNODE_QUERY_RANGE; z++) {
					NodeIndex query_origo = compound_origo + NodeIndex{ x,y,z };
					LIMAPOSITIONSYSTEM::applyBC<BoundaryCondition>(query_origo);

					// If the query node is NOT inside the box, which happens in some boundary conditions, we cannot continue, 
					// since the node wont exists, and thus compounds are not allowed to access it.
					if (!query_origo.isInBox(BOXGRID_N_NODES))
						continue;

					const int querynode_id = CompoundGrid::get1dIndex(query_origo);
					const Float3 querynode_pos = LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(query_origo);


					if (canCompoundInteractWithPoint<BoundaryCondition>(boundary_self, key_positions_self, querynode_pos)) {
						const bool success = nlist.addGridnode(querynode_id);
						if (!success) {
							sim_dev->signals->critical_error_encountered = true;
						}
					}
				}
			}
		}
	}
#endif

	// Push the new nlist
	if (compound_active) {
		//sim_dev->box->compound_neighborlists[compound_id] = nlist;
		sim_dev->compound_neighborlists[compound_id] = nlist;
	}
}
template __global__ void updateCompoundNlistsKernel<PeriodicBoundaryCondition>(SimulationDevice* sim_dev);
template __global__ void updateCompoundNlistsKernel<NoBoundaryCondition>(SimulationDevice* sim_dev);

const int nthreads_in_blockgridkernel = 128;
template <typename BoundaryCondition>
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

	__shared__ Float3 key_positions_buffer[nthreads_in_blockgridkernel * CompoundInteractionBoundary::k];
	__shared__ CompoundInteractionBoundary boundaries[nthreads_in_blockgridkernel];

	// Loop over all compounds in batches
	for (int offset = 0; offset < n_compounds; offset += blockDim.x) {
		const int compound_id = offset + threadIdx.x;
		__syncthreads();
		if (compound_id < n_compounds) {

			Float3* const positionsbegin = &key_positions_buffer[threadIdx.x * CompoundInteractionBoundary::k];
			getCompoundAbspositions<BoundaryCondition>(*sim_dev, compound_id, positionsbegin);
			boundaries[threadIdx.x] = sim_dev->box->compounds[compound_id].interaction_boundary;
		}
		__syncthreads();

		if (block_active) {
			for (int i = 0; i < blockDim.x; i++) {
				const int querycompound_id = i + offset;

				if (querycompound_id >= n_compounds) { break; }

				Float3* const positionsbegin = &key_positions_buffer[i * CompoundInteractionBoundary::k];
				if (canCompoundInteractWithPoint<BoundaryCondition>(boundaries[i], positionsbegin, block_abspos)) {
					if (!gridnode.addNearbyCompound(querycompound_id)) {
						sim_dev->signals->critical_error_encountered = true;
					}
				}
			}
		}
	}
	if (block_active) {
		CompoundGridNode* gridnode_global = sim_dev->compound_grid->getBlockPtr(block_id);
		gridnode_global->loadData(gridnode);
	}
}
template __global__ void updateBlockgridKernel<PeriodicBoundaryCondition>(SimulationDevice* sim_dev);
template __global__ void updateBlockgridKernel<NoBoundaryCondition>(SimulationDevice* sim_dev);






void NeighborLists::updateNlists(SimulationDevice* sim_dev, BoundaryConditionSelect bc_select, const BoxParams& boxparams, int& timing)
{
	const auto t0 = std::chrono::high_resolution_clock::now();

	// Technically we could only run if > 1, buuut running with any compounds lets us spot bugs easier.
	if (boxparams.n_compounds > 0) {
		const int n_blocks = boxparams.n_compounds / threads_in_compoundnlist_kernel + 1;
		LAUNCH_GENERIC_KERNEL(updateCompoundNlistsKernel, n_blocks, threads_in_compoundnlist_kernel, bc_select, sim_dev);
		//updateCompoundNlistsKernel<BoundaryCondition> << < n_blocks, threads_in_compoundnlist_kernel >> > (sim_dev);	// Must come before any other kernel()
	}

	cudaDeviceSynchronize();	// The above kernel overwrites the nlists, while the below fills ut the nlists present, so the above must be completed before progressing

//#ifdef ENABLE_SOLVENTS
	if (boxparams.n_solvents > 0){
		const int n_blocks = CompoundGrid::blocks_total / nthreads_in_blockgridkernel + 1;
		LAUNCH_GENERIC_KERNEL(updateBlockgridKernel, n_blocks, nthreads_in_blockgridkernel, bc_select, sim_dev);
		//updateBlockgridKernel<BoundaryCondition> <<<n_blocks, nthreads_in_blockgridkernel>>>(sim_dev);
	}
//#endif

	cudaDeviceSynchronize();
	const auto t1 = std::chrono::high_resolution_clock::now();
	timing += static_cast<int>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());
}