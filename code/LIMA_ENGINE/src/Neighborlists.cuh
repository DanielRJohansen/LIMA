#pragma once

#include "BoundaryCondition.cuh"
#include "LimaPositionSystem.cuh"
#include "Simulation.cuh"
#include "SimulationDevice.cuh"



namespace NeighborLists {
    void updateNlists(SimulationDevice*, BoundaryConditionSelect, const BoxParams&, cudaStream_t& s1, cudaStream_t& s2, cudaStream_t& s3);
};



























// Assumes the compound is active
template <typename BoundaryCondition>
__device__ void getCompoundAbspositions(const SimulationDevice& sim_dev, int compound_id, Float3* result, const NodeIndex& compoundOrigo, const CompoundInteractionBoundary& compoundInteractionBoundary)
{
	const Float3* const relPositions = &sim_dev.boxState->compoundsRelposNm[compound_id * MAX_COMPOUND_PARTICLES];

	for (int i = 0; i < CompoundInteractionBoundary::k; i++) {
		const int particle_index = compoundInteractionBoundary.key_particle_indices[i];
		const Float3 abspos = compoundOrigo.toFloat3() + relPositions[particle_index];
		result[i] = abspos;
	}
}

template <typename BoundaryCondition>
__device__ bool canCompoundsInteract(const CompoundInteractionBoundary& left, const CompoundInteractionBoundary& right, const Float3* const positionsLeft, const Float3* const positionsRight)
{
	for (int ileft = 0; ileft < CompoundInteractionBoundary::k; ileft++) {
		for (int iright = 0; iright < CompoundInteractionBoundary::k; iright++) {

            const float distSquared = LIMAPOSITIONSYSTEM::calcHyperDistSquaredNM<BoundaryCondition>(positionsLeft[ileft], positionsRight[iright]);
			const float max_dist = DeviceConstants::cutoffNM + left.radii[ileft] + right.radii[iright];

            if (distSquared < max_dist*max_dist)
				return true;
		}
	}

	return false;
}

template <typename BoundaryCondition>
__device__ bool canCompoundInteractWithPoint(const CompoundInteractionBoundary& boundary, const Float3* const positionsLeft, const Float3& point)
{
	for (int ileft = 0; ileft < CompoundInteractionBoundary::k; ileft++) {

        const float distSq = LIMAPOSITIONSYSTEM::calcHyperDistSquaredNM<BoundaryCondition>(positionsLeft[ileft], point);
		const float max_dist = DeviceConstants::cutoffNM + boundary.radii[ileft];

        if (distSq < max_dist*max_dist)
			return true;

	}

	return false;
}

// Returns false if an error occured
template <typename BoundaryCondition>
__device__ void addAllNearbyCompounds(const SimulationDevice& sim_dev, const Float3* const key_positions_others /*[n_compounds, k]*/,
	const Float3* const key_positions_self, int offset, int n_compounds, int compound_id, const CompoundInteractionBoundary& boundary_self,
	const CompoundInteractionBoundary* const boundaries_others,
	int n_bonded_compounds, const int* const bonded_compound_ids,
    const NodeIndex& myCompoundOrigo, const NodeIndex* const compoundOrigos,
    uint16_t& nNonbondedNeighbors, NlistUtil::IdAndRelshift* const nonbondedNeighbors, const uint8_t* const compoundsNParticles
	)
{
	// Now add all compounds nearby we are NOT bonded to. (They were added before this)
	for (int i = 0; i < blockDim.x; i++) {
		const int query_compound_id = offset + i;

		if (query_compound_id == n_compounds) { break; }

		if (query_compound_id == compound_id) { continue; }	// dont add self to self

		// Dont add bonded compounds to list again
		bool is_bonded_to_query = false;
		for (int j = 0; j < n_bonded_compounds; j++) {
			if (query_compound_id == bonded_compound_ids[j]) {
				is_bonded_to_query = true;
				break;
			}
		}
		if (is_bonded_to_query) { continue; }

		const Float3* const positionsbegin_other = &key_positions_others[i * CompoundInteractionBoundary::k];
		if (canCompoundsInteract<BoundaryCondition>(boundary_self, boundaries_others[i], key_positions_self, positionsbegin_other))
		{
			const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(myCompoundOrigo, compoundOrigos[i]);
			const Float3 relshift = LIMAPOSITIONSYSTEM_HACK::GetRelShiftFromOrigoShift_Float3(querycompound_hyperorigo, myCompoundOrigo);    
            nonbondedNeighbors[nNonbondedNeighbors++] = { static_cast<uint16_t>(query_compound_id), compoundsNParticles[i], relshift };
		}
	}
}



// This kernel creates a new nlist and pushes that to the kernel. Any other kernels that may
// interact with the neighborlist should come AFTER this kernel. Also, they should only run if this
// has run, and thus it is not allowed to comment out this kernel call.
const int threads_in_compoundnlist_kernel = 32;
template <typename BoundaryCondition>
__global__ void updateCompoundNlistsKernel(SimulationDevice* sim_dev) {

	const int n_compounds = sim_dev->boxparams.n_compounds;
	const int compound_id = blockIdx.x * blockDim.x + threadIdx.x;
	const bool compound_active = compound_id < n_compounds;
	const NodeIndex myCompoundOrigo = compound_active
		? sim_dev->boxState->compoundOrigos[compound_id]
		: NodeIndex{};

    uint16_t nNonbondedNeighborsTotal = 0;
	uint16_t nNonbondedNeighborsBatch = 0;
    NlistUtil::IdAndRelshift nonbondedNeighbors[threads_in_compoundnlist_kernel];

	const CompoundInteractionBoundary boundary_self = compound_active
		? sim_dev->boxConfig.compounds[compound_id].interaction_boundary
		: CompoundInteractionBoundary{};

	Float3 key_positions_self[CompoundInteractionBoundary::k];
	if (compound_active)
		getCompoundAbspositions<BoundaryCondition>(*sim_dev, compound_id, key_positions_self, myCompoundOrigo, boundary_self);



	int bonded_compound_ids[Compound::max_bonded_compounds];
	const int n_bonded_compounds = compound_active
		? sim_dev->boxConfig.compounds[compound_id].n_bonded_compounds
		: 0;
	for (int i = 0; i < n_bonded_compounds; i++) {
		bonded_compound_ids[i] = sim_dev->boxConfig.compounds[compound_id].bonded_compound_ids[i];
	}

	__shared__ Float3 key_positions_buffer[threads_in_compoundnlist_kernel * CompoundInteractionBoundary::k];
	__shared__ CompoundInteractionBoundary boundaries[threads_in_compoundnlist_kernel];
	__shared__ NodeIndex compoundOrigos[threads_in_compoundnlist_kernel];
    __shared__ uint8_t compoundNParticles[threads_in_compoundnlist_kernel];

	// Loop over all compounds and add all nearbys
	for (int offset = 0; offset < n_compounds; offset += blockDim.x) {

		// All threads help load a batch of compound_positions
		__syncthreads();

		auto block = cooperative_groups::this_thread_block();
		const int nElements = min(blockDim.x, n_compounds - offset);
		cooperative_groups::memcpy_async(block, compoundOrigos, &sim_dev->boxState->compoundOrigos[offset], sizeof(NodeIndex) * nElements);
		cooperative_groups::memcpy_async(block, boundaries, &sim_dev->compoundsInteractionBoundaryBuffer[offset], sizeof(CompoundInteractionBoundary) * nElements);
		cooperative_groups::memcpy_async(block, compoundNParticles, &sim_dev->nParticlesInCompoundsBuffer[offset], sizeof(uint8_t) * nElements);
		cooperative_groups::wait(block);

        const int query_compound_id = threadIdx.x + offset;
		if (query_compound_id < n_compounds) {
			Float3* const positionsbegin = &key_positions_buffer[threadIdx.x * CompoundInteractionBoundary::k];
			getCompoundAbspositions<BoundaryCondition>(*sim_dev, query_compound_id, positionsbegin, compoundOrigos[threadIdx.x], boundaries[threadIdx.x]);
        }
		__syncthreads();

		// All active-compound threads now loop through the batch
		if (compound_active) {
            addAllNearbyCompounds<BoundaryCondition>(*sim_dev, key_positions_buffer, key_positions_self, offset, n_compounds,
                compound_id, boundary_self, boundaries, n_bonded_compounds, bonded_compound_ids, myCompoundOrigo, compoundOrigos, nNonbondedNeighborsBatch, nonbondedNeighbors, compoundNParticles);

			// Push the neighbors found in this batch
			for (int i = 0; i < nNonbondedNeighborsBatch; i++) {
				sim_dev->nonbondedNeighborsBuffer[compound_id * NlistUtil::maxCompounds + nNonbondedNeighborsTotal + i] = nonbondedNeighbors[i];
			}
			nNonbondedNeighborsTotal += nNonbondedNeighborsBatch;
			nNonbondedNeighborsBatch = 0;

			// TODO only check when not pusjing
            if (nNonbondedNeighborsTotal >= NlistUtil::maxCompounds) {
                printf("\nFailed to insert compound neighbor id!\n");
            }

		}
    }

    if (compound_active)
        sim_dev->nNonbondedNeighborsBuffer[compound_id] = nNonbondedNeighborsTotal;
}

template <typename BoundaryCondition>
__global__ void updateCompoundGridnodes(SimulationDevice* sim_dev) {

    const int n_compounds = sim_dev->boxparams.n_compounds;
    const int compound_id = blockIdx.x * blockDim.x + threadIdx.x;
    const bool compound_active = compound_id < n_compounds;
    const NodeIndex myCompoundOrigo = compound_active
        ? sim_dev->boxState->compoundOrigos[compound_id]
        : NodeIndex{};

	const CompoundInteractionBoundary boundary_self = compound_active
		? sim_dev->boxConfig.compounds[compound_id].interaction_boundary
		: CompoundInteractionBoundary{};

    Float3 key_positions_self[CompoundInteractionBoundary::k];
    if (compound_active)
        getCompoundAbspositions<BoundaryCondition>(*sim_dev, compound_id, key_positions_self, myCompoundOrigo, boundary_self);

    

    int gridnode_ids[NeighborList::max_gridnodes];
    int n_gridnodes = 0;

	// Loop over the nearby gridnodes, and add them if they're within range
	if (compound_active)
	{

		const NodeIndex compound_origo = sim_dev->boxState->compoundOrigos[compound_id];

		for (int x = -GRIDNODE_QUERY_RANGE; x <= GRIDNODE_QUERY_RANGE; x++) {
			for (int y = -GRIDNODE_QUERY_RANGE; y <= GRIDNODE_QUERY_RANGE; y++) {
				for (int z = -GRIDNODE_QUERY_RANGE; z <= GRIDNODE_QUERY_RANGE; z++) {
					NodeIndex query_origo = compound_origo + NodeIndex{ x,y,z };
					BoundaryCondition::applyBC(query_origo);

					// If the query node is NOT inside the box, which happens in some boundary conditions, we cannot continue, 
					// since the node wont exists, and thus compounds are not allowed to access it.
					//printf("%d\n", DeviceConstants::boxSize.blocksPerDim);
					if (!query_origo.isInBox(DeviceConstants::boxSize.blocksPerDim))
						continue;

					const int querynode_id = BoxGrid::Get1dIndex(query_origo, DeviceConstants::boxSize.boxSizeNM_i);
					const Float3 querynode_pos = LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(query_origo);


					if (canCompoundInteractWithPoint<BoundaryCondition>(boundary_self, key_positions_self, querynode_pos)) {
                        gridnode_ids[n_gridnodes++] = querynode_id;
                        if (n_gridnodes > NeighborList::max_gridnodes) {
							sim_dev->signals->critical_error_encountered = true;
						}
					}
				}
			}
		}

        for (int i= 0; i < n_gridnodes; i++)
            sim_dev->compound_neighborlists[compound_id].gridnode_ids[i] = gridnode_ids[i];
        sim_dev->compound_neighborlists[compound_id].n_gridnodes = n_gridnodes;
	}
}



const int nthreads_in_blockgridkernel = 128;
template <typename BoundaryCondition>
__global__ void updateBlockgridKernel(SimulationDevice* sim_dev)
{
	const int block_id = blockIdx.x * blockDim.x + threadIdx.x;
	const bool block_active = block_id < BoxGrid::BlocksTotal(DeviceConstants::boxSize.blocksPerDim);
	const int n_compounds = sim_dev->boxparams.n_compounds;

	CompoundGridNode gridnode;

	const NodeIndex block_origo = block_active
		? BoxGrid::Get3dIndex(block_id, DeviceConstants::boxSize.boxSizeNM_i)
		: NodeIndex{};

	const Float3 block_abspos = LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(block_origo);

	__shared__ Float3 key_positions_buffer[nthreads_in_blockgridkernel * CompoundInteractionBoundary::k];
	__shared__ CompoundInteractionBoundary boundaries[nthreads_in_blockgridkernel];

	// Loop over all compounds in batches
	for (int offset = 0; offset < n_compounds; offset += blockDim.x) {
		const int compound_id = offset + threadIdx.x;
		__syncthreads();
		if (compound_id < n_compounds) {
			const NodeIndex compoundOrigo = sim_dev->boxState->compoundOrigos[compound_id];
			const CompoundInteractionBoundary boundary = sim_dev->compoundsInteractionBoundaryBuffer[compound_id];

			Float3* const positionsbegin = &key_positions_buffer[threadIdx.x * CompoundInteractionBoundary::k];
			getCompoundAbspositions<BoundaryCondition>(*sim_dev, compound_id, positionsbegin, compoundOrigo, boundary);
			boundaries[threadIdx.x] = sim_dev->boxConfig.compounds[compound_id].interaction_boundary;
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
		CompoundGridNode* gridnode_global = BoxGrid::GetNodePtr(sim_dev->compound_grid, block_id);
		gridnode_global->loadData(gridnode);
	}
}





template <typename BoundaryCondition>
void _updateNlists(SimulationDevice* sim_dev, const BoxParams& boxparams, cudaStream_t& s1, cudaStream_t& s2, cudaStream_t& s3)
{
	cudaDeviceSynchronize();
	if (boxparams.n_compounds > 0) {
		const int n_blocks = (boxparams.n_compounds + threads_in_compoundnlist_kernel -1) / threads_in_compoundnlist_kernel;
        updateCompoundNlistsKernel<BoundaryCondition><<<n_blocks, threads_in_compoundnlist_kernel, 0, s1>>>( sim_dev);
        LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: compounds");

        updateCompoundGridnodes<BoundaryCondition><<<(boxparams.n_compounds+31)/32, 32, 0, s2>>>(sim_dev);
        LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: gridnodes");
	}

	if (boxparams.n_solvents > 0) {
		const int n_blocks = BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)) / nthreads_in_blockgridkernel + 1;
        updateBlockgridKernel<BoundaryCondition> <<<n_blocks, nthreads_in_blockgridkernel, 0, s3>>>(sim_dev);
	}

	LIMA_UTILS::genericErrorCheck("Error during updateNlists: blockGrid");
}

void NeighborLists::updateNlists(SimulationDevice* sim_dev, BoundaryConditionSelect bc_select, const BoxParams& boxparams, cudaStream_t& s1, cudaStream_t& s2, cudaStream_t& s3)
{
	switch (bc_select) {
		
	case NoBC: 
        _updateNlists<NoBoundaryCondition>(sim_dev, boxparams, s1, s2, s3);
			break;
		
	case PBC:
        _updateNlists<PeriodicBoundaryCondition>(sim_dev, boxparams, s1, s2, s3);
			break;		
	default:
			throw std::runtime_error("Unsupported boundary condition in updateNlists");
	}
}
