#pragma once

#include "BoundaryCondition.cuh"
#include "LimaPositionSystem.cuh"
#include "Simulation.cuh"
#include "SimulationDevice.cuh"



namespace NeighborLists {
    void updateNlists(SimulationDevice*, BoundaryConditionSelect, const BoxParams&, cudaStream_t& s1, cudaStream_t& s2, cudaStream_t& s3, NeighborLists::Gridnode* const grid);

	struct CompoundInfo {
		Float3 keyPositions[CompoundInteractionBoundary::k];	// [nm] absolute pos
		float radii[CompoundInteractionBoundary::k];			// [nm]
		int compoundId = -1;
	};

    struct alignas(32) Gridnode {
        static const int maxCompoundsInNode = 16;

		CompoundInfo compoundInfos[maxCompoundsInNode];
		int nCompoundsInNode;
	};
    //static_assert(sizeof(Gridnode) == 256);
};







__global__ void PushCompoundsToGrid(const SimulationDevice* const simDev, NeighborLists::Gridnode* const grid, int nCompounds) {
	const int compoundId = blockIdx.x * blockDim.x + threadIdx.x;

	if (compoundId < nCompounds) {
		const NodeIndex compoundOrigo = simDev->boxState->compoundOrigos[compoundId];
		NeighborLists::CompoundInfo compoundInfo;
		compoundInfo.compoundId = compoundId;

		for (int i = 0; i < CompoundInteractionBoundary::k; i++) {
			const int particleLocalIndex = simDev->compoundsInteractionBoundaryBuffer[compoundId].key_particle_indices[i];
			compoundInfo.keyPositions[i] = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(compoundOrigo, simDev->boxState->compoundsRelposNm[compoundId * MAX_COMPOUND_PARTICLES + particleLocalIndex]);
			compoundInfo.radii[i] = simDev->boxConfig.compounds[compoundId].interaction_boundary.radii[i];
		}

		const int gridIndex = BoxGrid::Get1dIndex(compoundOrigo, DeviceConstants::boxSize.boxSizeNM_i); // CompoundOrigos are safely assumed always inside box
		const int localIndexInNode = atomicAdd(&grid[gridIndex].nCompoundsInNode, 1);
		grid[gridIndex].compoundInfos[localIndexInNode] = compoundInfo;

		if constexpr (!LIMA_PUSH) {
			if (localIndexInNode >= NeighborLists::Gridnode::maxCompoundsInNode) {
				printf("Too many compounds in node");
			}
		}
	}
}



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
__device__ bool canCompoundsInteract(const NeighborLists::CompoundInfo& left, const NeighborLists::CompoundInfo& right)
{
	for (int ileft = 0; ileft < CompoundInteractionBoundary::k; ileft++) {
		for (int iright = 0; iright < CompoundInteractionBoundary::k; iright++) {

            const float distSquared = LIMAPOSITIONSYSTEM::calcHyperDistSquaredNM<BoundaryCondition>(left.keyPositions[ileft], right.keyPositions[iright]);
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

template <typename BoundaryCondition>
__device__ bool canCompoundInteractWithPoint(const NeighborLists::CompoundInfo& queryCompound, const Float3& point)
{
	for (int ileft = 0; ileft < CompoundInteractionBoundary::k; ileft++) {
		const float distSq = LIMAPOSITIONSYSTEM::calcHyperDistSquaredNM<BoundaryCondition>(queryCompound.keyPositions[ileft], point);
		const float max_dist = DeviceConstants::cutoffNM + queryCompound.radii[ileft];

		if (distSq < max_dist * max_dist)
			return true;
	}
	return false;
}

template <typename BoundaryCondition>
__global__ void updateCompoundNlistsKernel(SimulationDevice* simDev, const NeighborLists::Gridnode* const grid, int nCompounds, int nodesPerDim) {

	const int compoundId = blockIdx.x * blockDim.x + threadIdx.x;
	const bool compoundActive = compoundId < nCompounds;

	if (!compoundActive)
		return;

	const NodeIndex myCompoundOrigo = compoundActive ? simDev->boxState->compoundOrigos[compoundId] : NodeIndex{};

	uint16_t nNonbondedNeighborsTotal = 0;
	NlistUtil::IdAndRelshift nonbondedNeighbors[NlistUtil::maxCompounds];

	const CompoundInteractionBoundary boundary_self = compoundActive ? simDev->compoundsInteractionBoundaryBuffer[compoundId] : CompoundInteractionBoundary{};

	Float3 key_positions_self[CompoundInteractionBoundary::k];
	if (compoundActive)
		getCompoundAbspositions<BoundaryCondition>(*simDev, compoundId, key_positions_self, myCompoundOrigo, boundary_self);

	// Load bonded compounds so we dont add them again
	int bondedCompoundIds[Compound::max_bonded_compounds];
	const int n_bonded_compounds = compoundActive ? simDev->boxConfig.compounds[compoundId].n_bonded_compounds : 0;
	for (int i = 0; i < n_bonded_compounds; i++) {
		bondedCompoundIds[i] = simDev->boxConfig.compounds[compoundId].bonded_compound_ids[i];
	}


	NeighborLists::CompoundInfo myCompoundInfo;
	myCompoundInfo.compoundId = compoundId;
	for (int i = 0; i < CompoundInteractionBoundary::k; i++) {
		myCompoundInfo.keyPositions[i] = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(myCompoundOrigo, simDev->boxState->compoundsRelposNm[compoundId * MAX_COMPOUND_PARTICLES + i]);
		myCompoundInfo.radii[i] = simDev->boxConfig.compounds[compoundId].interaction_boundary.radii[i];
	}


	// load one full x row at a time
	//Gridnode gridnodeBuffer[5]

	const int range = 2;
	for (int zOff = -range; zOff <= range; zOff++) {
		for (int yOff = -range; yOff <= range; yOff++) {
			for (int xOff = -range; xOff <= range; xOff++) {
				const NodeIndex queryNode = BoundaryCondition::applyBC(myCompoundOrigo + NodeIndex{ xOff, yOff, zOff }, nodesPerDim);

				const int queryNodeIndex = BoxGrid::Get1dIndex(queryNode, nodesPerDim);
				const NeighborLists::Gridnode& gridnode = grid[queryNodeIndex];

				for (int i = 0; i < gridnode.nCompoundsInNode; i++) {
					const NeighborLists::CompoundInfo& queryCompoundInfo = gridnode.compoundInfos[i];

					const int queryCompoundId = queryCompoundInfo.compoundId;

					if (queryCompoundId == compoundId) { continue; }	// dont add self to self
					// Dont add bonded compounds to list again
					bool is_bonded_to_query = false;
					for (int j = 0; j < n_bonded_compounds; j++) {
						if (queryCompoundId == bondedCompoundIds[j]) {
							is_bonded_to_query = true;
							break;
						}
					}
					if (is_bonded_to_query) { continue; }



					if (canCompoundsInteract<BoundaryCondition>(myCompoundInfo, queryCompoundInfo)) {
						const uint8_t queryCompoundNParticles = simDev->nParticlesInCompoundsBuffer[queryCompoundId];

						const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(myCompoundOrigo, simDev->boxState->compoundOrigos[queryCompoundId]);
						const Float3 relshift = LIMAPOSITIONSYSTEM_HACK::GetRelShiftFromOrigoShift_Float3(querycompound_hyperorigo, myCompoundOrigo);

						if constexpr (!LIMA_PUSH) {
							if (nNonbondedNeighborsTotal >= NlistUtil::maxCompounds) 
								printf("Too many nonbonded neighbors %d %d\n", nNonbondedNeighborsTotal, compoundId);							
						}

						nonbondedNeighbors[nNonbondedNeighborsTotal++] = { static_cast<uint16_t>(queryCompoundId), queryCompoundNParticles, relshift };
					}
				}
			}
		}
	}

	for (int i = 0; i < nNonbondedNeighborsTotal; i++) {
		simDev->nonbondedNeighborsBuffer[compoundId * NlistUtil::maxCompounds + i] = nonbondedNeighbors[i];
	}

	simDev->nNonbondedNeighborsBuffer[compoundId] = nNonbondedNeighborsTotal;
}




__device__ inline void Sort(NlistUtil::IdAndRelshift* const data, int nElements) { // Assuming that data is already in shared memory.
    int tid = threadIdx.x;
    for (int k = 2; k <= nElements; k <<= 1) {
        for (int j = k >> 1; j > 0; j >>= 1) {
            int ixj = tid ^ j;
            if (ixj > tid) {
                if ((tid & k) == 0) {
                    if (data[tid].id > data[ixj].id) {
                        // Swap data[tid] and data[ixj]
                        auto temp = data[tid];
                        data[tid] = data[ixj];
                        data[ixj] = temp;
                    }
                }
                else {
                    if (data[tid].id < data[ixj].id) {
                        // Swap data[tid] and data[ixj]
                        auto temp = data[tid];
                        data[tid] = data[ixj];
                        data[ixj] = temp;
                    }
                }
            }
            __syncthreads(); // Synchronize to ensure all threads complete this step before moving on
        }
    }
}

__global__ void SortNonbondedNeighborcompoundIds(SimulationDevice* const simDev) {
    __shared__ NlistUtil::IdAndRelshift neighbors[NlistUtil::maxCompounds];

    const int compoundId = blockIdx.x;

    auto block = cooperative_groups::this_thread_block();
    cooperative_groups::memcpy_async(block, neighbors, &simDev->nonbondedNeighborsBuffer[compoundId * NlistUtil::maxCompounds], sizeof(NlistUtil::IdAndRelshift) * NlistUtil::maxCompounds);
    cooperative_groups::wait(block);

    if (threadIdx.x >= simDev->nNonbondedNeighborsBuffer[compoundId])
        neighbors[threadIdx.x].id = UINT16_MAX;
    __syncthreads();

    Sort(neighbors, NlistUtil::maxCompounds);

    cooperative_groups::memcpy_async(block, &simDev->nonbondedNeighborsBuffer[compoundId * NlistUtil::maxCompounds], neighbors, sizeof(NlistUtil::IdAndRelshift) * NlistUtil::maxCompounds);
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
__global__ void updateBlockgridKernel(SimulationDevice* sim_dev, const NeighborLists::Gridnode* const grid, int nodesPerDim)
{
	const int block_id = blockIdx.x * blockDim.x + threadIdx.x;
	const bool block_active = block_id < BoxGrid::BlocksTotal(nodesPerDim);
	if (!block_active)
		return;

	CompoundGridNode compoundGridnode;

	const NodeIndex block_origo = block_active
		? BoxGrid::Get3dIndex(block_id, DeviceConstants::boxSize.boxSizeNM_i)
		: NodeIndex{};

	const Float3 block_abspos = LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(block_origo);


	const int range = 2;
	for (int zOff = -range; zOff <= range; zOff++) {
		for (int yOff = -range; yOff <= range; yOff++) {
			for (int xOff = -range; xOff <= range; xOff++) {
				const NodeIndex queryNode = BoundaryCondition::applyBC(block_origo + NodeIndex{ xOff, yOff, zOff }, nodesPerDim);

				const int queryNodeIndex = BoxGrid::Get1dIndex(queryNode, nodesPerDim);
				const NeighborLists::Gridnode& gridnode = grid[queryNodeIndex];

				for (int i = 0; i < gridnode.nCompoundsInNode; i++) {
					const NeighborLists::CompoundInfo& queryCompoundInfo = gridnode.compoundInfos[i];

					const int queryCompoundId = queryCompoundInfo.compoundId;

					// OPTIM: This should be treated as a point but as a box. This would lower the number of hits, making the force kernel faster
					if (canCompoundInteractWithPoint<BoundaryCondition>(queryCompoundInfo, block_abspos)) {
						if (!compoundGridnode.addNearbyCompound(queryCompoundId)) {
							printf("Adding too many compound to blockgrid\n");
						}
					}
				}
			}
		}
	}

	CompoundGridNode* gridnode_global = BoxGrid::GetNodePtr(sim_dev->compound_grid, block_id);
	gridnode_global->loadData(compoundGridnode);
}

__global__ void SortCompoundGridnodes(SimulationDevice* const simDev) {
	__shared__ uint16_t compoundIds[CompoundGridNode::max_nearby_compounds];

	const int nodeIndex = blockIdx.x;


	auto block = cooperative_groups::this_thread_block();
	cooperative_groups::memcpy_async(block, compoundIds, simDev->compound_grid[nodeIndex].compoundidsWithinLjCutoff, sizeof(uint16_t) * CompoundGridNode::max_nearby_compounds);
	cooperative_groups::wait(block);

	if (threadIdx.x >= simDev->compound_grid[nodeIndex].n_nearby_compounds)
		compoundIds[threadIdx.x] = UINT16_MAX;
	__syncthreads();

	LAL::Sort<uint16_t>(compoundIds, CompoundGridNode::max_nearby_compounds);

	cooperative_groups::memcpy_async(block, simDev->compound_grid[nodeIndex].compoundidsWithinLjCutoff, compoundIds, sizeof(uint16_t) * CompoundGridNode::max_nearby_compounds);
}



template <typename BoundaryCondition>
void _updateNlists(SimulationDevice* sim_dev, const BoxParams& boxparams, cudaStream_t& s1, cudaStream_t& s2, cudaStream_t& s3, NeighborLists::Gridnode* const grid)
{
	cudaDeviceSynchronize();
	if (boxparams.n_compounds > 0) {
		// Stream one
		cudaMemsetAsync(grid, 0, BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)) * sizeof(NeighborLists::Gridnode), s1);

		PushCompoundsToGrid<<<(boxparams.n_compounds + 31) / 32, 32, 0, s1 >> > (sim_dev, grid, boxparams.n_compounds);
		LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: PushCompoundsToGrid");

		cudaStreamSynchronize(s1);

		updateCompoundNlistsKernel<BoundaryCondition> << <(boxparams.n_compounds + 31) / 32, 32, 0, s1 >> > (sim_dev, grid, boxparams.n_compounds, boxparams.boxSize);
		LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: updateCompoundNlistsKernel");
		//cudaDeviceSynchronize();

        SortNonbondedNeighborcompoundIds<<<boxparams.n_compounds, NlistUtil::maxCompounds, 0, s1>>>(sim_dev);
		LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: SortNonbondedNeighborcompoundIds");

		// Stream 2
        updateCompoundGridnodes<BoundaryCondition><<<(boxparams.n_compounds+31)/32, 32, 0, s2>>>(sim_dev);
        LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: updateCompoundGridnodes");


		if (boxparams.n_solvents > 0) {
			// Stream 3
			const int n_blocks = BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)) / nthreads_in_blockgridkernel + 1;
			updateBlockgridKernel<BoundaryCondition> << <n_blocks, nthreads_in_blockgridkernel, 0, s3 >> > (sim_dev, grid, boxparams.boxSize);
			LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: updateBlockgridKernel");

			SortCompoundGridnodes << < BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), CompoundGridNode::max_nearby_compounds, 0, s3 >> > (sim_dev);
			LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: SortCompoundGridnodes");
		}
	}



	LIMA_UTILS::genericErrorCheck("Error during updateNlists: blockGrid");
}

void NeighborLists::updateNlists(SimulationDevice* sim_dev, BoundaryConditionSelect bc_select, const BoxParams& boxparams, cudaStream_t& s1, cudaStream_t& s2, cudaStream_t& s3, NeighborLists::Gridnode* const grid)
{
	switch (bc_select) {
		
	case NoBC: 
        _updateNlists<NoBoundaryCondition>(sim_dev, boxparams, s1, s2, s3, grid);
			break;
		
	case PBC:
        _updateNlists<PeriodicBoundaryCondition>(sim_dev, boxparams, s1, s2, s3, grid);
			break;		
	default:
			throw std::runtime_error("Unsupported boundary condition in updateNlists");
	}
}
