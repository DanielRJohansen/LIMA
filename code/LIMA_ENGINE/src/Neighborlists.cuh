#pragma once

#include "EngineBodies.cuh"
#include "BoundaryCondition.cuh"
#include "LimaPositionSystem.cuh"
#include "Simulation.cuh"
#include "SimulationDevice.cuh"



namespace NeighborList {

	static_assert(::MAX_COMPOUNDS <= UINT16_MAX, "Neighborlist cannot handle such large compound ids");

	// This value is only high because the electrostaticmanyparticles has the particles as compounds, when they should be tinymol. But pme doesnt yet support tinymol...
	static const int compoundsMaxNearbyCompounds = 256;	// TODO: We need to work on getting this number down!
	static const int gridnodesMaxNearbyCompounds = 128;

	struct IdAndRelshift {
		uint16_t id;
		uint16_t nParticles;
		Float3 relShift;
	};

	static const int compoundsMaxNearbyGridnodes = 128; // OPTIM too high

	struct CompoundInfo {
		Float3 keyPositions[CompoundInteractionBoundary::k];	// [nm] absolute pos
		float radii[CompoundInteractionBoundary::k];			// [nm]
        uint16_t compoundId = -1;
	};

    struct alignas(32) Gridnode {
        static const int maxCompoundsInNode = 16;

		CompoundInfo compoundInfos[maxCompoundsInNode];
		int nCompoundsInNode;
	};

	struct Buffers {
		int* compoundsNearbyGridnodes = nullptr;
		int* compoundsNNearbyGridnodes = nullptr;

		IdAndRelshift* compoundsNeighborNonbondedCompounds = nullptr;
		uint16_t* compoundsNNeighborNonbondedCompounds = nullptr;
		
		uint16_t* gridnodesNearbyCompounds = nullptr;
		int* gridnodesNNearbyCompounds = nullptr;

		__device__ constexpr Buffers() {}
		__host__ Buffers(int nCompounds, int nNodesTotal) {
			cudaMalloc(&compoundsNearbyGridnodes, nCompounds * compoundsMaxNearbyGridnodes * sizeof(int));
			cudaMalloc(&compoundsNNearbyGridnodes, nCompounds * sizeof(int));
			cudaMemset(compoundsNNearbyGridnodes, 0, nCompounds * sizeof(int));

			cudaMalloc(&compoundsNeighborNonbondedCompounds, nCompounds * compoundsMaxNearbyCompounds * sizeof(IdAndRelshift));
			cudaMalloc(&compoundsNNeighborNonbondedCompounds, nCompounds * sizeof(uint16_t));
			cudaMemset(compoundsNNeighborNonbondedCompounds, 0, nCompounds * sizeof(uint16_t));

			cudaMalloc(&gridnodesNearbyCompounds, sizeof(uint16_t) * gridnodesMaxNearbyCompounds * nNodesTotal);
			cudaMalloc(&gridnodesNNearbyCompounds, sizeof(int) * nNodesTotal);
			cudaMemset(gridnodesNNearbyCompounds, 0, sizeof(int) * nNodesTotal);
		}

		__host__ void Free() {
			cudaFree(compoundsNearbyGridnodes);
			cudaFree(compoundsNNearbyGridnodes);

			cudaFree(compoundsNeighborNonbondedCompounds);
			cudaFree(compoundsNNeighborNonbondedCompounds);

			cudaFree(gridnodesNearbyCompounds);
			cudaFree(gridnodesNNearbyCompounds);
		}
	};

	class Controller {
		Gridnode* grid = nullptr;
		Buffers buffers;

		template <typename BoundaryCondition>
		void _UpdateNlist(SimulationDevice* const simDev, const BoxParams&, std::array<cudaStream_t, 5>&);

	public:


		Controller(const BoxParams& boxParams) : buffers(boxParams.n_compounds, BoxGrid::BlocksTotal(boxParams.boxSize))
		{
			cudaMalloc(&grid, BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxParams.boxSize)) * sizeof(Gridnode));
			cudaMemset(grid, 0, sizeof(Gridnode) * BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxParams.boxSize)));

			cudaDeviceSynchronize();
		}
		~Controller() {
			buffers.Free();
		}

		void UpdateNlist(SimulationDevice* const simDev, const BoxParams&, BoundaryConditionSelect, std::array<cudaStream_t, 5>&);

		const Buffers& GetBuffers() const { return buffers; }
	};
};









__global__ void PushCompoundsToGrid(const SimulationDevice* const simDev, NeighborList::Gridnode* const grid, int nCompounds) {
	const int compoundId = blockIdx.x * blockDim.x + threadIdx.x;

	if (compoundId < nCompounds) {
		const NodeIndex compoundOrigo = simDev->boxState.compoundOrigos[compoundId];
		NeighborList::CompoundInfo compoundInfo;
		compoundInfo.compoundId = compoundId;

		for (int i = 0; i < CompoundInteractionBoundary::k; i++) {
			const int particleLocalIndex = simDev->compoundsInteractionBoundaryBuffer[compoundId].key_particle_indices[i];
			compoundInfo.keyPositions[i] = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(compoundOrigo, simDev->boxState.compoundsRelposNm[compoundId * MAX_COMPOUND_PARTICLES + particleLocalIndex]);
			compoundInfo.radii[i] = simDev->compoundsInteractionBoundaryBuffer[compoundId].radii[i];
		}

		const int gridIndex = BoxGrid::Get1dIndex(compoundOrigo, DeviceConstants::boxSize.boxSizeNM_i); // CompoundOrigos are safely assumed always inside box
		const int localIndexInNode = atomicAdd(&grid[gridIndex].nCompoundsInNode, 1);
		grid[gridIndex].compoundInfos[localIndexInNode] = compoundInfo;

		if constexpr (!LIMA_PUSH) {
			if (localIndexInNode >= NeighborList::Gridnode::maxCompoundsInNode) {
				printf("Too many compounds in node");
			}
		}
	}
}



// Assumes the compound is active
template <typename BoundaryCondition>
__device__ void getCompoundAbspositions(const SimulationDevice& simDev, int compound_id, Float3* result, const NodeIndex& compoundOrigo, const CompoundInteractionBoundary& compoundInteractionBoundary)
{
	const Float3* const relPositions = &simDev.boxState.compoundsRelposNm[compound_id * MAX_COMPOUND_PARTICLES];

	for (int i = 0; i < CompoundInteractionBoundary::k; i++) {
		const int particle_index = compoundInteractionBoundary.key_particle_indices[i];
		const Float3 abspos = compoundOrigo.toFloat3() + relPositions[particle_index];
		result[i] = abspos;
	}
}

template <typename BoundaryCondition>
__device__ bool canCompoundsInteract(const NeighborList::CompoundInfo& left, const NeighborList::CompoundInfo& right)
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
__device__ bool canCompoundInteractWithPoint(const NeighborList::CompoundInfo& queryCompound, const Float3& point)
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
__global__ void updateCompoundNlistsKernel(SimulationDevice* simDev, const NeighborList::Gridnode* const grid, int nCompounds, int nodesPerDim, NeighborList::Buffers nlistBuffers) {

    const bool compoundActive = blockIdx.x * blockDim.x + threadIdx.x < nCompounds;

	if (!compoundActive)
		return;

    const uint16_t compoundId = blockIdx.x * blockDim.x + threadIdx.x;
    const NodeIndex myCompoundOrigo = simDev->boxState.compoundOrigos[compoundId];

	uint16_t nNonbondedNeighborsTotal = 0;
	NeighborList::IdAndRelshift nonbondedNeighbors[NeighborList::compoundsMaxNearbyCompounds];

    const CompoundInteractionBoundary boundary_self = simDev->compoundsInteractionBoundaryBuffer[compoundId];

	Float3 key_positions_self[CompoundInteractionBoundary::k];
    getCompoundAbspositions<BoundaryCondition>(*simDev, compoundId, key_positions_self, myCompoundOrigo, boundary_self);

	// Load bonded compounds so we dont add them again
    uint16_t bondedCompoundIds[Compound::max_bonded_compounds];
    const int n_bonded_compounds = simDev->boxConfig.compounds[compoundId].n_bonded_compounds;
	for (int i = 0; i < n_bonded_compounds; i++) {
		bondedCompoundIds[i] = simDev->boxConfig.compounds[compoundId].bonded_compound_ids[i];
	}


	NeighborList::CompoundInfo myCompoundInfo;
	myCompoundInfo.compoundId = compoundId;
	for (int i = 0; i < CompoundInteractionBoundary::k; i++) {
		myCompoundInfo.keyPositions[i] = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(myCompoundOrigo, simDev->boxState.compoundsRelposNm[compoundId * MAX_COMPOUND_PARTICLES + i]);
		myCompoundInfo.radii[i] = simDev->compoundsInteractionBoundaryBuffer[compoundId].radii[i];
	}


	// load one full x row at a time
	//Gridnode gridnodeBuffer[5]

	const int range = 2;
	for (int zOff = -range; zOff <= range; zOff++) {
		for (int yOff = -range; yOff <= range; yOff++) {
			for (int xOff = -range; xOff <= range; xOff++) {
				const NodeIndex queryNode = BoundaryCondition::applyBC(myCompoundOrigo + NodeIndex{ xOff, yOff, zOff }, nodesPerDim);

				const int queryNodeIndex = BoxGrid::Get1dIndex(queryNode, nodesPerDim);
				const NeighborList::Gridnode& gridnode = grid[queryNodeIndex];

				for (int i = 0; i < gridnode.nCompoundsInNode; i++) {
					const NeighborList::CompoundInfo& queryCompoundInfo = gridnode.compoundInfos[i];

                    const uint16_t queryCompoundId = queryCompoundInfo.compoundId;

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

						const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(myCompoundOrigo, simDev->boxState.compoundOrigos[queryCompoundId]);
						const Float3 relshift = LIMAPOSITIONSYSTEM_HACK::GetRelShiftFromOrigoShift_Float3(querycompound_hyperorigo, myCompoundOrigo);

						if constexpr (!LIMA_PUSH) {
							if (nNonbondedNeighborsTotal >= NeighborList::compoundsMaxNearbyCompounds) 
								printf("Too many nonbonded neighbors %d %d\n", nNonbondedNeighborsTotal, compoundId);							
						}

						nonbondedNeighbors[nNonbondedNeighborsTotal++] = { static_cast<uint16_t>(queryCompoundId), queryCompoundNParticles, relshift };
					}
				}
			}
		}
	}

	for (int i = 0; i < nNonbondedNeighborsTotal; i++) {
		nlistBuffers.compoundsNeighborNonbondedCompounds[compoundId * NeighborList::compoundsMaxNearbyCompounds + i] = nonbondedNeighbors[i];
	}
	nlistBuffers.compoundsNNeighborNonbondedCompounds[compoundId] = nNonbondedNeighborsTotal;
}

__global__ void SortNonbondedNeighborcompoundIds(NeighborList::Buffers nlistBuffers) {
    __shared__ NeighborList::IdAndRelshift neighbors[NeighborList::compoundsMaxNearbyCompounds];

    const int compoundId = blockIdx.x;

    auto block = cooperative_groups::this_thread_block();
    cooperative_groups::memcpy_async(block, neighbors, &nlistBuffers.compoundsNeighborNonbondedCompounds[compoundId * NeighborList::compoundsMaxNearbyCompounds], sizeof(NeighborList::IdAndRelshift) * NeighborList::compoundsMaxNearbyCompounds);
    cooperative_groups::wait(block);

    if (threadIdx.x >= nlistBuffers.compoundsNNeighborNonbondedCompounds[compoundId])
        neighbors[threadIdx.x].id = UINT16_MAX;
    __syncthreads();

    LAL::Sort(neighbors, NeighborList::compoundsMaxNearbyCompounds, [](const NeighborList::IdAndRelshift& a){return a.id;});

    cooperative_groups::memcpy_async(block, &nlistBuffers.compoundsNeighborNonbondedCompounds[compoundId * NeighborList::compoundsMaxNearbyCompounds], neighbors, sizeof(NeighborList::IdAndRelshift) * NeighborList::compoundsMaxNearbyCompounds);
}














template <typename BoundaryCondition>
__global__ void UpdateCompoundsNeighborGridnodes(SimulationDevice* simDev, NeighborList::Buffers nlistBuffers) {

    const int n_compounds = simDev->boxparams.n_compounds;
    const int compound_id = blockIdx.x * blockDim.x + threadIdx.x;
    const bool compound_active = compound_id < n_compounds;
    const NodeIndex myCompoundOrigo = compound_active
        ? simDev->boxState.compoundOrigos[compound_id]
        : NodeIndex{};

	const CompoundInteractionBoundary boundary_self = compound_active
		? simDev->boxConfig.compounds[compound_id].interaction_boundary
		: CompoundInteractionBoundary{};

    Float3 key_positions_self[CompoundInteractionBoundary::k];
    if (compound_active)
        getCompoundAbspositions<BoundaryCondition>(*simDev, compound_id, key_positions_self, myCompoundOrigo, boundary_self);

    

    int compoundsNearbyGridnodes[NeighborList::compoundsMaxNearbyGridnodes];
    int nGridnodes = 0;

	// Loop over the nearby gridnodes, and add them if they're within range
	if (compound_active)
	{
		const NodeIndex compound_origo = simDev->boxState.compoundOrigos[compound_id];
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
						compoundsNearbyGridnodes[nGridnodes++] = querynode_id;
                        if (nGridnodes > NeighborList::compoundsMaxNearbyGridnodes) {
							simDev->signals->critical_error_encountered = true;
						}
					}
				}
			}
		}

        for (int i= 0; i < nGridnodes; i++)
			nlistBuffers.compoundsNearbyGridnodes[compound_id * NeighborList::compoundsMaxNearbyGridnodes + i] = compoundsNearbyGridnodes[i];
		nlistBuffers.compoundsNNearbyGridnodes[compound_id] = nGridnodes;
	}
}



const int nthreads_in_blockgridkernel = 128;
template <typename BoundaryCondition>
__global__ void updateBlockgridKernel(const NeighborList::Gridnode* const grid, int nodesPerDim, NeighborList::Buffers nlistBuffers)
{
	const int block_id = blockIdx.x * blockDim.x + threadIdx.x;
	const bool block_active = block_id < BoxGrid::BlocksTotal(nodesPerDim);
	if (!block_active)
		return;

	
	uint16_t nearbyCompoundIds[NeighborList::gridnodesMaxNearbyCompounds];
	int nNearbyCompounds = 0;

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
				const NeighborList::Gridnode& gridnode = grid[queryNodeIndex];

				for (int i = 0; i < gridnode.nCompoundsInNode; i++) {
					const NeighborList::CompoundInfo& queryCompoundInfo = gridnode.compoundInfos[i];

					const int queryCompoundId = queryCompoundInfo.compoundId;

					// OPTIM: This should be treated as a point but as a box. This would lower the number of hits, making the force kernel faster
					if (canCompoundInteractWithPoint<BoundaryCondition>(queryCompoundInfo, block_abspos)) {

						if constexpr (!LIMA_PUSH) {
							if (nNearbyCompounds >= NeighborList::gridnodesMaxNearbyCompounds) {
								printf("Too many compounds in blockgrid\n");
							}
						}
						nearbyCompoundIds[nNearbyCompounds++] = queryCompoundId;
					}
				}
			}
		}
	}

	for (int i = 0; i < nNearbyCompounds; i++) {
		nlistBuffers.gridnodesNearbyCompounds[block_id * NeighborList::gridnodesMaxNearbyCompounds + i] = nearbyCompoundIds[i];
	}
	nlistBuffers.gridnodesNNearbyCompounds[block_id] = nNearbyCompounds;
}

__global__ void SortCompoundGridnodes(NeighborList::Buffers nlistBuffers) {
	__shared__ uint16_t compoundIds[NeighborList::gridnodesMaxNearbyCompounds];

	const int nodeIndex = blockIdx.x;


	auto block = cooperative_groups::this_thread_block();
	cooperative_groups::memcpy_async(block, compoundIds, &nlistBuffers.gridnodesNearbyCompounds[nodeIndex * NeighborList::gridnodesMaxNearbyCompounds], sizeof(uint16_t) * NeighborList::gridnodesMaxNearbyCompounds);
	cooperative_groups::wait(block);

	if (threadIdx.x >= nlistBuffers.gridnodesNNearbyCompounds[nodeIndex])
		compoundIds[threadIdx.x] = UINT16_MAX;
	__syncthreads();

	static_assert(LAL::isPowerOf2(NeighborList::gridnodesMaxNearbyCompounds), "can only sort powers of 2");
    LAL::Sort(compoundIds, NeighborList::gridnodesMaxNearbyCompounds, [](const uint16_t& a) {return a;});

	cooperative_groups::memcpy_async(block, &nlistBuffers.gridnodesNearbyCompounds[nodeIndex * NeighborList::gridnodesMaxNearbyCompounds], compoundIds, sizeof(uint16_t) * NeighborList::gridnodesMaxNearbyCompounds);
}



template <typename BoundaryCondition>
void NeighborList::Controller::_UpdateNlist(SimulationDevice* simDev, const BoxParams& boxparams, std::array<cudaStream_t,5>& streams)
{
	cudaDeviceSynchronize();
	if (boxparams.n_compounds > 0) {
		// Stream one
		cudaMemsetAsync(grid, 0, BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)) * sizeof(NeighborList::Gridnode), streams[0]);

		PushCompoundsToGrid<<<(boxparams.n_compounds + 31) / 32, 32, 0, streams[0] >> > (simDev, grid, boxparams.n_compounds);
		LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: PushCompoundsToGrid");

		cudaStreamSynchronize(streams[0]);

		updateCompoundNlistsKernel<BoundaryCondition> << <(boxparams.n_compounds + 31) / 32, 32, 0, streams[0] >> > (simDev, grid, boxparams.n_compounds, boxparams.boxSize, buffers);
		LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: updateCompoundNlistsKernel");

        SortNonbondedNeighborcompoundIds<<<boxparams.n_compounds, NeighborList::compoundsMaxNearbyCompounds, 0, streams[0]>>>(buffers);
		LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: SortNonbondedNeighborcompoundIds");

		// Stream 2
        UpdateCompoundsNeighborGridnodes<BoundaryCondition><<<(boxparams.n_compounds+31)/32, 32, 0, streams[1] >>>(simDev, buffers);
        LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: UpdateCompoundsNeighborGridnodes");


		if (boxparams.n_solvents > 0) {
			// Stream 3
			const int n_blocks = BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)) / nthreads_in_blockgridkernel + 1;
			updateBlockgridKernel<BoundaryCondition> << <n_blocks, nthreads_in_blockgridkernel, 0, streams[2] >> > (grid, boxparams.boxSize, buffers);
			LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: updateBlockgridKernel");

			SortCompoundGridnodes << < BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxparams.boxSize)), NeighborList::gridnodesMaxNearbyCompounds, 0, streams[2] >> > (buffers);
			LIMA_UTILS::genericErrorCheckNoSync("Error during updateNlists: SortCompoundGridnodes");
		}
	}



	LIMA_UTILS::genericErrorCheck("Error during updateNlists: blockGrid");
}

void NeighborList::Controller::UpdateNlist(SimulationDevice* simDev, const BoxParams& boxparams, BoundaryConditionSelect bc_select, std::array<cudaStream_t,5>& streams)
{
	switch (bc_select) {
		
	case NoBC: 
        _UpdateNlist<NoBoundaryCondition>(simDev, boxparams, streams);
			break;
		
	case PBC:
        _UpdateNlist<PeriodicBoundaryCondition>(simDev, boxparams, streams);
			break;		
	default:
			throw std::runtime_error("Unsupported boundary condition in updateNlists");
	}
}
