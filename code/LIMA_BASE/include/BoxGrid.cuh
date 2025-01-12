#pragma once

#include "Constants.h"
#include "LimaTypes.cuh"
#include "Bodies.cuh"

// Highest concentration in smtv test is only 48 tinymols, but could be larger with more ions and less solvent
static const int MAX_PARTICLES_IN_BOXGRIDNODE = 64;


// blocks are notcentered 
struct SolventBlock {
	static constexpr int MAX_SOLVENTS_IN_BLOCK = MAX_PARTICLES_IN_BOXGRIDNODE * 2;

	__device__ __host__ void loadMeta(const SolventBlock& block) {
		n_solvents = block.n_solvents;
		if (n_solvents >= MAX_SOLVENTS_IN_BLOCK) {
			printf("Too many solvents in block!\n");
		}
	}
	__device__ __host__ void loadData(const SolventBlock& block) {
		rel_pos[threadIdx.x] = Coord{};	// temp
		if (threadIdx.x < n_solvents) {

			if (block.rel_pos[threadIdx.x] == Coord{ 0 }) {
				printf("Loading zeroes blockid %d nsol %d\n", blockIdx.x, n_solvents);
			}

			rel_pos[threadIdx.x] = block.rel_pos[threadIdx.x];
			ids[threadIdx.x] = block.ids[threadIdx.x];
			atomtypeIds[threadIdx.x] = block.atomtypeIds[threadIdx.x];
		}
	}

	__device__ static void Transfer(const SolventBlock& src, SolventBlock* const dst, 
		const int newIndex, const Coord& relposNext) {
		dst->rel_pos[newIndex] = relposNext;
		dst->ids[newIndex] = src.ids[threadIdx.x];
		dst->atomtypeIds[newIndex] = src.atomtypeIds[threadIdx.x];
	}

	__host__ bool addSolvent(const Coord& rel_position, uint32_t id, uint8_t atomtypeId) {
		if (n_solvents == MAX_SOLVENTS_IN_BLOCK) {
			return false;
		}
		ids[n_solvents] = id;
		atomtypeIds[n_solvents] = atomtypeId;
		rel_pos[n_solvents++] = rel_position;
		return true;
	}
	 
	int n_solvents = 0;
	Coord rel_pos[MAX_SOLVENTS_IN_BLOCK];	// Pos rel to lower left forward side of block, or floor() of pos
	uint32_t ids[MAX_SOLVENTS_IN_BLOCK];
	uint8_t atomtypeIds[MAX_SOLVENTS_IN_BLOCK];

	// Not sure this is the ideal place, as it is an interim and never transferred.. 
	ForceEnergy forceEnergiesCompoundinteractions[MAX_SOLVENTS_IN_BLOCK];
	ForceEnergy forceEnergiesTinymolinteractions[MAX_SOLVENTS_IN_BLOCK];
};






namespace BoxGrid {
	static const int blocksizeNM = 1;
	__device__ __host__ constexpr int NodesPerDim(int boxlenNM) { return boxlenNM; }
	__device__ __host__ constexpr int BlocksTotal(int blocksPerDim) { return blocksPerDim * blocksPerDim * blocksPerDim; }

	__device__ __host__ inline int Get1dIndex(const NodeIndex& index3d, int boxSizeNM) {
		const int bpd = NodesPerDim(boxSizeNM);
		return index3d.x + index3d.y * bpd + index3d.z * bpd * bpd;
	}

	__device__ __host__ inline NodeIndex Get3dIndexWithNNodes(int index1d, int npd) {
		int z = index1d / (npd * npd);
		index1d -= z * npd * npd;
		int y = index1d / npd;
		index1d -= y * npd;
		int x = index1d;
		return NodeIndex{ x, y, z };
	}
	__device__ __host__ inline NodeIndex Get3dIndex(int index1d, int boxlenNM) {
		return Get3dIndexWithNNodes(index1d, NodesPerDim(boxlenNM));
	}




	/// <summary>
	/// Create a BoxGrid of NodeTyoe. This function assumes the type NodeType does not need to be initialized
	/// </summary>
	template <typename NodeType>
	__host__ static NodeType* MallocOnDevice(int boxSizeNM) {
		const int blocksTotal = BoxGrid::NodesPerDim(boxSizeNM) * BoxGrid::NodesPerDim(boxSizeNM) * BoxGrid::NodesPerDim(boxSizeNM);

		NodeType* grid_dev;
		cudaMalloc(&grid_dev, sizeof(NodeType) * blocksTotal);
		cudaMemset(grid_dev, 0, sizeof(NodeType) * blocksTotal);
		return grid_dev;
	}

	// This function assumes the user has used PBC
	template <typename NodeType>
	__device__ __host__ inline NodeType* GetNodePtr(NodeType* grid, const int index1d) {
		return &grid[index1d];
	}

	namespace TinymolBlockAdjacency {
        static const int nNearbyBlocks = 32;

		struct BlockRef {
			int blockId = -1;
			Float3 relShift{};
		};

		// Returns a cudapointer to the data
        BlockRef* PrecomputeNeabyBlockIds(int boxlenNM, float ljCutoffNm);

		__device__ static const BlockRef* GetPtrToNearbyBlockids(int blockId, const BlockRef* const nearbyBlockIdsData) {
			return &nearbyBlockIdsData[blockId * nNearbyBlocks];
		}
	}
};


namespace SolventBlocksCircularQueue {
	static const int STEPS_PER_SOLVENTBLOCKTRANSFER = 5;	// If we go below 2, we might see issue in solventtransfers
	static const int SOLVENTBLOCK_TRANSFERSTEP = STEPS_PER_SOLVENTBLOCKTRANSFER - 1;
	static const int queue_len = STEPS_PER_SOLVENTBLOCKTRANSFER;


	static __host__ int nElementsTotal(int boxlenNM) {
		return BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxlenNM)) * queue_len;
	}

	static std::vector<SolventBlock> createQueue(int boxlenNM) {
		return std::vector<SolventBlock>(nElementsTotal(boxlenNM));
	}

	__device__ __host__ bool static isTransferStep(int64_t step) {
		return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == SOLVENTBLOCK_TRANSFERSTEP;
	}
	__device__ bool static isFirstStepAfterTransfer(int64_t step) {
		return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == 0;
	}

	// This function assumes the user has used PBC
	__device__ __host__ static SolventBlock* getBlockPtr(SolventBlock* queue, int blocksPerDim, const size_t index1d, const size_t step) {
		const size_t step_offset = (step % queue_len) * BoxGrid::BlocksTotal(blocksPerDim);
		return &queue[index1d + step_offset];
	}

	__host__ static SolventBlock& GetBlockRef(std::vector<SolventBlock>& queue, NodeIndex index3d, const int64_t step, int boxSizeNm) {
		if (index3d.x >= boxSizeNm || index3d.y >= boxSizeNm || index3d.z >= boxSizeNm
			|| index3d.x < 0 || index3d.y < 0 || index3d.z < 0) {
			throw std::runtime_error("Bad 3d index for blockptr\n");
		}
		return *getBlockPtr(queue.data(), boxSizeNm, BoxGrid::Get1dIndex(index3d, boxSizeNm), step);
	}
};

