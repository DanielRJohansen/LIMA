#pragma once

#include "Constants.h"
#include "LimaTypes.cuh"
#include "Bodies.cuh"

static const int MAX_PARTICLES_IN_BOXGRIDNODE = 64 + 32+32;


// blocks are notcentered 
struct SolventBlock {
	static constexpr int MAX_SOLVENTS_IN_BLOCK = MAX_PARTICLES_IN_BOXGRIDNODE;

	__device__ __host__ void loadMeta(const SolventBlock& block) {
		origo = block.origo; // Not necessary, is given by the blockIdx.x
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

	NodeIndex origo{};
	int n_solvents = 0;
	Coord rel_pos[MAX_SOLVENTS_IN_BLOCK];	// Pos rel to lower left forward side of block, or floor() of pos
	uint32_t ids[MAX_SOLVENTS_IN_BLOCK];
	uint8_t atomtypeIds[MAX_SOLVENTS_IN_BLOCK];

	// Not sure this is the ideal place, as it is an interim and never transferred.. 
	ForceEnergy forceEnergies[MAX_SOLVENTS_IN_BLOCK];
};






namespace BoxGrid {
	static const int blocksizeNM = 1;
	static const int64_t blocksizeLM = NANO_TO_LIMA_i * blocksizeNM;
	__device__ __host__ constexpr int NodesPerDim(int boxlenNM) { return boxlenNM / blocksizeNM; }
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
};


class SolventBlocksCircularQueue {
	static const int STEPS_PER_SOLVENTBLOCKTRANSFER = 5;	// If we go below 2, we might see issue in solventtransfers
	static const int SOLVENTBLOCK_TRANSFERSTEP = STEPS_PER_SOLVENTBLOCKTRANSFER - 1;

	// Please dont add other non-static vars to this class without checking movetodevice is not broken
	// {Step,z,y,x}
	SolventBlock* blocks = nullptr;
	int blocksInGrid = 0;
	bool has_allocated_data = false;
	bool is_on_device = false;


public:
	static const int queue_len = STEPS_PER_SOLVENTBLOCKTRANSFER;

	SolventBlocksCircularQueue() {};	// C

	__host__ static std::unique_ptr<SolventBlocksCircularQueue> createQueue(int boxlenNM) {
		auto queue = std::make_unique<SolventBlocksCircularQueue>();
		queue->allocateData(BoxGrid::NodesPerDim(boxlenNM));
		queue->initializeBlocks(boxlenNM);
		return queue;
	}

	//__host__ bool addSolventToGrid(const NodeIndex& nodeindex, const Coord& coord, uint32_t solvent_id, int64_t step, int boxSizeNM) {
	//	// TODO: Implement safety feature checking and failing if PBC is not met!
	//	return getBlockPtr(nodeindex, step, boxSizeNM)->addSolvent(coord, solvent_id);
	//}

	__host__ void allocateData(int blocksPerDim) {
		if (has_allocated_data) {
			throw std::runtime_error("BoxGridCircularQueue may not allocate data multiple times");
		}
		
		blocksInGrid = BoxGrid::BlocksTotal(blocksPerDim);
		const int blocksTotal = blocksInGrid * queue_len;
		blocks = new SolventBlock[blocksTotal]();
		has_allocated_data = true;
	}

	__host__ void initializeBlocks(int boxlenNM) {
		for (int64_t step = 0; step < queue_len; step++) {
			for (int i = 0; i < blocksInGrid; i++) {
				getBlockPtr(i, step)->origo = BoxGrid::Get3dIndex(i, boxlenNM);
			}
		}
	}

	__host__ SolventBlocksCircularQueue* CopyToDevice() const {
		SolventBlocksCircularQueue queueTemp = *this;
		const int blocksTotal = blocksInGrid * queue_len;
		queueTemp.blocks = GenericCopyToDevice(blocks, blocksTotal);

		return GenericCopyToDevice(&queueTemp, 1);
	}
	__host__ void CopyDataFromDevice(const SolventBlocksCircularQueue* const queue) const {
		const int blocksTotal = blocksInGrid * queue_len;
		cudaMemcpy(blocks, queue->blocks, sizeof(SolventBlock) * blocksTotal, cudaMemcpyDeviceToHost);
	}


	// This function assumes the user has used PBC
	__host__ SolventBlock* getBlockPtr(const NodeIndex& index3d, const int64_t step, int boxSizeNm) {
#if defined LIMASAFEMODE
		if (index3d.x >= BOXGRID_N_NODES || index3d.y >= BOXGRID_N_NODES || index3d.z >= BOXGRID_N_NODES
			|| index3d.x < 0 || index3d.y < 0 || index3d.z < 0) {
			throw std::runtime_error("Bad 3d index for blockptr\n");
		}
#endif

		return getBlockPtr(BoxGrid::Get1dIndex(index3d, boxSizeNm), step);
	}

	__device__ __host__ bool static isTransferStep(int64_t step) {
		return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == SOLVENTBLOCK_TRANSFERSTEP;
	}
	__device__ bool static isFirstStepAfterTransfer(int64_t step) {
		return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == 0;
	}

	// This function assumes the user has used PBC
	__device__ __host__ SolventBlock* getBlockPtr(const size_t index1d, const size_t step) {
		const size_t step_offset = (step % queue_len) * blocksInGrid;
		return &blocks[index1d + step_offset];
	}
};

