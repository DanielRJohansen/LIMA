#pragma once

#include "Constants.h"
#include "LimaTypes.cuh"
#include "Bodies.cuh"

// blocks are notcentered 
struct SolventBlock {
	static const int MAX_SOLVENTS_IN_BLOCK = 64 + 32;

	__device__ __host__ SolventBlock() {};
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
		}
	}

	__host__ bool addSolvent(const Coord& rel_position, uint32_t id) {
		if (n_solvents == MAX_SOLVENTS_IN_BLOCK) {
			return false;
		}
		ids[n_solvents] = id;
		rel_pos[n_solvents++] = rel_position;
		return true;
	}

	NodeIndex origo{};
	int n_solvents = 0;
	Coord rel_pos[MAX_SOLVENTS_IN_BLOCK];	// Pos rel to lower left forward side of block, or floor() of pos
	uint32_t ids[MAX_SOLVENTS_IN_BLOCK];
};






namespace BoxGrid {
	static const int blocksizeNM = 1;
	static const int64_t blocksizeLM = NANO_TO_LIMA_i * blocksizeNM;
	__device__ __host__ constexpr int NodesPerDim(int boxlenNM) { return boxlenNM / blocksizeNM; }

	const int blocksPerDim = NodesPerDim(static_cast<int>(BOX_LEN_NM));


	__device__ __host__ static int Get1dIndex(const NodeIndex& index3d, int boxSizeNM) {
		if (boxSizeNM != 7) printf("FAIL!");
		const int bpd = NodesPerDim(boxSizeNM);
		return index3d.x + index3d.y * bpd + index3d.z * bpd * bpd;
	}
	__device__ __host__ static NodeIndex Get3dIndex(int index1d, int boxlenNM) {
		if (boxlenNM != 7) printf("FAIL!");
		const int bpd = NodesPerDim(boxlenNM);
		int z = index1d / (bpd * bpd);
		index1d -= z * bpd * bpd;
		int y = index1d / bpd;
		index1d -= y * bpd;
		int x = index1d;
		return NodeIndex{ x, y, z };
	}
};


class SolventBlocksCircularQueue {
	static const int queue_len = STEPS_PER_SOLVENTBLOCKTRANSFER;


	// Please dont add other non-static vars to this class without checking movetodevice is not broken
	// {Step,z,y,x}
	SolventBlock* blocks = nullptr;
	bool has_allocated_data = false;
	bool is_on_device = false;


public:
	static const int boxlenNmTemp = (int)BOX_LEN_NM;	// TODO: TEMP!
	static const int blocks_per_grid = BoxGrid::NodesPerDim(boxlenNmTemp) * BoxGrid::NodesPerDim(boxlenNmTemp) * BoxGrid::NodesPerDim(boxlenNmTemp);
	static const int grid_bytesize = sizeof(SolventBlock) * blocks_per_grid;
	static const int blocks_total = blocks_per_grid * queue_len;
	static const int gridqueue_bytesize = sizeof(SolventBlock) * blocks_total;
	static const int first_step_prev = queue_len - 1;

	SolventBlocksCircularQueue() {};	// C

	__host__ static SolventBlocksCircularQueue* createQueue(int boxlenNM) {
		auto queue = new SolventBlocksCircularQueue();
		queue->allocateData();
		queue->initializeBlocks(boxlenNM);
		return queue;
	}

	__host__ bool addSolventToGrid(const SolventCoord& coord, uint32_t solvent_id, int step, int boxSizeNM) {
		// TODO: Implement safety feature checking and failing if PBC is not met!
		return getBlockPtr(coord.origo, step, boxSizeNM)->addSolvent(coord.rel_position, solvent_id);
	}

	__host__ void allocateData() {
		if (has_allocated_data) {
			throw std::runtime_error("BoxGridCircularQueue may not allocate data multiple times");
		}
		blocks = new SolventBlock[blocks_total]();
		has_allocated_data = true;
	}

	__host__ void initializeBlocks(int boxlenNM) {
		for (int step = 0; step < queue_len; step++) {
			for (int i = 0; i < blocks_per_grid; i++) {
				getBlockPtr(i, step)->origo = BoxGrid::Get3dIndex(i, boxlenNM);
			}
		}
	}


	__host__ SolventBlocksCircularQueue* moveToDevice() {
		blocks = genericMoveToDevice(blocks, blocks_total);

		is_on_device = true;
		return genericMoveToDevice(this, 1);
	}

	// Fuck me this seems overcomplicated
	__host__ SolventBlocksCircularQueue* copyToHost() {
		SolventBlocksCircularQueue* this_host = new SolventBlocksCircularQueue();
		this_host->allocateData();
		cudaMemcpy(this_host->blocks, blocks, gridqueue_bytesize, cudaMemcpyDeviceToHost);
		return this_host;
	}




	// This function assumes the user has used PBC
	__host__ SolventBlock* getBlockPtr(const NodeIndex& index3d, const int step, int boxSizeNm) {
#if defined LIMASAFEMODE
		if (index3d.x >= BOXGRID_N_NODES || index3d.y >= BOXGRID_N_NODES || index3d.z >= BOXGRID_N_NODES
			|| index3d.x < 0 || index3d.y < 0 || index3d.z < 0) {
			throw std::runtime_error("Bad 3d index for blockptr\n");
		}
#endif

		return getBlockPtr(BoxGrid::Get1dIndex(index3d, boxSizeNm), step);
	}

	__device__ __host__ bool static isTransferStep(int step) {
		return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == SOLVENTBLOCK_TRANSFERSTEP;
	}
	__device__ bool static isFirstStepAfterTransfer(int step) {
		return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == 0;
	}

	// This function assumes the user has used PBC
	__device__ __host__ SolventBlock* getBlockPtr(const size_t index1d, const size_t step) {
		const size_t step_offset = (step % queue_len) * blocks_per_grid;
		return &blocks[index1d + step_offset];
	}
};

