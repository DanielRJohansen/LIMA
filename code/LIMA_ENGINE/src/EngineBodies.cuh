#pragma once

#include "Bodies.cuh"





struct CompoundGridNode {
	__device__ __host__ bool addNearbyCompound(int16_t compound_id) {
		if (n_nearby_compounds >= max_nearby_compounds) {
			//throw std::runtime_error("Failed to add compound to CompoundGridNode\n");
			printf("Failed to add compound to CompoundGridNode\n");
			return false;
		}
		nearby_compound_ids[n_nearby_compounds++] = compound_id;
		return true;
	}

	__device__ void loadData(const CompoundGridNode& other) {
		n_nearby_compounds = other.n_nearby_compounds;
		for (int i = 0; i < n_nearby_compounds; i++) {
			nearby_compound_ids[i] = other.nearby_compound_ids[i];
		}
	}
	// Compounds that are near this specific node
	// A particle belonging to this node coord, can iterate through this list
	// to find all appropriate nearby compounds;	// This is insanely high
	static const int max_nearby_compounds = 64 + 16;
	int16_t nearby_compound_ids[max_nearby_compounds]{};	// MAX_COMPOUNDS HARD LIMIT
	int n_nearby_compounds = 0;
};


//template <typename NodeType>

// Class for signaling compound origo's and quickly searching nearby compounds using a coordinate on the grid
class CompoundGrid {
public:
	static const int blocks_total = BOXGRID_N_NODES * BOXGRID_N_NODES * BOXGRID_N_NODES;
	CompoundGridNode blocks[blocks_total];
	static const int first_step_prev = STEPS_PER_SOLVENTBLOCKTRANSFER - 1;


	// This function assumes the user has used PBC
	__host__ CompoundGridNode* getBlockPtr(const NodeIndex& index3d) {
		if (index3d.x >= BOXGRID_N_NODES || index3d.y >= BOXGRID_N_NODES || index3d.z >= BOXGRID_N_NODES
			|| index3d.x < 0 || index3d.y < 0 || index3d.z < 0) {
			throw std::runtime_error("Bad 3d index for blockptr\n");
		}
		return getBlockPtr(get1dIndex(index3d));
	}

	// This function assumes the user has used PBC
	__device__ __host__ CompoundGridNode* getBlockPtr(const int index1d) {
		return &blocks[index1d];
	}

	__device__ __host__ static int get1dIndex(const NodeIndex& index3d) {
		static const int bpd = BOXGRID_N_NODES;
		return index3d.x + index3d.y * bpd + index3d.z * bpd * bpd;
	}
	__device__ static NodeIndex get3dIndex(int index1d) {
		static const int bpd = BOXGRID_N_NODES;
		auto z = index1d / (bpd * bpd);
		index1d -= z * bpd * bpd;
		auto y = index1d / bpd;
		index1d -= y * bpd;
		auto x = index1d;
		return NodeIndex{ x, y, z };
	}
};

//class CompoundGrid : public BoxGrid<CompoundGridNode> {};







template <int size>
struct SolventTransferqueue {
	Coord rel_positions[size];
	uint32_t ids[size];
	int n_elements = 0;

	// Do NOT call on queue residing in global memory
	__device__ bool addElement(const Coord& pos, const Coord& pos_prev, uint32_t id) {
		if (n_elements >= size) {
			return false;
		}
		rel_positions[n_elements] = pos;
		ids[n_elements] = id;
		n_elements++;
		return true;
	}

	// Insert relative to thread calling.
	__device__ void fastInsert(const Coord& relpos, const int id) {
		rel_positions[threadIdx.x] = relpos;
		ids[threadIdx.x] = id;
	}
};

struct SolventBlockTransfermodule {
	// Only use directly (full plane contact) adjacent blocks
	static const int n_queues = 6;			// or, adjecent_solvent_blocks
	static const int max_queue_size = 32;	// Maybe this is a bit dangerous

	// Each queue will be owned solely by 1 adjecent solventblock
	SolventTransferqueue<max_queue_size> transfer_queues[n_queues];

	int n_remain = 0;

	/// <param name="transfer_direction">Relative to the originating block</param>
	__device__ static int getQueueIndex(const NodeIndex& transfer_direction) {
		// Fucking magic yo...
		// Only works if transfer_direction.len() == 1
		// First op leaves a space at index 3:
		//{-z, -y, -x, _ x, y, z}
		//const int tmp_index = transfer_direction.dot(NodeIndex{ 1, 2, 3 });
		int tmp_index = transfer_direction.x * 1 + transfer_direction.y * 2 + transfer_direction.z * 3;
		// Shift positive values left.
		tmp_index = tmp_index > 0 ? tmp_index + 2 : tmp_index + 3;
		return tmp_index;
	}
};

using STransferQueue = SolventTransferqueue<SolventBlockTransfermodule::max_queue_size>;
using SRemainQueue = SolventTransferqueue<SolventBlock::MAX_SOLVENTS_IN_BLOCK>;
































class NeighborList {
public:
	__device__ __host__ bool addCompound(uint16_t new_id) {
		if (n_compound_neighbors >= NEIGHBORLIST_MAX_COMPOUNDS) {
			printf("\nFailed to insert compound neighbor id %d!\n", new_id);
			return false;
			//throw std::runtime_error("Neighborlist overflow");
		}
		neighborcompound_ids[n_compound_neighbors++] = new_id;
		return true;
	}

	__device__ void loadMeta(NeighborList* nl_ptr) {	// Called from thread 0
		n_compound_neighbors = nl_ptr->n_compound_neighbors;
#ifdef ENABLE_SOLVENTS
		n_gridnodes = nl_ptr->n_gridnodes;
#endif
	}
	__device__ void loadData(NeighborList* nl_ptr) {
		//static_assert(MAX_COMPOUND_PARTICLES >= NEIGHBORLIST_MAX_COMPOUNDS, "nlist_loaddata broken: not enough threads");
		//if (threadIdx.x < n_compound_neighbors)			// DANGER Breaks when threads < mAX_COMPOUND_Ns
		//	neighborcompound_ids[threadIdx.x] = nl_ptr->neighborcompound_ids[threadIdx.x];

		static_assert(MAX_COMPOUND_PARTICLES < NEIGHBORLIST_MAX_COMPOUNDS, "No need to use a for loop then");
		for (int i = threadIdx.x; i < n_compound_neighbors; i += blockDim.x) {
			neighborcompound_ids[i] = nl_ptr->neighborcompound_ids[i];
		}

#ifdef ENABLE_SOLVENTS
		for (int i = threadIdx.x; i < n_gridnodes; i += blockDim.x) {
			gridnode_ids[i] = nl_ptr->gridnode_ids[i];
		}
#endif
	}

	// It is guaranteed that the first compoudns are bonded_compounds. How many, that is something the compound
	// itself keeps track of
	static_assert(MAX_COMPOUNDS < UINT16_MAX, "Neighborlist cannot handle such large compound ids");
	uint16_t neighborcompound_ids[NEIGHBORLIST_MAX_COMPOUNDS];
	int n_compound_neighbors = 0;

#ifdef ENABLE_SOLVENTS
	// returns false if an error occured
	__device__ __host__ bool addGridnode(uint16_t gridnode_id) {
		if (n_gridnodes >= max_gridnodes) {
			//throw std::runtime_error("No room for more nearby gridnodes"); }
			printf("No room for more nearby gridnodes");
			return false;
		}
		gridnode_ids[n_gridnodes++] = gridnode_id;
		return true;
	}

	static const int max_gridnodes = 64 + 4;	// Arbitrary value
	static_assert(SolventBlocksCircularQueue::blocks_total < UINT16_MAX, "Neighborlist cannot handle such large gridnode_ids");
	uint16_t gridnode_ids[max_gridnodes];
	int n_gridnodes = 0;
#endif

	int associated_id = -1;
};