#pragma once

#include "Bodies.cuh"
#include "BoxGrid.cuh"
#include "KernelConstants.cuh"

#include <cuda_runtime.h>


// Half the interaction between two distant set of particles
struct ForceAndPotential {
	Float3 forcePart;		// [1/C J/mol/nm] // force = myCharge * forcePart
	float potentialPart;	// [1/C J/mol] // potential = myCharge * potential
};

struct CompoundGridNode {
	__device__ __host__ bool addNearbyCompound(uint16_t compound_id) {
		if (n_nearby_compounds >= max_nearby_compounds) {
			printf("Failed to add compound to CompoundGridNode\n");
			return false;
		}
		compoundidsWithinLjCutoff[n_nearby_compounds++] = compound_id;
		return true;
	}

	__device__ void loadData(const CompoundGridNode& other) {
		n_nearby_compounds = other.n_nearby_compounds;
		for (int i = 0; i < n_nearby_compounds; i++) {
			compoundidsWithinLjCutoff[i] = other.compoundidsWithinLjCutoff[i];
		}
	}
	// Compounds that are near this specific node
	// A particle belonging to this node coord, can iterate through this list
	// to find all appropriate nearby compounds;	// This is insanely high
	static_assert(MAX_COMPOUNDS <= UINT16_MAX-1, "CompoundGridNode cannot handle such large compound ids");
	static const int max_nearby_compounds = 128;
	uint16_t compoundidsWithinLjCutoff[max_nearby_compounds]{};
	uint16_t compoundidsWithinShortRangeESCutoff[max_nearby_compounds]{};
	int n_nearby_compounds = 0;
};

// Extra functions that require access to kernel constants
namespace BoxGrid {

	// This function assumes the user has used PBC
	template <typename NodeType>
	__device__ NodeType* GetNodePtr(NodeType* grid, const NodeIndex& index3d) { // Dont like this function, it hides using constant mem...
		//if (index3d.x >= boxSize_device.boxSizeNM_i || index3d.y >= boxSize_device.boxSizeNM_i 
		//	|| index3d.z >= boxSize_device.boxSizeNM_i
		//	|| index3d.x < 0 || index3d.y < 0 || index3d.z < 0) {
		//	printf("Bad 3d index for blockptr %d %d %d\n", index3d.x, index3d.y, index3d.z);
		//	return nullptr;
		//}

		return GetNodePtr<NodeType>(grid, Get1dIndex(index3d, boxSize_device.boxSizeNM_i));
	}

	__device__ static NodeIndex Get3dIndex(int index1d) {
		const int bpd = NodesPerDim(boxSize_device.boxSizeNM_i);
		int z = index1d / (bpd * bpd);
		index1d -= z * bpd * bpd;
		int y = index1d / bpd;
		index1d -= y * bpd;
		int x = index1d;
		return NodeIndex{ x, y, z };
	}
};




template <int size>
struct SolventTransferqueue {
	Coord rel_positions[size];
	uint32_t ids[size];
	uint8_t atomtypeIds[size];
	int n_elements = 0;

	// Do NOT call on queue residing in global memory
	__device__ bool addElement(const Coord& pos, const Coord& pos_prev, uint32_t id) {
		if (n_elements >= size) {
			return false;
		}
		rel_positions[n_elements] = pos;
		ids[n_elements] = id;
		atomtypeIds[n_elements] = 0;
		n_elements++;
		return true;
	}

	// Insert relative to thread calling.
	__device__ void fastInsert(const Coord& relpos, const int id, uint8_t atomtypeId) {
		rel_positions[threadIdx.x] = relpos;
		ids[threadIdx.x] = id;
		atomtypeIds[threadIdx.x] = atomtypeId;
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
		if (nNonbondedNeighbors >= NEIGHBORLIST_MAX_COMPOUNDS) {
			printf("\nFailed to insert compound neighbor id %d!\n", new_id);
			return false;
			//throw std::runtime_error("Neighborlist overflow");
		}
		nonbondedNeighborcompoundIds[nNonbondedNeighbors++] = new_id;
		return true;
	}

	static_assert(MAX_COMPOUNDS <= UINT16_MAX, "Neighborlist cannot handle such large compound ids");
	uint16_t nonbondedNeighborcompoundIds[NEIGHBORLIST_MAX_COMPOUNDS];
	int nNonbondedNeighbors = 0;

#ifdef ENABLE_SOLVENTS
	// returns false if an error occured
	__device__ __host__ bool addGridnode(int gridnode_id) {
		if (n_gridnodes >= max_gridnodes) {
			printf("No room for more nearby gridnodes\n");
			return false;
		}
		gridnode_ids[n_gridnodes++] = gridnode_id;
		return true;
	}

	static const int max_gridnodes = 128;	// Arbitrary value
	int gridnode_ids[max_gridnodes];
	int n_gridnodes = 0;
#endif

	int associated_id = -1;
};