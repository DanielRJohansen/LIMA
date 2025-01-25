#pragma once

#include "Bodies.cuh"
#include "BoxGrid.cuh"
#include "KernelConstants.cuh"

#include <cuda_runtime.h>


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
	__device__ constexpr NodeType* GetNodePtr(NodeType* grid, const NodeIndex& index3d) { // Dont like this function, it hides using constant mem...
		//if (index3d.x >= DeviceConstants::boxSize.boxSizeNM_i || index3d.y >= DeviceConstants::boxSize.boxSizeNM_i 
		//	|| index3d.z >= DeviceConstants::boxSize.boxSizeNM_i
		//	|| index3d.x < 0 || index3d.y < 0 || index3d.z < 0) {
		//	printf("Bad 3d index for blockptr %d %d %d\n", index3d.x, index3d.y, index3d.z);
		//	return nullptr;
		//}

		return GetNodePtr<NodeType>(grid, Get1dIndex(index3d, DeviceConstants::boxSize.boxSizeNM_i));
	}

	__device__ constexpr static NodeIndex Get3dIndex(int index1d) {
		const int bpd = NodesPerDim(DeviceConstants::boxSize.boxSizeNM_i);
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
			return false; // TODO: LOsing atoms here!!!
		}
		rel_positions[n_elements] = pos;
		ids[n_elements] = id;
		atomtypeIds[n_elements] = 0;
		n_elements++;
		return true;
	}
};

struct SolventBlockTransfermodule {
	// Only use directly (full plane contact) adjacent blocks
	static const int n_queues = 6;			// or, adjecent_solvent_blocks
	static const int max_queue_size = 32;	// Maybe this is a bit dangerous

	// Each queue will be owned solely by 1 adjecent solventblock
	SolventTransferqueue<max_queue_size> transfer_queues[n_queues];

	int n_remain = 0; 
};

namespace SolventblockTransferUtils {
	__device__ constexpr NodeIndex indexToDirection[6] = {
		NodeIndex{0, 0, -1}, NodeIndex{0, 0, 1},
		NodeIndex{0, -1, 0}, NodeIndex{0, 1, 0},
		NodeIndex{-1, 0, 0}, NodeIndex{1, 0, 0}
	};
}

using STransferQueue = SolventTransferqueue<SolventBlockTransfermodule::max_queue_size>;
using SRemainQueue = SolventTransferqueue<SolventBlock::MAX_SOLVENTS_IN_BLOCK>;






namespace NlistUtil {
    static_assert(::MAX_COMPOUNDS <= UINT16_MAX, "Neighborlist cannot handle such large compound ids");

    static const int maxCompounds = 512;	// TODO: We need to work on getting this number down!
    struct IdAndRelshift {
        int id;
        Float3 relShift;
    };

    __device__ inline bool AddCompound(uint16_t new_id, const Float3& relshift, IdAndRelshift* const neighborIds, uint16_t& nNeighbors) {
        // OPTIM only check when not LIMA_PUSH
        if (nNeighbors >= maxCompounds) {
            printf("\nFailed to insert compound neighbor id %d!\n", new_id);
            return false;
            //throw std::runtime_error("Neighborlist overflow");
        }
        neighborIds[nNeighbors++] = { (int)new_id, relshift };
        return true;
    }
}




// OPTIM use alignas 128 here
class alignas(128) NeighborList {

public:


    //__device__ __host__ bool addCompound(uint16_t new_id, const Float3& relshift) {
    //	// OPTIM only check when not LIMA_PUSH
    //	if (nNonbondedNeighbors >= maxCompounds) {
    //		printf("\nFailed to insert compound neighbor id %d!\n", new_id);
    //		return false;
    //		//throw std::runtime_error("Neighborlist overflow");
    //	}
    //	nonbondedNeighborCompounds[nNonbondedNeighbors++] = { (int)new_id, relshift };
    //	return true;
    //}



    //IdAndRelshift nonbondedNeighborCompounds[maxCompounds];

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
	int nNonbondedNeighbors = 0;
};

