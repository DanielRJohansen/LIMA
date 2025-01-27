#pragma once

#include "Bodies.cuh"
#include "BoxGrid.cuh"
#include "KernelConstants.cuh"

#include <cuda_runtime.h>
//#include <cuda_fp8.h>


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



//class SmallShift {
//    __nv_fp8_e4m3 data[4];
//
//    SmallShift() {}
//    constexpr SmallShift(const Float3& shift) {
//        data[0] = __float2fp8(shift.x);
//        data[1] = __float2fp8(shift.y);
//        data[2] = __float2fp8(shift.z);
//    }
//
//    constexpr Float3 ToFloat3() const {
//        return Float3 { __fp82float(data[0]), __fp82float(data[1]), __fp82float(data[2]) }
//    }
//};





