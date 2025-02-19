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


class TinymolTransferModule {
public:

	static const int maxOutgoingBondgroups = 8;
	static const int maxOutgoingParticles = maxOutgoingBondgroups * 3;
	/*static const int maxIncomingBondgroups = 64;
	static const int maxIncomingParticles = 192;*/


	// 6 elements per block
	int* nIncomingParticles;
	int* nIncomingBondgroups;

	Coord* incomingPositions;
	uint32_t* incomingIds;
	uint8_t* incomingAtomtypeIds;
	uint8_t* incomingBondgroupIds;
	TinyMolParticleState* incomingStates;

	// TODO: This should just be an index into a variant of the bondgroup kept in constant memory
	BondgroupTinymol* incomingBondgroups; // 64 elements per block 
	int* incomingBondgroupsParticlesOffset;

	static TinymolTransferModule Create(int nBlocksTotal) {
		TinymolTransferModule transferModule;
		cudaMalloc(&transferModule.nIncomingParticles, sizeof(int) * 6 * nBlocksTotal);		
		cudaMalloc(&transferModule.nIncomingBondgroups, sizeof(int) * 6 * nBlocksTotal);
		cudaMemset(transferModule.nIncomingParticles, 0, sizeof(int) * 6 * nBlocksTotal);
		cudaMemset(transferModule.nIncomingBondgroups, 0, sizeof(int) * 6 * nBlocksTotal);

		cudaMalloc(&transferModule.incomingPositions, sizeof(Coord) * 6 * maxOutgoingParticles * nBlocksTotal);
		cudaMalloc(&transferModule.incomingIds, sizeof(uint32_t) * 6 * maxOutgoingParticles * nBlocksTotal);
		cudaMalloc(&transferModule.incomingAtomtypeIds, sizeof(uint8_t) * 6 * maxOutgoingParticles * nBlocksTotal);
		cudaMalloc(&transferModule.incomingBondgroupIds, sizeof(uint8_t) * 6 * maxOutgoingParticles * nBlocksTotal);
		cudaMalloc(&transferModule.incomingStates, sizeof(TinyMolParticleState) * 6 * maxOutgoingParticles * nBlocksTotal);

		cudaMalloc(&transferModule.incomingBondgroups, sizeof(BondgroupTinymol) * 6 * maxOutgoingBondgroups * nBlocksTotal);
		cudaMalloc(&transferModule.incomingBondgroupsParticlesOffset, sizeof(int) * 6 * maxOutgoingBondgroups  * nBlocksTotal);

		return transferModule;
	}
	void Free() const {
		cudaFree(nIncomingParticles);
		cudaFree(nIncomingBondgroups);

		cudaFree(incomingPositions);
		cudaFree(incomingIds);
		cudaFree(incomingAtomtypeIds);
		cudaFree(incomingBondgroupIds);

		cudaFree(incomingBondgroups);
		cudaFree(incomingBondgroupsParticlesOffset);
	}
};




//template <int nParticles, int nBondgroups>
//struct SolventTransferqueue {
//	Coord rel_positions[nParticles];
//	uint32_t ids[nParticles];
//	uint8_t atomtypeIds[nParticles];
//	int n_elements;
//
//	BondgroupTinymol bondgroups[nBondgroups];
//	int nBondgroups;
//
//	//// Do NOT call on queue residing in global memory
//	//__device__ bool addElement(const Coord& pos, const BondgroupTinymol& bondgroup, uint32_t id) {
//	//	if (n_elements >= size) {
//	//		return false; // TODO: LOsing atoms here!!!
//	//	}
//	//	rel_positions[n_elements] = pos;
//	//	ids[n_elements] = id;
//	//	atomtypeIds[n_elements] = 0;
//	//	n_elements++;
//	//	return true;
//	//}
//};
//
//struct SolventBlockTransfermodule {
//	// Only use directly (full plane contact) adjacent blocks
//	static const int n_queues = 6;			// or, adjecent_solvent_blocks
//	static const int max_queue_size = 32;	// Maybe this is a bit dangerous
//
//	// Each queue will be owned solely by 1 adjecent solventblock
//	SolventTransferqueue<max_queue_size> transfer_queues[n_queues];
//
//	int n_remain = 0; 
//};
//
//namespace SolventblockTransferUtils {
//	__device__ constexpr NodeIndex indexToDirection[6] = {
//		NodeIndex{0, 0, -1}, NodeIndex{0, 0, 1},
//		NodeIndex{0, -1, 0}, NodeIndex{0, 1, 0},
//		NodeIndex{-1, 0, 0}, NodeIndex{1, 0, 0}
//	};
//}
//const int MAX_BONDGROUP_TO_TRANSFER = 16;
//using STransferQueue = SolventTransferqueue<SolventBlockTransfermodule::max_queue_size, MAX_BONDGROUP_TO_TRANSFER>;
//using SRemainQueue = SolventTransferqueue<SolventBlock::MAX_SOLVENTS_IN_BLOCK, SolventBlock::maxBondgroups>;


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





