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





