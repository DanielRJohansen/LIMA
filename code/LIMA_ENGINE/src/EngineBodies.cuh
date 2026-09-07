#pragma once

#include "Bodies.cuh"
#include "BoxGrid.cuh"

#include <cuda_runtime.h>
//#include <cuda_fp8.h>


// Extra BoxGrid functions used by engine kernels.
namespace BoxGrid {

	// This function assumes the user has used PBC.
	template <typename NodeType>
	__device__ constexpr NodeType* GetNodePtr(NodeType* grid, const NodeIndex& index3d, const Int3& boxSize) {
		return GetNodePtr<NodeType>(grid, Get1dIndex(index3d, boxSize));
	}
};


// Todo move this impl to ParticleClusters.cuh
class PClusterTransfermodule {
public:
	static const int maxClustersPerBlock = 64;
	static const int blockLen = 1;	
	static const int maxOutgoingClusters = 8;

	// Set by GetPclusterPositions kernel
	Float3* meanPositionOfPClustersPerBlock = nullptr; // 1 value per pCluster per block
	int* idsOfPclustersInBlocks = nullptr; // 1 value per pclusters per block
	int* nPClustersPerBlock = nullptr;		// 1 value per block

	int* nIncomingClusters = nullptr;
	int* idsOfIncomingClusters = nullptr;
	Float3* meanpositionsOfIncomingClusters = nullptr;

	__host__ static PClusterTransfermodule Create(Int3 boxSize) {
		const int nBlocksTotal = boxSize.InnerProduct();
		PClusterTransfermodule transferModule;
		cudaMalloc(&transferModule.meanPositionOfPClustersPerBlock, sizeof(Float3) * maxClustersPerBlock * nBlocksTotal);		
		cudaMalloc(&transferModule.idsOfPclustersInBlocks, sizeof(int) * maxClustersPerBlock * nBlocksTotal);
		cudaMalloc(&transferModule.nPClustersPerBlock, sizeof(int) * nBlocksTotal);

		cudaMalloc(&transferModule.nIncomingClusters, sizeof(int) * 6 * maxOutgoingClusters * nBlocksTotal);
		cudaMalloc(&transferModule.idsOfIncomingClusters, sizeof(int) * 6 * maxOutgoingClusters * nBlocksTotal);
		cudaMalloc(&transferModule.meanpositionsOfIncomingClusters, sizeof(Float3) * 6 * maxOutgoingClusters * nBlocksTotal);
		transferModule.Reset(boxSize);

		cudaMemset(transferModule.idsOfIncomingClusters, 0, sizeof(int) * 6 * maxOutgoingClusters * nBlocksTotal);
		cudaMemset(transferModule.meanpositionsOfIncomingClusters, 0, sizeof(Float3) * 6 * maxOutgoingClusters * nBlocksTotal);

		return transferModule;
	}
	__host__ void Reset(Int3 boxSize) {
		const int nBlocksTotal = boxSize.InnerProduct();
		cudaMemset(nPClustersPerBlock, 0, sizeof(int) * nBlocksTotal);
		cudaMemset(nIncomingClusters, 0, sizeof(int) * 6 * maxOutgoingClusters * nBlocksTotal);
	}
	__host__ void Free() {
		cudaFree(meanPositionOfPClustersPerBlock);
		cudaFree(nPClustersPerBlock);
		cudaFree(idsOfPclustersInBlocks);
		cudaFree(nIncomingClusters);
		cudaFree(idsOfIncomingClusters);
		cudaFree(meanpositionsOfIncomingClusters);
	}
};
	

struct SuperClustersControl {
	//int* nSuperclustersInGrid;
	//int* pclusterIdsInSuperclusters;
	static const int maxClustersPerBlock = 16;

	SuperClusterMeta* scMeta = nullptr;
	SuperCluster* scData = nullptr;

	int* scIdsInBlocks = nullptr;
	int* nSuperclustersInBlocks = nullptr;

	__host__ SuperClustersControl (Int3 boxSize, int maxSuperclusters) {
		cudaMalloc(&scMeta, sizeof(SuperClusterMeta) * maxSuperclusters);
		cudaMalloc(&scData, sizeof(SuperCluster) * maxSuperclusters);
		
		cudaMalloc(&scIdsInBlocks, sizeof(int) * maxClustersPerBlock * boxSize.InnerProduct());
		cudaMalloc(&nSuperclustersInBlocks, sizeof(int) * boxSize.InnerProduct());

		Reset(boxSize);
	}
	__host__ void Reset(Int3 boxSize/*int nSuperclustersMax*/ /*The struct does not track this number itself*/) {
		//cudaMemset(scMeta, 0, sizeof(SuperClusterMeta) * nSuperclustersMax); // doesnt matter
		//cudaMemset(scData, 0, sizeof(SuperCluster) * nSuperclustersMax);
		//cudaMemset(nSuperclustersAtomic, 0, sizeof(int));
		cudaMemset(nSuperclustersInBlocks, 0, sizeof(int) * boxSize.InnerProduct());
	}
	__host__ void Free() {
		cudaFree(scMeta);
		cudaFree(scData);

		cudaFree(scIdsInBlocks);
		cudaFree(nSuperclustersInBlocks);
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





