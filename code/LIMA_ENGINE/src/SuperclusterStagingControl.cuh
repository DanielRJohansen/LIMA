#pragma once

#include "EngineBodies.cuh"

#include <cuda_runtime.h>

class SuperclusterStagingControl {
public:
	SuperclusterStagingControl() {}
	__host__ SuperclusterStagingControl(int nBlocks) {
		//const int nElements = _nBlocks + 1; 
		//printf("Bytesize %f MB\n", static_cast<float>(byteSize) / 1024.f / 1024.f);

		cudaMalloc(&nClustersPerBlock, sizeof(int) * (nBlocks + 1)); // 1 extra element allows is to see the sum at the final prefixsum index
		cudaMalloc(&nClustersPrefixSum, sizeof(int) * (nBlocks + 1));
		cudaMalloc(&scData, sizeof(SuperCluster) * nBlocks * SuperClustersControl::maxClustersPerBlock);
		cudaMalloc(&scMeta, sizeof(SuperClusterMeta) * nBlocks * SuperClustersControl::maxClustersPerBlock);

		cudaMemset(nClustersPerBlock, 0, sizeof(int) * (nBlocks + 1));
		cudaMemset(nClustersPrefixSum, 0, sizeof(int) * (nBlocks + 1));
	}

	__host__ void Free() {
		if (nClustersPerBlock == nullptr)
			return;

		cudaFree(nClustersPerBlock);
		cudaFree(nClustersPrefixSum);
		cudaFree(scData);
		cudaFree(scMeta);
	}

	int* nClustersPerBlock = nullptr;
	int* nClustersPrefixSum = nullptr;
	SuperCluster* scData = nullptr;
	SuperClusterMeta* scMeta = nullptr;
};
