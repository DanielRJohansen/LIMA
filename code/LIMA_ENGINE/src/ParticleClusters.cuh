#pragma once

#include "EngineBodies.cuh"
#include "Engine.cuh"
#include "DebugUtils.h"
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/scan.h>
#include "BoundaryCondition.cuh"
#include <numeric>





class SuperclusterStagingControl {
public:
	SuperclusterStagingControl() {}
	__host__ SuperclusterStagingControl(Int3 boxSize) {
		const int nBlocks = BoxGrid::BlocksTotal(boxSize);
		//const int nElements = _nBlocks + 1; 

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







// nBlocks = nPclusters/32
// blockdim = (32, 1, 1)
__global__ void GetPclusterPositions(PClusterTransfermodule transferModule, PersistentCluster* const pClustersData, const int nPclusters, Int3 boxSize) {
	const int pcId = blockIdx.x * blockDim.x + threadIdx.x;
	if (pcId >= nPclusters)
		return;

	//ParticleToCompoundOrSolventMapping mapping = mappings[particleId];

	// First ensure all particles in pcluster are same hyperpos
	for (int i = 1; i < PersistentCluster::maxParticles; i++) {
		if (pClustersData[pcId].pqd[i].Valid()) {
			PeriodicBoundaryCondition::applyHyperposNM(pClustersData[pcId].pqd[0].position, pClustersData[pcId].pqd[i].position);
		}
	}

	Float3 meanPos{};
	int count = 0;
	for (int i = 0; i < 4; i++) {
		if (pClustersData[pcId].pqd[i].Valid()) {
			meanPos += pClustersData[pcId].pqd[i].position;
			count++;
		}
	}
	meanPos = meanPos / static_cast<float>(count);

	for (int i = 0; i < count; i++) {
		if ((meanPos - pClustersData[pcId].pqd[i].position).len() > 1.f) {
			printf("meanpos %f %f %f mypos %f %f %f\n", meanPos.x, meanPos.y, meanPos.z, pClustersData[pcId].pqd[i].position.x, pClustersData[pcId].pqd[i].position.y, pClustersData[pcId].pqd[i].position.z);
		}
	}


	//const Float3 pos = pClustersData[pcId].pqd[0].position; // TODO: use actual mean pos?
	Float3 temp = meanPos;
	Float3 gridPosF = meanPos.Floor();
	PeriodicBoundaryCondition::applyBCNM(gridPosF);
	PeriodicBoundaryCondition::applyHyperposNM(gridPosF, meanPos);

	NodeIndex blockId = NodeIndex(gridPosF.x, gridPosF.y, gridPosF.z);
	if constexpr (INDEXING_CHECKS) {
		if (blockId.x < 0 || blockId.y < 0 || blockId.z < 0 || blockId.x >= boxSize.x || blockId.y >= boxSize.y || blockId.z >= boxSize.z)
			printf("Storing pc %d pos %f %f %f at block %d %d %d particleCount %d\n", pcId, meanPos.x, meanPos.y, meanPos.z, blockId.x, blockId.y, blockId.z, count);
	}
	const int blockIndex = BoxGrid::Get1dIndex(blockId, boxSize);
	int indexInBlock = atomicAdd(&transferModule.nPClustersPerBlock[blockIndex], 1);
	int index = blockIndex * PClusterTransfermodule::maxClustersPerBlock + indexInBlock;

	if constexpr (INDEXING_CHECKS){
		if (indexInBlock >= PClusterTransfermodule::maxClustersPerBlock) {
			printf("Too many pclusters in block %d %d %d. Count %d temppos %f %f %f\n", blockId.x, blockId.y, blockId.z, indexInBlock, temp.x, temp.y, temp.z);
		}
	}

	transferModule.idsOfPclustersInBlocks[index] = pcId;
	transferModule.meanPositionOfPClustersPerBlock[index] = meanPos; // TODO: THis is not even the actual meanposition is it?
}


// Since we use non-deterministic atomics to place pclusters, we need to sort the indices to make the result deterministic again
// blockDim = 32,1,1
__global__ void SortPClusterIndicesInBlocks(PClusterTransfermodule transferModule) {
	__shared__ int ids[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ Float3 positions[PClusterTransfermodule::maxClustersPerBlock];

	const int nPclustersToSort = transferModule.nPClustersPerBlock[blockIdx.x];

	if constexpr (INDEXING_CHECKS) {
		if (threadIdx.x == 0 && nPclustersToSort > PClusterTransfermodule::maxClustersPerBlock) {
			printf("Not allowed to sort %d pclusters\n", nPclustersToSort);
		}
	}

	for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i+=blockDim.x) {
		if (i < nPclustersToSort) {
			const int globalIndex = blockIdx.x * PClusterTransfermodule::maxClustersPerBlock + i;
			ids[i] = transferModule.idsOfPclustersInBlocks[globalIndex];
			positions[i] = transferModule.meanPositionOfPClustersPerBlock[globalIndex];
		}
		else {
			ids[i] = INT_MAX;
		}
	}
	__syncthreads();

	LAL::Sort<PClusterTransfermodule::maxClustersPerBlock>(ids, positions);

	for (int i = threadIdx.x; i < nPclustersToSort; i+=blockDim.x) {
		if (i < nPclustersToSort) {
			const int globalIndex = blockIdx.x * PClusterTransfermodule::maxClustersPerBlock + i;

			transferModule.idsOfPclustersInBlocks[globalIndex] = ids[i];
			transferModule.meanPositionOfPClustersPerBlock[globalIndex] = positions[i];
		}
	}
}

template <typename BoundaryCondition>
__global__ void ClusteringPretransferKernel(PClusterTransfermodule transferModule, Int3 boxSize) {	


	static const NodeIndex directions[6]{
	{1, 0, 0},
	{-1, 0, 0},
	{0, 1, 0},
	{0, -1, 0},
	{0, 0, 1},
	{0, 0, -1}
	};

	__shared__ int nPClusters;
	__shared__ int indexOfFirstCluster;
	__shared__ int directionIndexOfPCluster[PClusterTransfermodule::maxClustersPerBlock]; // -1 for stay

	if (threadIdx.x == 0) {
		nPClusters = transferModule.nPClustersPerBlock[blockIdx.x];
		indexOfFirstCluster = blockIdx.x * PClusterTransfermodule::maxClustersPerBlock;
	}
	for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i += blockDim.x) {
		directionIndexOfPCluster[i] = -1;
	}
	__syncthreads();

	// TODO FIx this part
	const Float3 blockCenter = BoxGrid::Get3dIndex(blockIdx.x, boxSize).toFloat3() + Float3{ 0.5f };

	for (int i = threadIdx.x; i < nPClusters; i += blockDim.x) {
		const int pcIndex = indexOfFirstCluster + i;
		Float3 absPos = transferModule.meanPositionOfPClustersPerBlock[pcIndex];
		PeriodicBoundaryCondition::applyHyperposNM(blockCenter, absPos);	// TODO: OPTIM: shoudn't be necessary if pClusters are placed correctly in blocks...

		

		const Float3 posRelativeToBlockCenter = absPos - blockCenter;
		const NodeIndex direction = LIMAPOSITIONSYSTEM::GetTransferDirection(posRelativeToBlockCenter);


		//if (std::abs(absPos.x - 2.23) < 0.00001 || (LAL::Fequal(blockCenter.z, 1.5) && absPos.z < 0.f)) {
		//	marked = true;
		//	printf("blockcenter %f %f %f abspos %f %f %f direction %d %d %d\n", blockCenter.x, blockCenter.y, blockCenter.z, absPos.x, absPos.y, absPos.z, direction.x, direction.y, direction.z);
		//}
		for (int directionIndex = 0; directionIndex < 6; directionIndex++) {
			if (direction == directions[directionIndex]) {
				directionIndexOfPCluster[i] = directionIndex;
				break;
			}
		}
	}
	__syncthreads();

	// First 6 threads are responsible for marking a direction
	__shared__ int nClustersThisDirection[6];
	__shared__ int clusterIdsThisDirectionRelativeToBlock[PClusterTransfermodule::maxOutgoingClusters * 6];
	if (threadIdx.x < 6) {
		int myCount = 0;
		if constexpr (INDEXING_CHECKS) {
			if (myCount >= PClusterTransfermodule::maxOutgoingClusters)
				printf("Too many pClusters in one direction");
		}
		for (int i = 0; i < nPClusters; i++) {
			if (directionIndexOfPCluster[i] == threadIdx.x) {				
				clusterIdsThisDirectionRelativeToBlock[threadIdx.x * PClusterTransfermodule::maxOutgoingClusters + myCount] = i;
				myCount++;
			}
		}	
		nClustersThisDirection[threadIdx.x] = myCount;
	}
	__syncthreads();



	// Now all threads loop over the direction, and if they have a particle, they push it directy to the incoming queue in global memory
	const NodeIndex blockOrigo = BoxGrid::Get3dIndex(blockIdx.x, boxSize);
	for (int directionIndex = 0; directionIndex < 6; directionIndex++) {
		const NodeIndex& direction = directions[directionIndex];		
		const NodeIndex targetBlock = BoundaryCondition::applyBC(blockOrigo + direction, boxSize);
		const int targetBlockId = BoxGrid::Get1dIndex(targetBlock, boxSize);
		if constexpr (INDEXING_CHECKS) {
			if (targetBlockId < 0 || targetBlockId >= BoxGrid::BlocksTotal(boxSize))
				printf("Target block %d %d %d was out of bounds. Id %d out of %d\n", targetBlock.x, targetBlock.y, targetBlock.z, targetBlockId, BoxGrid::BlocksTotal(boxSize));
		}

		//const Coord relposShift = Coord{ -direction.toFloat3() };

		
		// Write results directly to global mem
		if (threadIdx.x == 0) {
			transferModule.nIncomingClusters[targetBlockId * 6 + directionIndex] = nClustersThisDirection[directionIndex];
		}
		// Each thread takes 1 bondgroup
		if (threadIdx.x < nClustersThisDirection[directionIndex]) {
			const int clusterIndexRelativeToBlock = clusterIdsThisDirectionRelativeToBlock[directionIndex * PClusterTransfermodule::maxOutgoingClusters + threadIdx.x];
			const int clusterTargetGlobalIndex = targetBlockId * 6 * PClusterTransfermodule::maxOutgoingClusters
				+ directionIndex * PClusterTransfermodule::maxOutgoingClusters
				+ threadIdx.x;
			const int clusterSrcGlobalIndex = indexOfFirstCluster + clusterIndexRelativeToBlock;

			//printf("pcluster leaving to index %d %d %d\n", targetBlock.x, targetBlock.y, targetBlock.z);
			transferModule.idsOfIncomingClusters[clusterTargetGlobalIndex] = transferModule.idsOfPclustersInBlocks[clusterSrcGlobalIndex];
			transferModule.meanpositionsOfIncomingClusters[clusterTargetGlobalIndex] = transferModule.meanPositionOfPClustersPerBlock[clusterSrcGlobalIndex];
		}
	}
	__syncthreads();


	// Finally compress remainders. Reuse the direction-of-particle buffer as prefixsum buffer
	__shared__ int indexRelativeToBlockOfRemainingClusters[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int nClustersRemaining;	

	// TODO OPTIM: Find a not idiot way of doing this compression...
	if (threadIdx.x == 0) {
		int nextDestIndex = 0;
		for (int i = 0; i < nPClusters; i++) {
			if (directionIndexOfPCluster[i] == -1) {
				indexRelativeToBlockOfRemainingClusters[nextDestIndex++] = i;
			}
		}
		nClustersRemaining = nextDestIndex;
		transferModule.nPClustersPerBlock[blockIdx.x] = nClustersRemaining;
	}
	__syncthreads();


	// Now finally compress data in global mem
	for (int i = threadIdx.x; i < nClustersRemaining; i += blockDim.x) {
		//printf("pcluster remains block %d thread %d\n", blockIdx.x, threadIdx.x);
		const int srcIndex = indexOfFirstCluster + indexRelativeToBlockOfRemainingClusters[i];
		const int destIndex = indexOfFirstCluster + i;
		transferModule.idsOfPclustersInBlocks[destIndex] = transferModule.idsOfPclustersInBlocks[srcIndex];
		transferModule.meanPositionOfPClustersPerBlock[destIndex] = transferModule.meanPositionOfPClustersPerBlock[srcIndex];
	}
}


//uint32_t MakeKey(
//	uint32_t mask1, // first 4 bits
//	uint32_t mask2, // second 4 bits
//	float value,
//	float minValue,
//	float maxValue
//) {
//	float fraction = std::clamp((value - minValue) / (maxValue - minValue), 0.f, 1.f);
//	uint32_t fractionAsInt = fraction * (float)(1 << 24); // 24 bits for value
//	return (mask1 & 0xF) << 28 | (mask2 & 0xF) << 24 | (fractionAsInt & 0xFFFFFF);
//}

// Called with 32 threads
__global__ void ClusteringKernel(const PClusterTransfermodule transferModule, const PersistentCluster* const pClusters, SuperclusterStagingControl scStagingControl, 
	const PersistentClusterMeta* const persistentClusterMeta, Int3 boxSize)
{
	__shared__ Float3 meanPositionsOfPClusters[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int idsOfPclustersInBlock[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int nPclustersInBlock;

	__shared__ float sortKeys[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int sortIds[PClusterTransfermodule::maxClustersPerBlock]; // starts out as iota, tracks relative
	__shared__ int idsOfPclustersSorted[PClusterTransfermodule::maxClustersPerBlock];// Used for deterministic tie-breaking

	__shared__ int assignedScIds[PClusterTransfermodule::maxClustersPerBlock];
	//__shared__ int scOutStartIndex;

	const Float3 blockCenter = BoxGrid::Get3dIndex(blockIdx.x, boxSize).toFloat3() + Float3{ 0.5f };

	if (threadIdx.x == 0) {
		nPclustersInBlock = transferModule.nPClustersPerBlock[blockIdx.x];
	}

	for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i += blockDim.x) {
		sortIds[i] = i;
		meanPositionsOfPClusters[i] = Float3{ INFINITY,INFINITY, INFINITY };
		idsOfPclustersInBlock[i] = -1;
	}

	__syncthreads();

	// Load 
	for (int i = threadIdx.x; i < nPclustersInBlock; i += blockDim.x) {
		meanPositionsOfPClusters[i] = transferModule.meanPositionOfPClustersPerBlock[blockIdx.x * PClusterTransfermodule::maxClustersPerBlock + i];
		idsOfPclustersInBlock[i] = transferModule.idsOfPclustersInBlocks[blockIdx.x * PClusterTransfermodule::maxClustersPerBlock + i];
	}
	for (int dir = 0; dir < 6; dir++){
		const int incomingBaseIndex = (blockIdx.x * 6 + dir) * PClusterTransfermodule::maxOutgoingClusters;
		if (threadIdx.x < transferModule.nIncomingClusters[blockIdx.x * 6 + dir]) {
			const int indexInBlock = nPclustersInBlock + threadIdx.x;
			meanPositionsOfPClusters[indexInBlock] = transferModule.meanpositionsOfIncomingClusters[incomingBaseIndex + threadIdx.x];
			idsOfPclustersInBlock[indexInBlock] = transferModule.idsOfIncomingClusters[incomingBaseIndex + threadIdx.x];
		}
		__syncthreads();
		if (threadIdx.x == 0) {
			nPclustersInBlock += transferModule.nIncomingClusters[blockIdx.x * 6 + dir];
			if constexpr (INDEXING_CHECKS) {
				if (nPclustersInBlock > PClusterTransfermodule::maxClustersPerBlock)
					printf("Too many clusters in block after adding incoming. Block %d, count %d\n", blockIdx.x, nPclustersInBlock);
			}
		}
		__syncthreads();
	}
	// We've loaded positions from this a neighboring blocks, make sure there at the same hyperpos
	{
		for (int i = threadIdx.x; i < nPclustersInBlock; i += blockDim.x) {
			PeriodicBoundaryCondition::applyHyperposNM(blockCenter, meanPositionsOfPClusters[i]);
		}
		__syncthreads();
	}


	for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i+=blockDim.x) {
		idsOfPclustersSorted[i] = idsOfPclustersInBlock[i];
	}
	__syncthreads();





	const int bucketsPerDim = 4;
	// Sort along z 
	{
		for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i+=blockDim.x) {
			sortKeys[i] = meanPositionsOfPClusters[sortIds[i]].z;
		}
		__syncthreads();
		const int bucketSize = PClusterTransfermodule::maxClustersPerBlock;
		LAL::SortInBins<1, bucketSize>(sortKeys, idsOfPclustersSorted, sortIds);
		__syncthreads();
	}
	// Sort the 4 buckets along y
	{
		for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i += blockDim.x) {
			sortKeys[i] = meanPositionsOfPClusters[sortIds[i]].y;
		}
		__syncthreads();
		const int bucketSize = PClusterTransfermodule::maxClustersPerBlock / bucketsPerDim;		
		LAL::SortInBins<bucketsPerDim, bucketSize>(sortKeys, idsOfPclustersSorted, sortIds);
		__syncthreads();
	}
	// Sort the 4x4 buckets along x
	{
		for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i += blockDim.x) {
			sortKeys[i] = meanPositionsOfPClusters[sortIds[i]].x;
		}
		__syncthreads();
		const int bucketSize = PClusterTransfermodule::maxClustersPerBlock / (bucketsPerDim * bucketsPerDim);
		LAL::SortInBins<bucketsPerDim * bucketsPerDim, bucketSize>(sortKeys, idsOfPclustersSorted, sortIds);
		__syncthreads();
	}

	// Assign scIds
	__shared__ int nClustersToMake;
	if (threadIdx.x == 0) {
		int scId = 0;
		int scParticlesSum = 0;
		for (int i = 0; i < nPclustersInBlock; i++) {
			const int pcIdRelativeToBlock = sortIds[i];
			const int pcIdGlobal = idsOfPclustersInBlock[pcIdRelativeToBlock];
			const int nParticles = persistentClusterMeta[pcIdGlobal].nParticles;
			if (scParticlesSum + nParticles > 16) {
				scId++;
				scParticlesSum = 0;
			}
			assignedScIds[i] = scId;
			scParticlesSum += nParticles;
		}

		if (scParticlesSum > 0)
			scId++;
		nClustersToMake = scId;
		scStagingControl.nClustersPerBlock[blockIdx.x] = nClustersToMake;
		if constexpr (INDEXING_CHECKS) {
			if (nClustersToMake > SuperClustersControl::maxClustersPerBlock) {
				printf("Trying to make too many superclusters in block %d: %d (max %d)\n", blockIdx.x, nClustersToMake, SuperClustersControl::maxClustersPerBlock);
			}
		}		
	}
	__syncthreads();


	if (threadIdx.x < nClustersToMake){
		SuperClusterMeta scMeta{};
		SuperCluster sc{};


		for (int i = 0; i < nPclustersInBlock; i++) {
			if (assignedScIds[i] != threadIdx.x) {
				continue;
			}

			const int pcIdRelativeToBlock = sortIds[i];
			const int pcIdGlobal = idsOfPclustersInBlock[pcIdRelativeToBlock];
			const int nParticles = persistentClusterMeta[pcIdGlobal].nParticles;

			if (pcIdRelativeToBlock < 0 || pcIdRelativeToBlock >= PClusterTransfermodule::maxClustersPerBlock || pcIdGlobal < 0 || nParticles < 0) {
				printf("Invalid cluster id. Relative %d global %d\n", pcIdRelativeToBlock, pcIdGlobal);
			}

			for (int indexInPc = 0; indexInPc < nParticles; indexInPc++) {
				const int indexInSc = scMeta.nParticles + indexInPc;
				PData pData = pClusters[pcIdGlobal].pqd[indexInPc];
				PeriodicBoundaryCondition::applyHyperposNM(blockCenter, pData.position);
				sc.pData[indexInSc] = pData;
				scMeta._pclusterIds[indexInSc] = pcIdGlobal;
				scMeta.globalParticleIds[indexInSc] = persistentClusterMeta[pcIdGlobal].particleIdsGlobal[indexInPc];
				scMeta.indexInPcluster[indexInSc] = indexInPc;

				if (scMeta.nUniquePcIds == 0 || scMeta.uniquePclusterIds[scMeta.nUniquePcIds - 1] != pcIdGlobal) {
					scMeta.uniquePclusterIds[scMeta.nUniquePcIds] = pcIdGlobal;
					scMeta.nUniquePcIds++;
				}
			}
			scMeta.nParticles += nParticles;
		}
		for (int i = scMeta.nParticles; i < SuperCluster::maxParticles; i++) {
			sc.pData[i] = PData{};
			scMeta._pclusterIds[i] = -1;
			scMeta.globalParticleIds[i] = -1;
			scMeta.indexInPcluster[i] = -1;
		}

		const int stagingStartIndex = blockIdx.x * SuperClustersControl::maxClustersPerBlock;
		scStagingControl.scData[stagingStartIndex + threadIdx.x] = sc;
		scStagingControl.scMeta[stagingStartIndex + threadIdx.x] = scMeta;
	}
	__syncthreads();
}

// 1 cudablock per block, blockdim = 32,1,1
__global__ void CompressSuperclusters(SuperClustersControl scControl, const SuperclusterStagingControl scStagingControl) {
	const int srcBlockId = blockIdx.x;
	const int nClusters = scStagingControl.nClustersPerBlock[srcBlockId];

	const int srcStartIndex = srcBlockId * SuperClustersControl::maxClustersPerBlock;
	const int dstStartIndex = scStagingControl.nClustersPrefixSum[srcBlockId];

	if (nClusters > 0) {
		auto tb = cooperative_groups::this_thread_block();
		cooperative_groups::memcpy_async(tb, &scControl.scData[dstStartIndex], &scStagingControl.scData[srcStartIndex], sizeof(SuperCluster) * nClusters);
		cooperative_groups::memcpy_async(tb, &scControl.scMeta[dstStartIndex], &scStagingControl.scMeta[srcStartIndex], sizeof(SuperClusterMeta) * nClusters);
		cooperative_groups::wait(tb);
	}

	// The dstIndex is implicitly the SC id. Now assign that to the blocks for taskbuilding
	for (int indexInBlock = threadIdx.x; indexInBlock < nClusters; indexInBlock += blockDim.x) {
		const int scIdGlobal = dstStartIndex + indexInBlock;
		const int targetIndex = srcBlockId * SuperClustersControl::maxClustersPerBlock + indexInBlock;
		scControl.scIdsInBlocks[targetIndex] = scIdGlobal;
	}
	
	if (threadIdx.x == 0) {
		scControl.nSuperclustersInBlocks[srcBlockId] = nClusters;		
	}
}


//#pragma optimize("", off)
//
//struct PosAndId {
//	Float3 pos;
//	int id;
//};
//
//std::vector<PosAndId> Sort(std::vector<Float3> unsorted, std::vector<int> ids) {
//	assert(unsorted.size() == 64);
//	assert(ids.size() == 64);
//
//	std::vector<PosAndId> combined(unsorted.size());
//	for (size_t i = 0; i < unsorted.size(); ++i) {
//		combined[i] = { unsorted[i], ids[i] };
//	}
//
//	// Sort all along z
//	std::sort(combined.begin(), combined.end(), [](const PosAndId& a, const PosAndId& b) {
//		return a.pos.z < b.pos.z || (a.pos.z == b.pos.z && a.id < b.id);
//	});
//	// Now sort y in batches of 16
//	for (int i = 0; i < combined.size(); i += 16) {
//		std::sort(combined.begin() + i, combined.begin() + i + 16, [](const PosAndId& a, const PosAndId& b) {
//			return a.pos.y < b.pos.y || (a.pos.y == b.pos.y && a.id < b.id);
//		});
//	}
//	// Now sort x in batches of 4
//	for (int i = 0; i < combined.size(); i += 4) {
//		std::sort(combined.begin() + i, combined.begin() + i + 4, [](const PosAndId& a, const PosAndId& b) {
//			return a.pos.x < b.pos.x || (a.pos.x == b.pos.x && a.id < b.id);
//		});
//	}
//	return combined;
//}
//
//
//void SortAndCompare(std::vector<Float3> unsorted, std::vector<Float3> sorted, std::vector<int> idsUnsorted) {
//	int n = unsorted.size() / 64;
//
//	for (auto& id : idsUnsorted)
//		if (id == -1)
//			id = INT_MAX;
//
//	for (int i = 0; i < n; i++) {
//		std::vector<Float3> unsortedBatch(unsorted.begin() + i * 64, unsorted.begin() + (i + 1) * 64);
//		std::vector<Float3> sortedBatch(sorted.begin() + i * 64, sorted.begin() + (i + 1) * 64);
//		std::vector<int> idsUnsortedBatch(idsUnsorted.begin() + i * 64, idsUnsorted.begin() + (i + 1) * 64);
//		std::vector<PosAndId> sortedBatchCpu = Sort(unsortedBatch, idsUnsortedBatch);
//
//
//
//		for (int j = 0; j < 64; j++) {
//			if (sortedBatch[j] != sortedBatchCpu[j].pos) {
//				printf("Mismatch at batch %d index %d: GPU %f %f %f, CPU %f %f %f\n", i, j, sortedBatch[j].x, sortedBatch[j].y, sortedBatch[j].z, sortedBatchCpu[j].pos.x, sortedBatchCpu[j].pos.y, sortedBatchCpu[j].pos.z);
//			}
//		}
//	}
//}
//#pragma optimize("", on)


void Engine::RunClustering(bool getPclusters) {
	Int3 boxSize = simulation->box->boxparams.boxSize;
	const int nBlocks = BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxSize));
	const int nPclusters = simulation->box->persistentClusters.size();

	if (!superclusterStagingControl) {
		superclusterStagingControl = std::make_unique<SuperclusterStagingControl>(boxSize);
	}
	
	//auto pClusters = GenericCopyToHost(pClusterDevice, nPclusters);
	//DebugUtils::VerifyIdentical(pClusters, "PClustersBeforeClustering" + std::to_string(simulation->getStep()));
	//DebugUtils::VerifyIdentical(pClusterDevice, nPclusters, "PClustersBeforeClustering", simulation->getStep());

	if (getPclusters) {
		
		int nCudablocks = (nPclusters + 31) / 32;
		GetPclusterPositions<<<nCudablocks, 32>>>(
			*pclusterTransfermodule,
			pClusterDevice,
			nPclusters,
			simulation->box->boxparams.boxSize);
		LIMA_UTILS::genericErrorCheckNoSync("Error after GetPclusterPositions kernel");	

		//DebugUtils::VerifyIdentical(pclusterTransfermodule->meanPositionOfPClustersPerBlock)

		auto nClusters = GenericCopyToHost(pclusterTransfermodule->nPClustersPerBlock, nBlocks);

		SortPClusterIndicesInBlocks <<<nBlocks, 32>>> (*pclusterTransfermodule);
		LIMA_UTILS::genericErrorCheckNoSync("Error after SortPClusterIndicesInBlocks kernel");
	}

	ClusteringPretransferKernel<PeriodicBoundaryCondition>
		<<<nBlocks, 32 >>> (*pclusterTransfermodule, boxSize);
	LIMA_UTILS::genericErrorCheckNoSync("Error after ClusteringPretransferKernel");

	//DebugUtils::VerifyIdentical(pclusterTransfermodule->idsOfIncomingClusters, 6 * PClusterTransfermodule::maxOutgoingClusters * nBlocks, "idsOfIncomingCLusters", simulation->getStep());
	//DebugUtils::VerifyIdentical(pclusterTransfermodule->meanpositionsOfIncomingClusters, 6 * PClusterTransfermodule::maxOutgoingClusters * nBlocks, "meanpositionsOfIncomingCLusters", simulation->getStep());
	
	

	// This simply stages the SC's per block, need to compress after
	ClusteringKernel <<<nBlocks, 32>>> (*pclusterTransfermodule, pClusterDevice, *superclusterStagingControl, pClusterMetaDevice, boxSize);
	LIMA_UTILS::genericErrorCheckNoSync("Error after ClusteringKernel");
	//DebugUtils::VerifyIdentical(superclusterStagingControl->scData, nBlocks * SuperClustersControl::maxClustersPerBlock, "RunClustering_SCData", simulation->getStep());

	/*SortAndCompare(
		GenericCopyToHost(positionsUnsorted, nBlocks * PClusterTransfermodule::maxClustersPerBlock), 
		GenericCopyToHost(positionsSorted, nBlocks * PClusterTransfermodule::maxClustersPerBlock),
		GenericCopyToHost(idsUnsorted, nBlocks * PClusterTransfermodule::maxClustersPerBlock)
	);*/

	// Compute prefixsum buffer on nScPerBlock
	{
		const int nElements = nBlocks + 1; // Extra element for total sum at end
		//DebugUtils::VerifyIdentical(superclusterStagingControl->nClustersPerBlock, nElements, "nClustersPerBlock", simulation->getStep());
		thrust::exclusive_scan(thrust::device, superclusterStagingControl->nClustersPerBlock, superclusterStagingControl->nClustersPerBlock + nElements, superclusterStagingControl->nClustersPrefixSum);
		LIMA_UTILS::genericErrorCheckNoSync("Error after Prefixsum");
		nSuperclusters = GenericCopyToHost<int>(superclusterStagingControl->nClustersPrefixSum + nElements-1);
	}

	CompressSuperclusters<<<nBlocks, 32>>>(*superClustersControl, *superclusterStagingControl);

	pclusterTransfermodule->Reset(boxSize);


	//DebugUtils::VerifyIdentical(superClustersControl->scData, nSuperclusters, "RunClustering_SCData_Compressed", simulation->getStep());
}

void Engine::BootstrapClustering() {
	Int3 boxSize = simulation->box->boxparams.boxSize;
	const int nBlocks = BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxSize));;
	const Float3 boxSizeF = Float3(boxSize.x, boxSize.y, boxSize.z);
	const Box& box = *simulation->box;
	
	std::vector<Float3> meanPositionsofPclusters(nBlocks * PClusterTransfermodule::maxClustersPerBlock);	
	std::vector<int> idsOfPclustersInBlocks(nBlocks * PClusterTransfermodule::maxClustersPerBlock);
	std::vector<int> nPclustersPerBlock(nBlocks);

	for (int pcId = 0; pcId < box.persistentClusters.size(); pcId++) {

		Float3 pos = box.persistentClusters[pcId].pqd[0].position;
		NodeIndex targetBlock{ static_cast<int>(floorf(pos.x)), static_cast<int>(floorf(pos.y)), static_cast<int>(floorf(pos.z)) };
		BoundaryConditionPublic::applyBC(targetBlock, boxSize, BoundaryConditionSelect::PBC);
		BoundaryConditionPublic::applyHyperposNM(targetBlock.toFloat3(), pos, boxSizeF, BoundaryConditionSelect::PBC);

		int blockIndex = BoxGrid::Get1dIndex(targetBlock, BoxGrid::NodesPerDim(boxSize));

		int targetDataIndex = blockIndex * PClusterTransfermodule::maxClustersPerBlock + nPclustersPerBlock[blockIndex];
		meanPositionsofPclusters[targetDataIndex] = pos;
		idsOfPclustersInBlocks[targetDataIndex] = pcId;
		nPclustersPerBlock[blockIndex]++;
	}

	cudaMemcpy(pclusterTransfermodule->meanPositionOfPClustersPerBlock, meanPositionsofPclusters.data(), meanPositionsofPclusters.size() * sizeof(Float3), cudaMemcpyHostToDevice);
	cudaMemcpy(pclusterTransfermodule->idsOfPclustersInBlocks, idsOfPclustersInBlocks.data(), idsOfPclustersInBlocks.size() * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(pclusterTransfermodule->nPClustersPerBlock, nPclustersPerBlock.data(), nPclustersPerBlock.size() * sizeof(int), cudaMemcpyHostToDevice);

	LIMA_UTILS::genericErrorCheckNoSync("Error after uploading pCluster bootstrap data");

	RunClustering(false);	
}

