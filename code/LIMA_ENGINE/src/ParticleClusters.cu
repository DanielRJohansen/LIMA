// The particle clustering kernels and Engine::RunClustering. A translation unit of its own, compiled in parallel with Engine.cu

#include "SuperclusterStagingControl.cuh"

#include "EngineBodies.cuh"
#include "Engine.cuh"
#include "SimulationData.h"
#include "CubWrappers.h"
#include "BoundaryCondition.cuh"
#include "DeviceAlgorithms.cuh"
#include "Utilities.h"
#include "LimitTesting.cuh"
#include "PME.cuh"
#include "LimaPositionSystem.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/memcpy_async.h>
#include <numeric>
#include <limits>

// blockDim = (32, 1, 1)
__global__ void ApplyBoundaryCondition(PersistentCluster* const pClusters, int nPclusters, Float3 boxSize, Float3 boxSizeInv) {
	const int workId = blockIdx.x * blockDim.x + threadIdx.x;
	if (workId >= nPclusters) return;
	const int pcId = workId;


	PersistentCluster pcluster = pClusters[pcId];
	PeriodicBoundaryCondition::ApplyBC(pcluster.pqd[0].position, boxSize, boxSizeInv); // TODO: Templated BC

	for (int i = 1; i < PersistentCluster::maxParticles; i++) {
		PeriodicBoundaryCondition::ApplyHyperpos(pcluster.pqd[0].position, pcluster.pqd[i].position, boxSize, boxSizeInv);
	}

	pClusters[pcId] = pcluster;
}

// nBlocks = nPclusters/32
// blockdim = (32, 1, 1)
__global__ void GetPclusterPositions(PClusterTransfermodule transferModule, PersistentCluster* const pClustersData, int pclusterOffset, const int nPclusters, int gridOffset, Int3 boxSize, const Float3 boxSizeFloat, const Float3 boxSizeFloatInv) {
	const int workId = blockIdx.x * blockDim.x + threadIdx.x;
	if (workId >= nPclusters) return;
	const int pcId = pclusterOffset + workId;


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
			printf("meanpos %f %f %f mypos %f %f %f\n", meanPos.x, meanPos.y, meanPos.z, pClustersData[pcId].pqd[i].position.x, pClustersData[pcId].pqd[i].position.y, pClustersData[pcId].pqd[i].position.z);	// TODO: Put this in constexper
		}
	}


	//const Float3 pos = pClustersData[pcId].pqd[0].position; // TODO: use actual mean pos?
	Float3 gridPosF = meanPos.Floor();
	PeriodicBoundaryCondition::ApplyBC(gridPosF, boxSizeFloat, boxSizeFloatInv);
	PeriodicBoundaryCondition::ApplyHyperpos(gridPosF, meanPos, boxSizeFloat, boxSizeFloatInv);

	NodeIndex blockId = NodeIndex(gridPosF.x, gridPosF.y, gridPosF.z);
	if constexpr (INDEXING_CHECKS) {
		if (blockId.x < 0 || blockId.y < 0 || blockId.z < 0 || blockId.x >= boxSize.x || blockId.y >= boxSize.y || blockId.z >= boxSize.z)
			printf("Storing pc %d pos %f %f %f at block %d %d %d particleCount %d\n", pcId, meanPos.x, meanPos.y, meanPos.z, blockId.x, blockId.y, blockId.z, count);
	}
	const int blockIndex = gridOffset + BoxGrid::Get1dIndex(blockId, boxSize);
	int indexInBlock = atomicAdd(&transferModule.nPClustersPerBlock[blockIndex], 1);
	if (indexInBlock >= PClusterTransfermodule::maxClustersPerBlock) {
		// The counter keeps growing past capacity; every reader clamps it
		transferModule.overflow.Report(CapacityOverflow::ClusterOccupancy, blockIndex, indexInBlock + 1, PClusterTransfermodule::maxClustersPerBlock);
		return;
	}
	int index = blockIndex * PClusterTransfermodule::maxClustersPerBlock + indexInBlock;

	transferModule.idsOfPclustersInBlocks[index] = pcId;
	transferModule.meanPositionOfPClustersPerBlock[index] = meanPos; // TODO: THis is not even the actual meanposition is it?
}


// Since we use non-deterministic atomics to place pclusters, we need to sort the indices to make the result deterministic again
// blockDim = 32,1,1
__global__ void SortPClusterIndicesInBlocks(PClusterTransfermodule transferModule) {
	__shared__ int ids[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ Float3 positions[PClusterTransfermodule::maxClustersPerBlock];

	const int nPclustersToSort = min(transferModule.nPClustersPerBlock[blockIdx.x], PClusterTransfermodule::maxClustersPerBlock);

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
__global__ void ClusteringPretransferKernel(PClusterTransfermodule transferModule, Int3 boxSize, Float3 boxSizeFloat) {


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
		nPClusters = min(transferModule.nPClustersPerBlock[blockIdx.x], PClusterTransfermodule::maxClustersPerBlock);
		indexOfFirstCluster = blockIdx.x * PClusterTransfermodule::maxClustersPerBlock;
	}
	for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i += blockDim.x) {
		directionIndexOfPCluster[i] = -1;
	}
	__syncthreads();

	// TODO FIx this part
	const Float3 blockCenter = BoxGrid::Get3dIndex(blockIdx.x % boxSize.InnerProduct(), boxSize).toFloat3() + Float3{ 0.5f };

	for (int i = threadIdx.x; i < nPClusters; i += blockDim.x) {
		const int pcIndex = indexOfFirstCluster + i;
		Float3 absPos = transferModule.meanPositionOfPClustersPerBlock[pcIndex];
		PeriodicBoundaryCondition::applyHyperposNM(blockCenter, absPos, boxSizeFloat);	// TODO: OPTIM: shoudn't be necessary if pClusters are placed correctly in blocks...

		

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
		for (int i = 0; i < nPClusters; i++) {
			if (directionIndexOfPCluster[i] == threadIdx.x) {
				if (myCount < PClusterTransfermodule::maxOutgoingClusters)
					clusterIdsThisDirectionRelativeToBlock[threadIdx.x * PClusterTransfermodule::maxOutgoingClusters + myCount] = i;
				myCount++;
			}
		}
		if (myCount > PClusterTransfermodule::maxOutgoingClusters) {
			transferModule.overflow.Report(CapacityOverflow::ClusterTransfer, blockIdx.x, myCount, PClusterTransfermodule::maxOutgoingClusters);
			myCount = PClusterTransfermodule::maxOutgoingClusters;
		}
		nClustersThisDirection[threadIdx.x] = myCount;
	}
	__syncthreads();

	// Now all threads loop over the direction, and if they have a particle, they push it directy to the incoming queue in global memory
	const NodeIndex blockOrigo = BoxGrid::Get3dIndex(blockIdx.x % boxSize.InnerProduct(), boxSize);
	for (int directionIndex = 0; directionIndex < 6; directionIndex++) {
		const NodeIndex& direction = directions[directionIndex];		
		const NodeIndex targetBlock = BoundaryCondition::applyBC(blockOrigo + direction, boxSize);
		const int targetBlockId = (blockIdx.x / boxSize.InnerProduct()) * boxSize.InnerProduct() + BoxGrid::Get1dIndex(targetBlock, boxSize);
		if constexpr (INDEXING_CHECKS) {
			if (targetBlockId < 0 || targetBlockId % boxSize.InnerProduct() >= BoxGrid::BlocksTotal(boxSize))
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
	const PersistentClusterMeta* const persistentClusterMeta, Int3 boxSize, Float3 boxSizeFloat)
{
	__shared__ Float3 meanPositionsOfPClusters[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int idsOfPclustersInBlock[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int nPclustersInBlock;

	__shared__ float sortKeys[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int sortIds[PClusterTransfermodule::maxClustersPerBlock]; // starts out as iota, tracks relative
	__shared__ int idsOfPclustersSorted[PClusterTransfermodule::maxClustersPerBlock];// Used for deterministic tie-breaking

	__shared__ int assignedScIds[PClusterTransfermodule::maxClustersPerBlock];
	//__shared__ int scOutStartIndex;

	const Float3 blockCenter = BoxGrid::Get3dIndex(blockIdx.x % boxSize.InnerProduct(), boxSize).toFloat3() + Float3{ 0.5f };

	if (threadIdx.x == 0) {
		nPclustersInBlock = min(transferModule.nPClustersPerBlock[blockIdx.x], PClusterTransfermodule::maxClustersPerBlock);
	}

	for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i += blockDim.x) {
		sortIds[i] = i;
		meanPositionsOfPClusters[i] = Float3{ std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity() };
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
		const int nIncoming = transferModule.nIncomingClusters[blockIdx.x * 6 + dir];
		if (threadIdx.x < nIncoming && nPclustersInBlock + threadIdx.x < PClusterTransfermodule::maxClustersPerBlock) {
			const int indexInBlock = nPclustersInBlock + threadIdx.x;
			meanPositionsOfPClusters[indexInBlock] = transferModule.meanpositionsOfIncomingClusters[incomingBaseIndex + threadIdx.x];
			idsOfPclustersInBlock[indexInBlock] = transferModule.idsOfIncomingClusters[incomingBaseIndex + threadIdx.x];
		}
		__syncthreads();
		if (threadIdx.x == 0) {
			nPclustersInBlock += nIncoming;
			if (nPclustersInBlock > PClusterTransfermodule::maxClustersPerBlock) {
				transferModule.overflow.Report(CapacityOverflow::ClusterOccupancy, blockIdx.x, nPclustersInBlock, PClusterTransfermodule::maxClustersPerBlock);
				nPclustersInBlock = PClusterTransfermodule::maxClustersPerBlock;
			}
		}
		__syncthreads();
	}
	// We've loaded positions from this a neighboring blocks, make sure there at the same hyperpos
	{
		for (int i = threadIdx.x; i < nPclustersInBlock; i += blockDim.x) {
			PeriodicBoundaryCondition::applyHyperposNM(blockCenter, meanPositionsOfPClusters[i], boxSizeFloat);
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
		scMeta.simulationId = blockIdx.x / boxSize.InnerProduct();
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
				PeriodicBoundaryCondition::applyHyperposNM(blockCenter, pData.position, boxSizeFloat);
				sc.SetPdata(pData, indexInSc);
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
			sc.SetPdata(PData{}, i);
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


__global__ void SetSuperclusterSimulationId(SuperClusterMeta* metadata, int firstSupercluster, int count, int simulationId) {
	const int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index < count) metadata[firstSupercluster + index].simulationId = simulationId;
}

// blockDim = (16, 4, 1)
__global__ void MapPclusterParticlesToSlots(const SuperClusterMeta* const scMeta, int nSuperclusters, int* const pclusterParticleSlots) {
	const int scId = blockIdx.x * blockDim.y + threadIdx.y;
	if (scId >= nSuperclusters)
		return;
	const int pcId = scMeta[scId]._pclusterIds[threadIdx.x];
	if (pcId != -1)
		pclusterParticleSlots[pcId * PersistentCluster::maxParticles + scMeta[scId].indexInPcluster[threadIdx.x]] = scId * SuperCluster::maxParticles + threadIdx.x;
}

// blockDim = (16, 4, 1)
__global__ void LoadIntegrationStatesKernel(const SuperClusterMeta* const scMeta, int nSuperclusters, const PersistentclusterInterimState* const pcStates,
	const PersistentClusterMeta* const pcMeta, ParticleIntegrationState* const states) {
	const int scId = blockIdx.x * blockDim.y + threadIdx.y;
	if (scId >= nSuperclusters)
		return;
	const int pcId = scMeta[scId]._pclusterIds[threadIdx.x];
	const int index = scMeta[scId].indexInPcluster[threadIdx.x];
	states[scId * SuperCluster::maxParticles + threadIdx.x] = pcId == -1 ? ParticleIntegrationState{} : ParticleIntegrationState{
		pcStates[pcId].vels_prev[index], pcMeta[pcId].mass[index], pcStates[pcId].forces_prev[index], pcId * PersistentCluster::maxParticles + index };
}

// forcesMagnitudeSquared may be nullptr. blockDim = (16, 4, 1)
__global__ void StoreIntegrationStatesKernel(const SuperCluster* const superClusters, const SuperClusterMeta* const scMeta, int nSuperclusters,
	const ParticleIntegrationState* const states, PersistentCluster* const pclusters, PersistentclusterInterimState* const pcStates, float* const forcesMagnitudeSquared) {
	const int scId = blockIdx.x * blockDim.y + threadIdx.y;
	if (scId >= nSuperclusters)
		return;
	const ParticleIntegrationState state = states[scId * SuperCluster::maxParticles + threadIdx.x];
	if (state.mass == 0.f) // Unused slot
		return;
	const int pcId = state.pclusterParticle / PersistentCluster::maxParticles;
	const int index = state.pclusterParticle % PersistentCluster::maxParticles;
	pclusters[pcId].pqd[index].position = superClusters[scId].Position(threadIdx.x);
	pcStates[pcId].vels_prev[index] = state.velocity;
	pcStates[pcId].forces_prev[index] = state.forcePrev;
	if (forcesMagnitudeSquared)
		forcesMagnitudeSquared[scMeta[scId].globalParticleIds[threadIdx.x]] = state.forcePrev.lenSquared();
}

// Copies the positions and MD integration states from the supercluster slots back to the pclusters, and the force magnitudes to
// forcesMagnitudeSquareDevice. Anything reading the pclusters' positions or states must call this first
void Engine::StoreIntegrationStates(cudaStream_t stream) {
	if (!batch->integrationStatesLoaded || batch->nSuperclusters == 0)
		return;
	StoreIntegrationStatesKernel<<<(batch->nSuperclusters + 3) / 4, dim3(SuperCluster::maxParticles, 4, 1), 0, stream>>>(
		batch->superClustersControl->scData, batch->superClustersControl->scMeta, batch->nSuperclusters, batch->integrationStates.Get(),
		batch->pClusterDevice.Get(), batch->boxState.pclusterInterimStates,
		batch->forceMagnitudesInStates ? batch->forcesMagnitudeSquareDevice.Get() : nullptr);
	batch->forceMagnitudesInStates = false;
}

void Engine::LoadIntegrationStates(cudaStream_t stream) {
	batch->integrationStates.Expand(size_t(batch->nSuperclusters) * SuperCluster::maxParticles, 1.2);
	if (batch->nSuperclusters > 0)
		LoadIntegrationStatesKernel<<<(batch->nSuperclusters + 3) / 4, dim3(SuperCluster::maxParticles, 4, 1), 0, stream>>>(
			batch->superClustersControl->scMeta, batch->nSuperclusters, batch->boxState.pclusterInterimStates, batch->pClusterMetaDevice.Get(),
			batch->integrationStates.Get());
	batch->integrationStatesLoaded = true;
}

__global__ void MapBondgroupParticlesToSlots(const BondGroup::ParticleRef* const particles, int nParticles, const int* const pclusterParticleSlots,
	int* const bondgroupParticleSlots) {
	const int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index < nParticles)
		bondgroupParticleSlots[index] = pclusterParticleSlots[particles[index].pcid * PersistentCluster::maxParticles + particles[index].pid];
}

void Engine::RunClustering(cudaStream_t stream, bool getPclusters) {
	const Int3 boxSize = batch->boxSize;
	const int nBlocks = batch->nGridnodes;
	const int nPclusters = batch->nPclusters;
	const Float3 boxSizeF = Float3{ boxSize.x, boxSize.y, boxSize.z };
	const Float3 boxSizeFInv = boxSizeF.Inv();


	if (!batch->superclusterStagingControl) {
		batch->superclusterStagingControl = std::make_unique<SuperclusterStagingControl>(batch->nGridnodes);
	}

	// The integration states move with the particles to their new supercluster slots
	StoreIntegrationStates(stream);
	
	//auto pClusters = GenericCopyToHost(batch->pClusterDevice, nPclusters);
	//DebugUtils::VerifyIdentical(pClusters, "PClustersBeforeClustering" + std::to_string(batch->step));
	//DebugUtils::VerifyIdentical(batch->pClusterDevice, nPclusters, "PClustersBeforeClustering", batch->step);

	if (nPclusters > 0) ApplyBoundaryCondition << <(nPclusters + 31) / 32, 32, 0, stream >> > (
		batch->pClusterDevice.Get(), nPclusters, boxSizeF, boxSizeFInv);

	if (getPclusters && nPclusters > 0) {
		
		for (const auto& sim : batch->simulations) {
			if (!sim.device.active) continue;
			const int count = sim.device.pclusters.count;
			GetPclusterPositions<<<(count + 31) / 32, 32, 0, stream>>>(
				*batch->pclusterTransfermodule, batch->pClusterDevice.Get(), sim.device.pclusters.offset,
				count, sim.device.gridnodes.offset, boxSize, boxSizeF, boxSizeFInv);
		}
		LIMA_UTILS::genericErrorCheckNoSync("Error after GetPclusterPositions kernel");	

		//DebugUtils::VerifyIdentical(batch->pclusterTransfermodule->meanPositionOfPClustersPerBlock)

		SortPClusterIndicesInBlocks <<<nBlocks, 32, 0, stream>>> (*batch->pclusterTransfermodule);
		LIMA_UTILS::genericErrorCheckNoSync("Error after SortPClusterIndicesInBlocks kernel");
	}

	ClusteringPretransferKernel<PeriodicBoundaryCondition>
		<<<nBlocks, 32, 0, stream >>> (*batch->pclusterTransfermodule, boxSize, boxSizeF);
	LIMA_UTILS::genericErrorCheckNoSync("Error after ClusteringPretransferKernel");

	//DebugUtils::VerifyIdentical(batch->pclusterTransfermodule->idsOfIncomingClusters, 6 * PClusterTransfermodule::maxOutgoingClusters * nBlocks, "idsOfIncomingCLusters", batch->step);
	//DebugUtils::VerifyIdentical(batch->pclusterTransfermodule->meanpositionsOfIncomingClusters, 6 * PClusterTransfermodule::maxOutgoingClusters * nBlocks, "meanpositionsOfIncomingCLusters", batch->step);
	
	

	// This simply stages the SC's per block, need to compress after
	ClusteringKernel <<<nBlocks, 32, 0, stream>>> (*batch->pclusterTransfermodule, batch->pClusterDevice.Get(), *batch->superclusterStagingControl, batch->pClusterMetaDevice.Get(), boxSize, boxSizeF);
	LIMA_UTILS::genericErrorCheckNoSync("Error after ClusteringKernel");
	//DebugUtils::VerifyIdentical(batch->superclusterStagingControl->scData, nBlocks * SuperClustersControl::maxClustersPerBlock, "RunClustering_SCData", batch->step);

	/*SortAndCompare(
		GenericCopyToHost(positionsUnsorted, nBlocks * PClusterTransfermodule::maxClustersPerBlock), 
		GenericCopyToHost(positionsSorted, nBlocks * PClusterTransfermodule::maxClustersPerBlock),
		GenericCopyToHost(idsUnsorted, nBlocks * PClusterTransfermodule::maxClustersPerBlock)
	);*/

	// Compute prefixsum buffer on nScPerBlock
	{
		const int nElements = nBlocks + 1; // Extra element for total sum at end
		//DebugUtils::VerifyIdentical(batch->superclusterStagingControl->nClustersPerBlock, nElements, "nClustersPerBlock", batch->step);
		CubWrappers::ExclusiveScan(batch->superclusterStagingControl->nClustersPerBlock,
			batch->superclusterStagingControl->nClustersPerBlock + nElements, batch->superclusterStagingControl->nClustersPrefixSum, stream);
		LIMA_UTILS::genericErrorCheckNoSync("Error after Prefixsum");
		cudaMemcpyAsync(&batch->nSuperclusters, batch->superclusterStagingControl->nClustersPrefixSum + nElements - 1,
			sizeof(int), cudaMemcpyDeviceToHost, stream);
		// Capacity status rides on this existing sync. PME runs on its own stream, so its status may lag
		// by a step; FinalizeSimulation does a fully synchronized check.
		int overflow[2][CapacityOverflow::nValues]{};
		cudaMemcpyAsync(overflow[0], batch->pclusterTransfermodule->overflow.status, sizeof(overflow[0]), cudaMemcpyDeviceToHost, stream);
		const CapacityOverflow* pmeOverflow = batch->pmeController ? batch->pmeController->Overflow() : nullptr;
		if (pmeOverflow)
			cudaMemcpyAsync(overflow[1], pmeOverflow->status, sizeof(overflow[1]), cudaMemcpyDeviceToHost, stream);
		cudaStreamSynchronize(stream);
		CapacityOverflow::ThrowIfSet(overflow[0]);
		CapacityOverflow::ThrowIfSet(overflow[1]);
	}

	// Clusters are compressed in bin order, so each simulation remains a range.
	std::vector<int> boundaries(batch->simulations.size() * 2);
	for (size_t id = 0; id < batch->simulations.size(); ++id) {
		const auto& sim = batch->simulations[id];
		if (!sim.device.active) continue;
		cudaMemcpyAsync(&boundaries[2 * id], batch->superclusterStagingControl->nClustersPrefixSum + sim.device.gridnodes.offset,
			sizeof(int), cudaMemcpyDeviceToHost, stream);
		cudaMemcpyAsync(&boundaries[2 * id + 1], batch->superclusterStagingControl->nClustersPrefixSum + sim.device.gridnodes.offset + sim.device.gridnodes.count,
			sizeof(int), cudaMemcpyDeviceToHost, stream);
	}
	cudaStreamSynchronize(stream);
	for (size_t id = 0; id < batch->simulations.size(); ++id) {
		if (!batch->simulations[id].device.active) continue;
		batch->simulations[id].superclusters = {boundaries[2 * id], boundaries[2 * id + 1] - boundaries[2 * id]};
	}

	CompressSuperclusters<<<nBlocks, 32, 0, stream>>>(*batch->superClustersControl, *batch->superclusterStagingControl);
	for (int id = 0; id < batch->simulations.size(); ++id) {
		const auto range = batch->simulations[id].superclusters;
		if (range.count > 0)
			SetSuperclusterSimulationId<<<(range.count + 127) / 128, 128, 0, stream>>>(
				batch->superClustersControl->scMeta, range.offset, range.count, id);
	}

	batch->pclusterTransfermodule->Reset(batch->nGridnodes, stream);

	batch->pclusterParticleSlots.Expand(size_t(nPclusters) * PersistentCluster::maxParticles);
	cudaMemsetAsync(batch->pclusterParticleSlots.Get(), 0xFF, sizeof(int) * nPclusters * PersistentCluster::maxParticles, stream);
	if (batch->nSuperclusters > 0)
		MapPclusterParticlesToSlots<<<(batch->nSuperclusters + 3) / 4, dim3(SuperCluster::maxParticles, 4, 1), 0, stream>>>(
			batch->superClustersControl->scMeta, batch->nSuperclusters, batch->pclusterParticleSlots.Get());
	const int nBondgroupParticles = batch->nBondgroupParticles;
	batch->bondgroupParticleSlots.Expand(nBondgroupParticles);
	if (nBondgroupParticles > 0)
		MapBondgroupParticlesToSlots<<<(nBondgroupParticles + 127) / 128, 128, 0, stream>>>(batch->bondgroupParticles.Get(), nBondgroupParticles,
			batch->pclusterParticleSlots.Get(), batch->bondgroupParticleSlots.Get());
	LoadIntegrationStates(stream);


	//DebugUtils::VerifyIdentical(batch->superClustersControl->scData, batch->nSuperclusters, "RunClustering_SCData_Compressed", batch->step);
}

void Engine::BootstrapClustering(cudaStream_t stream) {
	const auto boxSize = batch->boxSize;
	const Float3 boxSizeF = NodeIndex(boxSize).toFloat3();
	const int nBlocks = batch->nGridnodes;
	std::vector<Float3> positions(nBlocks * PClusterTransfermodule::maxClustersPerBlock);
	std::vector<int> ids(nBlocks * PClusterTransfermodule::maxClustersPerBlock);
	std::vector<int> counts(nBlocks);
	for (const auto& sim : batch->simulations) {
		if (!sim.device.active) continue;
		const auto& box = *sim.simulation->box;
		for (int pc = 0; pc < sim.device.pclusters.count; ++pc) {
			const int pcId = pc + sim.device.pclusters.offset;
			Float3 pos = box.persistentClusters[pc].pqd[0].position;
			NodeIndex node{static_cast<int>(floorf(pos.x)), static_cast<int>(floorf(pos.y)), static_cast<int>(floorf(pos.z))};
			BoundaryConditionPublic::applyBC(node, boxSize, PBC);
			BoundaryConditionPublic::applyHyperposNM(node.toFloat3(), pos, boxSizeF, PBC);
			const int bin = sim.device.gridnodes.offset + BoxGrid::Get1dIndex(node, boxSize);
			if (counts[bin] >= PClusterTransfermodule::maxClustersPerBlock)
				throw std::runtime_error("Too many pclusters in bootstrap grid bin");
			const int dst = bin * PClusterTransfermodule::maxClustersPerBlock + counts[bin]++;
			positions[dst] = pos;
			ids[dst] = pcId;
		}
	}
	cudaMemcpyAsync(batch->pclusterTransfermodule->meanPositionOfPClustersPerBlock, positions.data(), sizeof(Float3) * positions.size(), cudaMemcpyHostToDevice, stream);
	cudaMemcpyAsync(batch->pclusterTransfermodule->idsOfPclustersInBlocks, ids.data(), sizeof(int) * ids.size(), cudaMemcpyHostToDevice, stream);
	cudaMemcpyAsync(batch->pclusterTransfermodule->nPClustersPerBlock, counts.data(), sizeof(int) * counts.size(), cudaMemcpyHostToDevice, stream);
	RunClustering(stream, false);
}

namespace EngineLimitTesting {
	void ClusterTransfer(int count) {
		constexpr Int3 boxSize{4, 4, 4};
		constexpr int blocks = 64;
		const int source = BoxGrid::Get1dIndex(Int3{1, 1, 1}, boxSize);
		const int target = BoxGrid::Get1dIndex(Int3{2, 1, 1}, boxSize);
		auto transfer = MakeTransferModule(blocks);
		std::vector<int> counts(blocks), ids(blocks * PClusterTransfermodule::maxClustersPerBlock, -1);
		std::vector<Float3> positions(ids.size());
		counts[source] = count;
		for (int i = 0; i < count; ++i) {
			const int index = source * PClusterTransfermodule::maxClustersPerBlock + i;
			ids[index] = i;
			positions[index] = Float3{2.1f, 1.5f, 1.5f};
		}
		LIMA_UTILS::genericErrorCheck(cudaMemcpy(transfer->nPClustersPerBlock, counts.data(), counts.size() * sizeof(int), cudaMemcpyHostToDevice));
		LIMA_UTILS::genericErrorCheck(cudaMemcpy(transfer->idsOfPclustersInBlocks, ids.data(), ids.size() * sizeof(int), cudaMemcpyHostToDevice));
		LIMA_UTILS::genericErrorCheck(cudaMemcpy(transfer->meanPositionOfPClustersPerBlock, positions.data(), positions.size() * sizeof(Float3), cudaMemcpyHostToDevice));
		// Sentinel slots catch writes crossing a logical queue even if they stay within the allocation.
		std::vector<int> incoming(blocks * 6 * PClusterTransfermodule::maxOutgoingClusters, -1);
		LIMA_UTILS::genericErrorCheck(cudaMemcpy(transfer->idsOfIncomingClusters, incoming.data(), incoming.size() * sizeof(int), cudaMemcpyHostToDevice));
		CheckCuda();
		ClusteringPretransferKernel<PeriodicBoundaryCondition><<<blocks, 32>>>(*transfer, boxSize, Float3{4.f});
		CheckCuda();
		const auto actual = GenericCopyToHost(transfer->idsOfIncomingClusters, incoming.size());
		const auto actualCounts = GenericCopyToHost(transfer->nIncomingClusters, blocks * 6);
		const auto remaining = GenericCopyToHost(transfer->nPClustersPerBlock, blocks);
		const int first = target * 6 * PClusterTransfermodule::maxOutgoingClusters;
		for (size_t i = 0; i < actual.size(); ++i) {
			const bool inQueue = i >= first && i < first + PClusterTransfermodule::maxOutgoingClusters;
			Require(inQueue || actual[i] == -1, "Cluster transfer overwrote a neighboring queue");
		}
		transfer->overflow.Check();
		Require(count <= PClusterTransfermodule::maxOutgoingClusters, "Cluster transfer accepted more than 8 outgoing clusters without rejecting the overflow");
		Require(remaining[source] == 0 && actualCounts[target * 6] == count, "Cluster transfer lost particles");
		for (int i = 0; i < count; ++i)
			Require(actual[first + i] == i, "Cluster transfer changed particle identities");
	}

	void ClusterOccupancy(int count) {
		// One cell starts with count-1 clusters and receives one. This exercises the runtime merge,
		// independently of bootstrap validation and GetPclusterPositions rebuilding membership.
		auto transfer = MakeTransferModule(1);
		const int resident = count - 1;
		LIMA_UTILS::genericErrorCheck(cudaMemcpy(transfer->nPClustersPerBlock, &resident, sizeof(int), cudaMemcpyHostToDevice));
		std::vector<int> ids(PClusterTransfermodule::maxClustersPerBlock);
		std::iota(ids.begin(), ids.end(), 0);
		std::vector<Float3> positions(ids.size(), Float3{0.5f});
		std::vector<PersistentCluster> clusters(count);
		std::vector<PersistentClusterMeta> metadata(count);
		for (int i = 0; i < count; ++i) {
			clusters[i].pqd[0] = PData{Float3{0.5f}, NBParams{0.f, 0.f, 0.f}};
			metadata[i].nParticles = 1;
			metadata[i].particleIdsGlobal[0] = i;
		}
		LIMA_UTILS::genericErrorCheck(cudaMemcpy(transfer->idsOfPclustersInBlocks, ids.data(), ids.size() * sizeof(int), cudaMemcpyHostToDevice));
		LIMA_UTILS::genericErrorCheck(cudaMemcpy(transfer->meanPositionOfPClustersPerBlock, positions.data(), positions.size() * sizeof(Float3), cudaMemcpyHostToDevice));
		const int incoming = 1;
		LIMA_UTILS::genericErrorCheck(cudaMemcpy(transfer->nIncomingClusters, &incoming, sizeof(int), cudaMemcpyHostToDevice));
		LIMA_UTILS::genericErrorCheck(cudaMemcpy(transfer->idsOfIncomingClusters, &resident, sizeof(int), cudaMemcpyHostToDevice));
		LIMA_UTILS::genericErrorCheck(cudaMemcpy(transfer->meanpositionsOfIncomingClusters, positions.data(), sizeof(Float3), cudaMemcpyHostToDevice));
		CudaBuffer<PersistentCluster> clustersDevice;
		CudaBuffer<PersistentClusterMeta> metadataDevice;
		clustersDevice.SetData(clusters);
		metadataDevice.SetData(metadata);
		std::unique_ptr<SuperclusterStagingControl, FreeDeviceMembers<SuperclusterStagingControl>> staging(new SuperclusterStagingControl(1));
		CheckCuda();
		ClusteringKernel<<<1, 32>>>(*transfer, clustersDevice.Get(), *staging, metadataDevice.Get(), Int3{1, 1, 1}, Float3{1.f});
		CheckCuda();
		transfer->overflow.Check();
		Require(count <= PClusterTransfermodule::maxClustersPerBlock, "Runtime clustering accepted more than 128 clusters without rejecting the overflow");
		const auto nClusters = GenericCopyToHost(staging->nClustersPerBlock, 1).front();
		Require(nClusters > 0 && nClusters <= SuperClustersControl::maxClustersPerBlock, "Invalid supercluster count");
		const auto output = GenericCopyToHost(staging->scMeta, nClusters);
		std::vector<int> actualIds;
		for (const auto& meta : output) {
			Require(meta.nParticles > 0 && meta.nParticles <= SuperCluster::maxParticles, "Invalid supercluster particle count");
			for (int i = 0; i < meta.nParticles; ++i) actualIds.push_back(meta.globalParticleIds[i]);
		}
		std::ranges::sort(actualIds);
		Require(actualIds.size() == count, "Runtime clustering lost or duplicated particles");
		for (int i = 0; i < count; ++i) Require(actualIds[i] == i, "Runtime clustering changed particle identities");
	}

}
