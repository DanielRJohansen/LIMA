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
__global__ void GetPclusterPositions(PClusterTransfermodule transferModule, const PersistentCluster* const pClustersData, const int nPclusters, Int3 boxSize) {
	const int pcId = blockIdx.x * blockDim.x + threadIdx.x;
	if (pcId >= nPclusters)
		return;

	//ParticleToCompoundOrSolventMapping mapping = mappings[particleId];


	const Float3 pos = pClustersData[pcId].pqd[0].position; // TODO: use actual mean pos?
	Float3 gridPosF = pos.round();
	PeriodicBoundaryCondition::applyBCNM(gridPosF);
	NodeIndex blockId = NodeIndex(gridPosF.x, gridPosF.y, gridPosF.z);
	const int blockIndex = BoxGrid::Get1dIndex(blockId, boxSize);
	int indexInBlock = atomicAdd(&transferModule.nPClustersPerBlock[blockIndex], 1);
	int index = blockIndex * PClusterTransfermodule::maxClustersPerBlock + indexInBlock;

	transferModule.idsOfPclustersInBlocks[index] = pcId;
	transferModule.meanPositionOfPClustersPerBlock[index] = pos; // TODO: THis is not even the actual meanposition is it?
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
	if (threadIdx.x < PClusterTransfermodule::maxClustersPerBlock) {
		directionIndexOfPCluster[threadIdx.x] = -1;
	}
	__syncthreads();

	// TODO FIx this part
	const Float3 blockOrigoF = BoxGrid::Get3dIndex(blockIdx.x, boxSize).toFloat3();

	for (int i = threadIdx.x; i < nPClusters; i += blockDim.x) {
		const int pcIndex = indexOfFirstCluster + i;
		Float3 absPos = transferModule.meanPositionOfPClustersPerBlock[pcIndex];
		PeriodicBoundaryCondition::applyHyperposNM(blockOrigoF, absPos);	// TODO: OPTIM: shoudn't be necessary if pClusters are placed correctly in blocks...

		const Float3 posRelativeToBlockCenter = absPos - (blockOrigoF + Float3{ 0.5f, 0.5f, 0.5f });
		const NodeIndex direction = LIMAPOSITIONSYSTEM::getTransferDirection(Coord{ posRelativeToBlockCenter });					//OPTIM here
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
				if constexpr (INDEXING_CHECKS) {
					if (myCount >= PClusterTransfermodule::maxOutgoingClusters)
						printf("Too many pClusters in one direction");
				}
				
				clusterIdsThisDirectionRelativeToBlock[threadIdx.x * PClusterTransfermodule::maxOutgoingClusters + myCount] = i;
				myCount++;
			}
		}
		nClustersThisDirection[threadIdx.x] = myCount;
	}
	__syncthreads();



	// Now all threads loop over the direction, and if they have a particle, they push it directy to the incoming queue in global memory
	for (int directionIndex = 0; directionIndex < 6; directionIndex++) {
		const NodeIndex& direction = directions[directionIndex];
		const NodeIndex blockOrigo = BoxGrid::Get3dIndex(blockIdx.x, boxSize);
		const NodeIndex targetBlock = BoundaryCondition::applyBC(blockOrigo + direction, boxSize);
		const int targetBlockId = BoxGrid::Get1dIndex(targetBlock, boxSize);
		if constexpr (INDEXING_CHECKS) {
			if (targetBlockId < 0 || targetBlockId > BoxGrid::BlocksTotal(boxSize))
				printf("Target block %d %d %d was out of bounds. Id %d out of %d\n", targetBlock.x, targetBlock.y, targetBlock.z, targetBlockId, BoxGrid::BlocksTotal(DeviceConstants::boxSize.blocksPerDim));
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
			//printf("Here: %d\n", directionIndexOfPCluster[i]);
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
__global__ void ClusteringKernel(const PClusterTransfermodule transferModule, const PersistentCluster* const pClusters, SuperclusterStagingControl scStagingControl)
{
	__shared__ Float3 meanPositionsOfPClusters[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int idsOfPclustersInBlock[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int nPclustersInBlock;

	__shared__ float sortKeys[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int sortIds[PClusterTransfermodule::maxClustersPerBlock]; // starts out as iota, tracks relative
	//__shared__ int scOutStartIndex;

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
		if (threadIdx.x == 0) {
			nPclustersInBlock += transferModule.nIncomingClusters[blockIdx.x * 6 + dir];
		}
	}
	__syncthreads();
	/*if (threadIdx.x == 0 && nPclustersInBlock > 0) {
		printf("pCID %d\n", idsOfPclustersInBlock[0]);
	}*/

	const int bucketsPerDim = 4;
	// Sort along z 
	{
		for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i+=blockDim.x) {
			sortKeys[i] = meanPositionsOfPClusters[i].z;
		}
		const int bucketSize = PClusterTransfermodule::maxClustersPerBlock;
		LAL::SortBins<1, bucketSize>(sortKeys, sortIds);
		__syncthreads();
	}


	// Sort the 4 buckets along y
	{
		for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i += blockDim.x) {
			sortKeys[i] = meanPositionsOfPClusters[sortIds[i]].y;
		}
		const int bucketSize = PClusterTransfermodule::maxClustersPerBlock / bucketsPerDim;
		LAL::SortBins<bucketsPerDim, bucketSize>(sortKeys, sortIds);
		__syncthreads();
	}


	// Sort the 4x4 buckets along x
	{
		for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i += blockDim.x) {
			sortKeys[i] = meanPositionsOfPClusters[sortIds[i]].x;
		}
		const int bucketSize = PClusterTransfermodule::maxClustersPerBlock / (bucketsPerDim * bucketsPerDim);
		LAL::SortBins<bucketsPerDim * bucketsPerDim, bucketSize>(sortKeys, sortIds);
		__syncthreads();
	}


	const int nClustersToMake = (nPclustersInBlock + 3) / 4;	
	if (threadIdx.x == 0) {
		scStagingControl.nClustersPerBlock[blockIdx.x] = nClustersToMake;
		if constexpr (INDEXING_CHECKS) {
			if (nClustersToMake >= SuperClustersControl::maxClustersPerBlock) {
				printf("Trying to make too many superclusters in block %d: %d (max %d)\n", blockIdx.x, nClustersToMake, SuperClustersControl::maxClustersPerBlock);
			}				
		}
	}

	
	if (threadIdx.x < nClustersToMake){
		SuperClusterMeta scMeta{};
		SuperCluster sc{};

		for (int pcId = 0; pcId < 4; pcId++) {
			const int srcIndex = threadIdx.x * 4 + pcId;
			if (srcIndex < nPclustersInBlock) {


				const int pcIdRelativeToBlock = sortIds[srcIndex];
				const int pcIdGlobal = idsOfPclustersInBlock[pcIdRelativeToBlock];
				//printf("rel %d global %d value %d\n", pcIdRelativeToBlock, pcIdGlobal, idsOfPclustersInBlock[pcIdGlobal]);
				scMeta.pclusterIds[pcId] = pcIdGlobal;
				//pClusters[pcIdGlobal].pqd[0].position.print('p');
				/*printf("posx %f block %d, idsat0 %d sortId %d, pcIdGlobal %d, pcIdRelativeToBlock %d\n", pClusters[pcIdGlobal].pqd[0].position.x, blockIdx.x, idsOfPclustersInBlock[0], pcIdRelativeToBlock, pcIdGlobal, pcIdRelativeToBlock);
				printf("Clusting %f %f %f\n", pClusters[pcIdGlobal].pqd[0].position.x, pClusters[pcIdGlobal].pqd[0].position.y, pClusters[pcIdGlobal].pqd[0].position.z);*/
				for (int particleId = 0; particleId < 4; particleId++) {
					sc.pData[pcId * 4 + particleId] = pClusters[pcIdGlobal].pqd[particleId];					
				}						
			}
			else {
				//scDataOut[blockIdx.x * SuperClusterGridData::maxSuperClustersPerBlock + i].constituentPClusterIds[j] = -1;
				scMeta.pclusterIds[pcId] = -1;
			}
		}

		const int stagingStartIndex = blockIdx.x * SuperClustersControl::maxClustersPerBlock;
		scStagingControl.scData[stagingStartIndex + threadIdx.x] = sc;
		scStagingControl.scMeta[stagingStartIndex + threadIdx.x] = scMeta;
	}
}

// 1 cudablock per block, blockdim = 32,1,1
__global__ void CompressSuperclusters(SuperClustersControl scControl, const SuperclusterStagingControl scStagingControl) {
	const int nClusters = scStagingControl.nClustersPerBlock[blockIdx.x];
	if (nClusters == 0)
		return;

	const int srcStartIndex = blockIdx.x * SuperClustersControl::maxClustersPerBlock;
	const int dstStartIndex = scStagingControl.nClustersPrefixSum[blockIdx.x];


	auto tb = cooperative_groups::this_thread_block();
	cooperative_groups::memcpy_async(tb, &scControl.scData[dstStartIndex], &scStagingControl.scData[srcStartIndex], sizeof(SuperCluster) * nClusters);
	cooperative_groups::memcpy_async(tb, &scControl.scMeta[dstStartIndex], &scStagingControl.scMeta[srcStartIndex], sizeof(SuperClusterMeta) * nClusters);
	cooperative_groups::wait(tb);
}





void Engine::RunClustering(bool getPclusters) {
	Int3 boxSize = simulation->box_host->boxparams.boxSize;
	const int nBlocks = BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxSize));


	if (!superclusterStagingControl) {
		superclusterStagingControl = std::make_unique<SuperclusterStagingControl>(boxSize);
	}


	if (getPclusters) {
		const int nPclusters = simulation->box_host->persistentClusters.size();
		int nCudablocks = (nPclusters + 31) / 32;
		GetPclusterPositions<<<nCudablocks, 32>>>(
			*pclusterTransfermodule,
			pClusterDevice,
			nPclusters,
			simulation->box_host->boxparams.boxSize);
		LIMA_UTILS::genericErrorCheckNoSync("Error after GetPclusterPositions kernel");	

		SortPClusterIndicesInBlocks <<<nBlocks, 32>>> (*pclusterTransfermodule);
		LIMA_UTILS::genericErrorCheckNoSync("Error after SortPClusterIndicesInBlocks kernel");
	}

	ClusteringPretransferKernel<PeriodicBoundaryCondition>
		<<<nBlocks, 32 >>> (*pclusterTransfermodule, boxSize);

	LIMA_UTILS::genericErrorCheckNoSync("Error after ClusteringPretransferKernel");

	// This simply stages the SC's per block, need to compress after
	ClusteringKernel <<<nBlocks, 32>>> (*pclusterTransfermodule, pClusterDevice, *superclusterStagingControl);
	LIMA_UTILS::genericErrorCheckNoSync("Error after ClusteringKernel");


	// Compute prefixsum buffer on nScPerBlock
	{
		const int nElements = nBlocks + 1; // Extra element for total sum at end
		thrust::exclusive_scan(thrust::device, superclusterStagingControl->nClustersPerBlock, superclusterStagingControl->nClustersPerBlock + nElements, superclusterStagingControl->nClustersPrefixSum);
		LIMA_UTILS::genericErrorCheckNoSync("Error after Prefixsum");
		nSuperclusters = GenericCopyToHost<int>(superclusterStagingControl->nClustersPrefixSum + nElements-1);

		std::vector<int> counts = GenericCopyToHost(superclusterStagingControl->nClustersPerBlock, nElements);
		int sum = std::accumulate(counts.begin(), counts.end(), 0);
		if (sum != nSuperclusters)
			throw std::runtime_error("Prefixsum mismatch in clustering");
	}

	CompressSuperclusters<<<nBlocks, 32>>>(*superClustersControl, *superclusterStagingControl);

	pclusterTransfermodule->Reset(boxSize);


	// temp
	//const int nSuperclusters = GenericCopyToHost(superClustersControl->nSuperclustersAtomic);
	/*std::vector<SuperClusterMeta> meta = GenericCopyToHost(superClustersControl->scMeta, nSuperclusters);
	DebugUtils::VerifyIdentical(meta, "SCMeta" + std::to_string(simulation->getStep()));*/
}

void Engine::BootstrapClustering() {
	Int3 boxSize = simulation->box_host->boxparams.boxSize;
	const int nBlocks = BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxSize));;
	const Float3 boxSizeF = Float3(boxSize.x, boxSize.y, boxSize.z);
	const Box& box = *simulation->box_host;
	
	std::vector<Float3> meanPositionsofPclusters(nBlocks * PClusterTransfermodule::maxClustersPerBlock);	
	std::vector<int> idsOfPclustersInBlocks(nBlocks * PClusterTransfermodule::maxClustersPerBlock);
	std::vector<int> nPclustersPerBlock(nBlocks);

	for (int pcId = 0; pcId < box.persistentClusters.size(); pcId++) {

		Float3 pos = box.persistentClusters[pcId].pqd[0].position;
		BoundaryConditionPublic::applyBCNM(pos, boxSizeF, BoundaryConditionSelect::PBC);	// TODO: OPTIM: Shouldnt be necessary, maybe just BC every step before resorting??

		NodeIndex targetBlock{ static_cast<int>(floorf(pos.x)), static_cast<int>(floorf(pos.y)), static_cast<int>(floorf(pos.z)) };
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

