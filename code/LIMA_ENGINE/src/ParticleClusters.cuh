#include "EngineBodies.cuh"
#include "Engine.cuh"


static const NodeIndex directions[6]{
	{1, 0, 0},
	{-1, 0, 0},
	{0, 1, 0},
	{0, -1, 0},
	{0, 0, 1},
	{0, 0, -1}
};



template <typename BoundaryCondition>
__global__ void SolventPretransferKernel(PClusterTransfermodule transferModule, Int3 boxSize) {
	

	__shared__ int nClusters;
	__shared__ int indexOfFirstCluster;
	__shared__ int directionIndexOfPCluster[PClusterTransfermodule::maxClustersPerBlock]; // -1 for stay

	if (threadIdx.x == 0) {
		nClusters = transferModule.nPClustersPerBlock[blockIdx.x];
		indexOfFirstCluster = blockIdx.x * PClusterTransfermodule::maxClustersPerBlock;
	}
	if (threadIdx.x < PClusterTransfermodule::maxClustersPerBlock) {
		directionIndexOfPCluster[threadIdx.x] = -1;
	}
	__syncthreads();

	// TODO FIx this part
	//for (int i = threadIdx.x; i < nParticlesInCluster; i += blockDim.x) {
	//	const int particleIndex = indexOfFirstCluster + i;
	//	const Float3 absPos = transferModule->positionsOfFirstParticleInPCluster[particleIndex];
	//	//const Float3 posRelativeToBlock = 
	//	// Account for blockLen > 1
	//	const NodeIndex direction = LIMAPOSITIONSYSTEM::getTransferDirection(posRelativeToBlock);					//OPTIM here
	//	for (int directionIndex = 0; directionIndex < 6; directionIndex++) {
	//		if (direction == directions[directionIndex]) {
	//			directionIndexOfPCluster[i] = directionIndex;
	//			break;
	//		}
	//	}
	//}
	//__syncthreads();

	// First 6 threads are responsible for marking a direction
	__shared__ int nClustersThisDirection[6];
	__shared__ int clusterIdsThisDirectionRelativeToBlock[PClusterTransfermodule::maxOutgoingClusters * 6];
	if (threadIdx.x < 6) {
		int myCount = 0;
		for (int i = 0; i < nClusters; i++) {
			if (directionIndexOfPCluster[i] == threadIdx.x) {
				if constexpr (INDEXING_CHECKS) {
					if (myCount >= PClusterTransfermodule::maxOutgoingClusters)
						printf("Too many pClusters in one direction");
				}
				
				clusterIdsThisDirectionRelativeToBlock[threadIdx.x * PClusterTransfermodule::maxOutgoingClusters + myCount] = i;
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
			if (targetBlockId < 0 || targetBlockId > BoxGrid::BlocksTotal(DeviceConstants::boxSize.blocksPerDim))
				printf("Target block was out of bounds");

			if (targetBlockId >= DeviceConstants::boxSize.blocksPerDim.x * DeviceConstants::boxSize.blocksPerDim.y * DeviceConstants::boxSize.blocksPerDim.z)
				printf("Target block was out of bounds");
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

			transferModule.idsOfIncomingClusters[clusterTargetGlobalIndex] = transferModule.idsOfPclustersInBlocks[clusterSrcGlobalIndex];
			transferModule.meanpositionsOfIncomingClusters[clusterTargetGlobalIndex] = transferModule.meanPositionOfPClusters[clusterSrcGlobalIndex];
		}
	}
	__syncthreads();


	// Finally compress remainders. Reuse the direction-of-particle buffer as prefixsum buffer
	__shared__ int indexRelativeToBlockOfRemainingClusters[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int nClustersRemaining;	
	/*for (int i = threadIdx.x; i < nClusters; i += blockDim.x) {
		bondgroupRemains[i] = directionIndexOfPCluster[i] == -1;
	}*/

	// TODO OPTIM: Find a not idiot way of doing this compression...
	if (threadIdx.x == 0) {
		int nextDestIndex = 0;
		for (int i = 0; i < nClusters; i++) {
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
		const int srcIndex = indexOfFirstCluster + indexRelativeToBlockOfRemainingClusters[i];
		const int destIndex = indexOfFirstCluster + i;
		transferModule.idsOfPclustersInBlocks[destIndex] = transferModule.idsOfPclustersInBlocks[srcIndex];
		transferModule.meanPositionOfPClusters[destIndex] = transferModule.meanPositionOfPClusters[srcIndex];
	}
}


//__global__ void ClusteringKernel(const PClusterTransfermodule transferModule, const Float3* const pclustersParticlesPositions, SuperClusterGridData scGridData,
//	SuperClusterMeta* scMetaOut, SuperCluster* scDataOut) 
//{
//	__shared__ Float3 meanPositionsOfPClusters[PClusterTransfermodule::maxClustersPerBlock];
//	__shared__ int idsOfPclustersInBlock[PClusterTransfermodule::maxClustersPerBlock];
//	__shared__ int nPclustersInBlock;
//
//	if (threadIdx.x == 0) {
//		nPclustersInBlock = transferModule.nPClustersPerBlock[blockIdx.x];	
//	}
//	__syncthreads();
//
//	// Load remaining
//	for (int i = threadIdx.x; i < nPclustersInBlock; i += blockDim.x) {
//		meanPositionsOfPClusters[i] = transferModule.meanPositionOfPClusters[blockIdx.x * PClusterTransfermodule::maxClustersPerBlock + i];
//		transferModule.idsOfPclustersInBlocks[i] = transferModule.idsOfPclustersInBlocks[blockIdx.x * PClusterTransfermodule::maxClustersPerBlock + i];
//	}
//
//	// Load incoming
//	for (int i = 0; i < 6; i++) {
//		if (threadIdx.x < PClusterTransfermodule::maxOutgoingClusters) {
//			size_t srcIndex = blockIdx.x * 6 * PClusterTransfermodule::maxOutgoingClusters + i * PClusterTransfermodule::maxOutgoingClusters + threadIdx.x;
//			meanPositionsOfPClusters[nPclustersInBlock + threadIdx.x] = transferModule.meanpositionsOfIncomingClusters[srcIndex];
//			idsOfPclustersInBlock[nPclustersInBlock + threadIdx.x] = transferModule.idsOfIncomingClusters[srcIndex];
//		}
//		__syncthreads();
//
//		if (threadIdx.x == 0) {
//			nPclustersInBlock += transferModule.nIncomingClusters[blockIdx.x * 6 + i];
//		}
//		__syncthreads();
//	}
//
//	// 
//
//
//
//}

//
//template <int nBins, int nValuesPerBin>
//__device__ Sort(int* keys, int* ids) {
//
//}



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
__global__ void ClusteringKernel(const PClusterTransfermodule transferModule, const PersistentCluster* const pClusters, const Float3* const pclustersParticlesPositions, SuperClustersControl scControl)
{
	__shared__ Float3 meanPositionsOfPClusters[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int idsOfPclustersInBlock[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int nPclustersInBlock;

	__shared__ float sortKeys[PClusterTransfermodule::maxClustersPerBlock];
	__shared__ int sortIds[PClusterTransfermodule::maxClustersPerBlock]; // starts out as iota, tracks relative
	__shared__ int scOutStartIndex;

	if (threadIdx.x == 0) {
		nPclustersInBlock = transferModule.nPClustersPerBlock[blockIdx.x];
	}

	for (int i = threadIdx.x; i < PClusterTransfermodule::maxClustersPerBlock; i += blockDim.x) {
		sortIds[i] = i;
		meanPositionsOfPClusters[i] = Float3{ INFINITY,INFINITY, INFINITY };
	}

	__syncthreads();

	// Load 
	for (int i = threadIdx.x; i < nPclustersInBlock; i += blockDim.x) {
		meanPositionsOfPClusters[i] = transferModule.meanPositionOfPClusters[blockIdx.x * PClusterTransfermodule::maxClustersPerBlock + i];
		idsOfPclustersInBlock[i] = transferModule.idsOfPclustersInBlocks[blockIdx.x * PClusterTransfermodule::maxClustersPerBlock + i];
	}
	__syncthreads();


	const int bucketsPerDim = 4;
	// Sort along z 
	{
		for (int i = threadIdx.x; i < nPclustersInBlock; i+=blockDim.x) {
			sortKeys[i] = meanPositionsOfPClusters[i].z;
		}
		const int bucketSize = PClusterTransfermodule::maxClustersPerBlock;
		LAL::SortBins<1, bucketSize>(sortKeys, sortIds);
		__syncthreads();
	}
	

	// Sort the 4 buckets along y
	{
		for (int i = threadIdx.x; i < nPclustersInBlock; i += blockDim.x) {
			sortKeys[i] = meanPositionsOfPClusters[sortIds[i]].y;
		}
		const int bucketSize = PClusterTransfermodule::maxClustersPerBlock / bucketsPerDim;
		LAL::SortBins<bucketsPerDim, bucketSize>(sortKeys, sortIds);
		__syncthreads();
	}


	// Sort the 4x4 buckets along x
	{
		for (int i = threadIdx.x; i < nPclustersInBlock; i += blockDim.x) {
			sortKeys[i] = meanPositionsOfPClusters[sortIds[i]].x;
		}
		const int bucketSize = PClusterTransfermodule::maxClustersPerBlock / (bucketsPerDim * bucketsPerDim);
		LAL::SortBins<bucketsPerDim * bucketsPerDim, bucketSize>(sortKeys, sortIds);
		__syncthreads();
	}


	const int nClustersToMake = (nPclustersInBlock + 3) / 4;
	if (threadIdx.x == 0) {
		scOutStartIndex = atomicAdd(scControl.nSuperclustersAtomic, nClustersToMake);
	}
	__syncthreads();

	if (threadIdx.x < nClustersToMake){
		SuperClusterMeta scMeta{};
		SuperCluster sc{};

		for (int pcId = 0; pcId < 4; pcId++) {
			const int srcIndex = threadIdx.x * 4 + pcId;
			if (srcIndex < nPclustersInBlock) {
				const int pcIdRelativeToBlock = sortIds[srcIndex];
				const int pcIdGlobal = idsOfPclustersInBlock[pcIdRelativeToBlock];
				scMeta.pclusterIds[pcId] = idsOfPclustersInBlock[pcIdGlobal];
				for (int particleId = 0; particleId < 4; particleId++) {
					sc.pData[pcId * 4 + particleId] = pClusters[pcIdGlobal].pqd[particleId];
				}						
			}
			else {
				//scDataOut[blockIdx.x * SuperClusterGridData::maxSuperClustersPerBlock + i].constituentPClusterIds[j] = -1;
				scMeta.pclusterIds[pcId] = -1;
			}
		}

		scControl.scData[scOutStartIndex + threadIdx.x] = sc;
		scControl.scMeta[scOutStartIndex + threadIdx.x] = scMeta;
	}
}





void Engine::BootstrapClustering() {
	



	Int3 boxSize = simulation->box_host->boxparams.boxSize;
	const int nBlocks = BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(boxSize));;
	
	SolventPretransferKernel<PeriodicBoundaryCondition>
		<<<nBlocks, 32>>>(pclusterTransfermodule, boxSize);


}