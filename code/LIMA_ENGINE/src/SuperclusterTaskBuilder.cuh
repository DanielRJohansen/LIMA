#include "Engine.cuh"
#include "EngineBodies.cuh"


class InteractionToken {
	uint32_t data{};
public:
	constexpr InteractionToken() {}
	constexpr InteractionToken(int queryId, bool useNointeractionMatrix) {
		data = ((uint32_t)queryId << 1) | (uint32_t)useNointeractionMatrix;
	}
	constexpr int GetQueryId() const { return (int)(data >> 1); }
	constexpr bool UseNointeractionMatrix() const { return data & 1; }
};


// Pass to GPU
struct TaskBuilderControlContents {
	static const int maxTasksPerSc = 256;

	std::array<float4, 4>* superclusterPositionSpheres;

	ParticlesBondedToParticle* particlesBondedToParticle;
	PclustersBondedToPcluster* pclustersBondedToPcluster;

	int* nInteractionsOwned;
	InteractionToken* interactionsOwned;
	int* nInteractionsNonowned;
	int* scIdsQueryNonowned;
	int* nResults;

	int* nResultsPrefixsum;
	int* nQueryBuffersPrefixsum;
};

// Keep on CPU
class TaskBuilderControl {

public:
	const int nSuperclustersUpperbound;
	TaskBuilderControlContents contents;

	TaskBuilderControl(const TaskBuilderControl&) = delete;
	TaskBuilderControl& operator=(const TaskBuilderControl&) = delete;
	TaskBuilderControl(int nSuperclustersUpperbound, const std::vector<ParticlesBondedToParticle>& particlesBondedToParticle,
		const std::vector<PclustersBondedToPcluster>& pclustersBondedToPcluster, cudaStream_t stream)
		: nSuperclustersUpperbound(nSuperclustersUpperbound)
	{
		contents.particlesBondedToParticle = GenericCopyToDevice(particlesBondedToParticle);
		contents.pclustersBondedToPcluster = GenericCopyToDevice(pclustersBondedToPcluster);

		cudaMalloc(&contents.superclusterPositionSpheres, sizeof(std::array<float4, 4>) * nSuperclustersUpperbound);
		cudaMalloc(&contents.nInteractionsOwned, sizeof(int) * (nSuperclustersUpperbound + 1));
		cudaMalloc(&contents.interactionsOwned, sizeof(InteractionToken) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		cudaMalloc(&contents.nInteractionsNonowned, sizeof(int) * nSuperclustersUpperbound);
		cudaMalloc(&contents.scIdsQueryNonowned, sizeof(int) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		cudaMalloc(&contents.nResults, sizeof(int) * (nSuperclustersUpperbound + 1));

		cudaMalloc(&contents.nResultsPrefixsum, sizeof(int) * (nSuperclustersUpperbound + 1));
		cudaMalloc(&contents.nQueryBuffersPrefixsum, sizeof(int) * (nSuperclustersUpperbound + 1));

		Reset(stream);
	}

	void Reset(cudaStream_t stream) {
		cudaMemsetAsync(contents.interactionsOwned, 0xFF,
			sizeof(InteractionToken) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound, stream);
		//cudaMemset(contents.scIdsQueryNonowned, 0xFF, sizeof(int) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		thrust::fill_n(thrust::cuda::par.on(stream), contents.scIdsQueryNonowned,
			TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound, INT_MAX);
	}

	~TaskBuilderControl() {
		cudaFree(contents.particlesBondedToParticle);
		cudaFree(contents.pclustersBondedToPcluster);
		cudaFree(contents.superclusterPositionSpheres);
		cudaFree(contents.nInteractionsOwned);
		cudaFree(contents.interactionsOwned);
		cudaFree(contents.nInteractionsNonowned);
		cudaFree(contents.scIdsQueryNonowned);
		cudaFree(contents.nResults);

		cudaFree(contents.nResultsPrefixsum);
		cudaFree(contents.nQueryBuffersPrefixsum);
	}
};


//constexpr std::array<int, 16> GetParticleIdsOfSuperCluster(const PersistentClusterMeta* const pClusterMeta, const SuperClusterMeta& scMeta) {
//	std::array<int, 16> particleIds{};
//	for (auto& e : particleIds) { e = -1; }
//	int cnt = 0;
//	for (auto pcId : scMeta.pclusterIds) {
//		if (pcId == -1) break;
//		for (int particleId : pClusterMeta[pcId].particleIdsGlobal) {
//			particleIds[cnt++] = particleId;
//		}
//	}
//	return particleIds;
//}

//__global__ void ComputeMeanposAndRadiiForEachPclusterInEachSuperclusterKernel(const SuperCluster* const superclusters, const SuperClusterMeta* const scMeta, std::array<float4, 16>* const out, int nSuperclusters) {
//	const int scId = blockIdx.x * blockDim.x + threadIdx.x;
//
//	if (scId >= nSuperclusters)
//		return;
//
//	std::array<float4, 16> positionSpheres{};
//
//	Float3 sum{};
//	int cnt = 0;
//	int positionClusterIndex = 0;
//	for (int i = 0; i < scMeta[scId].nParticles; i++){
//		const int pcId = scMeta[scId]._pclusterIds[i];
//		const PData& pdata = superclusters[scId].pData[i];
//
//		sum += pdata.position;
//		cnt++;
//
//		int nextPcId = i == 15 ? -1 : scMeta[scId]._pclusterIds[i + 1];
//		if (nextPcId != pcId) {
//			Float3 meanPos = sum * (1.f / static_cast<float>(cnt));
//			float radius = 0;
//			for (int ii = i -cnt+1; ii <= i; ii++) {
//				const PData& pData = superclusters[scId].pData[ii];
//				radius = std::max(radius, (pData.position - meanPos).len());
//			}
//			out[scId][positionClusterIndex] = float4{ meanPos.x, meanPos.y, meanPos.z, radius };
//			sum = {}; 
//			cnt = 0;
//			positionClusterIndex++;
//		}
//	}
//
//
//	for (int i = positionClusterIndex; i < 16; i++) {
//		out[scId][i] = float4{ 0,0,0,-1. };
//	}
//}

namespace InteractionSpheres {
	__device__ inline float4 Merge(const float4& a, const float4& b) {
		const Float3 newCenter{ (a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f, (a.z + b.z) * 0.5f};
		const float moveDist = (Float3(a) - Float3(newCenter)).len();
		const float newRadius = std::max(a.w, b.w) + moveDist;
		return float4{ newCenter.x, newCenter.y, newCenter.z, newRadius };
	}

	__device__ inline float SphereVolumeProxy(const float4& s) {// The 4/3pi scalar doesnt matter for comparisons
		return s.w * s.w * s.w;
	}

	__device__ inline float MergeCost(const float4& a, const float4& b) {
		const float4 merged = Merge(a, b);
		return SphereVolumeProxy(merged) - SphereVolumeProxy(a) - SphereVolumeProxy(b);
	}
}


__global__ void ComputeMeanposAndRadiiForEachPclusterInEachSuperclusterKernel(const SuperCluster* const superclusters, const SuperClusterMeta* const scMeta, std::array<float4, 4>* const out, int nSuperclusters) {
	const int scId = blockIdx.x * blockDim.x + threadIdx.x;
	if (scId >= nSuperclusters)
		return;

	if constexpr (INDEXING_CHECKS) {
		if (scMeta[scId].nParticles < 0 || scMeta[scId].nParticles > 16)
			printf("Illegal number of particles in supercluster %d: %d\n", scId, scMeta[scId].nParticles);
	}

	int sphereCount = 0;
	static const int nSpheresToPush = 4;
	std::array<float4, 16> spheres;


	Float3 sum{};
	int cnt = 0;

	for (int i = 0; i < scMeta[scId].nParticles; i++) {
		const int pcId = scMeta[scId]._pclusterIds[i];

		sum += superclusters[scId].Position(i);
		cnt++;

		int nextPcId = i == 15 ? -1 : scMeta[scId]._pclusterIds[i + 1];
		if (nextPcId != pcId) {
			Float3 meanPos = sum * (1.f / static_cast<float>(cnt));
			float radius = 0;
			for (int ii = i - cnt + 1; ii <= i; ii++) {
				//const PData& pData = ;

				radius = std::max(radius, (superclusters[scId].Position(ii) - meanPos).len());
			}
			//out[scId][positionClusterIndex] = float4{ meanPos.x, meanPos.y, meanPos.z, radius };
			spheres[sphereCount] = float4{ meanPos.x, meanPos.y, meanPos.z, radius };
			sum = {};
			cnt = 0;
			sphereCount++;
		}
	}

	// Reduce untill 4 spheres left
	while (sphereCount > 4) {
		int bestI = 0;
		int bestJ = 1;
		float bestCost = InteractionSpheres::MergeCost(spheres[0], spheres[1]);

		for (int i = 0; i < sphereCount; i++) {
			for (int j = i + 1; j < sphereCount; j++) {
				const float cost = InteractionSpheres::MergeCost(spheres[i], spheres[j]);
				if (cost < bestCost) {
					bestCost = cost;
					bestI = i;
					bestJ = j;
				}
			}
		}

		spheres[bestI] = InteractionSpheres::Merge(spheres[bestI], spheres[bestJ]);
		spheres[bestJ] = spheres[sphereCount - 1];
		sphereCount--;
	}

	for (int i = 0; i < nSpheresToPush; i++) {
		out[scId][i] = spheres[i];
	}
}

__device__ inline bool Warp_DoesSuperclustersInteractFine(const SuperCluster* const scData, int scId0, int scId1, float cutoffDistance, const Float3& boxSize, const Float3& boxSizeInv) {
	constexpr unsigned int mask = 0xFFFFFFFFu;

	const int lane = threadIdx.x & 31;
	const float cutoffDistanceSq = cutoffDistance * cutoffDistance;

	bool interacts = false;

	for (int pairId = lane; pairId < 16 * 16; pairId += 32) {
		const int i = pairId / 16;
		const int j = pairId % 16;

		const Float3 pos0 = scData[scId0].Position(i);
		Float3 pos1 = scData[scId1].Position(j);
		float eps0 = scData[scId0].epsilonSqrt[i];	// %TODO: OPTIM: pos being nan would mean eps is not needed
		float eps1 = scData[scId1].epsilonSqrt[j];


		if (eps0 != -1 && eps1 != -1) {//if (p0.Valid() && p1.Valid()) {
			PeriodicBoundaryCondition::ApplyHyperpos(pos0, pos1, boxSize, boxSizeInv);

			const Float3 delta = pos0 - pos1;
			const float distanceSq = delta.dot(delta);

			interacts |= distanceSq <= cutoffDistanceSq;
		}
	}

	return __any_sync(mask, interacts);
}

__device__ inline bool Warp_DoesSuperclustersInteract(const std::array<float4, 4>* const superclusterPositionSpheres, const SuperCluster* const scData, int scId0, int scId1, float cutoffDistance, const Float3& boxSize, const Float3& boxSizeInv) {
	constexpr unsigned int mask = 0xFFFFFFFFu;

	if (scId0 > scId1) {
		const int tmp = scId0;
		scId0 = scId1;
		scId1 = tmp;
	}

	const int lane = threadIdx.x & 31;

	bool coarseHit = false;

	for (int pairId = lane; pairId < 4 * 4; pairId += 32) {
		const int posSphereId0 = pairId / 4;
		const int posSphereId1 = pairId % 4;

		float4 p0 = superclusterPositionSpheres[scId0][posSphereId0];
		float4 p1 = superclusterPositionSpheres[scId1][posSphereId1];

		if (p0.w >= 0 && p1.w >= 0) {
			Float3 pos0 = Float3{ p0 };
			Float3 pos1 = Float3{ p1 };

			PeriodicBoundaryCondition::ApplyHyperpos(pos0, pos1, boxSize, boxSizeInv);
			const Float3 delta = pos0 - pos1;
			const float radiusSum = p0.w + p1.w;
			const float coarseCutoff = cutoffDistance + radiusSum;

			coarseHit |= delta.dot(delta) <= coarseCutoff * coarseCutoff;
		}
	}

	if (!__any_sync(mask, coarseHit)) {
		return false;
	}

	return Warp_DoesSuperclustersInteractFine(scData, scId0, scId1, cutoffDistance, boxSize, boxSizeInv);
}

__device__ inline bool Warp_ScAreBonded(const SuperClusterMeta& sc0, const SuperClusterMeta& sc1, const PclustersBondedToPcluster* const pclustersBondedToPcluster) {
	constexpr unsigned int mask = 0xFFFFFFFFu;

	const int lane = threadIdx.x & 31;
	const int nPairs = sc0.nUniquePcIds * sc1.nUniquePcIds;

	bool bonded = false;

	for (int pairId = lane; pairId < nPairs; pairId += 32) {
		const int i = pairId / sc1.nUniquePcIds;
		const int j = pairId % sc1.nUniquePcIds;

		bonded |= pclustersBondedToPcluster[sc0.uniquePclusterIds[i]].Contains(sc1.uniquePclusterIds[j]);
	}

	return __any_sync(mask, bonded);
}

__global__ void ReserveInteractions(SuperClustersControl scControl, Int3 boxSize, TaskBuilderControlContents tbContents, float cutoffNm, Float3 boxSizeFloat, Float3 boxSizeFloatInv) {

	const int gridOffset = (blockIdx.x / boxSize.InnerProduct()) * boxSize.InnerProduct();
	NodeIndex nodeIndex = BoxGrid::Get3dIndex(blockIdx.x % boxSize.InnerProduct(), boxSize);
	const int nodeId = gridOffset + BoxGrid::Get1dIndex(nodeIndex, boxSize);

	__shared__ int nOwnedInteractions;
	__shared__ int nNonownedInteractions;
	__shared__ int nNointeractionMatrices;

	if (blockIdx.y >= scControl.nSuperclustersInBlocks[nodeId]) {
		return;
	}

	const int lane = threadIdx.x & 31;
	const int warpId = threadIdx.x >> 5;
	const int nWarps = blockDim.x >> 5;

	const int scIndexInSelf = blockIdx.y;
	const int scId = scControl.scIdsInBlocks[nodeId * SuperClustersControl::maxClustersPerBlock + blockIdx.y];

	if (threadIdx.x == 0) {
		nOwnedInteractions = 0;
		nNonownedInteractions = 0;
		nNointeractionMatrices = 0;
	}
	__syncthreads();

	constexpr int nNeighborBlocks = 27;
	constexpr int nCandidates =
		nNeighborBlocks * SuperClustersControl::maxClustersPerBlock;

	for (int candidateId = warpId; candidateId < nCandidates; candidateId += nWarps) {
		const int neighborBlockId = candidateId / SuperClustersControl::maxClustersPerBlock;
		const int scIndexInQueryblock = candidateId % SuperClustersControl::maxClustersPerBlock;

		NodeIndex targetBlockRelative =
			BoxGrid::Get3dIndex(neighborBlockId, Int3(3, 3, 3)) - Int3(1, 1, 1);

		NodeIndex targetBlock =
			PeriodicBoundaryCondition::applyBC(nodeIndex + targetBlockRelative, boxSize);

		int targetIndex = gridOffset + BoxGrid::Get1dIndex(targetBlock, boxSize);

		const bool validQuery = scIndexInQueryblock < scControl.nSuperclustersInBlocks[targetIndex];

		int queryScId = -1;
		bool doesInteract = false;
		bool useNointeractionMatrix = false;

		if (validQuery) {
			queryScId = scControl.scIdsInBlocks[targetIndex * SuperClustersControl::maxClustersPerBlock + scIndexInQueryblock];

			doesInteract = Warp_DoesSuperclustersInteract(tbContents.superclusterPositionSpheres, scControl.scData, scId, queryScId, cutoffNm, boxSizeFloat, boxSizeFloatInv);

			if (doesInteract) {
				useNointeractionMatrix = scId == queryScId || Warp_ScAreBonded(scControl.scMeta[scId], scControl.scMeta[queryScId], tbContents.pclustersBondedToPcluster);
			}
		}

		if (lane == 0 && doesInteract) {
			if (scId <= queryScId) {
				int putIndex = atomicAdd(&nOwnedInteractions, 1);
				tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc + putIndex] = InteractionToken(queryScId, useNointeractionMatrix);

				if (useNointeractionMatrix) {
					atomicAdd(&nNointeractionMatrices, 1);
				}
			}
			else {
				int putIndex = atomicAdd(&nNonownedInteractions, 1);
				tbContents.scIdsQueryNonowned[scId * TaskBuilderControlContents::maxTasksPerSc + putIndex] = queryScId;
			}
		}
	}

	__syncthreads();

	if constexpr (INDEXING_CHECKS) {
		if (threadIdx.x == 0 && (nOwnedInteractions + nNonownedInteractions > TaskBuilderControlContents::maxTasksPerSc))
			printf("Too many interactions for scId %d: %d owned + %d nonowned\n", scId, nOwnedInteractions, nNonownedInteractions);
	}

	if (threadIdx.x == 0) {
		tbContents.nInteractionsOwned[scId] = nOwnedInteractions;
		tbContents.nInteractionsNonowned[scId] = nNonownedInteractions;

		tbContents.nResults[scId] = 1 + nNonownedInteractions; // 1 for the sum of all owned interactions, 1 result per unowned interaction
	}
}


__global__ void SortReserveInteractionsOutput(TaskBuilderControlContents tbContents) {
	const int scId = blockIdx.x;

	static_assert(sizeof(InteractionToken) == sizeof(uint32_t), "InteractionToken must be 32 bits");
	uint32_t* intTokenBufferRaw = reinterpret_cast<uint32_t*>(tbContents.interactionsOwned);
	LAL::Sort(&intTokenBufferRaw[scId * TaskBuilderControlContents::maxTasksPerSc], TaskBuilderControlContents::maxTasksPerSc, [](const uint32_t& token) {
		return token;
		});
	LAL::Sort(&tbContents.scIdsQueryNonowned[scId * TaskBuilderControlContents::maxTasksPerSc], TaskBuilderControlContents::maxTasksPerSc, [](const int& id) {
		return id;
		});
}

constexpr int IndexOfId(int* ids, int nIds, int idToFind) {
	for (int i = 0; i < nIds; i++) {
		if (ids[i] == idToFind)
			return i;
	}
	return -1;
}

__device__ inline int GetResultIndexOfQuery(TaskBuilderControlContents tbContents, int scId, int scIdQuery, int selfResultIndex) {
	if (scId == scIdQuery)
		return selfResultIndex;

	const int indexInQuery = IndexOfId(
		&tbContents.scIdsQueryNonowned[scIdQuery * TaskBuilderControlContents::maxTasksPerSc],
		tbContents.nInteractionsNonowned[scIdQuery],
		scId
	);

	if (indexInQuery == -1) {
		printf("Illegal query index for scId %d querying scId %d. Check %d ids\n", scId, scIdQuery, tbContents.nInteractionsNonowned[scIdQuery]);
		return -1;
	}

	const int queryHasOwnedResult = tbContents.nInteractionsOwned[scIdQuery] > 0 ? 1 : 0;
	return tbContents.nResultsPrefixsum[scIdQuery] + queryHasOwnedResult + indexInQuery;
}

__global__ void BuildTasks(
	TaskBuilderControlContents tbContents,
	SuperClusterMeta* const superClusterMeta,
	int nSuperclusters,
	ScScTask* const tasks,
	int* const idsOfQuerySuperclusters,
	int* const resultIndices
) {
	const int scId = blockIdx.x * blockDim.x + threadIdx.x;
	if (scId >= nSuperclusters)
		return;

	const int nOwnedInteractions = tbContents.nInteractionsOwned[scId];
	const int queryBufferStart = tbContents.nQueryBuffersPrefixsum[scId];
	const int selfResultIndex = tbContents.nResultsPrefixsum[scId];

	ScScTask task;
	task.startIndexInQueriesBuffers = queryBufferStart;
	task.nQueryScs = nOwnedInteractions;

	for (int i = 0; i < nOwnedInteractions; i++) {
		const InteractionToken token = tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc + i];
		const int scIdQuery = token.GetQueryId();

		idsOfQuerySuperclusters[queryBufferStart + i] = scIdQuery;
		resultIndices[queryBufferStart + i] = GetResultIndexOfQuery(tbContents, scId, scIdQuery, selfResultIndex);
	}

	tasks[scId] = task;

	superClusterMeta[scId].resultsStartIndex = tbContents.nResultsPrefixsum[scId];
	superClusterMeta[scId].nResults = tbContents.nResults[scId];
}



//// gridDim = (nSuperclusters, 1, 1)
//// blockDim = (16, 1, 1)
__global__ void BuildNointeractionMatricesKernel(
	const SuperClusterMeta* const superClusterMetas,
	const PersistentClusterMeta* const pClustersMeta,
	TaskBuilderControlContents tbContents,
	BoolMatrix16x16* const nointeractionMatrices,
	int nSuperclusters
) {
	const int scId = blockIdx.x;
	const int row = threadIdx.x;

	const int nOwnedInteractions = tbContents.nInteractionsOwned[scId];
	const int queryBufferStart = tbContents.nQueryBuffersPrefixsum[scId];

	for (int i = 0; i < nOwnedInteractions; i++) {
		const InteractionToken token = tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc + i];
		const int matrixIndex = queryBufferStart + i;

		uint16_t rowData = 0;

		if (token.UseNointeractionMatrix()) {
			const int scIdQuery = token.GetQueryId();
			const bool isSelfInteractionTask = scId == scIdQuery;

			for (int col = 0; col < 16; ++col) {
				const int pidSelf = superClusterMetas[scId].globalParticleIds[col];
				const int pidQuery = superClusterMetas[scIdQuery].globalParticleIds[row];

				if (pidSelf == -1 || pidQuery == -1)
					continue;

				bool noInteraction = tbContents.particlesBondedToParticle[pidSelf].Contains(pidQuery);

				if (isSelfInteractionTask && row == col)
					noInteraction = true;

				if (noInteraction)
					BoolMatrix16x16::SetValueInRow(col, rowData);
			}
		}

		nointeractionMatrices[matrixIndex].SetRow(row, rowData);
	}
}















bool Engine::MakeSuperClusterTasksGPU(cudaStream_t stream) {
	if (batch->nSuperclusters == 0)
		return true;

	const Int3 boxSize = batch->boxSize;
	const Float3 boxSizeF = NodeIndex(boxSize).toFloat3();

	const int nSuperclustersUpperbound = batch->nSuperclusters * 2;

	if (!batch->taskbuilderControl || batch->taskbuilderControl->nSuperclustersUpperbound < batch->nSuperclusters)
		batch->taskbuilderControl = std::make_unique<TaskBuilderControl>(
			nSuperclustersUpperbound, batch->particlesBondedToParticle, batch->pclustersBondedToPcluster, stream);

	ComputeMeanposAndRadiiForEachPclusterInEachSuperclusterKernel << <(batch->nSuperclusters + 31) / 32, 32, 0, stream >> > (
		batch->superClustersControl->scData,
		batch->superClustersControl->scMeta,
		batch->taskbuilderControl->contents.superclusterPositionSpheres,
		batch->nSuperclusters
		);
	cudaStreamSynchronize(stream);

	{
		batch->taskbuilderControl->Reset(stream);

		const uint32_t nGridnodes = batch->nGridnodes;

		ReserveInteractions << <
			dim3(nGridnodes, SuperClustersControl::maxClustersPerBlock, 1),
			256, 0, stream >> > (
				*batch->superClustersControl,
				boxSize,
				batch->taskbuilderControl->contents,
				batch->params.cutoff_nm,
				boxSizeF,
				Float3{ 1.0f } / boxSizeF
				);

		LIMA_UTILS::genericErrorCheck(stream, "ReserveInteractions");

		SortReserveInteractionsOutput << <batch->nSuperclusters, TaskBuilderControlContents::maxTasksPerSc, 0, stream >> > (
			batch->taskbuilderControl->contents
			);

		cudaStreamSynchronize(stream);
	}

	cudaMemsetAsync(batch->taskbuilderControl->contents.nResults + batch->nSuperclusters, 0, sizeof(int), stream);
	cudaMemsetAsync(batch->taskbuilderControl->contents.nInteractionsOwned + batch->nSuperclusters, 0, sizeof(int), stream);

	thrust::exclusive_scan(
		thrust::cuda::par.on(stream),
		batch->taskbuilderControl->contents.nResults,
		batch->taskbuilderControl->contents.nResults + batch->nSuperclusters + 1,
		batch->taskbuilderControl->contents.nResultsPrefixsum
	);

	thrust::exclusive_scan(
		thrust::cuda::par.on(stream),
		batch->taskbuilderControl->contents.nInteractionsOwned,
		batch->taskbuilderControl->contents.nInteractionsOwned + batch->nSuperclusters + 1,
		batch->taskbuilderControl->contents.nQueryBuffersPrefixsum
	);

	cudaMemcpyAsync(&batch->nResults, batch->taskbuilderControl->contents.nResultsPrefixsum + batch->nSuperclusters,
		sizeof(int), cudaMemcpyDeviceToHost, stream);
	int nTasks = batch->nSuperclusters;
	int nQueryBufferEntries = 0;
	cudaMemcpyAsync(&nQueryBufferEntries, batch->taskbuilderControl->contents.nQueryBuffersPrefixsum + batch->nSuperclusters,
		sizeof(int), cudaMemcpyDeviceToHost, stream);
	cudaStreamSynchronize(stream);

	batch->scscTasksDevice.Expand(nTasks, 1.2);
	batch->idsOfQuerySuperclustersDevice.Expand(nQueryBufferEntries, 1.2);
	batch->resultIndicesDevice.Expand(nQueryBufferEntries, 1.2);
	batch->noInteractionMatricesDevice.Expand(nQueryBufferEntries, 1.2);
	batch->scResultsDevice.Expand(batch->nResults, 1.2);

	BuildTasks << <(batch->nSuperclusters + 31) / 32, 32, 0, stream >> > (
		batch->taskbuilderControl->contents,
		batch->superClustersControl->scMeta,
		batch->nSuperclusters,
		batch->scscTasksDevice.Get(),
		batch->idsOfQuerySuperclustersDevice.Get(),
		batch->resultIndicesDevice.Get()
		);

	LIMA_UTILS::genericErrorCheck(stream, "BuildTasks");

	BuildNointeractionMatricesKernel << <batch->nSuperclusters, 16, 0, stream >> > (
		batch->superClustersControl->scMeta,
		batch->pClusterMetaDevice.Get(),
		batch->taskbuilderControl->contents,
		batch->noInteractionMatricesDevice.Get(),
		batch->nSuperclusters
		);

	LIMA_UTILS::genericErrorCheck(stream, "BuildNointeractionMatricesKernel");

	return true;
}
