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
	static const int maxTasksPerSc = 256;// Need this as a powerof2 to be sort-able.. 8 * 27; // 8 sc/node 3^3 nodes..

	std::array<float4, 4>* superclusterPositionSpheres;

	// Static, made at creation
	ParticlesBondedToParticle* particlesBondedToParticle;
	PclustersBondedToPcluster* pclustersBondedToPcluster;

	// Written in reserve pass
	int* nInteractionsOwned;				// 1 per sc
	InteractionToken* interactionsOwned;	// maxTasksPerSc per sc
	int* nInteractionsNonowned;				// 1 per sc
	int* scIdsQueryNonowned;				// maxTasksPerSc per sc1
	int* nResults;							// 1 per sc, is simply nInteractionsOwned+nInteractionsNonowned
	int* nTasksOwned;						// 1 per sc, is simply ceil(nInteractionsOwned / 2)
	int* nNointeractionmatricesOwned;		// 1 per sc

	// Written in second pass
	int* nTasksPrefixsum;
	int* nResultsPrefixsum;
	int* nNointeractionmatricesPrefixsum;

	// Written in buildTasks pass
	//ScScTask* tasks;
};

// Keep on CPU
class TaskBuilderControl {

public:
	const int nSuperclustersUpperbound;
	TaskBuilderControlContents contents;

	TaskBuilderControl(const TaskBuilderControl&) = delete;
	TaskBuilderControl& operator=(const TaskBuilderControl&) = delete;
	TaskBuilderControl(int nSuperclustersUpperbound, const std::vector<ParticlesBondedToParticle>& particlesBondedToParticle, const std::vector<PclustersBondedToPcluster>& pclustersBondedToPcluster)
		: nSuperclustersUpperbound(nSuperclustersUpperbound)
	{
		contents.particlesBondedToParticle = GenericCopyToDevice(particlesBondedToParticle);
		contents.pclustersBondedToPcluster = GenericCopyToDevice(pclustersBondedToPcluster);

		cudaMalloc(&contents.superclusterPositionSpheres, sizeof(std::array<float4, 4>) * nSuperclustersUpperbound);
		cudaMalloc(&contents.nInteractionsOwned, sizeof(int) * nSuperclustersUpperbound);
		cudaMalloc(&contents.interactionsOwned, sizeof(InteractionToken) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		cudaMalloc(&contents.nInteractionsNonowned, sizeof(int) * nSuperclustersUpperbound);
		cudaMalloc(&contents.scIdsQueryNonowned, sizeof(int) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		cudaMalloc(&contents.nResults, sizeof(int) * nSuperclustersUpperbound);
		cudaMalloc(&contents.nTasksOwned, sizeof(int) * nSuperclustersUpperbound);
		cudaMalloc(&contents.nNointeractionmatricesOwned, sizeof(int) * nSuperclustersUpperbound);

		cudaMalloc(&contents.nResultsPrefixsum, sizeof(int) * (nSuperclustersUpperbound + 1));
		cudaMalloc(&contents.nNointeractionmatricesPrefixsum, sizeof(int) * (nSuperclustersUpperbound + 1));
		cudaMalloc(&contents.nTasksPrefixsum, sizeof(int) * (nSuperclustersUpperbound + 1));

		const size_t totalMemUsageMB = (sizeof(std::array<float4, 4>) * nSuperclustersUpperbound
			+ sizeof(int) * nSuperclustersUpperbound
			+ sizeof(InteractionToken) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound
			+ sizeof(int) * nSuperclustersUpperbound
			+ sizeof(int) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound
			+ sizeof(int) * nSuperclustersUpperbound
			+ sizeof(int) * nSuperclustersUpperbound
			+ sizeof(int) * nSuperclustersUpperbound
			+ sizeof(int) * (nSuperclustersUpperbound + 1)
			+ sizeof(int) * (nSuperclustersUpperbound + 1)
			+ sizeof(int) * (nSuperclustersUpperbound + 1)) / (1024.0f * 1024.0f);

		//cudaMalloc(&contents.tasks, sizeof(ScScTask) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		//Reset(nSuperclustersUpperbound);
		Reset();
	}

	void Reset() {
		cudaMemset(contents.interactionsOwned, 0xFF, sizeof(InteractionToken) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		//cudaMemset(contents.scIdsQueryNonowned, 0xFF, sizeof(int) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		thrust::fill_n(thrust::device, contents.scIdsQueryNonowned, TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound, INT_MAX);
	}

	~TaskBuilderControl() {
		cudaFree(contents.superclusterPositionSpheres);
		cudaFree(contents.nInteractionsOwned);
		cudaFree(contents.interactionsOwned);
		cudaFree(contents.nInteractionsNonowned);
		cudaFree(contents.scIdsQueryNonowned);
		cudaFree(contents.nResults);
		cudaFree(contents.nTasksOwned);
		cudaFree(contents.nNointeractionmatricesOwned);

		cudaFree(contents.nResultsPrefixsum);
		cudaFree(contents.nNointeractionmatricesPrefixsum);
		cudaFree(contents.nTasksPrefixsum);

		//cudaFree(contents.tasks);
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
			PeriodicBoundaryCondition::applyHyperposNM(pos0, pos1);
			//PeriodicBoundaryCondition::ApplyHyperpos(pos0, pos1, boxSize, boxSizeInv);

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

			PeriodicBoundaryCondition::applyHyperposNM(pos0, pos1);
			//PeriodicBoundaryCondition::ApplyHyperpos(pos0, pos1, boxSize, boxSizeInv);
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

	NodeIndex nodeIndex = BoxGrid::Get3dIndex(blockIdx.x, boxSize);
	const int nodeId = BoxGrid::Get1dIndex(nodeIndex, boxSize);

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

		int targetIndex = BoxGrid::Get1dIndex(targetBlock, boxSize);

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
		tbContents.nNointeractionmatricesOwned[scId] = nNointeractionMatrices;

		const int nOwnedTasks = (nOwnedInteractions + 1) / 2;

		tbContents.nResults[scId] = nOwnedTasks + nNonownedInteractions;
		tbContents.nTasksOwned[scId] = nOwnedTasks;
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
	//int searchStartIndex = //
	int indexInQuery = scId != scIdQuery
		? IndexOfId(&tbContents.scIdsQueryNonowned[scIdQuery * TaskBuilderControlContents::maxTasksPerSc], tbContents.nInteractionsNonowned[scIdQuery], scId)
		: selfResultIndex - tbContents.nResultsPrefixsum[scId];

	if (indexInQuery == -1) {
		// This means that the query sc does not have an interaction with the sc of this task. This can only happen if the query sc has less interactions than the sc of this task, 
		// and thus we can be sure that the result index of the query sc is after the result index of this task, and thus we can safely set it to -1 to indicate that it should be ignored.
		printf("Illegal query index for scId %d querying scId %d. Check %d ids\n", scId, scIdQuery, tbContents.nInteractionsNonowned[scIdQuery]);
		//printf("Something went very wrong!\n");
		return -1;
	}
	else {
		const int queryNumOwnedResults = tbContents.nTasksOwned[scIdQuery];
		return tbContents.nResultsPrefixsum[scIdQuery] + queryNumOwnedResults + indexInQuery;
	}
}

__global__ void BuildTasks(TaskBuilderControlContents tbContents, const SuperCluster* const superClusters, SuperClusterMeta* const superClusterMeta, int nSuperclusters, ScScTask* const tasks, Float3 boxSize, Float3 boxSizeInv) {
	const int scId = blockIdx.x * blockDim.x + threadIdx.x;
	if (scId >= nSuperclusters)
		return;

	int nMatricesUsed = 0;

	for (int i = 0; i < tbContents.nInteractionsOwned[scId]; i += 2) {
		ScScTask task;
		task.scIds[0] = scId;
		task.scIds[1] = -1;
		task.scIds[2] = -1;
		task.resultIndices[0] = tbContents.nResultsPrefixsum[scId] + i / 2;
		task.resultIndices[1] = -1;
		task.resultIndices[2] = -1;
		task.nointeractionMatrixIndex[0] = -1;
		task.nointeractionMatrixIndex[1] = -1;

		{
			InteractionToken token = tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc + i];
			const int scIdQuery = token.GetQueryId();

			if (token.UseNointeractionMatrix()) {
				int nointeractionMatrixIdGlobal = tbContents.nNointeractionmatricesPrefixsum[scId] + nMatricesUsed;
				task.nointeractionMatrixIndex[0] = nointeractionMatrixIdGlobal;
				nMatricesUsed++;
			}

			task.scIds[1] = scIdQuery;
			task.resultIndices[1] = GetResultIndexOfQuery(tbContents, scId, scIdQuery, task.resultIndices[0]);
		}

		if (i + 1 < tbContents.nInteractionsOwned[scId]) {
			InteractionToken token = tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc + i + 1];
			const int scIdQuery = token.GetQueryId();

			if (token.UseNointeractionMatrix()) {
				int nointeractionMatrixIdGlobal = tbContents.nNointeractionmatricesPrefixsum[scId] + nMatricesUsed;
				task.nointeractionMatrixIndex[1] = nointeractionMatrixIdGlobal;
				nMatricesUsed++;
			}

			task.scIds[2] = scIdQuery;
			task.resultIndices[2] = GetResultIndexOfQuery(tbContents, scId, scIdQuery, task.resultIndices[0]);
		}

		/*const Float3 sc0Pos0 = superClusters[scId].Position(0);
		const Float3 sc1Pos0 = superClusters[scIdQuery].Position(0);
		task.sc1Translation = PeriodicBoundaryCondition::GetHyperposTranslation(sc0Pos0, sc1Pos0, boxSize, boxSizeInv);*/

		tasks[tbContents.nTasksPrefixsum[scId] + i / 2] = task;
	}


	superClusterMeta[scId].resultsStartIndex = tbContents.nResultsPrefixsum[scId];
	superClusterMeta[scId].nResults = tbContents.nResults[scId];
}



//// gridDim = (nSuperclusters, 1, 1)
//// blockDim = (16, 1, 1)
__global__ void BuildNointeractionMatricesKernel(const SuperClusterMeta* const superClusterMetas, const PersistentClusterMeta* const pClustersMeta,
	TaskBuilderControlContents tbContents, BoolMatrix16x16* const nointeractionMatrices, int nSuperclusters) {

	const int scId = blockIdx.x;

	//std::array<int, 16> particleIdsSelf = GetParticleIdsOfSuperCluster(pClustersMeta, superClusterMetas[scId]);
	int matrixCount = 0;

	for (int i = 0; i < tbContents.nInteractionsOwned[scId]; i++) {
		InteractionToken token = tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc + i];

		if (!token.UseNointeractionMatrix()) {
			continue;
		}

		const int scIdQuery = token.GetQueryId();

		//std::array<int, 16> particleIdsQuery = GetParticleIdsOfSuperCluster(pClustersMeta, superClusterMetas[scIdQuery]);
		const bool isSelfInteractionTask = scId == scIdQuery;

		const int row = threadIdx.x;
		uint16_t rowData = 0;
		for (int col = 0; col < 16; ++col) {
			int pidSelf = superClusterMetas[scId].globalParticleIds[row];
			int pidQuery = superClusterMetas[scIdQuery].globalParticleIds[col];
			if (pidSelf == -1 || pidQuery == -1)
				continue;

			bool noInteraction = tbContents.particlesBondedToParticle[pidSelf].Contains(pidQuery);
			if (isSelfInteractionTask && row == col) {
				noInteraction = true;
			}
			if (noInteraction)
				BoolMatrix16x16::SetValueInRow(col, rowData);
		}




		const int matrixIndex = tbContents.nNointeractionmatricesPrefixsum[scId] + matrixCount;
		nointeractionMatrices[matrixIndex].SetRow(row, rowData);
		matrixCount++;
	}
}















bool Engine::MakeSuperClusterTasksGPU() {
	if (nSuperclusters == 0)
		return true;

	const Box& box = *simulation->box;
	Int3 boxSize = box.boxparams.boxSize;
	Float3 boxSizeF = simulation->box->boxparams.BoxSizeFloat();

	const int nSuperclustersUpperbound = nSuperclusters * 2;// simulation->box->persistentClusters.size(); // a little pessimistic

	if (!taskbuilderControl)
		taskbuilderControl = std::make_unique<TaskBuilderControl>(nSuperclustersUpperbound, box.particlesBondedToParticle, box.pclustersBondedToPcluster);

	cudaDeviceSynchronize();

	ComputeMeanposAndRadiiForEachPclusterInEachSuperclusterKernel << <(nSuperclusters + 31) / 32, 32 >> > (superClustersControl->scData, superClustersControl->scMeta, taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters);
	cudaDeviceSynchronize();

	/*std::vector<std::array<float4, 4>> scps = GenericCopyToHost(taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters);
	std::vector<float> scpsFlat(reinterpret_cast<float*>(scps.data()), reinterpret_cast<float*>(scps.data()) + scps.size() * 4 * 4);
	DebugUtils::VerifyIdentical(scpsFlat, "scPositionsSpheres_" + std::to_string(simulation->getStep()));*/

	//auto iSpheres = GenericCopyToHost(taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters);
	//DebugUtils::VerifyIdentical(taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters, "SCPositionSpheres", simulation->getStep());

	{
		taskbuilderControl->Reset();
		cudaDeviceSynchronize();

		dim3 gridDim{ (uint32_t)boxSize.InnerProduct(), (uint32_t)SuperClustersControl::maxClustersPerBlock, 1u };
		dim3 blockDim{ 3 * 3 * 3 * SuperClustersControl::maxClustersPerBlock, 1, 1 };
		uint32_t nGridnodes = boxSize.InnerProduct();
		//ReserveInteractions << <gridDim, blockDim >> > (*superClustersControl, boxSize, taskbuilderControl->contents, simulation->simParams.cutoff_nm);
		ReserveInteractions << <
			dim3(nGridnodes, SuperClustersControl::maxClustersPerBlock, 1),
			256 >> > (*superClustersControl, boxSize, taskbuilderControl->contents, simulation->simParams.cutoff_nm, boxSizeF, Float3{ 1.0f } / boxSizeF);
		LIMA_UTILS::genericErrorCheck("ReserveInteractions");

		cudaDeviceSynchronize();
		SortReserveInteractionsOutput << <nSuperclusters, TaskBuilderControlContents::maxTasksPerSc >> > (taskbuilderControl->contents);
		cudaDeviceSynchronize();
		//DebugUtils::VerifyIdentical(taskbuilderControl->contents.nInteractionsOwned, nSuperclusters, "NInteractionsOwned", simulation->getStep());
	}

	/*auto ninteractionsOwnedHost = GenericCopyToHost(taskbuilderControl->contents.nInteractionsOwned, nSuperclusters);
	auto nSuperclusterPerNode = GenericCopyToHost(superClustersControl->nSuperclustersInBlocks, boxSize.InnerProduct());*/



	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nResults, taskbuilderControl->contents.nResults + nSuperclusters + 1, taskbuilderControl->contents.nResultsPrefixsum);
	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nTasksOwned, taskbuilderControl->contents.nTasksOwned + nSuperclusters + 1, taskbuilderControl->contents.nTasksPrefixsum);
	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nNointeractionmatricesOwned, taskbuilderControl->contents.nNointeractionmatricesOwned + nSuperclusters + 1, taskbuilderControl->contents.nNointeractionmatricesPrefixsum);
	cudaDeviceSynchronize();
	nResults = GenericCopyToHost(taskbuilderControl->contents.nResultsPrefixsum + nSuperclusters);
	nTasks = GenericCopyToHost(taskbuilderControl->contents.nTasksPrefixsum + nSuperclusters);
	const int nNointeractionMatrices = GenericCopyToHost(taskbuilderControl->contents.nNointeractionmatricesPrefixsum + nSuperclusters);

	scscTasksDevice.Expand(nTasks, 1.2);
	noInteractionMatricesDevice.Expand(nNointeractionMatrices, 1.2);
	scResultsDevice.Expand(nResults, 1.2);



	BuildTasks << <(nSuperclusters + 31) / 32, 32 >> > (taskbuilderControl->contents, superClustersControl->scData, superClustersControl->scMeta, nSuperclusters, scscTasksDevice.Get(), boxSizeF, boxSizeF.Inv());
	BuildNointeractionMatricesKernel << <nSuperclusters, 16 >> > (superClustersControl->scMeta, pClusterMetaDevice.Get(), taskbuilderControl->contents, noInteractionMatricesDevice.Get(), nSuperclusters);
	cudaDeviceSynchronize();

	//auto resCounts = GenericCopyToHost(taskbuilderControl->contents.nResults, nSuperclustersUpperbound);

	//DebugUtils::VerifyIdentical(scscTasksDevice.Get(), nTasks, "SCSCTasks", simulation->getStep());

	cudaDeviceSynchronize();

	//auto tasksHost = GenericCopyToHost(scscTasksDevice, nSuperclustersUpperbound * TaskBuilderControlContents::maxTasksPerSc);
	//auto 


	//taskbuilderControl->Reset();

	return true;
}