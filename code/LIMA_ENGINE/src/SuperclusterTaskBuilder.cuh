#include "Engine.cuh"
#include "EngineBodies.cuh"

#include <set>
#include <execution>
#include "EngineCore.h"
#include "Neighborlists.cuh"


class InteractionToken {
	uint32_t data{};
public:
	constexpr InteractionToken(){}
	constexpr InteractionToken(int queryId, bool useNointeractionMatrix) {
		data = (uint32_t)queryId & 0x7FFFFFFF | ((uint32_t)useNointeractionMatrix << 31);
	}
	constexpr int GetQueryId() const { return (int) (data & 0x7FFFFFFF); }
	constexpr bool UseNointeractionMatrix() const { return (data >> 31); }
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
	TaskBuilderControl(int nSuperclustersUpperbound, const std::vector<ParticlesBondedToParticle>& particlesBondedToParticle,const std::vector<PclustersBondedToPcluster>& pclustersBondedToPcluster) 
		: nSuperclustersUpperbound(nSuperclustersUpperbound)
	{
		contents.particlesBondedToParticle = GenericCopyToDevice(particlesBondedToParticle);
		contents.pclustersBondedToPcluster = GenericCopyToDevice(pclustersBondedToPcluster);

		cudaMalloc(&contents.superclusterPositionSpheres, sizeof(std::array<float4, 4>) * nSuperclustersUpperbound);
		cudaMalloc(&contents.nInteractionsOwned, sizeof(int) * nSuperclustersUpperbound);
		cudaMalloc(&contents.interactionsOwned, sizeof(int) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		cudaMalloc(&contents.nInteractionsNonowned, sizeof(int) * nSuperclustersUpperbound);
		cudaMalloc(&contents.scIdsQueryNonowned, sizeof(int) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		cudaMalloc(&contents.nResults, sizeof(int) * nSuperclustersUpperbound);
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
			+ sizeof(int) * (nSuperclustersUpperbound + 1)
			+ sizeof(int) * (nSuperclustersUpperbound + 1)
			+ sizeof(int) * (nSuperclustersUpperbound + 1)) / (1024.0f * 1024.0f);

		//cudaMalloc(&contents.tasks, sizeof(ScScTask) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		//Reset(nSuperclustersUpperbound);
		Reset();
	}

	void Reset() {
		//cudaMemset(contents.interactionsOwned, 0xFF, sizeof(InteractionToken) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound); 
		//cudaMemset(contents.scIdsQueryNonowned, 0xFF, sizeof(int) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
	}

	~TaskBuilderControl() {
		cudaFree(contents.superclusterPositionSpheres);
		cudaFree(contents.nInteractionsOwned);
		cudaFree(contents.interactionsOwned);
		cudaFree(contents.nInteractionsNonowned);
		cudaFree(contents.scIdsQueryNonowned);
		cudaFree(contents.nResults);
		cudaFree(contents.nNointeractionmatricesOwned);

		cudaFree(contents.nResultsPrefixsum);
		cudaFree(contents.nNointeractionmatricesPrefixsum);
		cudaFree(contents.nTasksPrefixsum);

		//cudaFree(contents.tasks);
	}
};


__host__ __device__ inline bool ScAreBonded(const SuperClusterMeta& sc0, const SuperClusterMeta& sc1, const PclustersBondedToPcluster* const pclustersBondedToPcluster) {
	for (int i = 0; i < sc0.nUniquePcIds; i++) {
		for (int j = 0; j < sc1.nUniquePcIds; j++) {
			if (pclustersBondedToPcluster[sc0.uniquePclusterIds[i]].Contains(sc1.uniquePclusterIds[j])) {
				return true;
			}
		}
	}
	return false;
}

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
		const PData& pdata = superclusters[scId].pData[i];

		sum += pdata.position;
		cnt++;

		int nextPcId = i == 15 ? -1 : scMeta[scId]._pclusterIds[i + 1];
		if (nextPcId != pcId) {
			Float3 meanPos = sum * (1.f / static_cast<float>(cnt));
			float radius = 0;
			for (int ii = i - cnt + 1; ii <= i; ii++) {
				const PData& pData = superclusters[scId].pData[ii];
				radius = std::max(radius, (pData.position - meanPos).len());
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

// Fine pass
__host__ __device__ inline bool _DoesSuperclustersInteract(const SuperCluster* const scData, int scId0, int scId1, float cutoffDistance) {
	for (int i = 0; i < 16; i++) {
		const PData& p0 = scData[scId0].pData[i];
		if (!p0.Valid())
			break;
		for (int j = 0; j < 16; j++) {
			const PData& p1 = scData[scId1].pData[j];
			if (!p1.Valid())
				break;
			Float3 pos0 = p0.position;
			Float3 pos1 = p1.position;
			PeriodicBoundaryCondition::applyHyperposNM(pos0, pos1);
			float distance = (pos0 - pos1).len();
			if (distance <= cutoffDistance) {
				return true;
			}
		}
	}
	return false;
}

// Coarse pass
__host__ __device__ inline bool DoesSuperclustersInteract(const std::array<float4, 4>* const superclusterPositionSpheres, const SuperCluster* const scData, int scId0, int scId1, float cutoffDistance, Float3 boxSize) {
	// Ensure deterministic by comparing smallerId with largerId
	if (scId0 > scId1) {
		std::swap(scId0, scId1);
	}

	for (int posSphereId0 = 0; posSphereId0 < 4; posSphereId0++) {
		float4 p0 = superclusterPositionSpheres[scId0][posSphereId0];
		if (p0.w < 0)
			break;
		for (int posSphereId1 = 0; posSphereId1 < 4; posSphereId1++) {			
			float4 p1 = superclusterPositionSpheres[scId1][posSphereId1];
			if (p1.w < 0)
				break;

			Float3 pos0 = Float3{ p0 };
			Float3 pos1 = Float3{ p1 };
			PeriodicBoundaryCondition::applyHyperposNM(pos0, pos1);

			//float distance = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos0, pos1, boxSize, BoundaryConditionSelect::PBC);
			float distance = (pos0 - pos1).len(); // OPTIM Compute in squared-space instead..
			float radiusSum = p0.w + p1.w;

			if (distance <= cutoffDistance + radiusSum) {	// optim use LenSq
//				return true;
				return _DoesSuperclustersInteract(scData, scId0, scId1, cutoffDistance);
			}
		}
	}
	//return true;
	return false;
}


// gridDim = (nGridnodes, SuperClustersControl::maxClustersPerBlock, 1)
// blockDim = (3^3 * SuperClustersControl::maxClustersPerBlock, 1, 1) 
//template <typename BoundaryCondition>
__global__ void ReserveInteractions(SuperClustersControl scControl, Int3 boxSize, TaskBuilderControlContents tbContents, float cutoffNm) {

	NodeIndex nodeIndex = BoxGrid::Get3dIndex(blockIdx.x, boxSize);
	const int nodeId = BoxGrid::Get1dIndex(nodeIndex, boxSize);

	__shared__ int nOwnedInteractions;
	__shared__ int nNonownedInteractions;
	__shared__ int nNointeractionMatrices;


	// we can use gridDIm.y for this dimension, but then we'd have to use atomicAdds to the global counters at the bottom of this kernel, and reset those between runs. Which isnt great
	if (blockIdx.y >= scControl.nSuperclustersInBlocks[nodeId]) {
		return; // All threads in the block return
	}

	Float3 boxSizeF{ boxSize.x, boxSize.y, boxSize.z };
	const int scIndexInSelf = blockIdx.y;
	const int scIndexInQueryblock = threadIdx.x % SuperClustersControl::maxClustersPerBlock;
	const int scId = scControl.scIdsInBlocks[nodeId * SuperClustersControl::maxClustersPerBlock + blockIdx.y];

	//const bool isControlThread = threadIdx.x == 0 && threadIdx.y == 0;
	if (threadIdx.x == 0) {
		nOwnedInteractions = 0;
		nNonownedInteractions = 0;
		nNointeractionMatrices = 0;
	}
	if (threadIdx.x < TaskBuilderControlContents::maxTasksPerSc) {
		tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc + threadIdx.x] = InteractionToken(INT_MAX, true);
		tbContents.scIdsQueryNonowned[scId * TaskBuilderControlContents::maxTasksPerSc + threadIdx.x] = INT_MAX;
	}
	__syncthreads();

	NodeIndex targetBlockRelative = BoxGrid::Get3dIndex(threadIdx.x / SuperClustersControl::maxClustersPerBlock, Int3(3, 3, 3)) - Int3(1,1,1);
	NodeIndex targetBlock = PeriodicBoundaryCondition::applyBC(nodeIndex + targetBlockRelative, boxSize);
	int targetIndex = BoxGrid::Get1dIndex(targetBlock, boxSize);
	const int queryScId = scControl.scIdsInBlocks[targetIndex * SuperClustersControl::maxClustersPerBlock + scIndexInQueryblock];
	const bool validQuery = scIndexInQueryblock < scControl.nSuperclustersInBlocks[targetIndex];

	// TODO: Rethink the block-dimensions of this kernel, so more threads collaborate on the DoesSuperclustersInteract check
	// Beware that we might have to atomicadd the count in global memory instead then..
	if (validQuery && DoesSuperclustersInteract(tbContents.superclusterPositionSpheres, scControl.scData, scId, queryScId, cutoffNm, boxSizeF)) {
		const bool useNointeractionMatrix = scId == queryScId || ScAreBonded(scControl.scMeta[scId], scControl.scMeta[queryScId], tbContents.pclustersBondedToPcluster);
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
	__syncthreads();

	if constexpr (INDEXING_CHECKS) {
		if (threadIdx.x == 0 && (nOwnedInteractions + nNonownedInteractions > TaskBuilderControlContents::maxTasksPerSc))
			printf("Too many interactions for scId %d: %d owned + %d nonowned\n", scId, nOwnedInteractions, nNonownedInteractions);
	}
	 
	// TODO: Move this part into a separate kernel, so this kernel focuses on the DoesSuperclustersInteract hot-path
	// Now sort the elements to obtain deterministic results
	if (threadIdx.x < TaskBuilderControlContents::maxTasksPerSc) {
		LAL::Sort(&tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc], TaskBuilderControlContents::maxTasksPerSc, [](const InteractionToken& token) {
			return token.GetQueryId();
			});
		LAL::Sort(&tbContents.scIdsQueryNonowned[scId * TaskBuilderControlContents::maxTasksPerSc], TaskBuilderControlContents::maxTasksPerSc, [](const int& id) {
			return id;
			});
	}
	__syncthreads();
	//

	if (threadIdx.x == 0) {
		tbContents.nInteractionsOwned[scId] = nOwnedInteractions;
		tbContents.nInteractionsNonowned[scId] = nNonownedInteractions;
		tbContents.nNointeractionmatricesOwned[scId] = nNointeractionMatrices;
		tbContents.nResults[scId] = nOwnedInteractions + nNonownedInteractions;
	}
}


constexpr int IndexOfId(int* ids, int nIds, int idToFind) {
	for (int i = 0; i < nIds; i++) {
		if (ids[i] == idToFind)
			return i;
	}
	return -1;
}

__global__ void BuildTasks(TaskBuilderControlContents tbContents, SuperClusterMeta* const superClusterMeta, int nSuperclusters, ScScTask* const tasks) {
	const int scId = blockIdx.x * blockDim.x + threadIdx.x;
	if (scId >= nSuperclusters)
		return;
	
	int nMatricesUsed = 0;

	for (int i = 0; i < tbContents.nInteractionsOwned[scId]; i++) {
		ScScTask task;
		InteractionToken token = tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc + i];
		const int scIdQuery = token.GetQueryId();

		if (token.UseNointeractionMatrix()) {
			int nointeractionMatrixIdGlobal = tbContents.nNointeractionmatricesPrefixsum[scId] + nMatricesUsed;
			task.nointeractionMatrixIndex = nointeractionMatrixIdGlobal;
			nMatricesUsed++;
		}

		task.scIds[0] = scId;
		task.scIds[1] = scIdQuery;

		task.resultIndices[0] = tbContents.nResultsPrefixsum[scId] + i;
		{
			//int searchStartIndex = //
			int indexInQuery = scId != scIdQuery
				? IndexOfId(&tbContents.scIdsQueryNonowned[scIdQuery * TaskBuilderControlContents::maxTasksPerSc], tbContents.nInteractionsNonowned[token.GetQueryId()], scId)
				: task.resultIndices[0];

			 if (indexInQuery == -1) {
				 // This means that the query sc does not have an interaction with the sc of this task. This can only happen if the query sc has less interactions than the sc of this task, and thus we can be sure that the result index of the query sc is after the result index of this task, and thus we can safely set it to -1 to indicate that it should be ignored.
				 task.resultIndices[1] = -1;
				 printf("Something went very wrong!\n");
			 }
			 else {
				 const int queryNumOwnedTasks = tbContents.nInteractionsOwned[scIdQuery];
				 task.resultIndices[1] = tbContents.nResultsPrefixsum[scIdQuery] + queryNumOwnedTasks + indexInQuery;
			 }
		}

		tasks[tbContents.nTasksPrefixsum[scId] + i] = task;
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

	const Box& box = *simulation->box_host;
	Int3 boxSize = box.boxparams.boxSize;
	Float3 boxSizeF = simulation->box_host->boxparams.BoxSizeFloat();

	const int nSuperclustersUpperbound = nSuperclusters * 2;// simulation->box_host->persistentClusters.size(); // a little pessimistic

	if (!taskbuilderControl)
		taskbuilderControl = std::make_unique<TaskBuilderControl>(nSuperclustersUpperbound, box.particlesBondedToParticle, box.pclustersBondedToPcluster);

	cudaDeviceSynchronize();

	ComputeMeanposAndRadiiForEachPclusterInEachSuperclusterKernel<<<(nSuperclusters + 31) / 32, 32 >>>(superClustersControl->scData, superClustersControl->scMeta, taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters);
	cudaDeviceSynchronize();

	/*std::vector<std::array<float4, 4>> scps = GenericCopyToHost(taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters);
	std::vector<float> scpsFlat(reinterpret_cast<float*>(scps.data()), reinterpret_cast<float*>(scps.data()) + scps.size() * 4 * 4);
	DebugUtils::VerifyIdentical(scpsFlat, "scPositionsSpheres_" + std::to_string(simulation->getStep()));*/

	//auto iSpheres = GenericCopyToHost(taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters);
	//DebugUtils::VerifyIdentical(taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters, "SCPositionSpheres", simulation->getStep());

	{
		dim3 gridDim{ (uint32_t)boxSize.InnerProduct(), (uint32_t)SuperClustersControl::maxClustersPerBlock, 1u };
		dim3 blockDim{ 3 * 3 * 3 * SuperClustersControl::maxClustersPerBlock, 1, 1 };
		ReserveInteractions << <gridDim, blockDim >> > (*superClustersControl, boxSize, taskbuilderControl->contents, simulation->simparams_host.cutoff_nm);
		LIMA_UTILS::genericErrorCheck("ReserveInteractions");

		//DebugUtils::VerifyIdentical(taskbuilderControl->contents.nInteractionsOwned, nSuperclusters, "NInteractionsOwned", simulation->getStep());
	}

	/*auto ninteractionsOwnedHost = GenericCopyToHost(taskbuilderControl->contents.nInteractionsOwned, nSuperclusters);
	auto nSuperclusterPerNode = GenericCopyToHost(superClustersControl->nSuperclustersInBlocks, boxSize.InnerProduct());*/



	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nResults, taskbuilderControl->contents.nResults + nSuperclusters + 1, taskbuilderControl->contents.nResultsPrefixsum);
	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nInteractionsOwned, taskbuilderControl->contents.nInteractionsOwned + nSuperclusters + 1, taskbuilderControl->contents.nTasksPrefixsum);
	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nNointeractionmatricesOwned, taskbuilderControl->contents.nNointeractionmatricesOwned + nSuperclusters + 1, taskbuilderControl->contents.nNointeractionmatricesPrefixsum);
	cudaDeviceSynchronize();
	nResults = GenericCopyToHost(taskbuilderControl->contents.nResultsPrefixsum + nSuperclusters);
	nTasks = GenericCopyToHost(taskbuilderControl->contents.nTasksPrefixsum + nSuperclusters);

	scscTasksDevice.Expand(nTasks, 1.2);
	noInteractionMatricesDevice.Expand(nTasks, 1.2);
	scResultsDevice.Expand(nResults, 1.2);



	BuildTasks << <(nSuperclusters + 31) / 32, 32 >> > (taskbuilderControl->contents, superClustersControl->scMeta, nSuperclusters, scscTasksDevice.Get());
	BuildNointeractionMatricesKernel << <nSuperclusters , 16 >> >(superClustersControl->scMeta, pClusterMetaDevice, taskbuilderControl->contents, noInteractionMatricesDevice.Get(), nSuperclusters);
	cudaDeviceSynchronize();

	//auto resCounts = GenericCopyToHost(taskbuilderControl->contents.nResults, nSuperclustersUpperbound);
	
	//DebugUtils::VerifyIdentical(scscTasksDevice.Get(), nTasks, "SCSCTasks", simulation->getStep());
	
	cudaDeviceSynchronize();

	//auto tasksHost = GenericCopyToHost(scscTasksDevice, nSuperclustersUpperbound * TaskBuilderControlContents::maxTasksPerSc);
	//auto 


	//taskbuilderControl->Reset();

	return true;
}
