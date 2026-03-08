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
	static const int maxTasksPerSc = 8*27; // 8 sc/node 3^3 nodes..

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
	TaskBuilderControlContents contents;

	TaskBuilderControl(const TaskBuilderControl&) = delete;
	TaskBuilderControl& operator=(const TaskBuilderControl&) = delete;
	TaskBuilderControl(int nSuperclustersUpperbound, const std::vector<ParticlesBondedToParticle>& particlesBondedToParticle,const std::vector<PclustersBondedToPcluster>& pclustersBondedToPcluster) {
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

		//cudaMalloc(&contents.tasks, sizeof(ScScTask) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		//Reset(nSuperclustersUpperbound);
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
	for (int i = 0; i < SuperCluster::nPclusters; i++) {
		if (sc0.pclusterIds[i] == -1)
			break;
		for (int j = 0; j < SuperCluster::nPclusters; j++) {
			if (sc1.pclusterIds[j] == -1)
				break;
			if (pclustersBondedToPcluster[sc0.pclusterIds[i]].Contains(sc1.pclusterIds[j])) {
				return true;
			}
		}
	}
	return false;
}

constexpr std::array<int, 16> GetParticleIdsOfSuperCluster(const PersistentClusterMeta* const pClusterMeta, const SuperClusterMeta& scMeta) {
	std::array<int, 16> particleIds{};
	for (auto& e : particleIds) { e = -1; }
	int cnt = 0;
	for (auto pcId : scMeta.pclusterIds) {
		if (pcId == -1) break;
		for (int particleId : pClusterMeta[pcId].particleIdsGlobal) {
			particleIds[cnt++] = particleId;
		}
	}
	return particleIds;
}

struct ReservedTask {
	int queryScId;
	int resultIndexRelativeSelf = -1;
	int resultIndexRelativeQuery = -1;
	int nointeractionMatrixIndexRelative = -1;
};

template <typename T>
std::vector<size_t> ExlusivePrefixsum(const std::vector<T>& counts) {
	static_assert(std::is_integral<T>::value, "ExlusivePrefixsum only supports integral types");
	std::vector<size_t> prefixsum(counts.size());
	//std::exclusive_scan(std::execution::par, counts.begin(), counts.end(), prefixsum.begin(), 0);
	std::exclusive_scan(counts.begin(), counts.end(), prefixsum.begin(), size_t{ 0 });
	return prefixsum;
}
template <typename T>
std::vector<size_t> ExlusivePrefixsum(const std::vector<std::vector<T>>& sizes) {
	std::vector<size_t> prefixsum(sizes.size());
	std::transform_exclusive_scan(std::execution::par, sizes.begin(), sizes.end(), prefixsum.begin(), size_t{ 0 }, std::plus<>{},
		[](const std::vector<T>& v) { return v.size(); }
	);
	return prefixsum;
}

std::vector<BoolMatrix16x16> BuildNointeractionMatrices(const std::vector<SuperClusterMeta>& superClusterMetas, const std::vector<int>& nBondedmatricesReserved, const std::vector<size_t>& nBondedMatricesPrefixsum,
	const std::vector<PersistentClusterMeta>& pClustersMeta, const std::vector<std::vector<ReservedTask>>& workPerSc, const std::vector<ParticlesBondedToParticle>& particlesBondedToParticle) {
	const size_t numNointeractionMatricesTotal = nBondedMatricesPrefixsum.back() + nBondedmatricesReserved.back();
	std::vector<BoolMatrix16x16> nointeractionMatrices(numNointeractionMatricesTotal);

	// Build all the nointeractionMatrices
	for (int scId = 0; scId < superClusterMetas.size(); ++scId) {
		std::array<int, 16> particleIdsSelf = GetParticleIdsOfSuperCluster(pClustersMeta.data(), superClusterMetas[scId]);
		for (int i = 0; i < workPerSc[scId].size(); i++) {
			if (workPerSc[scId][i].nointeractionMatrixIndexRelative == -1)
				continue;

			const int queryScId = workPerSc[scId][i].queryScId;
			BoolMatrix16x16 nointeractionMatrix{};

			std::array<int, 16> particleIdsQuery = GetParticleIdsOfSuperCluster(pClustersMeta.data(), superClusterMetas[queryScId]);
			const bool isSelfInteractionTask = scId == queryScId;

			for (int col = 0; col < 16; ++col) {
				for (int row = 0; row < 16; ++row) {
					if (particleIdsSelf[row] == -1)
						continue;
					int pidSelf = particleIdsSelf[row];
					int pidQuery = particleIdsQuery[col];
					if (pidSelf == pidQuery && pidSelf == 0)
						int a = 0;

					bool noInteraction = particlesBondedToParticle[particleIdsSelf[row]].Contains(particleIdsQuery[col]);
					if (isSelfInteractionTask && row == col) {
						noInteraction = true;
					}
					//		noInteraction = true;
					nointeractionMatrix.Set(row, col, noInteraction);
				}
			}


			const int matrixIndex = workPerSc[scId][i].nointeractionMatrixIndexRelative + nBondedMatricesPrefixsum[scId];
			nointeractionMatrices[matrixIndex] = nointeractionMatrix;
		}
	}

	return nointeractionMatrices;
}

//float MinDistanceBetweenPclustersInSupercluster(const SuperCluster& sc0, const SuperCluster& sc1, const Float3& boxSize) {
//	float minDist = FLT_MAX;
//	for (int pcid0 = 0; pcid0 < 4; pcid0++) {
//		if (!sc0.pData->Valid())
//			continue;
//		for (int pcid1 = 0; pcid1 < 4; pcid1++) {
//			if (!sc1.pData->Valid())
//				continue;
//			const float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(sc0.pData[pcid0].position, sc1.pData[pcid1].position, boxSize, BoundaryConditionSelect::PBC);
//			if (dist < minDist) {
//				minDist = dist;
//			}
//		}
//	}
//	return minDist;
//}

std::vector<std::array<float4, 4>> ComputeMeanposAndRadiiForEachPclusterInEachSupercluster(const std::vector<SuperCluster>& superclusters) {
	std::vector<std::array<float4, 4>> out(superclusters.size());

	// Debugging
	float maxRadius = 0;
	float maxIntraScDistance = 0;

	for (int scId = 0; scId < superclusters.size(); scId++) {
		for (int pcid = 0; pcid < 4; pcid++) {
			Float3 sum{};
			int cnt = 0;
			for (int pid = 0; pid < 4; pid++) {
				const PData& pData = superclusters[scId].pData[pcid * 4 + pid];
				if (pData.Valid()) {
					sum += pData.position;
					cnt++;
				}
			}

			const Float3 meanPos = sum * (1.0f / static_cast<float>(cnt));
			float radius = 0;
			for (int pid = 0; pid < cnt; pid++) {
				const PData& pData = superclusters[scId].pData[pcid * 4 + pid];
				radius = std::max(radius, (pData.position - meanPos).len());
			}
			out[scId][pcid] = float4{ meanPos.x, meanPos.y, meanPos.z, radius };

			// Debug
			maxRadius = std::max(maxRadius, radius);
			if (pcid != 0)
				maxIntraScDistance = std::max(maxIntraScDistance, (meanPos - Float3{ out[scId][pcid - 1] }).len());
			if (radius > .8f || maxIntraScDistance > 1.2f)
				int a = 0;
			//
		}
	}

	return out;
}
__global__ void ComputeMeanposAndRadiiForEachPclusterInEachSuperclusterKernel(const SuperCluster* const superclusters, std::array<float4, 4>* const out, int nSuperclusters) {
	const int scId = blockIdx.x * blockDim.x + threadIdx.x;

	if (scId >= nSuperclusters)
		return;

	for (int pcid = 0; pcid < 4; pcid++) { // optim: 1 thread per pcid
		Float3 sum{};
		int cnt = 0;
		for (int pid = 0; pid < 4; pid++) {
			const PData& pData = superclusters[scId].pData[pcid * 4 + pid];
			if (pData.Valid()) {
				sum += pData.position;
				cnt++;
			}
		}

		const Float3 meanPos = sum * (1.0f / static_cast<float>(cnt));
		float radius = 0;
		for (int pid = 0; pid < cnt; pid++) {
			const PData& pData = superclusters[scId].pData[pcid * 4 + pid];
			radius = std::max(radius, (pData.position - meanPos).len());
		}

		out[scId][pcid] = float4{ meanPos.x, meanPos.y, meanPos.z, radius };
	}
}



__host__ __device__ inline bool DoesSuperclustersInteract(const std::array<float4, 4>* const superclusterPositionSpheres, int scId0, int scId1, float cutoffDistance, Float3 boxSize) {
	// Ensure deterministic by comparing smallerId with largerId
	if (scId0 > scId1) {
		std::swap(scId0, scId1);
	}

	for (int pcid0 = 0; pcid0 < 4; pcid0++) {
		for (int pcid1 = 0; pcid1 < 4; pcid1++) {
			float4 p0 = superclusterPositionSpheres[scId0][pcid0];
			float4 p1 = superclusterPositionSpheres[scId1][pcid1];

			Float3 pos0 = Float3{ p0 };
			Float3 pos1 = Float3{ p1 };
			PeriodicBoundaryCondition::applyHyperposNM(pos0, pos1);

			//float distance = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos0, pos1, boxSize, BoundaryConditionSelect::PBC);
			float distance = (pos0 - pos1).len(); // OPTIM Compute in squared-space instead..
			float radiusSum = p0.w + p1.w;

			if (distance + radiusSum <= cutoffDistance) {	// optim use LenSq
				return true;
			}
		}
	}
	//return true;
	return false;
}




// gridDim = (nGridnodes, 1, 1)
// blockDim = (SuperClustersControl::maxClustersPerBlock, 1, 1) // THis is a silly dimension..
//template <typename BoundaryCondition>
__global__ void ReserveInteractions(SuperClustersControl scControl, Int3 boxSize, TaskBuilderControlContents tbContents, float cutoffNm) {
	
	NodeIndex nodeIndex = BoxGrid::Get3dIndex(blockIdx.x, boxSize);
	const int nodeId = BoxGrid::Get1dIndex(nodeIndex, boxSize);
	Float3 boxSizeF{ boxSize.x, boxSize.y, boxSize.z };


	if (threadIdx.x < scControl.nSuperclustersInBlocks[nodeId]) {
		const int scId = scControl.scIdsInBlocks[nodeId * SuperClustersControl::maxClustersPerBlock + threadIdx.x];
		int nOwnedInteractions = 0;
		int nNonownedInteractions = 0;
		int nNointeractionMatrices = 0;

		for (int zOff = -1; zOff <= 1; zOff++) {
			for (int yOff = -1; yOff <= 1; yOff++) {
				for (int xOff = -1; xOff <= 1; xOff++) {
					NodeIndex targetBlock = PeriodicBoundaryCondition::applyBC(nodeIndex + NodeIndex(xOff, yOff, zOff), boxSize);	// TODO: dont hardcode BC
					int targetIndex = BoxGrid::Get1dIndex(targetBlock, boxSize);


					for (int ii = 0; ii < scControl.nSuperclustersInBlocks[targetIndex]; ii++) {
						const int queryScId = scControl.scIdsInBlocks[targetIndex * SuperClustersControl::maxClustersPerBlock + ii];
						// TODO IMPORTANT: This algo requieres that scIdsInBlocks are sorted. But i think thats an implicit consequence of how ids are made already.. Maybe add a check tho in safe mode.

						if (DoesSuperclustersInteract(tbContents.superclusterPositionSpheres, scId, queryScId, cutoffNm, boxSizeF)) {
							const bool useNointeractionMatrix = scId == queryScId || ScAreBonded(scControl.scMeta[scId], scControl.scMeta[queryScId], tbContents.pclustersBondedToPcluster);

							//printf("Interaction %d %d nodeId %d\n", scId, queryScId, nodeId);

							if (scId <= queryScId) {
								tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc + nOwnedInteractions] = InteractionToken(queryScId, useNointeractionMatrix);
								nOwnedInteractions++;
								nNointeractionMatrices += useNointeractionMatrix ? 1 : 0;
							}
							else {
								tbContents.scIdsQueryNonowned[scId * TaskBuilderControlContents::maxTasksPerSc + nNonownedInteractions] = queryScId;
								nNonownedInteractions++;
							}
						}
					}
				}
			}
		}

		if constexpr (INDEXING_CHECKS) {
			if (nOwnedInteractions + nNonownedInteractions > TaskBuilderControlContents::maxTasksPerSc)
				printf("Too many interactions for scId %d: %d owned + %d nonowned\n", scId, nOwnedInteractions, nNonownedInteractions);
		}

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

	const int scId = blockIdx.x * blockDim.x + threadIdx.x;
	if (scId >= nSuperclusters)
		return;

	std::array<int, 16> particleIdsSelf = GetParticleIdsOfSuperCluster(pClustersMeta, superClusterMetas[scId]);
	int matrixCount = 0;

	for (int i = 0; i < tbContents.nInteractionsOwned[scId]; i++) {
		InteractionToken token = tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc + i];

		if (!token.UseNointeractionMatrix()) {
			continue;
		}

		const int scIdQuery = token.GetQueryId();

		BoolMatrix16x16 nointeractionMatrix{};
		nointeractionMatrix.Clear();

		std::array<int, 16> particleIdsQuery = GetParticleIdsOfSuperCluster(pClustersMeta, superClusterMetas[scIdQuery]);
		const bool isSelfInteractionTask = scId == scIdQuery;

		for (int col = 0; col < 16; ++col) {
			for (int row = 0; row < 16; ++row) {
				if (particleIdsSelf[row] == -1)
					continue;
				int pidSelf = particleIdsSelf[row];
				int pidQuery = particleIdsQuery[col];
				if (pidSelf == pidQuery && pidSelf == 0)
					int a = 0;

				bool noInteraction = tbContents.particlesBondedToParticle[particleIdsSelf[row]].Contains(particleIdsQuery[col]);
				if (isSelfInteractionTask && row == col) {
					noInteraction = true;
				}
				//		noInteraction = true;
				nointeractionMatrix.Set(row, col, noInteraction);
			}
		}


		const int matrixIndex = tbContents.nNointeractionmatricesPrefixsum[scId] + matrixCount;
		nointeractionMatrices[matrixIndex] = nointeractionMatrix;
		matrixCount++;
	}
}



















bool Engine::MakeSuperClusterTasksGPU() {
	if (nSuperclusters == 0)
		return true;

	const Box& box = *simulation->box_host;
	Int3 boxSize = box.boxparams.boxSize;
	Float3 boxSizeF = simulation->box_host->boxparams.BoxSizeFloat();

	const int nSuperclustersUpperbound = simulation->box_host->persistentClusters.size(); // a little pessimistic

	if (!taskbuilderControl)
		taskbuilderControl = std::make_unique<TaskBuilderControl>(nSuperclustersUpperbound, box.particlesBondedToParticle, box.pclustersBondedToPcluster);
	if (!scscTasksDevice) {
		int maxTasks = TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound;
		cudaMalloc(&scscTasksDevice, sizeof(ScScTask) * maxTasks); // THis is probably too many...
		cudaMalloc(&noInteractionMatricesDevice, sizeof(BoolMatrix16x16) * maxTasks);
		cudaMalloc(&scResultsDevice, sizeof(SCResult) * maxTasks * 2);
	}


	ComputeMeanposAndRadiiForEachPclusterInEachSuperclusterKernel<<<(nSuperclusters + 31) / 32, 32 >>>(superClustersControl->scData, taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters);


	//auto iSpheres = GenericCopyToHost(taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters);

	ReserveInteractions<< <boxSize.InnerProduct(), SuperClustersControl::maxClustersPerBlock >> > (*superClustersControl, boxSize, taskbuilderControl->contents, simulation->simparams_host.cutoff_nm);


	/*auto ninteractionsOwnedHost = GenericCopyToHost(taskbuilderControl->contents.nInteractionsOwned, nSuperclusters);
	auto nSuperclusterPerNode = GenericCopyToHost(superClustersControl->nSuperclustersInBlocks, boxSize.InnerProduct());*/



	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nResults, taskbuilderControl->contents.nResults + nSuperclusters + 1, taskbuilderControl->contents.nResultsPrefixsum);
	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nInteractionsOwned, taskbuilderControl->contents.nInteractionsOwned + nSuperclusters + 1, taskbuilderControl->contents.nTasksPrefixsum);
	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nNointeractionmatricesOwned, taskbuilderControl->contents.nNointeractionmatricesOwned + nSuperclusters + 1, taskbuilderControl->contents.nNointeractionmatricesPrefixsum);


	BuildTasks << <(nSuperclusters + 31) / 32, 32 >> > (taskbuilderControl->contents, superClustersControl->scMeta, nSuperclusters, scscTasksDevice);
	BuildNointeractionMatricesKernel << <(nSuperclusters + 31) / 32, 32 >> >(superClustersControl->scMeta, pClusterMetaDevice, taskbuilderControl->contents, noInteractionMatricesDevice, nSuperclusters);


	//auto resCounts = GenericCopyToHost(taskbuilderControl->contents.nResults, nSuperclustersUpperbound);

	nResults = GenericCopyToHost(taskbuilderControl->contents.nResultsPrefixsum + nSuperclusters);
	nTasks = GenericCopyToHost(taskbuilderControl->contents.nTasksPrefixsum + nSuperclusters);

	//auto tasksHost = GenericCopyToHost(scscTasksDevice, nSuperclustersUpperbound * TaskBuilderControlContents::maxTasksPerSc);
	//auto 

	return true;










	const std::vector<PersistentClusterMeta>& pClustersMeta = simulation->box_host->persistentClustersMetadata;
	const std::vector<SuperCluster> superClusters = GenericCopyToHost(superClustersControl->scData, nSuperclusters);
	std::vector<SuperClusterMeta> superClusterMetas = GenericCopyToHost(superClustersControl->scMeta, nSuperclusters);

	const std::vector<std::array<float4, 4>> superclusterPositionSpheres = ComputeMeanposAndRadiiForEachPclusterInEachSupercluster(superClusters);

	std::vector<std::vector<ReservedTask>> workPerSc(superClusters.size());
	std::vector<int> nResultsReserved(superClusters.size(), 0);
	std::vector<int> nBondedmatricesReserved(superClusters.size(), 0);

	for (int scId = 0; scId < superClusterMetas.size(); ++scId) {
		for (int queryScId = scId; queryScId < superClusterMetas.size(); ++queryScId) {

			if (DoesSuperclustersInteract(superclusterPositionSpheres.data(), scId, queryScId, simulation->simparams_host.cutoff_nm, boxSizeF)) {
				const bool useNointeractionMatrix = scId == queryScId || ScAreBonded(superClusterMetas[scId], superClusterMetas[queryScId], box.pclustersBondedToPcluster.data());

				workPerSc[scId].emplace_back(ReservedTask{
					queryScId,
					nResultsReserved[scId],
					scId != queryScId ? nResultsReserved[queryScId] : nResultsReserved[queryScId],
					useNointeractionMatrix ? nBondedmatricesReserved[scId] : -1
					});

				nResultsReserved[scId]++;
				if (scId != queryScId)
					nResultsReserved[queryScId]++;
				if (useNointeractionMatrix) {
					nBondedmatricesReserved[scId]++;
				}
			}
		}
	}

	// Make prefixsums
	const std::vector<size_t> nResultsPrefixsum = ExlusivePrefixsum(nResultsReserved);
	const std::vector<size_t> nBondedMatricesPrefixsum = ExlusivePrefixsum(nBondedmatricesReserved);
	const std::vector<size_t> nTasksPrefixsum = ExlusivePrefixsum(workPerSc);

	const size_t numTasksTotal = nTasksPrefixsum.back() + workPerSc.back().size();
	std::vector<ScScTask> tasks(numTasksTotal);

	// Build all the tasks and update the scMeta
	for (int scId = 0; scId < superClusterMetas.size(); ++scId) {
		for (int i = 0; i < workPerSc[scId].size(); i++) {
			//const bool bondedTask = workPerSc[scId][i].areBonded;
			ScScTask task;
			task.nointeractionMatrixIndex = workPerSc[scId][i].nointeractionMatrixIndexRelative != -1 ? workPerSc[scId][i].nointeractionMatrixIndexRelative + nBondedMatricesPrefixsum[scId] : -1;
			task.scIds[0] = scId;
			task.scIds[1] = workPerSc[scId][i].queryScId;
			task.resultIndices[0] = workPerSc[scId][i].resultIndexRelativeSelf + nResultsPrefixsum[scId];
			task.resultIndices[1] = scId != workPerSc[scId][i].queryScId ? (workPerSc[scId][i].resultIndexRelativeQuery + nResultsPrefixsum[workPerSc[scId][i].queryScId]) : -1;
			//task.resultIndices[1] = (workPerSc[scId][i].resultIndexRelativeQuery + nResultsPrefixsum[workPerSc[scId][i].queryScId]);
			tasks[nTasksPrefixsum[scId] + i] = task;
		}

		superClusterMetas[scId].resultsStartIndex = nResultsPrefixsum[scId];
		superClusterMetas[scId].nResults = nResultsReserved[scId];
	}

	//const size_t numNointeractionMatricesTotal = nBondedMatricesPrefixsum.back() + nBondedmatricesReserved.back();
	//nointeractionMatrices.resize(numNointeractionMatricesTotal);



	// Build all the nointeractionMatrices
	const std::vector<BoolMatrix16x16> nointeractionMatrices = BuildNointeractionMatrices(superClusterMetas, nBondedmatricesReserved, nBondedMatricesPrefixsum, pClustersMeta, workPerSc, box.particlesBondedToParticle);
	//for (const auto& mat : nointeractionMatrices) {
	//	mat.Print();
	//}
	//{
	//	std::vector<std::set<int>> expectedLjInteractions(16);
	//	for (int row = 0; row < 16; row++) {
	//		for (int col = 0; col < 16; col++) {
	//			if (col == 8 && row == 8)
	//				int aa = 0;
	//			auto _row = nointeractionMatrices[0].GetRow(row);
	//			if (!nointeractionMatrices[0].Get(_row, col)) {
	//				int pid0 = superClusterMetas[0].particlesIds[row];
	//				int pid1 = superClusterMetas[0].particlesIds[col];
	//				expectedLjInteractions[pid0].insert(pid1);
	//			}
	//		}
	//	}
	//	for (int pid = 0; pid < 16; pid++) {
	//		for (auto& interactPid : expectedLjInteractions[pid]) {
	//			printf("%d ", interactPid);
	//		}
	//		printf("\n");
	//	}
	//}


	//DebugUtils::VerifyIdentical(tasks, "ScScTasks" + std::to_string(simulation->getStep()));
	//DebugUtils::VerifyIdentical


	// Push back to device
	cudaMemcpy(superClustersControl->scMeta, superClusterMetas.data(), superClusterMetas.size() * sizeof(SuperClusterMeta), cudaMemcpyHostToDevice);
	cudaFree(scscTasksDevice);
	cudaFree(noInteractionMatricesDevice);
	scscTasksDevice = GenericCopyToDevice(tasks);
	noInteractionMatricesDevice = GenericCopyToDevice(nointeractionMatrices);

	nResults = nResultsPrefixsum.back() + nResultsReserved.back();
	cudaFree(scResultsDevice);
	cudaMalloc(&scResultsDevice, sizeof(SCResult) * nResults);
	cudaMemset(scResultsDevice, 0, sizeof(SCResult) * nResults);

	nTasks = numTasksTotal;
	if (simulation->getStep() == 787) {
		int a = 0;
	}
	return true;
}





bool Engine::MakeSuperClusterTasksCPU() {
	if (nSuperclusters == 0)
		return true;

	const Box& box = *simulation->box_host;
	Float3 boxSizeF = simulation->box_host->boxparams.BoxSizeFloat();

	const std::vector<PersistentClusterMeta>& pClustersMeta = simulation->box_host->persistentClustersMetadata;
	const std::vector<SuperCluster> superClusters = GenericCopyToHost(superClustersControl->scData, nSuperclusters);
	std::vector<SuperClusterMeta> superClusterMetas = GenericCopyToHost(superClustersControl->scMeta, nSuperclusters);

	const std::vector<std::array<float4, 4>> superclusterPositionSpheres = ComputeMeanposAndRadiiForEachPclusterInEachSupercluster(superClusters);

	std::vector<std::vector<ReservedTask>> workPerSc(superClusters.size());
	std::vector<int> nResultsReserved(superClusters.size(), 0);
	std::vector<int> nBondedmatricesReserved(superClusters.size(), 0);

	for (int scId = 0; scId < superClusterMetas.size(); ++scId) {
		for (int queryScId = scId; queryScId < superClusterMetas.size(); ++queryScId) {

			if (DoesSuperclustersInteract(superclusterPositionSpheres.data(), scId, queryScId, simulation->simparams_host.cutoff_nm, boxSizeF)) {
				const bool useNointeractionMatrix = scId == queryScId || ScAreBonded(superClusterMetas[scId], superClusterMetas[queryScId], box.pclustersBondedToPcluster.data());

				workPerSc[scId].emplace_back(ReservedTask{
					queryScId,
					nResultsReserved[scId],
					scId != queryScId ? nResultsReserved[queryScId] : nResultsReserved[queryScId],
					useNointeractionMatrix ? nBondedmatricesReserved[scId] : -1
					});

				nResultsReserved[scId]++;
				if (scId != queryScId)
					nResultsReserved[queryScId]++;
				if (useNointeractionMatrix) {
					nBondedmatricesReserved[scId]++;
				}
			}
		}
	}

	// Make prefixsums
	const std::vector<size_t> nResultsPrefixsum = ExlusivePrefixsum(nResultsReserved);
	const std::vector<size_t> nBondedMatricesPrefixsum = ExlusivePrefixsum(nBondedmatricesReserved);
	const std::vector<size_t> nTasksPrefixsum = ExlusivePrefixsum(workPerSc);

	const size_t numTasksTotal = nTasksPrefixsum.back() + workPerSc.back().size();
	std::vector<ScScTask> tasks(numTasksTotal);

	// Build all the tasks and update the scMeta
	for (int scId = 0; scId < superClusterMetas.size(); ++scId) {
		for (int i = 0; i < workPerSc[scId].size(); i++) {
			//const bool bondedTask = workPerSc[scId][i].areBonded;
			ScScTask task;
			task.nointeractionMatrixIndex = workPerSc[scId][i].nointeractionMatrixIndexRelative != -1 ? workPerSc[scId][i].nointeractionMatrixIndexRelative + nBondedMatricesPrefixsum[scId] : -1;
			task.scIds[0] = scId;
			task.scIds[1] = workPerSc[scId][i].queryScId;
			task.resultIndices[0] = workPerSc[scId][i].resultIndexRelativeSelf + nResultsPrefixsum[scId];
			task.resultIndices[1] = scId != workPerSc[scId][i].queryScId ? (workPerSc[scId][i].resultIndexRelativeQuery + nResultsPrefixsum[workPerSc[scId][i].queryScId]) : -1;
			//task.resultIndices[1] = (workPerSc[scId][i].resultIndexRelativeQuery + nResultsPrefixsum[workPerSc[scId][i].queryScId]);
			tasks[nTasksPrefixsum[scId] + i] = task;
		}

		superClusterMetas[scId].resultsStartIndex = nResultsPrefixsum[scId];
		superClusterMetas[scId].nResults = nResultsReserved[scId];
	}

	//const size_t numNointeractionMatricesTotal = nBondedMatricesPrefixsum.back() + nBondedmatricesReserved.back();
	//nointeractionMatrices.resize(numNointeractionMatricesTotal);



	// Build all the nointeractionMatrices
	const std::vector<BoolMatrix16x16> nointeractionMatrices = BuildNointeractionMatrices(superClusterMetas, nBondedmatricesReserved, nBondedMatricesPrefixsum, pClustersMeta, workPerSc, box.particlesBondedToParticle);
	//for (const auto& mat : nointeractionMatrices) {
	//	mat.Print();
	//}
	//{
	//	std::vector<std::set<int>> expectedLjInteractions(16);
	//	for (int row = 0; row < 16; row++) {
	//		for (int col = 0; col < 16; col++) {
	//			if (col == 8 && row == 8)
	//				int aa = 0;
	//			auto _row = nointeractionMatrices[0].GetRow(row);
	//			if (!nointeractionMatrices[0].Get(_row, col)) {
	//				int pid0 = superClusterMetas[0].particlesIds[row];
	//				int pid1 = superClusterMetas[0].particlesIds[col];
	//				expectedLjInteractions[pid0].insert(pid1);
	//			}
	//		}
	//	}
	//	for (int pid = 0; pid < 16; pid++) {
	//		for (auto& interactPid : expectedLjInteractions[pid]) {
	//			printf("%d ", interactPid);
	//		}
	//		printf("\n");
	//	}
	//}


	//DebugUtils::VerifyIdentical(tasks, "ScScTasks" + std::to_string(simulation->getStep()));
	//DebugUtils::VerifyIdentical


	// Push back to device
	cudaMemcpy(superClustersControl->scMeta, superClusterMetas.data(), superClusterMetas.size() * sizeof(SuperClusterMeta), cudaMemcpyHostToDevice);
	cudaFree(scscTasksDevice);
	cudaFree(noInteractionMatricesDevice);
	scscTasksDevice = GenericCopyToDevice(tasks);
	noInteractionMatricesDevice = GenericCopyToDevice(nointeractionMatrices);

	nResults = nResultsPrefixsum.back() + nResultsReserved.back();
	cudaFree(scResultsDevice);
	cudaMalloc(&scResultsDevice, sizeof(SCResult) * nResults);
	cudaMemset(scResultsDevice, 0, sizeof(SCResult) * nResults);

	nTasks = numTasksTotal;
	if (simulation->getStep() == 787) {
		int a = 0;
	}
	return true;
}



