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
	const int nSuperclustersUpperbound;

public:
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

__global__ void ComputeMeanposAndRadiiForEachPclusterInEachSuperclusterKernel(const SuperCluster* const superclusters, std::array<float4, 4>* const out, int nSuperclusters) {
	const int scId = blockIdx.x * blockDim.x + threadIdx.x;

	if (scId >= nSuperclusters)
		return;

	for (int pcid = 0; pcid < 4; pcid++) { // optim: 1 thread per pcid
		Float3 sum{};
		int cnt = 0;
		for (int pid = 0; pid < 4; pid++) {
			//const PData& pData = superclusters[scId].pData[pcid * 4 + pid];
			Float3 position = superclusters[scId].positions[pcid * 4 + pid];
			float epsilonSqrt = superclusters[scId].ljParams[pcid * 4 + pid].epsilonSqrt;
			if (epsilonSqrt != -1.f) {
				sum += position;
				cnt++;
			}
		}

		const Float3 meanPos = sum * (1.0f / static_cast<float>(cnt));
		float radius = 0;
		for (int pid = 0; pid < cnt; pid++) {
			//const PData& pData = superclusters[scId].pData[pcid * 4 + pid];
			Float3 position = superclusters[scId].positions[pcid * 4 + pid];
			radius = std::max(radius, (position - meanPos).len());
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

	if (validQuery && DoesSuperclustersInteract(tbContents.superclusterPositionSpheres, scId, queryScId, cutoffNm, boxSizeF)) {
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



// gridDim = (nSuperclusters, 1, 1)
// blockDim = (16, 1, 1)
__global__ void BuildNointeractionMatricesKernel(const SuperClusterMeta* const superClusterMetas, const PersistentClusterMeta* const pClustersMeta, 
	TaskBuilderControlContents tbContents, BoolMatrix16x16* const nointeractionMatrices, int nSuperclusters) {

	const int scId = blockIdx.x;

	std::array<int, 16> particleIdsSelf = GetParticleIdsOfSuperCluster(pClustersMeta, superClusterMetas[scId]);
	int matrixCount = 0;

	for (int i = 0; i < tbContents.nInteractionsOwned[scId]; i++) {
		InteractionToken token = tbContents.interactionsOwned[scId * TaskBuilderControlContents::maxTasksPerSc + i];

		if (!token.UseNointeractionMatrix()) {
			continue;
		}

		const int scIdQuery = token.GetQueryId();

		std::array<int, 16> particleIdsQuery = GetParticleIdsOfSuperCluster(pClustersMeta, superClusterMetas[scIdQuery]);
		const bool isSelfInteractionTask = scId == scIdQuery;

		const int row = threadIdx.x;
		uint16_t rowData = 0;
		for (int col = 0; col < 16; ++col) {
			int pidSelf = particleIdsSelf[row];
			int pidQuery = particleIdsQuery[col];
			if (pidSelf == -1 || pidQuery == -1)
				continue;

			bool noInteraction = tbContents.particlesBondedToParticle[particleIdsSelf[row]].Contains(particleIdsQuery[col]);
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

__global__ void TransposeNointeractionMatrices(const BoolMatrix16x16* const in, BoolMatrix16x16* const out) {
	__shared__ BoolMatrix16x16 matrixIn;
	__shared__ BoolMatrix16x16 matrixOut;

	if (threadIdx.x == 0) {
		matrixIn = in[blockIdx.x];
	}
	__syncthreads();

	uint16_t column = matrixIn.GetColumn(threadIdx.x);
	matrixOut.SetRow(threadIdx.x, column);
	__syncthreads();
	if (threadIdx.x == 0) {
		out[blockIdx.x] = matrixOut;
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

	ComputeMeanposAndRadiiForEachPclusterInEachSuperclusterKernel<<<(nSuperclusters + 31) / 32, 32 >>>(superClustersControl->scData, taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters);
	cudaDeviceSynchronize();

	//auto iSpheres = GenericCopyToHost(taskbuilderControl->contents.superclusterPositionSpheres, nSuperclusters);

	{
		dim3 gridDim{ (uint32_t)boxSize.InnerProduct(), (uint32_t)SuperClustersControl::maxClustersPerBlock, 1u };
		dim3 blockDim{ 3 * 3 * 3 * SuperClustersControl::maxClustersPerBlock, 1, 1 };
		ReserveInteractions << <gridDim, blockDim >> > (*superClustersControl, boxSize, taskbuilderControl->contents, simulation->simparams_host.cutoff_nm);
		LIMA_UTILS::genericErrorCheck("ReserveInteractions");
	}

	/*auto ninteractionsOwnedHost = GenericCopyToHost(taskbuilderControl->contents.nInteractionsOwned, nSuperclusters);
	auto nSuperclusterPerNode = GenericCopyToHost(superClustersControl->nSuperclustersInBlocks, boxSize.InnerProduct());*/



	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nResults, taskbuilderControl->contents.nResults + nSuperclusters + 1, taskbuilderControl->contents.nResultsPrefixsum);
	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nInteractionsOwned, taskbuilderControl->contents.nInteractionsOwned + nSuperclusters + 1, taskbuilderControl->contents.nTasksPrefixsum);
	thrust::exclusive_scan(thrust::device, taskbuilderControl->contents.nNointeractionmatricesOwned, taskbuilderControl->contents.nNointeractionmatricesOwned + nSuperclusters + 1, taskbuilderControl->contents.nNointeractionmatricesPrefixsum);
	cudaDeviceSynchronize();
	nResults = GenericCopyToHost(taskbuilderControl->contents.nResultsPrefixsum + nSuperclusters);
	nTasks = GenericCopyToHost(taskbuilderControl->contents.nTasksPrefixsum + nSuperclusters);
	const int nNointeractionMatrices = GenericCopyToHost(taskbuilderControl->contents.nNointeractionmatricesPrefixsum + nSuperclusters);


	scscTasksDevice.Expand(nTasks, 1.2);
	noInteractionMatricesDevice.Expand(nNointeractionMatrices, 1.2);
	noInteractionMatricesTransposedDevice.Expand(nNointeractionMatrices, 1.2);
	scResultsDevice.Expand(nResults, 1.2);



	BuildTasks << <(nSuperclusters + 31) / 32, 32 >> > (taskbuilderControl->contents, superClustersControl->scMeta, nSuperclusters, scscTasksDevice.Get());
	BuildNointeractionMatricesKernel << <nSuperclusters , 16 >> >(superClustersControl->scMeta, pClusterMetaDevice, taskbuilderControl->contents, noInteractionMatricesDevice.Get(), nSuperclusters);
	TransposeNointeractionMatrices<<<nNointeractionMatrices, 16 >>>(noInteractionMatricesDevice.Get(), noInteractionMatricesTransposedDevice.Get());
	cudaDeviceSynchronize();

	//auto resCounts = GenericCopyToHost(taskbuilderControl->contents.nResults, nSuperclustersUpperbound);

	
	cudaDeviceSynchronize();

	//auto tasksHost = GenericCopyToHost(scscTasksDevice, nSuperclustersUpperbound * TaskBuilderControlContents::maxTasksPerSc);
	//auto 


	//taskbuilderControl->Reset();

	return true;
}
