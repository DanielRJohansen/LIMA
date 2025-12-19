#include "Engine.cuh"
#include <set>
#include <execution>
#include "EngineCore.h"



bool ScAreBonded(const SuperClusterMeta& sc0, const SuperClusterMeta& sc1, const std::vector<std::set<int>>& bondedToPclusters) {
	for (auto pclusterId0 : sc0.pclusterIds) {
		for (auto pclusterId1 : sc1.pclusterIds) {
			if (bondedToPclusters[pclusterId0].contains(pclusterId1)) {
				return true;
			}
		}
	}
	return false;
}

bool MakeSuperClustersCPU(const std::vector<PersistentCluster>& pClusters, const BoxParams& params, const SimParams& simparams, const std::vector<SuperCluster>& superClusters, const std::vector<std::set<int>>& bondedToPclusters, /*1 set per pCluster*/
	std::vector<SuperClusterMeta>& superClusterMetasInOut, std::vector<ScScTask>& tasksOut) {	


	std::vector<Float3> scMeanPos(superClusters.size());
	std::transform(
		std::execution::par,
		superClusters.begin(),
		superClusters.end(),
		scMeanPos.begin(),
		[](const SuperCluster& sc) {
			Float3 sum{};
			int cnt = 0;

			for (const PData& p : sc.pData) {
				if (p.Valid()) {
					sum += p.position;
					++cnt;
				}
			}

			return sum * (1.0f / static_cast<float>(cnt));
		}
	);
	//for (int i = 0; i < superClusters.size(); i++) {
	//	const auto& sc = superClusters[i];
	//	Float3 meanPos{};
	//	int cnt = 0;
	//	for (const PData& pData : sc.pData) {
	//		if (pData.Valid()) {
	//			meanPos += pData.position;
	//			cnt++;
	//		}
	//	}
	//	meanPos *= 1.f/static_cast<float>(cnt);
	//	scMeanPos[i] = meanPos;
	//}



	const int maxExpectedTasksPerSc = 64; // seems high...

	const size_t maxTasks = superClusters.size() * maxExpectedTasksPerSc /2; // tasks are not distributed symmetrically

	std::atomic<size_t> taskCount = 0;
	std::atomic<size_t> bondedMatricesCount = 0;

	std::vector<ScScTask> scScTasks;
	scScTasks.resize(maxTasks);
	Float3 boxSizeF = params.BoxSizeFloat();


	std::vector<std::vector<std::pair<int, bool>>> workPerSc(superClusters.size());
	std::vector<int> nResultsReserved(superClusters.size());
	
	for (int scId = 0; scId < superClusterMetasInOut.size(); ++scId) {
		//std::vector<int> neighborScIds;
		//std::vector<bool> isBonded;



		for (int queryScId = scId; queryScId < superClusterMetasInOut.size(); ++queryScId) {			
			
			const float hyperDist = LIMAPOSITIONSYSTEM::calcHyperDistNM(scMeanPos[scId], scMeanPos[queryScId], boxSizeF, BoundaryConditionSelect::PBC);
			if (hyperDist < simparams.cutoff_nm) {
				workPerSc[scId].emplace_back(queryScId, ScAreBonded(superClusterMetasInOut[scId], superClusterMetasInOut[queryScId], bondedToPclusters));
				nResultsReserved[scId]++;
				if (scId != queryScId) {
					nResultsReserved[queryScId]++;
				}
			}
		}

		const int nBondedTasks = std::count(isBonded.begin(), isBonded.end(), true);

		const size_t startTaskIndex = taskCount.fetch_add(neighborScIds.size());
		const size_t startBondedMatrixIndex = bondedMatricesCount.fetch_add(nBondedTasks);

		for (int i = 0; i < neighborScIds.size(); i++) {
			ScScTask task;
			task.nointeractionMatrixIndex = startBondedMatrixIndex + i;
			task.scIds[0] = scId;
			task.scIds[1] = neighborScIds[i];
		}
	}




	

	// Dummy implementation: create tasks for each pair of SuperClusters
	for (size_t i = 0; i < superClusters.size(); ++i) {
		for (size_t j = i + 1; j < superClusters.size(); ++j) {
			ScScTask task;
			task.scIds[0] = static_cast<int>(i);
			task.scIds[1] = static_cast<int>(j);
			task.resultIndices[0] = 0; // Placeholder
			task.resultIndices[1] = 0; // Placeholder
			task.nointeractionMatrixIndex = 0; // Placeholder
			scScTasks.push_back(task);
		}
	}
	//return { superClusters, superClusterMetas, scScTasks };

	return true;
}