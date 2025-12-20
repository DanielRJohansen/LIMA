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

std::array<int, 16> GetParticleIdsOfSuperCluster(const std::vector<PersistentClusterMeta>& pClusterMeta, const SuperClusterMeta& scMeta) {
	std::array<int, 16> particleIds{};
	int cnt = 0;
	for (auto pcId : scMeta.pclusterIds) {
		for (int particleId : pClusterMeta[pcId].particleIdsGlobal) {
			particleIds[cnt++] = particleId;
		}
	}
	return particleIds;
}

struct ReservedTask {
	int queryScId;
	bool areBonded;
	int resultIndexRelativeSelf=-1;
	int resultIndexRelativeQuery=-1;
	int nointeractionMatrixIndexRelative = -1;
};

bool Engine::MakeSuperClusterTasksCPU(const std::vector<PersistentCluster>& pClusters, const std::vector< PersistentClusterMeta>& pClustersMeta,
	const std::vector<SuperCluster>& superClusters,
	std::vector<SuperClusterMeta>& superClusterMetasInOut, std::vector<ScScTask>& tasksOut, std::vector<BoolMatrix16x16>& nointeractionMatrices) {	


	const Box& box = *simulation->box_host;


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



	
	Float3 boxSizeF = simulation->box_host->boxparams.BoxSizeFloat();


	std::vector<std::vector<ReservedTask>> workPerSc(superClusters.size());
	std::vector<int> nResultsReserved(superClusters.size(), 0);
	std::vector<int> nBondedmatricesReserved(superClusters.size(), 0);

	for (int scId = 0; scId < superClusterMetasInOut.size(); ++scId) {
		for (int queryScId = scId; queryScId < superClusterMetasInOut.size(); ++queryScId) {			
			
			const float hyperDist = LIMAPOSITIONSYSTEM::calcHyperDistNM(scMeanPos[scId], scMeanPos[queryScId], boxSizeF, BoundaryConditionSelect::PBC);
			if (hyperDist < simulation->simparams_host.cutoff_nm) {
				const bool bonded = ScAreBonded(superClusterMetasInOut[scId], superClusterMetasInOut[queryScId], box.pclusterBondedToPcluster);
				workPerSc[scId].emplace_back(ReservedTask{
					queryScId,
					bonded,
					nResultsReserved[scId],
					scId != queryScId ? nResultsReserved[queryScId] : -1,
					bonded ? nBondedmatricesReserved[scId] : -1
				});

				nResultsReserved[scId]++;
				if (scId != queryScId) {
					nResultsReserved[queryScId]++;
				}
				if (bonded) {
					nBondedmatricesReserved[scId]++;
				}
			}
		}
	}

	// Make prefixsums
	std::vector<size_t> nResultsPrefixsum(superClusters.size());
	std::vector<size_t> nBondedMatricesPrefixsum(superClusters.size());
	std::vector<size_t> nTasksPrefixsum(superClusters.size());
	std::exclusive_scan(std::execution::par, nResultsReserved.begin(), nResultsReserved.end(), nResultsPrefixsum.begin(), 0);
	std::exclusive_scan(std::execution::par, nBondedmatricesReserved.begin(), nBondedmatricesReserved.end(), nBondedMatricesPrefixsum.begin(), 0);	
	//std::exclusive_scan(std::execution::par, workPerSc.begin(), workPerSc.end(), nTasksPrefixsum.begin(), 0);
	std::transform_exclusive_scan(std::execution::par, workPerSc.begin(), workPerSc.end(), nTasksPrefixsum.begin(), size_t{ 0 }, std::plus<>{},
		[](const std::vector<ReservedTask>& v) { return v.size(); }
	);

	/*std::exclusive_scan(std::execution::par, workPerSc.begin(), workPerSc.end(), nTasksPrefixsum.begin(), size_t{ 0 },
		[](size_t a, const std::vector<ReservedTask>& b) { return a + b.size(); });*/

	const size_t numTasksTotal = nTasksPrefixsum.back() + workPerSc.back().size();
	tasksOut.resize(numTasksTotal);
	std::vector<ScScTask> scScTasks(numTasksTotal);

	// Build all the tasks and update the scMeta
	for (int scId = 0; scId < superClusterMetasInOut.size(); ++scId) {
		for (int i = 0; i < workPerSc[scId].size(); i++) {
			const bool bondedTask = workPerSc[scId][i].areBonded;
			ScScTask task;
			task.nointeractionMatrixIndex = bondedTask ? workPerSc[scId][i].nointeractionMatrixIndexRelative + nBondedMatricesPrefixsum[scId] : -1;
			task.scIds[0] = scId;
			task.scIds[1] = workPerSc[scId][i].queryScId;
			task.resultIndices[0] = workPerSc[scId][i].resultIndexRelativeSelf + nResultsPrefixsum[scId];
			task.resultIndices[1] = scId != workPerSc[scId][i].queryScId ? (workPerSc[scId][i].resultIndexRelativeQuery + nResultsPrefixsum[workPerSc[scId][i].queryScId]) : -1;
			tasksOut[nTasksPrefixsum[scId] + i] = task;
		}

		superClusterMetasInOut[scId].resultsStartIndex = nResultsPrefixsum[scId];
		superClusterMetasInOut[scId].nResults = nResultsReserved[scId];
	}
	
	const size_t numNointeractionMatricesTotal = nBondedMatricesPrefixsum.back() + nBondedmatricesReserved.back();
	nointeractionMatrices.resize(numNointeractionMatricesTotal);
	


	// Build all the nointeractionMatrices
	for (int scId = 0; scId < superClusterMetasInOut.size(); ++scId) {
		std::array<int, 16> particleIdsSelf = GetParticleIdsOfSuperCluster(pClustersMeta, superClusterMetasInOut[scId]);
		for (int i = 0; i < workPerSc[scId].size(); i++) {
			if (!workPerSc[scId][i].areBonded)
				continue;

			const int queryScId = workPerSc[scId][i].queryScId;
			BoolMatrix16x16 nointeractionMatrix{};

			std::array<int, 16> particleIdsQuery = GetParticleIdsOfSuperCluster(pClustersMeta, superClusterMetasInOut[queryScId]);


			for (int col = 0; col < 16; ++col) {
				for (int row = 0; row < 16; ++row) {
					const bool noInteraction = box.particleBondedToParticle[particleIdsSelf[row]].contains(particleIdsQuery[col]);
					nointeractionMatrix.Set(row, col, noInteraction); 
				}
			}


			const int matrixIndex = workPerSc[scId][i].nointeractionMatrixIndexRelative + nBondedMatricesPrefixsum[scId];
			nointeractionMatrices[matrixIndex] = nointeractionMatrix;
		}
	}


	return true;
}