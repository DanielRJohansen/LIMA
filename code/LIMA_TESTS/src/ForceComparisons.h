#pragma once

#include "TestUtils.h"
#include "Filehandling.h"
#include "BoxImageBuilder.h"




namespace ForceComparisons {
	using namespace TestUtils;

	std::vector<Float3> GetForces(const Simulation& simulation, int64_t step) {
		const auto& boxImage = simulation.boxImage;
		std::vector<Float3> forces(boxImage->grofile.atoms.size());
		for (int clusterId = 0; clusterId < simulation.box->persistentClusters.size(); clusterId++) {
			for (int particleId = 0; particleId < PersistentCluster::maxParticles; particleId++) {
				const int globalId = simulation.box->persistentClustersMetadata[clusterId].particleIdsGlobal[particleId];
				if (globalId >= 0)
					forces[globalId] = simulation.forceBuffer->GetDatapointAtStep(clusterId, particleId, step) / KILO;
			}
		}
		return forces;
	}

	std::function<LimaUnittestResult()> DoAllForceComparisons(Environment& environment, EnvMode envmode) {
		const std::array<std::string, 4> directories{ "PoolNoES", "Singlebond", "Anglebond", "Dihedralbond" };
		std::vector<SimulationHandle> handles;
		handles.reserve(directories.size());
		for (const auto& directory : directories) {
			const fs::path workDir = HeavyTestsDir() / "CompareWithOtherMdEngines/Forcecomparison1step" / directory;
			SimulationJob job;
			job.workDir = workDir;
			job.groPath = workDir / "conf.gro";
			job.topPath = workDir / "topol.top";
			job.simParamsPath = fs::exists(workDir / "sim_params.txt")
				? workDir / "sim_params.txt" : workDir.parent_path() / "sim_params.txt";
			handles.push_back(environment.Submit(std::move(job)));
		}
		return [handles = std::move(handles), directories, envmode]() mutable {
		for (std::size_t index = 0; index < handles.size(); index++) {
			auto completed = handles[index].Get();
			const fs::path workDir = HeavyTestsDir() / "CompareWithOtherMdEngines/Forcecomparison1step" / directories[index];
			const auto reference = FileUtils::ReadCsvAsVectorOfFloat3(workDir / "forces.csv");
			const auto actual = GetForces(*completed.simulation, 0);
			if (actual.size() != reference.size())
				return LimaUnittestResult{ false, directories[index] + " force count mismatch", envmode == Full };
			for (std::size_t force = 0; force < actual.size(); force++) {
				if ((actual[force] - reference[force]).len() / reference[force].len() > 1e-3f)
					return LimaUnittestResult{ false, directories[index] + " force mismatch", envmode == Full };
			}
		}
		return LimaUnittestResult{ true, "Success", envmode == Full };
		};
	}
}
