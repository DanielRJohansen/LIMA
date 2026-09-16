#pragma once

#include "TestUtils.h"
#include "Filehandling.h"
#include "BoxImageBuilder.h"




namespace ForceComparisons {
	using namespace TestUtils;

	TestRoutine DoAllForceComparisons(Environment& environment, EnvMode envmode) {
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
		for (std::size_t index = 0; index < handles.size(); index++) {
			auto completed = co_await std::move(handles[index]);
			const fs::path workDir = HeavyTestsDir() / "CompareWithOtherMdEngines/Forcecomparison1step" / directories[index];
			const auto reference = FileUtils::ReadCsvAsVectorOfFloat3(workDir / "forces.csv");
			const auto actual = SimAnalysis::GetForces(*completed.simulation, 0);
			if (actual.size() != reference.size())
				co_return LimaUnittestResult{ false, directories[index] + " force count mismatch", envmode == Full };
			for (std::size_t force = 0; force < actual.size(); force++) {
				if ((actual[force] - reference[force]).len() / reference[force].len() > 1e-3f)
					co_return LimaUnittestResult{ false, directories[index] + " force mismatch", envmode == Full };
			}
		}
		co_return LimaUnittestResult{ true, "Success", envmode == Full };
	}
}
