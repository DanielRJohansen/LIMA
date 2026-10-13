#include "Tests.h"
#include "Format.h"

#include "TestUtils.h"
#include "Filehandling.h"
#include "BoxImageBuilder.h"




namespace ForceComparisons {
	using namespace TestUtils;

	TestRoutine DoAllForceComparisons(Environment& environment, EnvMode envmode) {
		// TODO: Can't these all use automatedtestsdir? 
		const fs::path heavyDir = HeavyTestsDir() / "CompareWithOtherMdEngines/Forcecomparison1step";
		const std::array<fs::path, 6> workDirs{ heavyDir / "PoolNoES", heavyDir / "Singlebond", heavyDir / "Anglebond", heavyDir / "Dihedralbond", heavyDir / "Cmapbond",
			AutomatedTestsDir() / "ForceComparisons/DihedralBackbone" };
		std::vector<SimulationHandle> handles;
		handles.reserve(workDirs.size());
		for (const fs::path& workDir : workDirs) {
			SimulationJob job;
			job.workDir = workDir;
			job.grofile.emplace(workDir / "conf.gro");
			job.topfile.emplace(workDir / "topol.top");
			job.simParams.emplace(fs::exists(workDir / "sim_params.txt")
				? workDir / "sim_params.txt" : workDir.parent_path() / "sim_params.txt");
			handles.push_back(environment.Submit(std::move(job)));
		}
		for (std::size_t index = 0; index < handles.size(); index++) {
			auto completed = co_await std::move(handles[index]);
			const fs::path& workDir = workDirs[index];
			const std::string name = workDir.filename().string();
			const auto reference = FileUtils::ReadCsvAsVectorOfFloat3(workDir / "forces.csv");
			const auto actual = SimAnalysis::GetForces(*completed.simulation, 0);
			if (actual.size() != reference.size())
				co_return LimaUnittestResult::Fail(name + " force count mismatch", envmode == Full);
			for (std::size_t force = 0; force < actual.size(); force++) {
				if ((actual[force] - reference[force]).len() / reference[force].len() > 1e-3f)
					co_return LimaUnittestResult::Fail(Lima::Format("{} force mismatch at atom {}: ({:.2f} {:.2f} {:.2f}) vs GROMACS ({:.2f} {:.2f} {:.2f})",
						name, force, actual[force].x, actual[force].y, actual[force].z,
						reference[force].x, reference[force].y, reference[force].z), envmode == Full);
			}
		}
		co_return LimaUnittestResult{ true, "Success", envmode == Full };
	}
}
