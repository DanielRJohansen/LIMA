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
		// Each atom must match to 1e-3, or with rmsTolerance the forces as a whole, as PME differs per atom between the engines by up
		// to ~0.4%. A reference with rmsTolerance may cover only the first atoms, the others must cover all
		struct Comparison {
			fs::path workDir;
			float rmsTolerance = 0.f;
		};
		// VillinAmber: the villin protein of the Salad benchmark (AMBER99SB) with its charges zeroed, so bonded terms and LJ are
		// compared without PME's differences. Reference from GROMACS 2024.2, one step, 5 nm cubic box.
		// VillinBox: the whole Salad villin system, in its rectangular, non-integer 4.30 x 3.29 x 3.52 nm box, with PME. Its water is
		// settled, for which LIMA generates bonds that GROMACS does not have, so only the protein's forces are compared. Reference from
		// GROMACS 2024.2 with fine PME (order 6, spacing 0.08 nm, ewald-rtol 1e-6)
		const std::array<Comparison, 8> comparisons{ { { heavyDir / "PoolNoES" }, { heavyDir / "Singlebond" }, { heavyDir / "Anglebond" },
			{ heavyDir / "Dihedralbond" }, { heavyDir / "Cmapbond" }, { AutomatedTestsDir() / "ForceComparisons/DihedralBackbone" },
			{ AutomatedTestsDir() / "ForceComparisons/VillinAmber" }, { AutomatedTestsDir() / "ForceComparisons/VillinBox", 1e-3f } } };
		std::vector<SimulationHandle> handles;
		handles.reserve(comparisons.size());
		for (const auto& [workDir, rmsTolerance] : comparisons) {
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
			const auto& [workDir, rmsTolerance] = comparisons[index];
			const std::string name = workDir.filename().string();
			const auto reference = FileUtils::ReadCsvAsVectorOfFloat3(workDir / "forces.csv");
			const auto actual = SimAnalysis::GetForces(*completed.simulation, 0);
			if (rmsTolerance > 0.f ? actual.size() < reference.size() : actual.size() != reference.size())
				co_return LimaUnittestResult::Fail(name + " force count mismatch", envmode == Full);
			if (rmsTolerance > 0.f) {
				double errorSq = 0.0, referenceSq = 0.0;
				for (std::size_t force = 0; force < reference.size(); force++) {
					errorSq += (actual[force] - reference[force]).lenSquared();
					referenceSq += reference[force].lenSquared();
				}
				const double rmsError = std::sqrt(errorSq / referenceSq);
				if (rmsError > rmsTolerance)
					co_return LimaUnittestResult::Fail(Lima::Format("{} RMS force error {:.2e} vs GROMACS, max {:.0e}", name, rmsError, rmsTolerance), envmode == Full);
				continue;
			}
			for (std::size_t force = 0; force < reference.size(); force++) {
				if ((actual[force] - reference[force]).len() / reference[force].len() > 1e-3f)
					co_return LimaUnittestResult::Fail(Lima::Format("{} force mismatch at atom {}: ({:.2f} {:.2f} {:.2f}) vs GROMACS ({:.2f} {:.2f} {:.2f})",
						name, force, actual[force].x, actual[force].y, actual[force].z,
						reference[force].x, reference[force].y, reference[force].z), envmode == Full);
			}
		}
		co_return LimaUnittestResult{ true, "Success", envmode == Full };
	}
}
