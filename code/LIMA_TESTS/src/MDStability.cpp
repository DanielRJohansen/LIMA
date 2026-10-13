#include "Tests.h"

#include "TestUtils.h"
#include "Environment.h"
#include "Printer.h"
#include "LimaTypes.cuh"
#include "Programs.h"

#include <string>


namespace TestMDStability {
	using namespace TestUtils;

	SimulationJob MakeEnergyMinJob(const fs::path& workDir, EnvMode envmode) {
		SimulationJob emJob;
		emJob.workDir = workDir;
		emJob.simParams = SimParams::BasicEMSimParams();
		emJob.mode = envmode;
		return emJob;
	}

	TestRoutine LoadEnergyMinAndRunBasicSimulation(
		Environment& environment, EnvMode envmode, std::string folderName, std::string testName) {
		const fs::path workDir = AutomatedTestsDir() / folderName;
		auto minimized = co_await environment.Submit(MakeEnergyMinJob(workDir, envmode));

		SimulationJob mdJob;
		mdJob.workDir = workDir;
		mdJob.simParams.emplace(workDir / "sim_params.txt");
		mdJob.initialSimulation = std::move(minimized.simulation);
		mdJob.mode = envmode;
		mdJob.postprocess = SimAnalysis::AnalyzeEnergy;
		auto completed = co_await environment.Submit(std::move(mdJob));
		if (!completed.analysis)
			co_return LimaUnittestResult::Fail("Environment returned no analysis", envmode == Full);

		co_return EvaluateVcTest(testName, { completed.analysis->variance_coefficient }, { completed.analysis->energy_gradient }, envmode == Full);
	}

	// constraints = h-bonds on T4 lysozyme in water: the constrained distances hold, and the energy is conserved, compared with the
	// same run with flexible bonds
	TestRoutine TestConstraints(Environment& environment, EnvMode envmode) {
		const fs::path workDir = AutomatedTestsDir() / "T4Lysozyme";
		std::array<SimulationHandle, 2> emHandles;
		for (int constrained = 0; constrained < 2; constrained++) {
			SimulationJob emJob = MakeEnergyMinJob(workDir, envmode);
			// The box is built by the minimization, so the constraints must be set there too
			emJob.simParams->constraints = constrained ? ConstraintsSelect::HBonds : ConstraintsSelect::None;
			emHandles[constrained] = environment.Submit(std::move(emJob));
		}
		std::array<SimulationHandle, 2> mdHandles;
		for (int constrained = 0; constrained < 2; constrained++) {
			SimulationJob mdJob;
			mdJob.workDir = workDir;
			mdJob.simParams.emplace(workDir / "sim_params.txt");
			mdJob.simParams->constraints = constrained ? ConstraintsSelect::HBonds : ConstraintsSelect::None;
			mdJob.initialSimulation = std::move((co_await std::move(emHandles[constrained])).simulation);
			mdJob.mode = envmode;
			mdJob.postprocess = SimAnalysis::AnalyzeEnergy;
			mdHandles[constrained] = environment.Submit(std::move(mdJob));
		}
		const auto flexible = co_await std::move(mdHandles[0]);
		const auto constrained = co_await std::move(mdHandles[1]);
		if (!flexible.analysis || !constrained.analysis)
			co_return LimaUnittestResult::Fail("Environment returned no analysis", envmode == Full);

		const Box& box = *constrained.simulation->box;
		if (box.pclusterConstraints.empty())
			co_return LimaUnittestResult::Fail("The box has no constraints", envmode == Full);
		const Float3 boxSize = box.boxparams.BoxSizeFloat();
		float maxError = 0.f;
		int64_t nConstraints = 0;
		for (size_t pc = 0; pc < box.pclusterConstraints.size(); pc++) {
			const PclusterConstraints& constraints = box.pclusterConstraints[pc];
			for (int c = 0; c < constraints.nConstraints; c++) {
				Float3 r = box.persistentClusters[pc].pqd[constraints.particles[c] & 3].position - box.persistentClusters[pc].pqd[constraints.particles[c] >> 2].position;
				r = r - boxSize * Float3{ std::round(r.x / boxSize.x), std::round(r.y / boxSize.y), std::round(r.z / boxSize.z) };
				maxError = std::max(maxError, std::abs(r.len() - constraints.distance[c]) / constraints.distance[c]);
				nConstraints++;
			}
		}
		const float vcFlexible = flexible.analysis->variance_coefficient;
		const float vcConstrained = constrained.analysis->variance_coefficient;
		if (maxError >= 1e-4f || !std::isfinite(vcConstrained) || vcConstrained >= 2.f * vcFlexible)
			co_return LimaUnittestResult::Fail(Lima::Format("{} constraints, max error {:.1e}, VC {:.3e} (flexible {:.3e})",
				nConstraints, maxError, vcConstrained, vcFlexible), envmode == Full);
		co_return LimaUnittestResult{ true, OneLineMessage("max error {:.1e}, VC {:.2e} / flexible {:.2e}",
			maxError, vcConstrained, vcFlexible), envmode == Full };
	}

	TestRoutine TestDeterministic(Environment& environment, EnvMode envmode) {
		const fs::path workDir = AutomatedTestsDir() / "T4Lysozyme";
		constexpr int nRuns = 2;
		std::array<SimulationHandle, nRuns> emHandles;
		for (auto& handle : emHandles)
			handle = environment.Submit(MakeEnergyMinJob(workDir, envmode));

		std::array<SimulationHandle, nRuns> mdHandles;
		for (int run = 0; run < nRuns; run++) {
			auto minimized = co_await std::move(emHandles[run]);
			SimulationJob mdJob;
			mdJob.workDir = workDir;
			mdJob.simParams.emplace(workDir / "sim_params.txt");
			mdJob.initialSimulation = std::move(minimized.simulation);
			mdJob.mode = envmode;
			mdJob.postprocess = SimAnalysis::AnalyzeEnergy;
			mdHandles[run] = environment.Submit(std::move(mdJob));
		}

		std::optional<float> referenceVc;
		std::optional<float> referenceGradient;
		for (auto& handle : mdHandles) {
			auto completed = co_await std::move(handle);
			const float vc = completed.analysis->variance_coefficient;
			const float gradient = completed.analysis->energy_gradient;
			if (referenceVc && (vc != *referenceVc || gradient != *referenceGradient))
				co_return LimaUnittestResult::Fail("Simulation results were not deterministic", envmode == Full);
			referenceVc = vc;
			referenceGradient = gradient;
		}
		co_return LimaUnittestResult{ true, "Success", envmode == Full };
	}

	bool doMoleculeTranslationTest(std::string foldername) {
		//auto env = TestUtils::basicSetup(foldername, SimulationParams{100, 10000});

		//const float vel = EngineUtils::tempToVelocity(300, 0.012f);	// [m/s] <=> [lm/ls]
		//const float dt = env.getSimparamRef()->dt;

		//const auto dir = Float3{ -1.f, -1.f, -1.f }.norm();	// Inv since we move the pos_prev

		//auto& coordarray_prev_ptr = env.getCoordarrayRef("prev");

		//// Too lazy to figure out which coords are vacant.. So i give same velocity to all possible
		//// particles, even inactive ones. If everything works, it should be no problem
		//for (auto& coords : coordarray_prev_ptr) {
		//	for (auto& pos : coords.rel_positions) {
		//		pos += Coord{ dir * vel * dt };
		//	}
		//}

		//env.run();

		//return TestUtils::verifyStability(env, 0.08);
		return true;
	}









	// Tests proposed by GPT 3.5, lol
	/*
		Atom Shift Analysis Test: This test verifies that the positions of atoms are accurately updated according to the integration scheme being used.

		Velocity Correlation Test: This test checks the correctness of the velocity distribution, ensuring that it follows the expected Maxwell-Boltzmann distribution.

		Force Field Validation Test: This test assesses the accuracy of the interatomic force field being used, comparing predicted structural and energetic properties to experimental data.

		Temperature Equilibration Test: This test confirms that the system reaches a stable temperature and that the temperature remains constant over time.

		Trajectory Analysis Test: This test analyzes the trajectory of the system to verify that it is physically reasonable and corresponds to expected behavior.

		Pressure Control Test: This test ensures that the system is correctly maintained at the desired pressure and that fluctuations are within expected limits.

		Energy Conservation Test: This test verifies that the total energy of the system is conserved over time, which is essential for accurate and stable simulations.

		Structural Stability Test: This test assesses the stability of the system by measuring fluctuations in bond lengths, angles, and other structural parameters.

		Hydrogen Bond Analysis Test: This test analyzes the presence and behavior of hydrogen bonds within the system, which can be important for understanding molecular interactions.

		Diffusion Coefficient Test: This test measures the rate of diffusion of particles in the system, which is related to properties such as viscosity and can provide insights into the behavior of the system.

		Radial Distribution Function Test: This test calculates the radial distribution function of the system, which provides information on the distribution of particles around each other.

		Potential Energy Surface Test: This test examines the potential energy surface of the system to ensure that it is consistent with expected behavior and accurately reflects the interatomic interactions.

		Solvation Free Energy Test: This test calculates the solvation free energy of a solute molecule in a solvent, which can provide insights into the energetics of molecular interactions.

		Conformational Sampling Test: This test examines the ability of the simulation to explore the conformational space of the system, which is important for understanding the behavior of flexible molecules.

		Viscosity Calculation Test: This test measures the viscosity of the system, which is related to diffusion and can be important for understanding the behavior of complex fluids.
	*/
}
