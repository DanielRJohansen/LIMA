#pragma once

#include "TestUtils.h"
#include "Environment.h"
#include "Printer.h"
#include "LimaTypes.cuh"
#include "Programs.h"

#include <future> // TODO: REmove
#include <string>


namespace TestMDStability {
	using namespace TestUtils;

	using PendingSimulation = std::shared_ptr<std::future<SimulationResult>>;

	// Schedule energy minimization now, then automatically submit the MD job as soon
	// as minimization finishes. std::async is used only as a test-side continuation:
	// its thread waits and submits jobs, while Environment still performs all costly
	// preprocessing and simulation work. The shared_ptr makes the move-only future
	// usable from LimaUnittest's copyable std::function callback.
	static PendingSimulation EnergyMinAndRun(
		Environment& environment, const fs::path& workDir, EnvMode envmode) {
		SimParams emParams;
		emParams.em_variant = true;
		emParams.dt = 1.5f * FEMTO_TO_NANO;
		emParams.em_force_tolerance = 100.f;
		emParams.data_logging_interval = 50;
		emParams.enable_electrostatics = true;
		emParams.n_steps = 20000;
		emParams.bc_select = BoundaryConditionSelect::PBC;

		SimulationJob emJob;
		emJob.workDir = workDir;
		emJob.groPath = workDir / "molecule/conf.gro";
		emJob.topPath = workDir / "molecule/topol.top";
		emJob.simParams = std::move(emParams);
		emJob.mode = envmode;
		auto emHandle = environment.Submit(std::move(emJob));

		return std::make_shared<std::future<SimulationResult>>(std::async(std::launch::async,
			[emHandle = std::move(emHandle), &environment, workDir, envmode]() mutable {
				auto minimized = emHandle.Get();

				SimulationJob mdJob;
				mdJob.workDir = workDir;
				mdJob.simParamsPath = workDir / "sim_params.txt";
				mdJob.initialSimulation = std::move(minimized.simulation);
				mdJob.mode = envmode;
				mdJob.analyze = true;
				return environment.Submit(std::move(mdJob)).Get();
			}));
	}

	static std::function<LimaUnittestResult()> LoadEnergyMinAndRunBasicSimulation(
		Environment& environment, std::string folderName, EnvMode envmode, std::string testName)
	{
		const fs::path workDir = AutomatedTestsDir() / folderName;

		// This starts both stages without running them on the TestManager thread:
		// EnergyMinAndRun submits EM immediately, and submits MD when EM completes.
		auto simulation = EnergyMinAndRun(environment, workDir, envmode);

		// TestManager invokes this callback later, in scheduled test order. By then the
		// simulation is normally complete; get() only waits if Environment is still busy.
		// All result validation deliberately stays on the TestManager thread.
		return [simulation = std::move(simulation), testName = std::move(testName), envmode]() mutable {
			auto completed = simulation->get();
			if (!completed.analysis)
				return LimaUnittestResult{ false, "Environment returned no analysis", envmode == Full };

			const auto evaluation = evaluateTest(testName,
				{ completed.analysis->variance_coefficient }, { completed.analysis->energy_gradient });
			return LimaUnittestResult{ evaluation.first, evaluation.second, envmode == Full };
		};
	}

	static std::function<LimaUnittestResult()> TestDeterministic(Environment& environment, EnvMode envmode) {
		const fs::path workDir = AutomatedTestsDir() / "T4Lysozyme";
		std::vector<PendingSimulation> simulations;
		for (int run = 0; run < 2; run++)
			simulations.push_back(EnergyMinAndRun(environment, workDir, envmode));

		return [simulations = std::move(simulations), envmode]() mutable {
		std::optional<float> referenceVc;
		std::optional<float> referenceGradient;
		for (auto& simulation : simulations) {
			auto completed = simulation->get();
			const float vc = completed.analysis->variance_coefficient;
			const float gradient = completed.analysis->energy_gradient;
			if (referenceVc && (vc != *referenceVc || gradient != *referenceGradient))
				return LimaUnittestResult{ false, "Simulation results were not deterministic", envmode == Full };
			referenceVc = vc;
			referenceGradient = gradient;
		}
		return LimaUnittestResult{ true, "Success", envmode == Full };
		};
	}

	static bool doMoleculeTranslationTest(std::string foldername) {
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
