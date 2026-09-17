#pragma once

#include "TestUtils.h"
#include "PhysicsUtils.cuh"
#include <format>

namespace ForceCorrectness {
	using namespace TestUtils;

	const fs::path TestsDir() { return AutomatedTestsDir(); }

	// Construct only the lightweight job description here. Environment performs all
	// file parsing, molecule construction, and GPU work on its bounded worker threads.
	SimulationJob MakeJob(const fs::path& workDir, EnvMode envmode, bool analyze = false) {
		SimulationJob job;
		job.workDir = workDir;
		job.mode = envmode;
		if (analyze)
			job.postprocess = SimAnalysis::AnalyzeEnergy;
		return job;
	}

	TestRoutine FinishStabilityTest(Environment& environment, std::string name, EnvMode envmode,
		SimulationJob job) {
		auto completed = co_await environment.Submit(std::move(job));
		const auto evaluation = evaluateTest(name,
			{ completed.analysis->variance_coefficient }, { completed.analysis->energy_gradient });
		co_return LimaUnittestResult{ evaluation.first, evaluation.second, envmode == Full };
	}

	
	TestRoutine doPoolBenchmark(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "Pool";
		// Test assumes two carbon particles in conf.gro.
		constexpr float mass = 12.011000f / 1000.f; // [kg/mol]
		constexpr float temperature = 400.f;
		const float velocity = PhysicsUtils::tempToVelocity(temperature, mass); // [m/s] == [nm/ns]
		auto job = MakeJob(workDir, envmode, true);
		job.simParams = SimParams{};
		job.simParams->enable_electrostatics = false;
		job.simParams->n_steps = LIMA_UTILS::roundUp(3000000 / static_cast<int>(velocity), 100);
		job.configureSimulation = [velocity](Simulation& simulation) {
			simulation.box->pclusterInterimStates[0].vels_prev[0] = Float3{ velocity, 0.f, 0.f };
			simulation.box->pclusterInterimStates[1].vels_prev[0] = Float3{ -velocity, 0.f, 0.f };
			};
		return FinishStabilityTest(environment, "doPoolBenchmark", envmode, std::move(job));
	}
	TestRoutine SinglebondForceAndPotentialSanityCheck(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "Singlebond";
		constexpr float expectedB0 = 0.133499995f; // [nm]
		constexpr float bondLengthError = 0.02f;   // (r-r0) [nm]

		// Displace the second carbon from equilibrium along x.
		auto job = MakeJob(workDir, envmode);
		job.preprocess = [](GroFile& grofile, TopologyFile&, SimParams& params) {
			params.n_steps = 1; params.data_logging_interval = 1;
			grofile.atoms[1].position = grofile.atoms[0].position + Float3{ expectedB0 + bondLengthError, 0.f, 0.f };
		};
		auto completed = co_await environment.Submit(std::move(job));
		auto& simulation = *completed.simulation;
		const auto parameters = simulation.box->bondgroups[0].singlebonds[0].params;

		// Analytic harmonic-bond force and potential.
		const double kB = parameters.kb / 2.; // [J/mol/nm^2]
		const Float3 expectedForce = Float3{ 1, 0, 0 } * 2.f * kB * bondLengthError; // [J/mol/nm]
		const float expectedPotential = kB * bondLengthError * bondLengthError;      // [J/mol]

		// Potential energy is split between the two particles, so sum both entries.
		const Float3 actualForce = simulation.box->pclusterInterimStates[0].forces_prev[0];
		const float actualPotential = simulation.potE_buffer->GetDatapoint(0, 0, 0)
			+ simulation.potE_buffer->GetDatapoint(0, 1, 0);
		const float forceError = (actualForce - expectedForce).len() / expectedForce.len();
		if (forceError >= 0.0001f)
			co_return LimaUnittestResult{ false, std::format("Force error {:.2e}", forceError), envmode == Full };
		const float potentialError = std::abs(actualPotential - expectedPotential) / expectedPotential;
		co_return LimaUnittestResult{ potentialError < 0.0001f,
			std::format("Potential error {:.2e}", potentialError), envmode == Full };
	}

	TestRoutine SinglebondOscillationTest(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "Singlebond";
		constexpr float expectedB0 = 0.133499995f; // [nm]
		constexpr float bondLengthError = 0.04f;   // (r-r0) [nm]

		// Simulate 1000 fs and record every step so oscillations can be counted.
		auto job = MakeJob(workDir, envmode);
		job.preprocess = [](GroFile& grofile, TopologyFile&, SimParams& params) {
			params.dt = 1.f * FEMTO_TO_NANO; params.n_steps = 1000; params.data_logging_interval = 1;
			grofile.atoms[1].position.x = grofile.atoms[0].position.x + expectedB0 + bondLengthError;
		};
		auto completed = co_await environment.Submit(std::move(job));
		auto& simulation = *completed.simulation;
		const auto parameters = simulation.box->bondgroups[0].singlebonds[0].params;
		const float massA = simulation.box->persistentClustersMetadata[0].mass[0];
		const float massB = simulation.box->persistentClustersMetadata[0].mass[1];
		const double reducedMass = massA * massB / (massA + massB); // [kg/mol]
		const double springConstant = parameters.kb / NANO / NANO;  // [J/(mol m^2)]
		const double expectedFrequency = std::sqrt(springConstant / reducedMass) / (2.f * PI) * FEMTO; // [1/fs]

		std::vector<float> lengths(simulation.simParams.n_steps); // [nm]
		for (int i = 0; i < simulation.simParams.n_steps; i++)
			lengths[i] = (simulation.traj_buffer->GetDatapoint(0, 0, i) - simulation.traj_buffer->GetDatapoint(0, 1, i)).len();
		const float elapsed = simulation.simParams.dt * simulation.simParams.n_steps * NANO_TO_FEMTO; // [fs]
		const float actualFrequency = static_cast<float>(SimAnalysis::CountOscillations(lengths)) / elapsed; // [1/fs]
		const float error = std::abs(actualFrequency - expectedFrequency) / expectedFrequency;
		co_return LimaUnittestResult{ error < 1e-2f,
			std::format("freq: {:.2e} / {:.2e} [1/fs]", actualFrequency, expectedFrequency), envmode == Full };
	}

	struct ExpectedForceEnergy { Float3 force{}; float potential = 0.f; };

	TestRoutine UreyBradleyForceAndPotentialSanityCheck(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "Anglebond";
		auto expected = std::make_shared<ExpectedForceEnergy>();
		auto job = MakeJob(workDir, envmode);
		job.preprocess = [](GroFile&, TopologyFile&, SimParams& params) { params.n_steps = 1; params.data_logging_interval = 1; };
		job.configureSimulation = [expected](Simulation& simulation) {
			// This used to require creating the simulation twice: once to discover
			// force-field parameters and again after adjusting the coordinates.
			// configure() runs after construction and gives direct access to both.
			auto& box = *simulation.box;
			const auto single = box.bondgroups[0].singlebonds[0].params;
			const auto angle = box.bondgroups[0].anglebonds[0].params;
			constexpr float angleError = 0.1f; // [rad]

			// First equilibrate both single bonds, then introduce only the angle error.
			auto& cluster = box.persistentClusters[0];
			cluster.pqd[1].position = Float3{};
			cluster.pqd[0].position = Float3{ single.b0, 0.f, 0.f };
			cluster.pqd[2].position = Float3::rodriguesRotatation(Float3{ single.b0, 0.f, 0.f },
				Float3{ 0.f, 1.f, 0.f }, -(angle.theta0 + angleError));
			const Float3 p0 = cluster.pqd[0].position, p1 = cluster.pqd[1].position, p2 = cluster.pqd[2].position;

			// Angular component.
			const Float3 angleForce = Float3{ 0.f, 0.f, 1.f }
				* (angle.kTheta * angleError / (p0 - p1).len()); // [J/mol/nm]
			const float anglePotential = angle.kTheta * angleError * angleError * .5f; // [J/mol]

			// Urey-Bradley 1-3 distance component.
			const float ubError = (p0 - p2).len() - angle.ub0; // [nm]
			const Float3 ubForce = (p0 - p2).norm() * -angle.kUB * ubError; // [J/mol/nm]
			const float ubPotential = angle.kUB * ubError * ubError * .5f;  // [J/mol]

			expected->force = angleForce + ubForce;
			expected->potential = anglePotential + ubPotential;
		};
		auto completed = co_await environment.Submit(std::move(job));
		auto& simulation = *completed.simulation;
		// Potential energy is distributed over all three atoms.
		const Float3 actualForce = simulation.box->pclusterInterimStates[0].forces_prev[0]; // [J/mol/nm]
		const float actualPotential = simulation.potE_buffer->GetDatapoint(0, 0, 0)
			+ simulation.potE_buffer->GetDatapoint(0, 1, 0) + simulation.potE_buffer->GetDatapoint(0, 2, 0);
		const float forceError = (actualForce - expected->force).len() / expected->force.len();
		const float potentialError = std::abs(actualPotential - expected->potential) / expected->potential;
		co_return LimaUnittestResult{ forceError < 0.0001f && potentialError < 0.0001f,
			std::format("Force error {:.2e}, potential error {:.2e}", forceError, potentialError), envmode == Full };
	}

	TestRoutine PairbondForceAndPotentialSanityCheck(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "Pairbond";
		auto expected = std::make_shared<ExpectedForceEnergy>();
		auto particleId = std::make_shared<int>(0);
		auto job = MakeJob(workDir, envmode);
		job.preprocess = [](GroFile&, TopologyFile&, SimParams& params) { params.n_steps = 1; params.data_logging_interval = 1; };
		job.configureSimulation = [expected, particleId](Simulation& simulation) {
			auto& group = simulation.box->bondgroups[0];

			// Isolate the pairbond by disabling the other bonded interactions.
			group.nSinglebonds = 0;
			group.nDihedralbonds = 0;
			const auto& pair = group.pairbonds[0];
			const int p0 = group.particles[pair.atom_indexes[0]].pid;
			const int p1 = group.particles[pair.atom_indexes[1]].pid;
			*particleId = p0;
			const Float3 diff = simulation.box->persistentClusters[0].pqd[p1].position
				- simulation.box->persistentClusters[0].pqd[p0].position;
			// Analytic Lennard-Jones 1-4 force and the per-particle potential.
			const float s = std::pow(pair.params.sigma / diff.len(), 6.f);
			const float forceScalar = 24.f * pair.params.epsilon * s
				/ diff.lenSquared() * (1.f - 2.f * s);
			expected->force = diff * forceScalar; // [J/mol/nm]
			expected->potential = 4.f * pair.params.epsilon * s * (s - 1.f) * .5f; // [J/mol]
		};
		auto completed = co_await environment.Submit(std::move(job));
		auto& simulation = *completed.simulation;
		const Float3 actualForce = simulation.forceBuffer->GetDatapoint(0, *particleId, 0);
		const float actualPotential = simulation.potE_buffer->GetDatapoint(0, *particleId, 0);
		const float forceError = (actualForce - expected->force).len() / expected->force.len();
		const float potentialError = std::abs(actualPotential - expected->potential) / expected->potential;
		co_return LimaUnittestResult{ forceError < 0.0001f && potentialError < 0.0001f,
			std::format("Force error {:.2e}, potential error {:.2e}", forceError, potentialError), envmode == Full };
	}



	TestRoutine doPoolCompSolBenchmark(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "PoolCompSol";
		constexpr float mass = 12.011000f / 1000.f; // [kg/mol]
		const float velocity = PhysicsUtils::tempToVelocity(400.f, mass); // [m/s] == [nm/ns]
		auto job = MakeJob(workDir, envmode, true);
		job.preprocess = [velocity](GroFile&, TopologyFile&, SimParams& params) {
			params.data_logging_interval = 1;
			params.n_steps = LIMA_UTILS::roundUp(6000000 / static_cast<int>(velocity), 100);
		};
		job.configureSimulation = [velocity](Simulation& simulation) {
			simulation.box->pclusterInterimStates[0].vels_prev[0] = Float3{ velocity, 0.f, 0.f };
		};
		return FinishStabilityTest(environment, "doPoolCompSolBenchmark", envmode, std::move(job));
	}

	TestRoutine doSinglebondBenchmark(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "Singlebond";
		auto job = MakeJob(workDir, envmode, true);
		job.preprocess = [](GroFile& grofile, TopologyFile&, SimParams& params) {
			params.data_logging_interval = 1; params.n_steps = 5000;
			constexpr float equilibriumLength = .1335f; // [nm]
			constexpr float bondLengthError = .02f;     // (r-r0) [nm]
			grofile.atoms[1].position.x += equilibriumLength + bondLengthError;
		};
		return FinishStabilityTest(environment, "doSinglebondBenchmark", envmode, std::move(job));
	}

	TestRoutine doAnglebondBenchmark(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "Anglebond";
		auto job = MakeJob(workDir, envmode, true);
		job.preprocess = [](GroFile& grofile, TopologyFile&, SimParams&) {
			constexpr float relaxedAngle = 1.8849f; // [rad]
			constexpr float angleError = .5f;       // (theta-theta0) [rad]
			grofile.atoms[2].position.rotateAroundOrigo(
				Float3{ 0.f, relaxedAngle + angleError, 0.f });
			// Keep every atom comfortably inside the test box.
			for (auto& atom : grofile.atoms) atom.position += Float3{ 3.f };
		};
		return FinishStabilityTest(environment, "doAnglebondBenchmark", envmode, std::move(job));
	}

	TestRoutine doDihedralbondBenchmark(Environment& environment, EnvMode envmode) {
		return LoadAndRunBasicSimulation(environment, envmode, "Dihedralbond", "doDihedralbondBenchmark");
	}

	TestRoutine doImproperDihedralBenchmark(Environment& environment, EnvMode envmode) {
		const fs::path workDir = TestsDir() / "Improperbond";
		auto job = MakeJob(workDir, envmode, true);
		job.preprocess = [](GroFile& grofile, TopologyFile& topfile, SimParams&) {
			const auto ids = topfile.GetMoleculeType().improperdihedralbonds[0].ids;
			// Translate the three vectors so atom i is the origin.
			const Float3 i = grofile.atoms[ids[0]].position;
			const Float3 j = grofile.atoms[ids[1]].position - i;
			const Float3 k = grofile.atoms[ids[2]].position - i;
			const Float3 l = grofile.atoms[ids[3]].position - i;
			// Rotate atom l about an axis in the i-j-k plane, introducing a
			// controlled improper-dihedral error of 0.4 rad.
			const Float3 axis = (j.cross(k).norm().cross(l.norm())).norm();
			const Float3 point = l / l.len();
			grofile.atoms[ids[3]].position += (Float3::rodriguesRotatation(point, axis, .4f) - point) * l.len();
		};
		return FinishStabilityTest(environment, "doImproperDihedralBenchmark", envmode, std::move(job));
	}
}

namespace StressTesting {}

namespace VerletintegrationTesting {
	using namespace TestUtils;

	// Apply a constant electric force and verify the resulting kinetic energy.
	TestRoutine TestIntegration(Environment& environment, EnvMode envmode) {
		const fs::path workDir = AutomatedTestsDir() / "Pool";
		constexpr float fieldStrength = .5f; // [V/nm]
		auto job = ForceCorrectness::MakeJob(workDir, envmode, true);
		job.simParams = SimParams{};
		job.simParams->n_steps = 1000;
		job.simParams->enable_electrostatics = true;
		job.simParams->data_logging_interval = 1;
		job.simParams->snf_select.insert(HorizontalChargeField);
		job.preprocess = [](GroFile& grofile, TopologyFile& topfile, SimParams&) {
			grofile.atoms.pop_back(); topfile.GetMoleculeType().atoms.pop_back();
		};
		job.configureSimulation = [](Simulation& simulation) {
			simulation.box->uniformElectricField = UniformElectricField{ Float3{ 1.f, 0.f, 0.f }, fieldStrength };
		};
		auto completed = co_await environment.Submit(std::move(job));
		const auto& simulation = *completed.simulation;
		const float charge = simulation.box->persistentClusters[0].pqd[0].params.charge * KILO; // [C/mol]
		const float mass = simulation.box->persistentClustersMetadata[0].mass[0]; // [kg/mol]
		const double elapsed = simulation.simParams.dt
			* static_cast<double>(simulation.simParams.n_steps) * NANO; // [s]
		const float expectedVelocity = charge * fieldStrength / NANO * elapsed / mass; // [m/s]
		const float expectedEnergy = PhysicsUtils::calcKineticEnergy(expectedVelocity, mass); // [J/mol]
		const float actualEnergy = completed.analysis->kin_energy.back();
		co_return LimaUnittestResult{ std::abs(actualEnergy - expectedEnergy) / expectedEnergy < .01f,
			std::format("Expected KE: {:.2e} Actual KE: {:.2e}", expectedEnergy, actualEnergy), envmode == Full };
	}
}
