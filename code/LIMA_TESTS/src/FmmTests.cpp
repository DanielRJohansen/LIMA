#include "Tests.h"
#include "Engine.cuh"
#include "BatchCompatibility.h"
#include "PhysicsUtils.cuh"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>
#include <random>

namespace FmmTests {
namespace {
	constexpr double pi = std::numbers::pi;
	void Require(bool condition, const char* message) {
		if (!condition) throw std::runtime_error(message);
	}

	std::unique_ptr<Simulation> MakeSimulation(int count, int order = 6, Float3 boxSize = Float3{4.f}, bool charged = false) {
		SimParams params;
		params.longrange_method = LongRangeElectrostaticsMethod::FMM;
		params.fmm_order = order;
		params.fmm_leaf_particles = 8;
		params.n_steps = 1;
		params.dt = 1e-8f;
		params.data_logging_interval = 1;
		params.stepsPerNlistupdate = 2;
		auto sim = std::make_unique<Simulation>(params, std::make_unique<Box>(boxSize));
		auto& box = *sim->box;
		box.boxparams.totalParticles = count;
		box.boxparams.degreesOfFreedom = count * 3;
		box.persistentClusters.resize(count);
		box.persistentClustersMetadata.resize(count);
		box.pclusterInterimStates.resize(count);
		box.particlesBondedToParticle.resize(count);
		box.pclustersBondedToPcluster.resize(count);
		std::mt19937 random(42);
		std::uniform_real_distribution<float> uniform(.08f, .92f);
		const int side = int(std::ceil(std::cbrt(double(count))));
		for (int i = 0; i < count; ++i) {
			const Float3 position{(i % side + uniform(random)) * boxSize.x / side,
				((i / side) % side + uniform(random)) * boxSize.y / side,
				(i / (side * side) + uniform(random)) * boxSize.z / side};
			box.persistentClusters[i].pqd[0] = PData{position, NBParams{0.f, 0.f, (charged || i % 2 == 0 ? 1.f : -1.f) * elementaryChargeToKiloCoulombPerMole}};
			auto& meta = box.persistentClustersMetadata[i];
			meta.particleIdsGlobal[0] = i;
			meta.mass[0] = .012f;
			++meta.nParticles;
		}
		sim->PrepareDataBuffers();
		return sim;
	}

	struct Reference {
		std::vector<std::array<double, 3>> forces;
		double energy = 0.;
	};

	// Independent, double precision periodic Ewald sum, including a neutralizing background.
	// A stride greater than one evaluates sampled target forces; its energy is undefined.
	Reference Ewald(const Simulation& sim, int targetStride = 1) {
		const auto& box = *sim.box;
		const int count = box.boxparams.totalParticles;
		const Float3 size = NodeIndex(box.boxparams.boxSize).toFloat3();
		const double volume = double(size.x) * size.y * size.z;
		const double alpha = 2. / std::min({size.x, size.y, size.z});
		const double coulomb = PhysicsUtils::modifiedCoulombConstant;
		Reference result{std::vector<std::array<double, 3>>(count)};
		const auto particle = [&](int i) -> const PData& { return box.persistentClusters[i].pqd[0]; };
		double totalCharge = 0.;
		for (int i = 0; i < count; ++i) totalCharge += particle(i).params.charge;
		for (int i = 0; i < count; i += targetStride) {
			const auto& a = particle(i);
			const double qi = a.params.charge;
			result.energy -= coulomb * qi * qi * alpha / std::sqrt(pi);
			for (int j = 0; j < count; ++j) {
				const auto& b = particle(j);
				const double product = coulomb * qi * b.params.charge;
				for (int z = -4; z <= 4; ++z) for (int y = -4; y <= 4; ++y) for (int x = -4; x <= 4; ++x) {
					if (i == j && x == 0 && y == 0 && z == 0) continue;
					const double r[3]{double(a.position.x) - b.position.x + x * size.x,
						double(a.position.y) - b.position.y + y * size.y, double(a.position.z) - b.position.z + z * size.z};
					const double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2], distance = std::sqrt(r2);
					if (alpha * distance > 7.) continue;
					const double erfc = std::erfc(alpha * distance);
					result.energy += .5 * product * erfc / distance;
					const double scalar = product * (erfc / distance + 2. * alpha / std::sqrt(pi) * std::exp(-alpha * alpha * r2)) / r2;
					for (int d = 0; d < 3; ++d) result.forces[i][d] += scalar * r[d];
				}
				if (i != j && box.particlesBondedToParticle[i].Contains(j)) {
					double r[3];
					for (int d = 0; d < 3; ++d) {
						r[d] = double(a.position[d]) - b.position[d];
						r[d] -= std::round(r[d] / size[d]) * size[d];
					}
					const double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2], distance = std::sqrt(r2);
					result.energy -= .5 * product / distance;
					for (int d = 0; d < 3; ++d) result.forces[i][d] -= product * r[d] / (r2 * distance);
				}
			}
		}
		result.energy -= pi * coulomb * totalCharge * totalCharge / (2. * alpha * alpha * volume);
		const double maxWave = 12. * alpha;
		const int nx = int(std::ceil(maxWave * size.x / (2. * pi)));
		const int ny = int(std::ceil(maxWave * size.y / (2. * pi)));
		const int nz = int(std::ceil(maxWave * size.z / (2. * pi)));
		for (int z = -nz; z <= nz; ++z) for (int y = -ny; y <= ny; ++y) for (int x = -nx; x <= nx; ++x) {
			if (x == 0 && y == 0 && z == 0) continue;
			const double k[3]{2. * pi * x / size.x, 2. * pi * y / size.y, 2. * pi * z / size.z};
			const double k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
			if (k2 > maxWave * maxWave) continue;
			const double weight = coulomb * 4. * pi / volume * std::exp(-k2 / (4. * alpha * alpha)) / k2;
			double real = 0., imaginary = 0.;
			for (int i = 0; i < count; ++i) {
				const auto& a = particle(i);
				const double phase = k[0] * a.position.x + k[1] * a.position.y + k[2] * a.position.z;
				real += a.params.charge * std::cos(phase);
				imaginary -= a.params.charge * std::sin(phase);
			}
			result.energy += .5 * weight * (real * real + imaginary * imaginary);
			for (int i = 0; i < count; i += targetStride) {
				const auto& a = particle(i);
				const double phase = k[0] * a.position.x + k[1] * a.position.y + k[2] * a.position.z;
				const double scalar = weight * a.params.charge * (std::sin(phase) * real + std::cos(phase) * imaginary);
				for (int d = 0; d < 3; ++d) result.forces[i][d] += scalar * k[d];
			}
		}
		if (targetStride > 1) result.energy = NAN;
		return result;
	}

	void Run(Simulation& sim) {
		Engine engine({&sim});
		while (!engine.IsFinished()) engine.step();
		engine.terminateSimulation();
	}

	double CheckReference(Simulation& sim, double forceTolerance, double energyTolerance) {
		const Reference reference = Ewald(sim);
		Run(sim);
		double forceNorm = 0., errorNorm = 0., energy = 0.;
		for (int i = 0; i < sim.box->boxparams.totalParticles; ++i) {
			const auto actual = sim.forceBuffer->GetDatapoint(i, 0, 0);
			energy += sim.potE_buffer->GetDatapoint(i, 0, 0);
			for (int d = 0; d < 3; ++d) {
				forceNorm += reference.forces[i][d] * reference.forces[i][d];
				errorNorm += std::pow(actual[d] - reference.forces[i][d], 2);
			}
		}
		const double forceError = std::sqrt(errorNorm / std::max(forceNorm, 1e-20));
		const double energyError = std::abs(energy - reference.energy) / std::max(std::abs(reference.energy), 1e-20);
		if (sim.simParams.longrange_method.value == LongRangeElectrostaticsMethod::FMM)
			std::cout << "FMM order " << sim.simParams.fmm_order.value;
		else std::cout << "PME";
		std::cout << ", N=" << sim.box->boxparams.totalParticles
			<< ": force RMS " << forceError << ", energy error " << energyError << '\n';
		Require(forceError < forceTolerance, "FMM forces differ from independent Ewald reference");
		Require(energyError < energyTolerance, "FMM energy differs from independent Ewald reference");
		return forceError;
	}
}

void RunAll() {
	{
		const auto path = std::filesystem::temp_directory_path() / "lima-fmm-simparams.txt";
		SimParams params;
		params.longrange_method = LongRangeElectrostaticsMethod::FMM;
		params.fmm_order = 8;
		params.DumpToFile(path);
		const SimParams parsed(path);
		Require(parsed.longrange_method.value == LongRangeElectrostaticsMethod::FMM && parsed.fmm_order == 8, "FMM parameters failed roundtrip");
		std::filesystem::remove(path);
		{
			std::ofstream file(path);
			file << "longrange_method=unknown\n";
		}
		bool rejected = false;
		try { const SimParams invalid(path); }
		catch (const std::runtime_error&) { rejected = true; }
		std::filesystem::remove(path);
		Require(rejected, "Unknown electrostatics method was accepted");
	}
	std::array<double, 5> errors{};
	for (int order = 4; order <= 8; ++order) {
		auto simulation = MakeSimulation(128, order);
		errors[order - 4] = CheckReference(*simulation, order == 4 ? .004 : .001, order == 4 ? .004 : .0005);
	}
	for (size_t i = 1; i < errors.size(); ++i)
		Require(errors[i] < errors[i - 1], "Increasing FMM order did not improve forces");
	{
		// More than sixteen parents exercises multiple spatial tiles in M2L.
		auto simulation = MakeSimulation(512, 6, Float3{8.f});
		simulation->simParams.fmm_leaf_particles = 4;
		CheckReference(*simulation, .001, .0005);
	}
	for (bool charged : {false, true}) {
		// Exercise the five-level production tree and non-binary coordinate scales.
		// Sampling targets keeps the independent periodic Ewald calculation bounded.
		auto simulation = MakeSimulation(4098, 6, charged ? Float3{22.f, 19.f, 23.f} : Float3{22.f}, charged);
		simulation->simParams.fmm_leaf_particles = 1;
		constexpr int targetStride = 128;
		const Reference reference = Ewald(*simulation, targetStride);
		Run(*simulation);
		double forceNorm = 0., errorNorm = 0.;
		for (int i = 0; i < simulation->box->boxparams.totalParticles; i += targetStride) {
			const auto actual = simulation->forceBuffer->GetDatapoint(i, 0, 0);
			for (int d = 0; d < 3; ++d) {
				forceNorm += reference.forces[i][d] * reference.forces[i][d];
				errorNorm += std::pow(actual[d] - reference.forces[i][d], 2);
			}
		}
		const double error = std::sqrt(errorNorm / forceNorm);
		std::cout << "FMM depth-5 " << (charged ? "charged rectangular " : "") << "sampled force RMS: " << error << '\n';
		Require(error < .001, "Deep FMM tree differs from independent Ewald reference");
	}
	{
		auto simulation = MakeSimulation(128);
		simulation->simParams.fmm_leaf_particles = SimParams{}.fmm_leaf_particles.value;
		CheckReference(*simulation, .004, .0005);
	}
	{
		auto simulation = MakeSimulation(65, 6, Float3{4.f, 5.f, 6.f}, true);
		CheckReference(*simulation, .004, .0005);
	}
	{
		auto simulation = MakeSimulation(32);
		auto& box = *simulation->box;
		box.persistentClusters[0].pqd[0].position = Float3{.01f, 2.f, 2.f};
		box.persistentClusters[1].pqd[0].position = Float3{3.99f, 2.f, 2.f};
		box.particlesBondedToParticle[0] = ParticlesBondedToParticle::Create({1, 31});
		box.particlesBondedToParticle[1] = ParticlesBondedToParticle::Create({0});
		box.particlesBondedToParticle[31] = ParticlesBondedToParticle::Create({0});
		CheckReference(*simulation, .004, .0005);
	}
	{
		auto simulation = MakeSimulation(32);
		simulation->simParams.em_variant = true;
		simulation->simParams.em_force_tolerance = 0.f;
		CheckReference(*simulation, .004, .0005);
	}
	for (int exponent : {-20, 10}) {
		// Normalized tensor inputs must handle charges outside FP16's range.
		auto simulation = MakeSimulation(32);
		simulation->simParams.em_variant = true;
		simulation->simParams.em_force_tolerance = 0.f;
		for (auto& cluster : simulation->box->persistentClusters)
			cluster.pqd[0].params.charge = std::ldexp(cluster.pqd[0].params.charge, exponent);
		CheckReference(*simulation, .004, .0005);
	}
	for (int smallCount : {17, 257}) for (bool em : {false, true}) {
		auto a = MakeSimulation(smallCount), b = MakeSimulation(512), refA = MakeSimulation(smallCount), refB = MakeSimulation(512);
		for (auto* simulation : {a.get(), refA.get(), b.get(), refB.get()}) {
			simulation->simParams.fmm_leaf_particles = 4;
			simulation->simParams.em_variant = em;
			simulation->simParams.em_force_tolerance = 0.f;
			simulation->simParams.n_steps = simulation->box->boxparams.totalParticles == smallCount ? 2 : 5;
			simulation->PrepareDataBuffers();
		}
		Run(*refA); Run(*refB);
		Engine engine({a.get(), b.get()});
		while (!engine.IsFinished()) engine.step();
		engine.terminateSimulation();
		Require(a->forceBuffer->GetBuffer() == refA->forceBuffer->GetBuffer(), "FMM batch changed first member forces");
		Require(b->forceBuffer->GetBuffer() == refB->forceBuffer->GetBuffer(), "FMM batch retirement changed survivor forces");
		Require(a->potE_buffer->GetBuffer() == refA->potE_buffer->GetBuffer(), "FMM batch changed first member energy");
		Require(b->potE_buffer->GetBuffer() == refB->potE_buffer->GetBuffer(), "FMM batch retirement changed survivor energy");
	}
	{
		// Force-only steps must produce the same trajectory as steps computing energy.
		auto everyStep = MakeSimulation(512), sparseLogging = MakeSimulation(512);
		for (auto* simulation : {everyStep.get(), sparseLogging.get()}) {
			simulation->simParams.fmm_leaf_particles = 4;
			simulation->simParams.n_steps = 5;
			simulation->simParams.dt = 1e-6f;
		}
		sparseLogging->simParams.data_logging_interval = 2;
		everyStep->PrepareDataBuffers(); sparseLogging->PrepareDataBuffers();
		Run(*everyStep); Run(*sparseLogging);
		for (int step = 0; step < 5; step += 2) for (int i = 0; i < 512; ++i) {
			Require(everyStep->forceBuffer->GetDatapoint(i, 0, step) == sparseLogging->forceBuffer->GetDatapoint(i, 0, step / 2),
				"FMM forces depend on potential-energy logging cadence");
			Require(everyStep->potE_buffer->GetDatapoint(i, 0, step) == sparseLogging->potE_buffer->GetDatapoint(i, 0, step / 2),
				"FMM energies depend on potential-energy logging cadence");
		}
	}
	{
		auto a = MakeSimulation(8), b = MakeSimulation(8);
		b->simParams.longrange_method = LongRangeElectrostaticsMethod::PME;
		Require(EngineBatch::FindIncompatibility(*a, *b) == "longrange_method", "Scheduler mixed PME and FMM");
	}
	{
		auto simulation = MakeSimulation(128);
		simulation->simParams.longrange_method = LongRangeElectrostaticsMethod::PME;
		CheckReference(*simulation, .005, .002);
	}
	{
		auto simulation = MakeSimulation(32);
		simulation->simParams.enable_electrostatics = false;
		Run(*simulation);
		for (const auto& force : simulation->forceBuffer->GetBuffer())
			Require(force.lenSquared() == 0.f, "Disabling electrostatics left FMM forces");
	}
	for (int order : {3, 9}) {
		auto simulation = MakeSimulation(8, order);
		bool rejected = false;
		try { Engine engine({simulation.get()}); }
		catch (const std::invalid_argument&) { rejected = true; }
		Require(rejected, "Engine accepted unsupported FMM order");
	}
	{
		auto empty = MakeSimulation(8), active = MakeSimulation(17), reference = MakeSimulation(17);
		empty->simParams.n_steps = 0;
		Run(*reference);
		Engine engine({empty.get(), active.get()});
		while (!engine.IsFinished()) engine.step();
		engine.terminateSimulation();
		Require(empty->getStep() == 0 && active->forceBuffer->GetBuffer() == reference->forceBuffer->GetBuffer(),
			"Zero-step FMM member affected active simulation");
	}
	{
		auto first = MakeSimulation(32), second = MakeSimulation(32);
		second->simParams.cutoff_nm = .8f;
		Run(*first); Run(*second);
		Require(first->forceBuffer->GetBuffer() == second->forceBuffer->GetBuffer(), "FMM Coulomb forces depend on LJ cutoff");
	}
	{
		auto simulation = MakeSimulation(64, 6, Float3{8.f});
		simulation->simParams.n_steps = 250;
		simulation->simParams.dt = 1e-6f;
		for (int i = 0; i < 64; ++i) {
			simulation->box->persistentClusters[i].pqd[0].params.charge =
				(i % 2 ? -1.f : 1.f) * elementaryChargeToKiloCoulombPerMole;
			simulation->box->pclusterInterimStates[i].vels_prev[0] = Float3{float(i % 2 ? -120 : 120), 80.f, -40.f};
		}
		// Cross a tree-cell boundary during the trajectory.
		simulation->box->persistentClusters[0].pqd[0].position.x = 3.99f;
		simulation->PrepareDataBuffers();
		Run(*simulation);
		double initial = 0., scale = 0., maxDrift = 0.;
		for (int step = 0; step < 250; ++step) {
			double potential = 0., kinetic = 0.;
			for (int i = 0; i < 64; ++i) {
				potential += simulation->potE_buffer->GetDatapoint(i, 0, step);
				const double velocity = simulation->vel_buffer->GetDatapoint(i, 0, step);
				kinetic += .5 * .012 * velocity * velocity;
			}
			if (step == 0) { initial = potential + kinetic; scale = std::abs(potential) + kinetic; }
			maxDrift = std::max(maxDrift, std::abs(potential + kinetic - initial) / scale);
		}
		std::cout << "FMM NVE maximum relative energy drift: " << maxDrift << '\n';
		Require(maxDrift < .001, "FMM NVE energy drift exceeds tolerance");
	}
	std::cout << "FMM tests passed\n";
}
}
