#include "Tests.h"
#include "Format.h"

#include "Engine.cuh"
#include "RenderDataPipe.h"
#include <stdexcept>
#include "Lipids.h"
#include "Workflow.h"
#include "MDFiles.h"
#include "BoxImageBuilder.h"
#include "BoxBuilder.cuh"
#include "SimulationBuilder.h"
#include "Analyzer.h"

#include <chrono>
#include <fstream>
#include <map>
#include <numeric>

namespace EngineBatchTests {
	void Require(bool value, const char* message) {
		if (!value) throw std::runtime_error(message);
	}

	std::unique_ptr<Simulation> MakeSimulation(int particles, int steps, float dt, bool electrostatics, bool em, int loggingInterval) {
		SimParams params;
		params.nsteps = steps;
		params.dt = dt;
		params.enable_electrostatics = electrostatics;
		params.em_variant = em;
		params.emtol = 0.f;
		params.nstlist = 3;
		params.nstenergy = loggingInterval;
		params.nsttcouple = 10;
		params.tcoupl = em ? TcouplSelect::No : TcouplSelect::VRescale;
		auto sim = std::make_unique<Simulation>(params, std::make_unique<Box>(Float3{4.f}));
		auto& box = *sim->box;
		box.boxparams.totalParticles = particles;
		box.boxparams.degreesOfFreedom = particles * 3;
		const int nPcs = (particles + 3) / 4;
		box.persistentClusters.resize(nPcs);
		box.persistentClustersMetadata.resize(nPcs);
		box.pclusterInterimStates.resize(nPcs);
		box.particlesBondedToParticle.resize(particles);
		box.pclustersBondedToPcluster.resize(nPcs);
		for (int id = 0; id < particles; ++id) {
			const int pc = id / 4, lane = id % 4;
			box.persistentClusters[pc].pqd[lane] = PData{
				Float3{0.5f + 0.3f * (id % 3), 0.5f + 0.3f * (id / 3), 0.5f},
				NBParams{0.1f, 0.1f, electrostatics ? (id % 2 ? -0.01f : 0.01f) : 0.f}};
			auto& meta = box.persistentClustersMetadata[pc];
			meta.particleIdsGlobal[lane] = id;
			meta.mass[lane] = 12.f;
			++meta.nParticles;
			box.pclusterInterimStates[pc].vels_prev[lane] = Float3{0.01f, 0.02f, 0.f};
		}
		// Exercise rebasing of topology, exclusions and per-particle gather references.
		BondGroup group;
		group.nParticles = 2;
		group.nSinglebonds = 1;
		box.bondgroups.groups.push_back(group);
		box.bondgroups.particles = {{0, 0}, {1, 0}};
		box.bondgroups.singlebonds.emplace_back(std::array<uint8_t, 2>{0, 1}, SingleBond::Parameters{0.28f, 100.f});
		box.persistentClustersMetadata[0].bondgroupReferences[0].Add({0, 0});
		box.persistentClustersMetadata[1].bondgroupReferences[0].Add({0, 1});
		box.particlesBondedToParticle[0] = ParticlesBondedToParticle::Create(std::set<int>{4});
		box.particlesBondedToParticle[4] = ParticlesBondedToParticle::Create(std::set<int>{0});
		box.pclustersBondedToPcluster[0] = PclustersBondedToPcluster::Create(std::set<int>{1});
		box.pclustersBondedToPcluster[1] = PclustersBondedToPcluster::Create(std::set<int>{0});
		sim->PrepareDataBuffers();
		return sim;
	}

	void Run(Simulation& sim) {
		Engine engine({&sim});
		while (!engine.IsFinished()) engine.step();
		engine.terminateSimulation();
	}

	double PotentialEnergyAtScale(bool electrostatics, double volumeScale) {
		auto simulation = MakeSimulation(29, 1, 1e-8f, electrostatics, false, 1);
		simulation->simParams.tcoupl = TcouplSelect::No;
		const float coordinateScale = static_cast<float>(std::cbrt(volumeScale));
		for (size_t pc = 0; pc < simulation->box->persistentClusters.size(); ++pc)
			for (int lane = 0; lane < simulation->box->persistentClustersMetadata[pc].nParticles; ++lane)
				simulation->box->persistentClusters[pc].pqd[lane].position *= coordinateScale;
		simulation->box->boxparams.SetBoxLength(simulation->box->boxparams.BoxSizeFloat() * coordinateScale);
		Run(*simulation);
		return std::accumulate(simulation->potE_buffer->GetBuffer().begin(), simulation->potE_buffer->GetBuffer().end(), 0.0);
	}

	double MeasuredPotentialPressure(bool electrostatics) {
		auto simulation = MakeSimulation(29, 1, 1e-8f, electrostatics, false, 1);
		simulation->simParams.tcoupl = TcouplSelect::No;
		simulation->simParams.nsttcouple = 1;
		simulation->simParams.pcoupl = PcouplSelect::Berendsen;
		simulation->simParams.nstpcouple = 1;
		simulation->simParams.tau_p = 1.f;
		simulation->simParams.compressibility = 1e-12f;
		const Float3 length = simulation->box->boxparams.BoxSizeFloat();
		const double volume = double(length.x) * length.y * length.z;
		const int64_t degreesOfFreedom = simulation->box->boxparams.degreesOfFreedom;
		Run(*simulation);
		Require(simulation->pressure_buffer.size() == 1 && simulation->temperature_buffer.size() == 1,
			"Pressure finite-difference test did not get one pressure and temperature");
		constexpr double gasConstant = BOLTZMANNCONSTANT * AVOGADROSNUMBER;
		constexpr double pressureUnit = 1e22 / AVOGADROSNUMBER;
		const double kineticEnergy = 0.5 * degreesOfFreedom * gasConstant * simulation->temperature_buffer.front();
		return simulation->pressure_buffer.front() / pressureUnit - 2.0 * kineticEnergy / (3.0 * volume);
	}

	void CheckPressureAgainstFiniteDifference(bool electrostatics) {
		constexpr double epsilon = 0.005;
		constexpr double volume = 4.0 * 4.0 * 4.0;
		const double lowerEnergy = PotentialEnergyAtScale(electrostatics, 1.0 - epsilon);
		const double upperEnergy = PotentialEnergyAtScale(electrostatics, 1.0 + epsilon);
		const double expected = -(upperEnergy - lowerEnergy) / (2.0 * epsilon * volume);
		const double actual = MeasuredPotentialPressure(electrostatics);
		const double relativeError = std::abs(actual - expected) / (std::max)(1.0, std::abs(expected));
		if (relativeError > 0.05)
			throw std::runtime_error(Lima::Format("{} pressure virial differs from finite-difference energy: {} vs {} ({:.1f}%)",
				electrostatics ? "PME" : "Non-PME", actual, expected, relativeError * 100.0));
	}

	void Compare(Simulation& a, Simulation& b) {
		Require(a.getStep() == b.getStep(), "Batch changed simulation step count");
		for (size_t pc = 0; pc < a.box->persistentClusters.size(); ++pc) {
			for (int lane = 0; lane < 4; ++lane) {
				Require(a.box->persistentClusters[pc].pqd[lane].position == b.box->persistentClusters[pc].pqd[lane].position,
					"Batch changed particle coordinates");
				Require(a.box->pclusterInterimStates[pc].vels_prev[lane] == b.box->pclusterInterimStates[pc].vels_prev[lane],
					"Batch changed particle velocity");
				Require(a.box->pclusterInterimStates[pc].forces_prev[lane] == b.box->pclusterInterimStates[pc].forces_prev[lane],
					"Batch changed particle forces");
			}
		}
		Require(a.traj_buffer->GetBuffer() == b.traj_buffer->GetBuffer(), "Batch changed logged coordinates");
		if (a.potE_buffer->GetBuffer() != b.potE_buffer->GetBuffer()) {
			const auto& actual = a.potE_buffer->GetBuffer();
			const auto& expected = b.potE_buffer->GetBuffer();
			for (size_t i = 0; i < std::min(actual.size(), expected.size()); ++i)
				if (actual[i] != expected[i]) throw std::runtime_error(Lima::Format("Batch changed logged energy at {}: {} vs {}", i, actual[i], expected[i]));
			throw std::runtime_error("Batch changed logged energy buffer length");
		}
		Require(a.forceBuffer->GetBuffer() == b.forceBuffer->GetBuffer(), "Batch changed logged forces");
		Require(a.temperature_buffer == b.temperature_buffer, "Batch mixed temperature reductions");
		Require(a.pressure_buffer == b.pressure_buffer, "Batch changed pressure reductions");
	}

	// The synthetic systems above are smaller than the cutoff and barely move, so re-clustering never changes their
	// neighbor lists. A real solvated system is needed to see a rebuild that perturbs the survivors.
	std::unique_ptr<Simulation> MakeMetsol(int steps, bool em) {
		const auto fixtures = FileUtils::GetLimaDir() / "tests" / "clitests" / "fixtures";
		SimParams params = em ? SimParams::BasicEMSimParams(0.f) : SimParams{};
		params.nsteps = steps;
		params.nstenergy = 20;
		params.nsttcouple = 50;
		params.tcoupl = em ? TcouplSelect::No : TcouplSelect::VRescale;
		GroFile gro{ fixtures / "metsol.gro" };
		TopologyFile top{ fixtures / "metsol.top" };
		auto image = LIMA_MOLECULEBUILD::buildMolecules(gro, top, params);
		auto sim = std::make_unique<Simulation>(params, BoxBuilder::BuildBox(params, *image, EnvMode::Headless));
		sim->PrepareDataBuffers();
		if (em) Engine::PrepareEnergyMinimization(*sim);
		return sim;
	}

	// A member retiring on a step that is not an nlist update step (nsteps, StopSimulation as used for maxRunTime)
	// must leave the survivor bitwise identical to a standalone run
	void RetirementOffNlistStep() {
		constexpr int survivorSteps = 300;
		constexpr int retireStep = 137; // Not a multiple of nstlist
		for (bool em : {false, true}) {
			for (bool stop : {false, true}) {
				const auto label = Lima::Format("EM {} {}", em, stop ? "StopSimulation" : "nsteps");
				auto retiring = MakeMetsol(stop ? survivorSteps : retireStep, em);
				auto survivor = MakeMetsol(survivorSteps, em);
				auto referenceRetiring = MakeMetsol(retireStep, em);
				auto referenceSurvivor = MakeMetsol(survivorSteps, em);
				Require(retireStep % survivor->simParams.nstlist != 0, "Retirement step is an nlist update step");
				Run(*referenceRetiring);
				Run(*referenceSurvivor);
				Engine engine({ retiring.get(), survivor.get() });
				while (!engine.IsFinished()) {
					engine.step();
					if (stop && retiring->getStep() == retireStep && !engine.GetRunStatus(0).simulation_finished)
						engine.StopSimulation(0);
				}
				engine.terminateSimulation();
				try { Compare(*survivor, *referenceSurvivor); }
				catch (const std::exception& error) { throw std::runtime_error(Lima::Format("{} survivor: {}", label, error.what())); }
				Require(retiring->getStep() == retireStep, "Retired member ran past its retirement step");
				for (size_t pc = 0; pc < retiring->box->persistentClusters.size(); ++pc)
					for (int lane = 0; lane < 4; ++lane)
						if (retiring->box->persistentClusters[pc].pqd[lane].position != referenceRetiring->box->persistentClusters[pc].pqd[lane].position)
							throw std::runtime_error(Lima::Format("{} retired member: final coordinates differ", label));
			}
		}
	}

	void RunAll() {
		RetirementOffNlistStep();
		const Int3 awkwardGrid = BoxParams::GridSizeForLength({2.8f, 3.8f, 4.9f});
		Require(awkwardGrid == Int3{2, 3, 4} && !GridSizeRequiresRebuild(awkwardGrid, awkwardGrid),
			"A freshly selected engine grid requested another rebuild");
		Require(GridSizeRequiresRebuild({2, 3, 4}, {3, 3, 4}), "Engine grid hysteresis retained a substantially undersized grid");
		CheckPressureAgainstFiniteDifference(false);
		CheckPressureAgainstFiniteDifference(true);
		{
			auto noPressure = MakeSimulation(5, 1, 0.00001f, false);
			noPressure->simParams.pcoupl = PcouplSelect::No;
			noPressure->simParams.nstpcouple = -1;
			Run(*noPressure);
		}
		{
			auto first = MakeSimulation(29, 20, 0.000002f, true);
			auto second = MakeSimulation(29, 20, 0.000002f, true);
			for (Simulation* simulation : {first.get(), second.get()}) {
				simulation->simParams.pcoupl = PcouplSelect::Berendsen;
				simulation->simParams.nstpcouple = 2;
				simulation->simParams.tau_p = 0.2f;
				simulation->simParams.compressibility = 4.5e-5f;
				simulation->simParams.constraints = ConstraintsSelect::HBonds;
				auto& box = *simulation->box;
				box.builtConstraints = ConstraintsSelect::HBonds;
				box.pclusterConstraints.resize(box.persistentClusters.size());
				box.pclusterConstraints[0].nConstraints = 1;
				box.pclusterConstraints[0].particles[0] = 0 | (1 << 2);
				box.pclusterConstraints[0].distance[0] = 0.3f;
			}
			const Float3 initialLength = first->box->boxparams.BoxSizeFloat();
			Run(*first);
			Run(*second);
			Require(first->pressure_buffer.size() == 10, "Barostat did not measure pressure at its own cadence");
			Require(std::ranges::all_of(first->pressure_buffer, [](float pressure) { return std::isfinite(pressure); }),
				"Barostat produced non-finite pressure");
			Require(first->box->boxparams.BoxSizeFloat() != initialLength, "Barostat did not scale the box");
			Require(first->box->boxparams.BoxSizeFloat() == second->box->boxparams.BoxSizeFloat(), "Barostat box scaling was not deterministic");
			Compare(*first, *second);
		}
		{
			auto simulation = MakeSimulation(5, 1, 0.00001f, false);
			RenderDataPipe renderDataPipe;
			Engine engine({ simulation.get() }, EngineRunMode::Simulation, { &renderDataPipe });
			while (!engine.IsFinished()) engine.step();
			engine.terminateSimulation();
			Require(renderDataPipe.GetState() == RenderDataPipe::State::Stopped,
				"Terminated engine left its render pipe active");
		}
		{
			auto first = MakeSimulation(29, 11, 0.00001f, true);
			auto middle = MakeSimulation(5, 3, 0.00002f, true);
			auto last = MakeSimulation(29, 13, 0.00003f, true);
			auto referenceFirst = MakeSimulation(29, 11, 0.00001f, true);
			auto referenceMiddle = MakeSimulation(5, 3, 0.00002f, true);
			auto referenceLast = MakeSimulation(29, 13, 0.00003f, true);
			Run(*referenceFirst);
			Run(*referenceMiddle);
			Run(*referenceLast);
			auto stepwiseFirst = MakeSimulation(29, 11, 0.00001f, true);
			Engine stepwiseEngine({stepwiseFirst.get()});
			Engine engine({first.get(), middle.get(), last.get()});
			while (!engine.IsFinished()) {
				engine.step();
				if (!stepwiseEngine.IsFinished()) stepwiseEngine.step();
				const int step = first->getStep();
				if (step > 0) {
					engine.CopySimulationToHost(0);
					stepwiseEngine.CopySimulationToHost(0);
					for (size_t pc = 0; pc < first->box->pclusterInterimStates.size(); ++pc)
						for (int lane = 0; lane < 4; ++lane)
							if (first->box->pclusterInterimStates[pc].vels_prev[lane] != stepwiseFirst->box->pclusterInterimStates[pc].vels_prev[lane]
								|| first->box->pclusterInterimStates[pc].forces_prev[lane] != stepwiseFirst->box->pclusterInterimStates[pc].forces_prev[lane])
								throw std::runtime_error(Lima::Format("First survivor state diverged at step {}, pc {}, lane {}", step, pc, lane));
					const auto state = GenericCopyToHost(engine.OffloadPclusterState(0).Get(), first->box->persistentClusters.size());
					const auto expected = GenericCopyToHost(stepwiseEngine.OffloadPclusterState(0).Get(), first->box->persistentClusters.size());
					for (size_t pc = 0; pc < state.size(); ++pc)
						for (int lane = 0; lane < 4; ++lane)
							if (state[pc].pqd[lane].Valid() && state[pc].pqd[lane].position != expected[pc].pqd[lane].position)
								throw std::runtime_error(Lima::Format("First survivor diverged at step {}, pc {}, lane {}", step, pc, lane));
				}
			}
			engine.terminateSimulation();
			try { Compare(*first, *referenceFirst); }
			catch (const std::exception& error) { throw std::runtime_error(std::string("First survivor: ") + error.what()); }
			try { Compare(*middle, *referenceMiddle); }
			catch (const std::exception& error) { throw std::runtime_error(std::string("Retired middle: ") + error.what()); }
			try { Compare(*last, *referenceLast); }
			catch (const std::exception& error) { throw std::runtime_error(std::string("Last survivor: ") + error.what()); }
		}
		for (bool electrostatics : {false, true}) {
			for (bool em : {false, true}) {
				auto a = MakeSimulation(5, 5, 0.00001f, electrostatics, em);
				auto b = MakeSimulation(29, 13, 0.00002f, electrostatics, em);
				auto referenceA = MakeSimulation(5, 5, 0.00001f, electrostatics, em);
				auto referenceB = MakeSimulation(29, 13, 0.00002f, electrostatics, em);
				Run(*referenceA);
				Run(*referenceB);
				Engine engine({a.get(), b.get()});
				while (!engine.IsFinished()) {
					engine.step();
					if (engine.GetRunStatus(0).simulation_finished) {
						try { Compare(*a, *referenceA); }
						catch (const std::exception& error) { throw std::runtime_error(Lima::Format("Batch electrostatics {} EM {} step {}: {}", electrostatics, em, b->getStep(), error.what())); }
					}
				}
				engine.terminateSimulation();
				engine.terminateSimulation(); // Finalization must be idempotent.
				engine.step(); // A completed batch must not advance.
				try { Compare(*a, *referenceA); Compare(*b, *referenceB); }
				catch (const std::exception& error) { throw std::runtime_error(Lima::Format("Batch electrostatics {} EM {} final: {}", electrostatics, em, error.what())); }
				Require(a->getStep() == 5 && b->getStep() == 13, "Batch ran past a step limit");
				const auto stopped = GenericCopyToHost(engine.OffloadPclusterState(0).Get(), a->box->persistentClusters.size());
				for (size_t pc = 0; pc < stopped.size(); ++pc)
					for (int lane = 0; lane < 4; ++lane)
						Require(stopped[pc].pqd[lane].position == referenceA->box->persistentClusters[pc].pqd[lane].position,
							"Retired simulation was still integrated on the GPU");
			}
		}
		for (int steps : {0, 1, 9, 10, 11}) {
			auto sim = MakeSimulation(5, steps, 0.00001f, false);
			Run(*sim);
			Require(sim->getStep() == steps, "Size-one batch mishandled step limit");
			if (steps > 0) {
				const int lastEntry = (steps - 1) / 2;
				Require(sim->traj_buffer->GetDatapoint(0, 0, lastEntry).lenSquared() > 0.f, "Final logging sample missing");
			}
		}
		auto disabledLogging = MakeSimulation(5, 7, 0.00001f, false, false, 0);
		Run(*disabledLogging);
		Require(disabledLogging->getStep() == 7, "Disabled logging broke completion");
		auto early = MakeSimulation(5, 200, 0.00001f, false, true);
		early->simParams.emtol = 1e20f;
		Run(*early);
		Require(early->getStep() == 1, "EM early stopping failed");
		// The live editor intentionally uses zero steps and disabled logging.
		auto interactive = MakeSimulation(5, 0, 0.00001f, false, false, 0);
		{
			Engine engine({interactive.get()}, EngineRunMode::Interactive);
			engine.step();
			interactive->simParams.em_variant = true;
			interactive->simParams.emtol = 1e20f;
			engine.step();
			interactive->simParams.em_variant = false;
			engine.step();
			Require(interactive->getStep() == 3 && !engine.IsFinished(), "Live editing stopped at an EM or step limit");
			engine.terminateSimulation();
		}
		{
			auto pressureInteractive = MakeSimulation(5, 0, 0.00001f, false, false, 0);
			pressureInteractive->simParams.pcoupl = PcouplSelect::Berendsen;
			bool rejected = false;
			try { Engine engine({pressureInteractive.get()}, EngineRunMode::Interactive); }
			catch (const std::invalid_argument&) { rejected = true; }
			Require(rejected, "Interactive pressure coupling was not rejected at construction");
		}
		// A zero-step member must not acquire a PME slot or advance another member.
		{
			auto zero = MakeSimulation(5, 0, 0.00001f, true);
			auto member = MakeSimulation(29, 6, 0.00002f, true);
			auto reference = MakeSimulation(29, 6, 0.00002f, true);
			Run(*reference);
			Engine engine({zero.get(), member.get()});
			while (!engine.IsFinished()) engine.step();
			try { Compare(*member, *reference); }
			catch (const std::exception& error) { throw std::runtime_error(std::string("Zero-step member: ") + error.what()); }
			Require(zero->getStep() == 0, "Zero-step member advanced");
		}
		bool rejected = false;
		try { Engine engine({}); } catch (const std::invalid_argument&) { rejected = true; }
		Require(rejected, "Empty batch accepted");
		auto incompatibleA = MakeSimulation(5, 2, 0.00001f, false);
		auto incompatibleB = MakeSimulation(5, 2, 0.00001f, false);
		incompatibleB->simParams.cutoff_nm = 1.f;
		rejected = false;
		try { Engine engine({incompatibleA.get(), incompatibleB.get()}); } catch (const std::invalid_argument&) { rejected = true; }
		Require(rejected, "Incompatible batch accepted");
		std::cout << "Engine batch regression tests passed\n";
	}

	// Water self-diffusion from 280 to 373 K: 32 small water boxes, simulated together as one batch, against experiment.
	// Prints how long the same simulations would take one at a time
	void ShowcaseMultisim() {
		const std::filesystem::path workDir = TestUtils::HeavyTestsDir() / "etc" / "showcase_multisim";
		Environment& environment = Environment::Get();
		constexpr int nTemperatures = 28;

		// One minimized box of SPC/E water, the starting point of every temperature
		SimulationJob minimization;
		minimization.workDir = workDir;
		GroFile emptyBox;
		emptyBox.box_size = Float3{ 6.f };
		minimization.grofile = std::move(emptyBox);
		minimization.topfile.emplace();
		minimization.simParams = SimParams::BasicEMSimParams(500.f);
		minimization.mode = EnvMode::Headless;
		minimization.preprocess = [](GroFile& coordinates, TopologyFile& topology, SimParams&) {
			SimulationBuilder::SolvateGrofile(coordinates, topology);
		};
		const MolecularSystem water = environment.Submit(std::move(minimization)).Get().FinalSystem();

		SimParams production;
		production.nsteps = 12500;	// 25 ps
		production.tcoupl = TcouplSelect::VRescale;
		production.self_diffusion = true;	// Each run writes its self_diffusion.csv
		// The thermostat changes velocities by at most 0.1% per measurement, so it needs frequent measurements to hold
		// each box at its temperature as kinetic energy flows into potential energy early in the run
		production.nsttcouple = 20;
		production.nstenergy = 100;
		const auto MakeJob = [&](float temperature, uint64_t steps, EnvMode mode) {
			SimParams params = production;
			params.ref_t = temperature;
			params.initial_temperature = temperature;
			params.nsteps = steps;
			return Programs::MakeSimulationJob(workDir / Lima::Format("{:.0f}K", temperature), water, params, mode);
		};

		// A short run alone gives the cost of simulating the temperatures one at a time
		SimulationJob reference = MakeJob(298.f, 2500, EnvMode::Headless);
		reference.simParams->self_diffusion = false;
		reference.mustRunAlone = true;
		const double aloneSecondsPerStep = environment.Submit(std::move(reference)).Get().engineTime.count() / 2500.0;

		const auto start = std::chrono::steady_clock::now();
		std::vector<SimulationHandle> handles;
		for (int i = 0; i < nTemperatures; ++i)
			handles.push_back(environment.Submit(MakeJob(280.f + 3.f * i, production.nsteps, EnvMode::Full)));
		std::vector<SimulationResult> results;
		for (auto& handle : handles)
			results.push_back(handle.Get());
		const double wallSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

		const std::filesystem::path csvPath = workDir / "water_self_diffusion.csv";
		std::ofstream csv(csvPath);
		csv << "ref_t_K,temperature_K,D_m2s\n";
		std::map<int, double> batchSeconds;	// [batchId] engine time of its slowest member
		for (const auto& result : results) {
			if (!result.simulation || result.simulation->getStep() != production.nsteps)
				throw std::runtime_error("A showcase simulation did not run fully");
			std::ifstream runCsv(workDir / Lima::Format("{:.0f}K", result.simulation->simParams.ref_t.value) / "self_diffusion.csv");
			std::string header, row;
			if (!std::getline(runCsv, header) || !std::getline(runCsv, row))
				throw std::runtime_error("A showcase simulation wrote no self-diffusion");
			csv << row << '\n';
			double& seconds = batchSeconds[result.execution.batchId];
			seconds = (std::max)(seconds, result.engineTime.count());
		}
		csv.close();

		double batchedSeconds = 0.0;
		for (const auto& [batchId, seconds] : batchSeconds) batchedSeconds += seconds;
		const double aloneSeconds = aloneSecondsPerStep * production.nsteps * nTemperatures;
		const size_t nAtoms = water.coordinates.atoms.size();
		std::cout << Lima::Format("\n{} simulations of {} atoms, {:.0f} ps each, in {} batch(es)\n", nTemperatures, nAtoms,
			production.nsteps * production.dt * NANO_TO_PICO, batchSeconds.size());
		std::cout << Lima::Format("Engine time: {:.1f} s batched, {:.1f} s one at a time ({:.1f}x faster). Wall time {:.1f} s\n",
			batchedSeconds, aloneSeconds, aloneSeconds / batchedSeconds, wallSeconds);
		std::cout << "Self-diffusion written to " << csvPath.string() << '\n';
		const auto script = FileUtils::GetLimaDir() / "dev" / "PyTools" / "WaterSelfDiffusion.py";
		if (std::system(Lima::Format("python \"{}\" \"{}\" --show", script.string(), csvPath.string()).c_str()) != 0)
			throw std::runtime_error("Matplotlib failed to plot the self-diffusion");
	}

	static std::filesystem::path LargeEmDir() { return TestUtils::HeavyTestsDir() / "etc" / "showcase_large_em"; }

	// Builds the system of ShowcaseLargeEM: a satellite tobacco mosaic virus (STMV) above a POPC membrane, solvated in a
	// 30 nm box, written to conf.gro and topol.top. The virus keeps its CHARMM parameters, but gets the lipids' 1-4 scaling
	// (see TODO.md), so the system is for showing speed only
	void PrepareLargeEM() {
		const std::filesystem::path workDir = LargeEmDir();
		const std::filesystem::path stmvDir = TestUtils::HeavyTestsDir() / "benchmarking" / "stmv";
		constexpr float boxSize = 30.f;	// [nm]

		MolecularSystem system;
		system.coordinates.box_size = Float3{ boxSize };
		SimulationBuilder::CreateMembrane(system.coordinates, system.topology,
			Lipids::Composition{ { { "POPC", 100. } }, workDir }.Load(), MembraneGeometry::Plane{ 6.f }, 1);
		// The virus without its water and ions, above the membrane
		const GroFile stmvCoordinates{ stmvDir / "conf.gro" };
		const TopologyFile stmvTopology{ stmvDir / "topol.top" };
		SimulationBuilder::InsertSystem(system.coordinates, system.topology, stmvCoordinates, stmvTopology,
			Float3{ boxSize / 2.f, boxSize / 2.f, 18.5f }, [](const TopologyFile::Moleculetype& moleculetype) {
				return std::ranges::any_of(moleculetype.atoms, [](const TopologyFile::AtomsEntry& atom) {
					return atom.residue != "TIP3" && atom.residue != "MG" && atom.residue != "CLA"; });
			});
		SimulationBuilder::SolvateGrofile(system.coordinates, system.topology);
		system.coordinates.printToFile(workDir / "conf.gro");
		system.topology.printToFile(workDir / "topol.top");
		std::cout << Lima::Format("Prepared {} atoms in {}\n", system.coordinates.atoms.size(), workDir.string());
	}

	// Energy minimization with PME of the few million atoms prepared by PrepareLargeEM, shown live
	void ShowcaseLargeEM() {
		const std::filesystem::path workDir = LargeEmDir();
		if (!std::filesystem::exists(workDir / "topol.top"))
			throw std::runtime_error("Run agenttesting --prepare-large-em first");
		SimulationJob job;
		job.workDir = workDir;
		job.grofile.emplace(workDir / "conf.gro");
		job.topfile.emplace(workDir / "topol.top");
		SimParams params = SimParams::BasicEMSimParams(1000.f);
		params.enable_electrostatics = true;
		job.simParams = params;
		job.mode = EnvMode::Full;
		job.mustRunAlone = true;

		const auto start = std::chrono::steady_clock::now();
		auto result = Environment::Get().Submit(std::move(job)).Get();
		const double totalSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
		if (!result.simulation)
			throw std::runtime_error("The large EM showcase produced no simulation");

		const auto& maxForces = result.simulation->maxForceBuffer;
		const double engineSeconds = result.engineTime.count();
		std::cout << Lima::Format("\n{} atoms minimized in {} steps: max force {:.0f} -> {:.0f} kJ/mol/nm\n",
			result.simulation->boxImage->grofile.atoms.size(), result.simulation->getStep(),
			maxForces.empty() ? 0.f : maxForces.front().second, maxForces.empty() ? 0.f : maxForces.back().second);
		std::cout << Lima::Format("Engine {:.1f} s ({:.1f} ms/step), total {:.1f} s including setup\n",
			engineSeconds, 1000.0 * engineSeconds / (std::max)(int64_t{ 1 }, result.simulation->getStep()), totalSeconds);
	}
}
