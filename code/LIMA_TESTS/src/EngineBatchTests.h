#pragma once

#include "Engine.cuh"
#include <stdexcept>

namespace EngineBatchTests {
	inline void Require(bool value, const char* message) {
		if (!value) throw std::runtime_error(message);
	}

	inline std::unique_ptr<Simulation> MakeSimulation(int particles, int steps, float dt, bool electrostatics, bool em = false, int loggingInterval = 2) {
		SimParams params;
		params.n_steps = steps;
		params.dt = dt;
		params.enable_electrostatics = electrostatics;
		params.em_variant = em;
		params.em_force_tolerance = 0.f;
		params.stepsPerNlistupdate = 3;
		params.data_logging_interval = loggingInterval;
		params.steps_per_temperature_measurement = 10;
		params.apply_thermostat = !em;
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
		box.particlesBondedToParticle[0] = ParticlesBondedToParticle::Create({4});
		box.particlesBondedToParticle[4] = ParticlesBondedToParticle::Create({0});
		box.pclustersBondedToPcluster[0] = PclustersBondedToPcluster::Create({1});
		box.pclustersBondedToPcluster[1] = PclustersBondedToPcluster::Create({0});
		sim->PrepareDataBuffers();
		return sim;
	}

	inline void Run(Simulation& sim) {
		Engine engine({&sim});
		while (!engine.IsFinished()) engine.step();
		engine.terminateSimulation();
	}

	inline void Compare(Simulation& a, Simulation& b) {
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
				if (actual[i] != expected[i]) throw std::runtime_error(std::format("Batch changed logged energy at {}: {} vs {}", i, actual[i], expected[i]));
			throw std::runtime_error("Batch changed logged energy buffer length");
		}
		Require(a.forceBuffer->GetBuffer() == b.forceBuffer->GetBuffer(), "Batch changed logged forces");
		Require(a.temperature_buffer == b.temperature_buffer, "Batch mixed temperature reductions");
	}

	inline void RunAll() {
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
								throw std::runtime_error(std::format("First survivor state diverged at step {}, pc {}, lane {}", step, pc, lane));
					const auto state = GenericCopyToHost(engine.OffloadPclusterState(0).Get(), first->box->persistentClusters.size());
					const auto expected = GenericCopyToHost(stepwiseEngine.OffloadPclusterState(0).Get(), first->box->persistentClusters.size());
					for (size_t pc = 0; pc < state.size(); ++pc)
						for (int lane = 0; lane < 4; ++lane)
							if (state[pc].pqd[lane].Valid() && state[pc].pqd[lane].position != expected[pc].pqd[lane].position)
								throw std::runtime_error(std::format("First survivor diverged at step {}, pc {}, lane {}", step, pc, lane));
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
						catch (const std::exception& error) { throw std::runtime_error(std::format("Batch electrostatics {} EM {} step {}: {}", electrostatics, em, b->getStep(), error.what())); }
					}
				}
				engine.terminateSimulation();
				engine.terminateSimulation(); // Finalization must be idempotent.
				engine.step(); // A completed batch must not advance.
				try { Compare(*a, *referenceA); Compare(*b, *referenceB); }
				catch (const std::exception& error) { throw std::runtime_error(std::format("Batch electrostatics {} EM {} final: {}", electrostatics, em, error.what())); }
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
		early->simParams.em_force_tolerance = 1e20f;
		Run(*early);
		Require(early->getStep() == 1, "EM early stopping failed");
		// The live editor intentionally uses zero steps and disabled logging.
		auto interactive = MakeSimulation(5, 0, 0.00001f, false, false, 0);
		{
			Engine engine({interactive.get()}, EngineRunMode::Interactive);
			engine.step();
			interactive->simParams.em_variant = true;
			interactive->simParams.em_force_tolerance = 1e20f;
			engine.step();
			interactive->simParams.em_variant = false;
			engine.step();
			Require(interactive->getStep() == 3 && !engine.IsFinished(), "Live editing stopped at an EM or step limit");
			engine.terminateSimulation();
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
}
