#include "BatchData.cuh"
#include "BatchCompatibility.h"
#include "EngineBodies.cuh"
#include "PhysicsUtils.cuh"
#include <algorithm>
#include <cmath>
#include "Format.h"
#include <type_traits>
#include <unordered_map>

namespace EngineBatch {
	template<typename T>
	void Append(std::vector<T>& dst, const std::vector<T>& src) {
		dst.insert(dst.end(), src.begin(), src.end());
	}

	int CheckedCount(size_t count) {
		if (count > INT_MAX) throw std::invalid_argument("Engine batch exceeds 32-bit work indices");
		return static_cast<int>(count);
	}

	bool FinitePositive(float value) { return value > 0.f && std::isfinite(value); }

	void ValidateParticles(const Box& box) {
		for (size_t pc = 0; pc < box.persistentClusters.size(); ++pc) {
			const auto& meta = box.persistentClustersMetadata[pc];
			for (int lane = 0; lane < meta.nParticles; ++lane) {
				const Float3 pos = box.persistentClusters[pc].pqd[lane].position;
				if (!std::isfinite(pos.x) || !std::isfinite(pos.y) || !std::isfinite(pos.z))
					throw std::invalid_argument(Lima::Format("Particle {} position must be finite", meta.particleIdsGlobal[lane]));
				if (!FinitePositive(meta.mass[lane]))
					throw std::invalid_argument(Lima::Format("Particle {} mass must be finite and positive", meta.particleIdsGlobal[lane]));
			}
		}
	}

	struct ParticleRef { int pc; int lane; };

	// Calls onPair(a, b, distance) for every pair of particles closer than maxDistance.
	// Pairs listed in particlesBondedToParticle are skipped unless includeBonded is set.
	template <typename OnPair>
	void ForEachClosePair(const Box& box, BoundaryConditionSelect bc, float maxDistance, bool includeBonded, OnPair&& onPair) {
		const Int3 size = box.boxparams.boxSize;
		const Int3 nCells{ std::max(1, static_cast<int>(size.x / maxDistance)), std::max(1, static_cast<int>(size.y / maxDistance)), std::max(1, static_cast<int>(size.z / maxDistance)) };
		const Float3 cellSize{ size.x / static_cast<float>(nCells.x), size.y / static_cast<float>(nCells.y), size.z / static_cast<float>(nCells.z) };
		const auto Wrap = [](int value, int n) { return ((value % n) + n) % n; };
		// Aliasing between keys (outside the box with NoBC) only adds candidate pairs; the distance test stays exact
		const auto Key = [&](int x, int y, int z) { return (static_cast<int64_t>(x) * nCells.y + y) * nCells.z + z; };
		struct Entry { Float3 pos; ParticleRef ref; };
		std::unordered_map<int64_t, std::vector<Entry>> cells;
		for (int pc = 0; pc < static_cast<int>(box.persistentClusters.size()); ++pc) {
			const auto& meta = box.persistentClustersMetadata[pc];
			for (int lane = 0; lane < meta.nParticles; ++lane) {
				Float3 pos = box.persistentClusters[pc].pqd[lane].position;
				if (bc == PBC) {
					pos.x -= std::floor(pos.x / size.x) * size.x;
					pos.y -= std::floor(pos.y / size.y) * size.y;
					pos.z -= std::floor(pos.z / size.z) * size.z;
				}
				const int id = meta.particleIdsGlobal[lane];
				const int cx = static_cast<int>(std::floor(pos.x / cellSize.x)), cy = static_cast<int>(std::floor(pos.y / cellSize.y)), cz = static_cast<int>(std::floor(pos.z / cellSize.z));
				for (int dx = -1; dx <= 1; ++dx) for (int dy = -1; dy <= 1; ++dy) for (int dz = -1; dz <= 1; ++dz) {
					const int nx = bc == PBC ? Wrap(cx + dx, nCells.x) : cx + dx;
					const int ny = bc == PBC ? Wrap(cy + dy, nCells.y) : cy + dy;
					const int nz = bc == PBC ? Wrap(cz + dz, nCells.z) : cz + dz;
					const auto found = cells.find(Key(nx, ny, nz));
					if (found == cells.end()) continue;
					for (const Entry& other : found->second) {
						Float3 diff = pos - other.pos;
						if (bc == PBC) {
							diff.x -= std::round(diff.x / size.x) * size.x;
							diff.y -= std::round(diff.y / size.y) * size.y;
							diff.z -= std::round(diff.z / size.z) * size.z;
						}
						if (diff.lenSquared() >= maxDistance * maxDistance) continue;
						const int otherId = box.persistentClustersMetadata[other.ref.pc].particleIdsGlobal[other.ref.lane];
						if (!includeBonded && id >= 0 && id < static_cast<int>(box.particlesBondedToParticle.size()) && box.particlesBondedToParticle[id].Contains(otherId)) continue;
						onPair(other.ref, ParticleRef{ pc, lane }, diff.len());
					}
				}
				cells[Key(cx, cy, cz)].push_back({ pos, { pc, lane } });
			}
		}
	}

	// MD has no way to resolve coincident nonbonded particles: the pair force is singular. EM is
	// expected to push such particles apart, so it skips this check.
	void ValidateNoOverlap(const Box& box, BoundaryConditionSelect bc) {
		constexpr float minDistance = 0.005f; // [nm]
		ForEachClosePair(box, bc, minDistance, false, [&](ParticleRef a, ParticleRef b, float distance) {
			throw std::invalid_argument(Lima::Format("Particle overlap: particles {} and {} are {:.4f} nm apart (minimum {} nm); run energy minimization first",
				box.persistentClustersMetadata[a.pc].particleIdsGlobal[a.lane], box.persistentClustersMetadata[b.pc].particleIdsGlobal[b.lane], distance, minDistance));
		});
	}

	// Exactly coincident particles have no defined force direction, so EM would see zero force and
	// report convergence. Nudge them apart deterministically so the minimizer has a gradient to follow.
	void SeparateCoincidentParticles(Box& box, BoundaryConditionSelect bc) {
		constexpr float coincident = 1e-4f;	// [nm]
		constexpr float nudge = 0.01f;		// [nm]
		for (int pass = 0; pass < 8; ++pass) {
			std::vector<ParticleRef> moves;
			ForEachClosePair(box, bc, coincident, true, [&](ParticleRef, ParticleRef b, float) { moves.push_back(b); });
			if (moves.empty()) return;
			for (const ParticleRef ref : moves) {
				// Golden-angle direction per particle id, so particles in a coincident cluster all move differently
				const int id = box.persistentClustersMetadata[ref.pc].particleIdsGlobal[ref.lane] + pass * 7919;
				const float z = 1.f - 2.f * std::fmod(id * 0.618034f + 0.5f, 1.f);
				const float r = std::sqrt(std::max(0.f, 1.f - z * z));
				const float phi = id * 2.399963f;
				auto& position = box.persistentClusters[ref.pc].pqd[ref.lane].position;
				position = position + Float3{ r * std::cos(phi), r * std::sin(phi), z } * nudge;
			}
		}
		throw std::invalid_argument("Could not separate coincident particles for energy minimization");
	}

	void Validate(const std::vector<Simulation*>& simulations) {
		if (simulations.empty()) throw std::invalid_argument("Engine requires a nonempty batch");
		for (size_t i = 0; i < simulations.size(); ++i) {
			const auto* sim = simulations[i];
			if (!sim || !sim->box) throw std::invalid_argument("Null simulation in batch");
			if (std::find(simulations.begin(), simulations.begin() + i, sim) != simulations.begin() + i)
				throw std::invalid_argument("Duplicate simulation in batch");
			if (sim->finished || sim->getStep() != 0)
				throw std::invalid_argument("Engine requires fresh simulations at step zero");
			const auto& params = sim->simParams;
			if (params.stepsPerNlistupdate <= 0 || params.steps_per_temperature_measurement <= 0 || params.data_logging_interval < 0)
				throw std::invalid_argument("Invalid engine measurement or update interval");
			if (!FinitePositive(params.dt))
				throw std::invalid_argument("Timestep must be finite and positive");
			if (!FinitePositive(params.cutoff_nm))
				throw std::invalid_argument("Cutoff must be finite and positive");
			if (params.n_steps.value > INT64_MAX)
				throw std::invalid_argument("Simulation step count exceeds engine limit");
			const auto& box = *sim->box;
			const auto dim = box.boxparams.boxSize;
			if (dim.x <= 0 || dim.y <= 0 || dim.z <= 0 || dim.x >= 1024 || dim.y >= 1024 || dim.z >= 1024)
				throw std::invalid_argument("Unsupported engine box dimensions");
			if (box.persistentClusters.size() != box.persistentClustersMetadata.size() || box.persistentClusters.size() != box.pclusterInterimStates.size()
				|| box.particlesBondedToParticle.size() != box.boxparams.totalParticles || box.pclustersBondedToPcluster.size() != box.persistentClusters.size())
				throw std::invalid_argument("Incomplete simulation particle metadata");
			if (params.data_logging_interval > 0 && (!sim->traj_buffer || !sim->potE_buffer || !sim->vel_buffer || !sim->forceBuffer))
				throw std::invalid_argument("Simulation logging buffers have not been prepared");
			if (sim->box->persistentClusters.empty()) throw std::invalid_argument("Cannot simulate an empty box");
			ValidateParticles(box);
			if (!params.em_variant) ValidateNoOverlap(box, params.bc_select);
			if (auto mismatch = FindIncompatibility(*simulations.front(), *sim))
				throw std::invalid_argument("Incompatible batch parameter: " + std::string(*mismatch));
		}
	}

	void Pack(EngineBatchData& batch, const std::vector<Simulation*>& simulations, const std::vector<bool>* active) {
		if (!active) Validate(simulations);
		else if (active->size() != simulations.size()) throw std::invalid_argument("Invalid active batch layout");
		batch.params = simulations.front()->simParams;
		batch.boxSize = simulations.front()->box->boxparams.boxSize;
		batch.ewaldKappa = PhysicsUtils::CalcEwaldkappa(batch.params.cutoff_nm);
		std::vector<PersistentCluster> pclusters;
		std::vector<PersistentClusterMeta> metadata;
		std::vector<PersistentclusterInterimState> states;
		BondGroups bonds;
		for (size_t simulationId = 0; simulationId < simulations.size(); ++simulationId) {
			auto* simulation = simulations[simulationId];
			if (active && !(*active)[simulationId]) {
				EngineSimulationData retired;
				retired.simulation = simulation;
				retired.device.active = false;
				retired.finalized = true;
				batch.simulations.push_back(retired);
				continue;
			}
			if (!active && simulation->simParams.em_variant) SeparateCoincidentParticles(*simulation->box, simulation->simParams.bc_select);
			const auto& box = *simulation->box;
			EngineSimulationData sim;
			sim.simulation = simulation;
			auto& layout = sim.device;
			layout.particles = { batch.nParticles, box.boxparams.totalParticles };
			layout.pclusters = { CheckedCount(pclusters.size()), CheckedCount(box.persistentClusters.size()) };
			layout.bondgroups = { CheckedCount(bonds.groups.size()), CheckedCount(box.bondgroups.size()) };
			layout.gridnodes = { batch.nGridnodes, batch.boxSize.InnerProduct() };
			layout.logOffset = pclusters.size() * PersistentCluster::maxParticles * DatabuffersDeviceController::nStepsInBuffer;
			layout.dt = simulation->simParams.dt;
			layout.uniformElectricField = box.uniformElectricField;
			layout.active = simulation->simParams.n_steps != 0;
			sim.runstatus.simulation_finished = !layout.active;
			batch.nParticles = CheckedCount(size_t(batch.nParticles) + layout.particles.count);
			batch.nGridnodes = CheckedCount(size_t(batch.nGridnodes) + layout.gridnodes.count);
			CheckedCount(size_t(batch.nGridnodes) * PClusterTransfermodule::maxClustersPerBlock);
			Append(pclusters, box.persistentClusters);
			Append(states, box.pclusterInterimStates);
			for (auto meta : box.persistentClustersMetadata) {
				for (int lane = 0; lane < PersistentCluster::maxParticles; ++lane) {
					if (meta.particleIdsGlobal[lane] >= 0) meta.particleIdsGlobal[lane] += layout.particles.offset;
					auto& refs = meta.bondgroupReferences[lane];
					for (int i = 0; i < refs.nBondgroupApperances; ++i) {
						refs.bondgroupApperances[i].indexInForceEnergiesBondgroups += CheckedCount(bonds.particles.size());
						refs.bondgroupApperances[i].bondgroupId += layout.bondgroups.offset;
					}
				}
				metadata.push_back(meta);
			}
			for (auto group : box.bondgroups.groups) {
				group.indexOfFirstParticle += CheckedCount(bonds.particles.size());
				group.indexOfFirstSinglebond += CheckedCount(bonds.singlebonds.size());
				group.indexOfFirstPairbond += CheckedCount(bonds.pairbonds.size());
				group.indexOfFirstAnglebond += CheckedCount(bonds.anglebonds.size());
				group.indexOfFirstDihedralbond += CheckedCount(bonds.dihedralbonds.size());
				group.indexOfFirstImproperdihedralbond += CheckedCount(bonds.improperdihedralbonds.size());
				bonds.groups.push_back(group);
			}
			for (auto ref : box.bondgroups.particles) {
				ref.pcid += layout.pclusters.offset;
				bonds.particles.push_back(ref);
			}
			Append(bonds.singlebonds, box.bondgroups.singlebonds);
			Append(bonds.pairbonds, box.bondgroups.pairbonds);
			Append(bonds.anglebonds, box.bondgroups.anglebonds);
			Append(bonds.dihedralbonds, box.bondgroups.dihedralbonds);
			Append(bonds.improperdihedralbonds, box.bondgroups.improperdihedralbonds);
			for (auto neighbors : box.particlesBondedToParticle) {
				neighbors.AddOffset(layout.particles.offset);
				batch.particlesBondedToParticle.push_back(neighbors);
			}
			for (auto neighbors : box.pclustersBondedToPcluster) {
				neighbors.AddOffset(layout.pclusters.offset);
				batch.pclustersBondedToPcluster.push_back(neighbors);
			}
			batch.simulations.push_back(sim);
		}
		batch.nPclusters = CheckedCount(pclusters.size());
		batch.nBondgroups = CheckedCount(bonds.groups.size());
		CheckedCount(pclusters.size() * PersistentCluster::maxParticles);
		batch.pClusterDevice.SetData(pclusters);
		batch.pClusterMetaDevice.SetData(metadata);
		batch.boxState.pclusterInterimStates = GenericCopyToDevice(states);
		batch.bondgroupDescriptors.SetData(bonds.groups);
		batch.bondgroupParticles.SetData(bonds.particles);
		batch.bondgroupSinglebonds.SetData(bonds.singlebonds);
		batch.bondgroupPairbonds.SetData(bonds.pairbonds);
		batch.bondgroupAnglebonds.SetData(bonds.anglebonds);
		batch.bondgroupDihedralbonds.SetData(bonds.dihedralbonds);
		batch.bondgroupImproperdihedralbonds.SetData(bonds.improperdihedralbonds);
		batch.forceEnergyInterims = std::make_unique<ForceEnergyInterims>(CheckedCount(bonds.particles.size()), batch.nParticles, batch.nPclusters);
		batch.forcesMagnitudeSquareDevice.Expand(batch.nParticles);
		cudaMemset(batch.forcesMagnitudeSquareDevice.Get(), 0, sizeof(float) * batch.nParticles);
		{ // Allocated for all batches, since interactive engines can switch to EM at any step. The preconditioner is uploaded by the Engine
			const size_t nSlots = size_t(batch.nPclusters) * PersistentCluster::maxParticles;
			batch.emParticles.Expand(nSlots);
			cudaMemset(batch.emParticles.Get(), 0, sizeof(EM::ParticleState) * nSlots);
			batch.emForces.Expand(nSlots);
			batch.emPreconditionedForce.Expand(nSlots);
			batch.emBlocksDone.Expand(batch.simulations.size());
			cudaMemset(batch.emBlocksDone.Get(), 0, sizeof(unsigned int) * batch.simulations.size());
			batch.emStates.Expand(batch.simulations.size());
			cudaMemset(batch.emStates.Get(), 0, sizeof(EM::SimState) * batch.simulations.size());
			batch.emStatesHost.resize(batch.simulations.size());
		}
	}
}
