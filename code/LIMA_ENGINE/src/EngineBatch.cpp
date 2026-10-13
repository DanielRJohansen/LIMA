#include "BatchData.cuh"
#include "BatchCompatibility.h"
#include "EngineBodies.cuh"
#include "PhysicsUtils.cuh"
#include <algorithm>
#include <array>
#include <exception>
#include "ParallelFor.h"
#include <cmath>
#include "Format.h"
#include <type_traits>
#include <unordered_map>
#include <cstdio>
#include <cstdlib>

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
	// For each particle, in (pcluster, lane) order, its pairs with earlier particles are reported, so b is always the later one.
	template <typename OnPair>
	void ForEachClosePair(const Box& box, BoundaryConditionSelect bc, float maxDistance, bool includeBonded, OnPair&& onPair) {
		const Int3 size = box.boxparams.boxSize;
		struct Entry { Float3 pos; ParticleRef ref; int id; };
		std::vector<Entry> entries;
		entries.reserve(box.persistentClusters.size() * PersistentCluster::maxParticles);
		for (int pc = 0; pc < static_cast<int>(box.persistentClusters.size()); ++pc) {
			const auto& meta = box.persistentClustersMetadata[pc];
			for (int lane = 0; lane < meta.nParticles; ++lane) {
				Float3 pos = box.persistentClusters[pc].pqd[lane].position;
				if (bc == PBC) {
					pos.x -= std::floor(pos.x / size.x) * size.x;
					pos.y -= std::floor(pos.y / size.y) * size.y;
					pos.z -= std::floor(pos.z / size.z) * size.z;
				}
				entries.push_back({ pos, { pc, lane }, meta.particleIdsGlobal[lane] });
			}
		}
		if (entries.empty()) return;

		// A flat grid of about one particle per cell. Cells are never smaller than maxDistance, so every close pair is
		// in the same or adjacent cells. Cells this small as a hash map were the cost of starting a large batch
		Float3 lo{}, hi{ static_cast<float>(size.x), static_cast<float>(size.y), static_cast<float>(size.z) };
		if (bc != PBC) {
			lo = hi = entries.front().pos;
			for (const Entry& entry : entries) {
				lo = Float3{ std::min(lo.x, entry.pos.x), std::min(lo.y, entry.pos.y), std::min(lo.z, entry.pos.z) };
				hi = Float3{ std::max(hi.x, entry.pos.x), std::max(hi.y, entry.pos.y), std::max(hi.z, entry.pos.z) };
			}
		}
		const Float3 extent{ std::max(hi.x - lo.x, maxDistance), std::max(hi.y - lo.y, maxDistance), std::max(hi.z - lo.z, maxDistance) };
		const float cell = std::max(maxDistance, std::cbrt(extent.x * extent.y * extent.z / static_cast<float>(entries.size())));
		const Int3 nCells{ std::max(1, static_cast<int>(extent.x / cell)), std::max(1, static_cast<int>(extent.y / cell)), std::max(1, static_cast<int>(extent.z / cell)) };
		const auto CellIndex = [&](float value, float origin, float length, int n) {
			return std::clamp(static_cast<int>(std::floor((value - origin) / length * static_cast<float>(n))), 0, n - 1);
		};
		std::vector<Int3> cellOf(entries.size());
		std::vector<int> cellStart(static_cast<size_t>(nCells.x) * nCells.y * nCells.z + 1, 0);
		const auto Flat = [&](Int3 c) { return (static_cast<size_t>(c.x) * nCells.y + c.y) * nCells.z + c.z; };
		for (size_t i = 0; i < entries.size(); ++i) {
			cellOf[i] = Int3{ CellIndex(entries[i].pos.x, lo.x, extent.x, nCells.x), CellIndex(entries[i].pos.y, lo.y, extent.y, nCells.y),
				CellIndex(entries[i].pos.z, lo.z, extent.z, nCells.z) };
			++cellStart[Flat(cellOf[i]) + 1];
		}
		for (size_t c = 1; c < cellStart.size(); ++c) cellStart[c] += cellStart[c - 1];
		std::vector<int> particlesInCells(entries.size());	// In particle order within each cell
		{
			std::vector<int> fill(cellStart.begin(), cellStart.end() - 1);
			for (size_t i = 0; i < entries.size(); ++i) particlesInCells[fill[Flat(cellOf[i])]++] = static_cast<int>(i);
		}

		for (size_t i = 0; i < entries.size(); ++i) {
			const Entry& entry = entries[i];
			std::array<size_t, 27> neighborCells;
			int nNeighborCells = 0;
			for (int dx = -1; dx <= 1; ++dx) for (int dy = -1; dy <= 1; ++dy) for (int dz = -1; dz <= 1; ++dz) {
				Int3 c{ cellOf[i].x + dx, cellOf[i].y + dy, cellOf[i].z + dz };
				if (bc == PBC) c = Int3{ (c.x + nCells.x) % nCells.x, (c.y + nCells.y) % nCells.y, (c.z + nCells.z) % nCells.z };
				else if (c.x < 0 || c.y < 0 || c.z < 0 || c.x >= nCells.x || c.y >= nCells.y || c.z >= nCells.z) continue;
				const size_t flat = Flat(c);
				// A box only a few cells wide reaches the same cell from both sides
				if (std::find(neighborCells.begin(), neighborCells.begin() + nNeighborCells, flat) == neighborCells.begin() + nNeighborCells)
					neighborCells[nNeighborCells++] = flat;
			}
			for (int n = 0; n < nNeighborCells; ++n) {
				for (int k = cellStart[neighborCells[n]]; k < cellStart[neighborCells[n] + 1]; ++k) {
					const int j = particlesInCells[k];
					if (j >= static_cast<int>(i)) break;	// Later particles report this pair themselves
					const Entry& other = entries[j];
					Float3 diff = entry.pos - other.pos;
					if (bc == PBC) {
						diff.x -= std::round(diff.x / size.x) * size.x;
						diff.y -= std::round(diff.y / size.y) * size.y;
						diff.z -= std::round(diff.z / size.z) * size.z;
					}
					if (diff.lenSquared() >= maxDistance * maxDistance) continue;
					if (!includeBonded && entry.id >= 0 && entry.id < static_cast<int>(box.particlesBondedToParticle.size()) && box.particlesBondedToParticle[entry.id].Contains(other.id)) continue;
					onPair(other.ref, entry.ref, diff.len());
				}
			}
		}
	}

	// Calls f(i) for every simulation in parallel. If several fail, the lowest index's error is rethrown, so the
	// message does not depend on thread timing
	template <typename F>
	void ForEachSimulationInParallel(size_t n, F&& f) {
		std::vector<std::exception_ptr> errors(n);
		ParallelUtils::ParallelFor(n, [&](size_t i) {
			try { f(i); }
			catch (...) { errors[i] = std::current_exception(); }
		});
		for (const auto& error : errors)
			if (error) std::rethrow_exception(error);
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
			if (auto mismatch = FindIncompatibility(*simulations.front(), *sim))
				throw std::invalid_argument("Incompatible batch parameter: " + std::string(*mismatch));
		}
		// The per-particle checks rely on the metadata checked above, and are by far the most expensive
		ForEachSimulationInParallel(simulations.size(), [&](size_t i) {
			const auto& box = *simulations[i]->box;
			ValidateParticles(box);
			if (!simulations[i]->simParams.em_variant) ValidateNoOverlap(box, simulations[i]->simParams.bc_select);
		});
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
		std::vector<PclusterConstraints> constraints;	// Stays empty unless a simulation has constraints
		BondGroups bonds;
		if (!active)
			ForEachSimulationInParallel(simulations.size(), [&](size_t i) {
				if (simulations[i]->simParams.em_variant) SeparateCoincidentParticles(*simulations[i]->box, simulations[i]->simParams.bc_select);
			});
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
			// The box's constraints and degrees of freedom are set when it is built, so MD must use the constraints it was built with
			if (!sim.simulation->simParams.em_variant && sim.simulation->simParams.constraints.value != box.builtConstraints)
				throw std::invalid_argument("The constraints setting differs from the one the box was built with. Use the same, eg. also "
					"in the energy minimization it continues from");
			// Energy minimization ignores the constraints, and uses the bonds they constrain
			if (!box.pclusterConstraints.empty() && !sim.simulation->simParams.em_variant) {
				constraints.resize(pclusters.size());
				Append(constraints, box.pclusterConstraints);
			}
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
				group.indexOfFirstCmapbond += CheckedCount(bonds.cmapbonds.size());
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
			{
				// Boxes usually share a forcefield, so their CMAP grids are deduplicated
				std::vector<int> gridIndices(box.bondgroups.cmapGrids.size());
				for (size_t g = 0; g < gridIndices.size(); g++) {
					const auto found = std::find(bonds.cmapGrids.begin(), bonds.cmapGrids.end(), box.bondgroups.cmapGrids[g]);
					gridIndices[g] = CheckedCount(found - bonds.cmapGrids.begin());
					if (found == bonds.cmapGrids.end())
						bonds.cmapGrids.push_back(box.bondgroups.cmapGrids[g]);
				}
				for (auto cmap : box.bondgroups.cmapbonds) {
					cmap.params.gridIndex = gridIndices.at(cmap.params.gridIndex);
					bonds.cmapbonds.push_back(cmap);
				}
			}
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
		if (!constraints.empty())
			constraints.resize(pclusters.size());
		batch.pclusterConstraintsDevice.SetData(constraints);
		batch.hasConstraints = !constraints.empty();
		batch.pClusterMetaDevice.SetData(metadata);
		batch.boxState.pclusterInterimStates = GenericCopyToDevice(states);
		batch.bondgroupDescriptors.SetData(bonds.groups);
		batch.bondgroupParticles.SetData(bonds.particles);
		batch.nBondgroupParticles = CheckedCount(bonds.particles.size());
		// One group owns the primary result in each particle slot. Only other appearances need gathered output.
		std::vector<int> extraResultIndices(bonds.particles.size(), -1);
		size_t nExtraResults = 0;
		for (const auto& meta : metadata) {
			for (const auto& refs : meta.bondgroupReferences) {
				for (int i = 1; i < refs.nBondgroupApperances; i++)
					extraResultIndices[refs.bondgroupApperances[i].indexInForceEnergiesBondgroups] = CheckedCount(nExtraResults++);
			}
		}
		batch.bondgroupExtraResultIndices.SetData(extraResultIndices);
		batch.extraBondForceResults.Expand(nExtraResults);
		batch.bondgroupSinglebonds.SetData(bonds.singlebonds);
		batch.bondgroupPairbonds.SetData(bonds.pairbonds);
		batch.bondgroupAnglebonds.SetData(bonds.anglebonds);
		batch.bondgroupDihedralbonds.SetData(bonds.dihedralbonds);
		batch.bondgroupImproperdihedralbonds.SetData(bonds.improperdihedralbonds);
		batch.bondgroupCmapbonds.SetData(bonds.cmapbonds);
		{
			batch.cmapGridSpacing = bonds.cmapGrids.empty() ? 0 : bonds.cmapGrids.front().gridSpacing;
			std::vector<float4> points;
			for (const CmapGrid& grid : bonds.cmapGrids) {
				if (grid.gridSpacing != batch.cmapGridSpacing)
					throw std::runtime_error("CMAP grids of different sizes in one batch are not supported");
				Append(points, grid.points);
			}
			batch.cmapGridPoints.SetData(points);
		}
		if (std::getenv("LIMA_BOND_LAYOUT_STATS")) {
			const auto PrintTable = [](const char* name, const auto& host, const auto& device) {
				std::printf("Bond table %s: %zu bonds, %zu parameters, %d-byte indices, %zu -> %zu bytes\n",
					name, host.size(), device.ParameterCount(), device.Get().parameterIdBytes,
					host.size() * sizeof(typename std::decay_t<decltype(host)>::value_type), device.Bytes());
			};
			PrintTable("single", bonds.singlebonds, batch.bondgroupSinglebonds);
			PrintTable("pair", bonds.pairbonds, batch.bondgroupPairbonds);
			PrintTable("angle", bonds.anglebonds, batch.bondgroupAnglebonds);
			PrintTable("dihedral", bonds.dihedralbonds, batch.bondgroupDihedralbonds);
			PrintTable("improper", bonds.improperdihedralbonds, batch.bondgroupImproperdihedralbonds);
			std::printf("Bond table cmap: %zu bonds, %zu grids of %d^2 points\n", bonds.cmapbonds.size(), bonds.cmapGrids.size(), batch.cmapGridSpacing);
		}
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
