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
			if (params.apply_thermostat && !FinitePositive(params.tau_t))
				throw std::invalid_argument("tau_t must be finite and positive");
			if (!FinitePositive(params.cutoff_nm))
				throw std::invalid_argument("Cutoff must be finite and positive");
			if (!(FinitePositive(params.verlet_buffer_tolerance) || params.verlet_buffer_tolerance == -1.f) || !std::isfinite(params.rlist.value))
				throw std::invalid_argument("verlet_buffer_tolerance must be positive, or -1 to use rlist");
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

	// The pair list buffer of an MD simulation, as GROMACS' Verlet buffer (calc_verletbuf.cpp, Pall & Hess 2013): the smallest buffer
	// for which the pairs that move from beyond the list radius into the cutoff before the next list update cause at most
	// verlet_buffer_tolerance energy drift. Each atom is displaced as a free particle at the reference temperature, and each atom
	// type pair's potential is expanded to third order around the cutoff. LIMA's potentials are not shifted, so unlike GROMACS
	// the potential at the cutoff is included too. A constrained atom moves less, as it rotates around its heaviest constraint
	// partner, and the lists hold 4x4 blocks with any pair in range, so also pairs beyond the list radius, which GROMACS credits
	// with surfaceFraction [nm]
	float VerletBuffer(const Simulation& simulation) {
		const SimParams& params = simulation.simParams;
		const float cutoff = params.cutoff_nm;
		if (params.verlet_buffer_tolerance < 0.f)
			return params.rlist == 0.f ? 0.05f : std::max(params.rlist.value, cutoff) - cutoff;

		const Box& box = *simulation.box;
		// partnerMass and constraintLength are of the heaviest constraint partner, 0 when unconstrained
		struct Atomtype { float mass, partnerMass, constraintLength, sigmaHalf, epsilonSqrt, charge; int64_t count; };
		std::vector<Atomtype> types;
		std::unordered_map<uint64_t, std::vector<int>> typesByKey;	// Hash of the parameters to candidate types
		double kineticEnergy = 0.0;	// [J/mol]
		int64_t nAtoms = 0;
		for (size_t pc = 0; pc < box.persistentClusters.size(); pc++) {
			const auto& meta = box.persistentClustersMetadata[pc];
			for (int p = 0; p < meta.nParticles; p++) {
				const NBParams& nb = box.persistentClusters[pc].pqd[p].params;
				Atomtype type{ meta.mass[p], 0.f, 0.f, nb.sigmaHalf, nb.epsilonSqrt, nb.charge, 1 };
				// MD constrains the box's constraints whatever the batch is, see Pack
				if (!box.pclusterConstraints.empty()) {
					const PclusterConstraints& constraints = box.pclusterConstraints[pc];
					for (int c = 0; c < constraints.nConstraints; c++) {
						const int a = constraints.particles[c] & 3, b = constraints.particles[c] >> 2;
						const int partner = a == p ? b : b == p ? a : -1;
						if (partner >= 0 && meta.mass[partner] > type.partnerMass) {
							type.partnerMass = meta.mass[partner];
							type.constraintLength = constraints.distance[c];
						}
					}
				}
				const uint64_t key = std::hash<float>{}(type.mass) ^ std::hash<float>{}(type.partnerMass) * 3 ^ std::hash<float>{}(type.sigmaHalf) * 7
					^ std::hash<float>{}(type.epsilonSqrt) * 13 ^ std::hash<float>{}(type.charge) * 17;
				auto& candidates = typesByKey[key];
				const auto same = std::ranges::find_if(candidates, [&](int t) { return types[t].mass == type.mass && types[t].partnerMass == type.partnerMass
					&& types[t].constraintLength == type.constraintLength && types[t].sigmaHalf == type.sigmaHalf
					&& types[t].epsilonSqrt == type.epsilonSqrt && types[t].charge == type.charge; });
				if (same == candidates.end()) {
					candidates.push_back(static_cast<int>(types.size()));
					types.push_back(type);
				}
				else
					types[*same].count++;
				const Float3 velocity = box.pclusterInterimStates[pc].vels_prev[p];	// [nm/ns] = [m/s]
				kineticEnergy += 0.5 * meta.mass[p] * velocity.lenSquared();
				nAtoms++;
			}
		}

		// As GROMACS: the reference temperature, or without a thermostat the starting temperature, and without one either a 10% buffer
		constexpr double gasConstant = BOLTZMANNCONSTANT * AVOGADROSNUMBER;	// [J/mol/K]
		const double temperature = params.apply_thermostat ? params.ref_t.value
			: params.initial_temperature > 0.f ? params.initial_temperature.value
			: kineticEnergy * 2.0 / (3.0 * nAtoms * gasConstant);
		if (temperature <= 0.0)
			return 0.1f * cutoff;

		// Over the list's lifetime, the steps after the one it was built for [ns]
		const double elapsed = (params.stepsPerNlistupdate - 1) * double(params.dt);
		if (elapsed <= 0.0)
			return 0.f;
		const double driftTarget = nAtoms * params.verlet_buffer_tolerance * 1e6 * params.stepsPerNlistupdate * double(params.dt);	// [J/mol]
		const Float3 boxSize = box.boxparams.BoxSizeFloat();
		const double volume = double(boxSize.x) * boxSize.y * boxSize.z;	// [nm^3]

		// The variance of an atom's displacement along one axis over the list's lifetime [nm^2]. A constrained atom moves with the
		// center of mass of it and its partner, and rotates around it on 2 axes, at most uniformly over the sphere of its arm.
		// GROMACS fits that saturation (constrained_atom_sigma2), here it is approximated by 1 - exp(-angle variance)
		const double kTt2 = gasConstant * temperature * elapsed * elapsed;	// [J/mol ns^2] = [kg/mol nm^2]
		const auto displacementVariance = [&](const Atomtype& type) {
			if (type.partnerMass == 0.f)
				return kTt2 / type.mass;
			const double totalMass = double(type.mass) + type.partnerMass;
			const double massFraction = type.partnerMass / totalMass;
			const double arm = type.constraintLength * massFraction;
			const double rotation = kTt2 * massFraction / type.mass;	// Per rotational degree of freedom, as arc length
			return kTt2 / totalMass + 2.0 / 3.0 * arm * arm * (1.0 - std::exp(-rotation / (arm * arm)));
		};
		const double kappa = PhysicsUtils::CalcEwaldkappa(cutoff);

		// Per atom type pair: the variance of the displacement along one axis [nm^2], and the potential and its first 3 derivatives at the cutoff
		struct PairTerms { double pairs, variance, v0, d1, d2, d3; };
		std::vector<PairTerms> pairTerms;
		for (size_t i = 0; i < types.size(); i++) {
			for (size_t j = i; j < types.size(); j++) {
				const Atomtype& a = types[i];
				const Atomtype& b = types[j];
				const double sigma = double(a.sigmaHalf) + b.sigmaHalf;
				const double epsilon = double(a.epsilonSqrt) * b.epsilonSqrt;	// [J/mol]
				const double chargeTerm = double(PhysicsUtils::modifiedCoulombConstant) * a.charge * b.charge;	// [J/mol nm]
				const auto potential = [&](double r) {	// [J/mol]
					const double s6 = std::pow(sigma / r, 6);
					return 4.0 * epsilon * s6 * (s6 - 1.0) + (params.enable_electrostatics ? chargeTerm * std::erfc(kappa * r) / r : 0.0);
				};
				// Central differences, the potential is smooth around the cutoff
				constexpr double h = 1e-3;	// [nm]
				const double vm2 = potential(cutoff - 2 * h), vm1 = potential(cutoff - h), v0 = potential(cutoff), vp1 = potential(cutoff + h), vp2 = potential(cutoff + 2 * h);
				const double pairs = i == j ? 0.5 * double(a.count) * (a.count - 1) : double(a.count) * b.count;
				if (pairs == 0.0)
					continue;
				pairTerms.push_back({ pairs, displacementVariance(a) + displacementVariance(b), v0,
					(vp1 - vm1) / (2 * h), (vp1 - 2 * v0 + vm1) / (h * h), (vp2 - 2 * vp1 + 2 * vm1 - vm2) / (2 * h * h * h) });
			}
		}

		// GROMACS' surface_frac for 4 atom clusters, assumed a regular tetrahedron with the average particle distance as edge: the
		// share of an atom's pairs at the list radius that are not in the list, as another pair of their blocks is within it
		const double particleDistance = std::cbrt(volume / nAtoms);	// [nm]
		const auto surfaceFraction = [&](double listRadius) {
			if (listRadius < 0.5 * particleDistance)
				return 1.0;
			const double d = particleDistance / listRadius;
			const double d2 = d * d;
			const double areaRelative = 1.0 + 1.0 / PI * (6.0 * std::acos(1.0 / std::sqrt(3.0)) * d
				+ std::sqrt(3.0) * d2 * (1.0 + 5.0 / 18.0 * d2 + 7.0 / 45.0 * d2 * d2 + 83.0 / 756.0 * d2 * d2 * d2));
			return areaRelative / 4.0;
		};

		// The energy of the pairs within the cutoff that the list misses, as GROMACS' energyDriftAtomPair, plus the 0th order term
		const auto drift = [&](double buffer) {
			const double blockFactor = surfaceFraction(cutoff + buffer) * surfaceFraction(cutoff + buffer);
			double total = 0.0;
			for (const PairTerms& t : pairTerms) {
				const double s = std::sqrt(t.variance), s2 = t.variance, b2 = buffer * buffer;
				const double cExp = std::exp(-b2 / (2 * s2)) / std::sqrt(2 * PI);
				const double cErfc = 0.5 * std::erfc(buffer / (std::sqrt(2.0) * s));
				const double pot0 = t.v0 * (s * cExp - buffer * cErfc);
				const double pot1 = -t.d1 / 2 * ((b2 + s2) * cErfc - buffer * s * cExp);
				const double pot2 = t.d2 / 6 * (s * (b2 + 2 * s2) * cExp - buffer * (b2 + 3 * s2) * cErfc);
				const double pot3 = -t.d3 / 24 * ((b2 * b2 + 6 * b2 * s2 + 3 * s2 * s2) * cErfc - buffer * s * (b2 + 5 * s2) * cExp);
				// The pair density at the surface. GROMACS uses (rlist + s)^2 as the effective r^2
				const double surface = 4 * PI * (cutoff + buffer + s) * (cutoff + buffer + s) / volume;
				total += t.pairs * surface * blockFactor * std::abs(pot0 + pot1 + pot2 + pot3);
			}
			return total;
		};

		// The drift falls with the buffer. Bisect to 1e-4 nm, then round up to GROMACS' 1e-3 nm resolution
		double lo = 0.0, hi = 0.1;
		while (drift(hi) > driftTarget) {
			hi *= 2;
			if (hi > cutoff)
				throw std::invalid_argument(Lima::Format("The pair list buffer for verlet_buffer_tolerance={} would exceed the cutoff. "
					"Lower nstlist or dt, or raise the tolerance", params.verlet_buffer_tolerance.value));
		}
		while (hi - lo > 1e-4) {
			const double mid = 0.5 * (lo + hi);
			(drift(mid) > driftTarget ? lo : hi) = mid;
		}
		return static_cast<float>(std::ceil(hi * 1000.0) / 1000.0);
	}

	void Pack(EngineBatchData& batch, const std::vector<Simulation*>& simulations, const std::vector<bool>* active) {
		if (!active) Validate(simulations);
		else if (active->size() != simulations.size()) throw std::invalid_argument("Invalid active batch layout");
		batch.params = simulations.front()->simParams;
		batch.boxSize = simulations.front()->box->boxparams.boxSize;
		batch.ewaldKappa = PhysicsUtils::CalcEwaldkappa(batch.params.cutoff_nm);
		// The largest of the simulations' buffers, as the batch shares one list radius. Only set at the start, so the lists of a
		// simulation do not change when other members retire
		if (!active) {
			batch.listBuffer = 0.f;
			for (const Simulation* simulation : simulations)
				if (!simulation->simParams.em_variant)
					batch.listBuffer = std::max(batch.listBuffer, VerletBuffer(*simulation));
		}
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
			sim.thermostat = VRescale(static_cast<uint64_t>(simulation->simParams.thermostat_seed.value));
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
