#pragma once
// Included only by Engine.cu, the engine is a single compilation unit.

#include "EnergyMinimizationTypes.h"
#include "BoundaryCondition.cuh"
#include "BatchLayout.cuh"
#include "Bodies.cuh"

#include <cub/block/block_reduce.cuh>

namespace EM {
	// Solves I w = torque for a 3x3 inertia tensor, regularized for linear and single-atom bodies
	__device__ inline Float3 SolveRotation(float I[3][3], const Float3& torque) {
		const float regularization = 1e-4f; // [nm^2]
		for (int a = 0; a < 3; a++) I[a][a] += regularization;
		const auto Determinant = [](const float M[3][3]) {
			return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
				+ M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
		};
		const float det = Determinant(I);
		const float t[3] = { torque.x, torque.y, torque.z };
		float w[3];
		for (int a = 0; a < 3; a++) {
			// Cramer's rule
			float M[3][3];
			for (int r = 0; r < 3; r++)
				for (int c = 0; c < 3; c++)
					M[r][c] = c == a ? t[r] : I[r][c];
			w[a] = Determinant(M) / det;
		}
		return Float3{ w[0], w[1], w[2] };
	}

	// One thread per pcluster. G = F/K, plus (1/rigidStiffness) * P_rigid F for pclusters that are whole molecules,
	// where P_rigid F is the rigid motion (translation + rotation) closest to F in the least squares sense
	template <typename BoundaryCondition>
	__global__ void PreconditionKernel(const Config config, const PersistentCluster* const pclusters, const PersistentClusterMeta* const pcMeta,
		const Float3* const forces, const float* const inverseStiffness, const uint8_t* const wholeMolecule,
		Float3* const preconditionedForce, int nPclusters, Float3 boxSize) {
		const int pcId = blockIdx.x * blockDim.x + threadIdx.x;
		if (pcId >= nPclusters) return;
		constexpr int n = PersistentCluster::maxParticles;

		Float3 laneForces[n], positions[n];
		bool valid[n];
		int nValid = 0;
		int first = -1;
		for (int lane = 0; lane < n; lane++) {
			valid[lane] = pcMeta[pcId].particleIdsGlobal[lane] != -1;
			if (!valid[lane]) continue;
			laneForces[lane] = forces[pcId * n + lane];
			positions[lane] = pclusters[pcId].pqd[lane].position;
			if (first == -1) first = lane;
			nValid++;
		}
		if (nValid == 0) return;

		Float3 translation{}, rotation{}, centroid{};
		const bool rigid = wholeMolecule[pcId];
		if (rigid) {
			for (int lane = 0; lane < n; lane++) {
				if (!valid[lane]) continue;
				BoundaryCondition::applyHyperposNM(positions[first], positions[lane], boxSize);
				centroid += positions[lane];
				translation += laneForces[lane];
			}
			centroid = centroid / float(nValid);
			translation = translation / float(nValid);

			// Unit-mass inertia tensor about the centroid
			Float3 torque{};
			float I[3][3] = {};
			for (int lane = 0; lane < n; lane++) {
				if (!valid[lane]) continue;
				const Float3 r = positions[lane] - centroid;
				torque += r.cross(laneForces[lane]);
				const float rr = r.dot(r);
				const float rv[3] = { r.x, r.y, r.z };
				for (int a = 0; a < 3; a++)
					for (int b = 0; b < 3; b++)
						I[a][b] += (a == b ? rr : 0.f) - rv[a] * rv[b];
			}
			rotation = SolveRotation(I, torque);
		}

		for (int lane = 0; lane < n; lane++) {
			if (!valid[lane]) continue;
			const int slot = pcId * n + lane;
			Float3 g = laneForces[lane] * inverseStiffness[slot];
			if (rigid)
				g += (translation + rotation.cross(positions[lane] - centroid)) * (1.f / config.rigidStiffness);
			preconditionedForce[slot] = g;
		}
	}

	struct CombineSums {
		__device__ Sums operator()(const Sums& a, const Sums& b) const {
			Sums sum;
			sum.forceDotVelocity = a.forceDotVelocity + b.forceDotVelocity;
			sum.velocityDotVelocity = a.velocityDotVelocity + b.velocityDotVelocity;
			sum.gDotVelocity = a.gDotVelocity + b.gDotVelocity;
			sum.gDotG = a.gDotG + b.gDotG;
			sum.maxForceSq = fmaxf(a.maxForceSq, b.maxForceSq);
			sum.nParticles = a.nParticles + b.nParticles;
			return sum;
		}
	};

	// FIRE 2.0 decisions, from the sums over the whole simulation
	__device__ inline void Decide(const Config& config, const Sums& sums, SimState& state) {
		state.maxForce = sqrtf(sums.maxForceSq);
		if (state.iteration == 0) {
			state.dt = config.dtStart;
			state.alpha = config.alphaStart;
		}

		double velocityDotVelocity = sums.velocityDotVelocity;
		double gDotVelocity = sums.gDotVelocity;
		if (sums.forceDotVelocity > 0.) {
			state.nPositive++;
			if (state.nPositive > config.nDelay) {
				state.dt = fminf(state.dt * config.dtGrow, config.dtMax);
				state.alpha = fmaxf(state.alpha * config.alphaShrink, config.alphaMin);
			}
			state.reset = 0;
		}
		else {
			state.nPositive = 0;
			if (state.iteration >= config.nDelay) {
				state.dt = fmaxf(state.dt * config.dtShrink, config.dtMin);
				state.alpha = config.alphaStart;
			}
			state.reset = 1;
			velocityDotVelocity = 0.;
			gDotVelocity = 0.;
		}

		// |v| after the kick v += dt * G, used to mix v towards G
		const double dt = state.dt;
		const double velocityNormSq = velocityDotVelocity + 2. * dt * gDotVelocity + dt * dt * sums.gDotG;
		state.velocityScale = 1.f - state.alpha;
		state.forceScale = sums.gDotG > 0. ? static_cast<float>(state.alpha * sqrt(velocityNormSq / sums.gDotG)) : 0.f;
		state.iteration++;
	}

	// gridDim = (ceil(maxSlotsPerSimulation / BlockSize), nSimulations). Each block stores its partial sums, and the last block
	// to finish for a simulation sums the partials in a fixed order and makes the FIRE decisions
	template <int BlockSize>
	__global__ void ReduceAndDecideKernel(const Config config, const IntegrationSimulationData* const simulations, const PersistentClusterMeta* const pcMeta,
		const Float3* const forces, const Float3* const preconditionedForce, const ParticleState* const particles,
		Sums* const blockSums, unsigned int* const nBlocksDone, SimState* const states) {
		using BlockReduce = cub::BlockReduce<Sums, BlockSize>;
		__shared__ typename BlockReduce::TempStorage tempStorage;
		__shared__ bool isLastBlock;

		const int simulationId = blockIdx.y;
		const BatchRange pclusters = simulations[simulationId].pclusters;
		const int localSlot = blockIdx.x * BlockSize + threadIdx.x;

		Sums values{};
		if (localSlot < pclusters.count * PersistentCluster::maxParticles) {
			const int slot = pclusters.offset * PersistentCluster::maxParticles + localSlot;
			if (pcMeta[slot / PersistentCluster::maxParticles].particleIdsGlobal[slot % PersistentCluster::maxParticles] != -1) {
				const Float3 force = forces[slot];
				const Float3 g = preconditionedForce[slot];
				const Float3 velocity = particles[slot].velocity;
				values.forceDotVelocity = force.dot(velocity);
				values.velocityDotVelocity = velocity.dot(velocity);
				values.gDotVelocity = g.dot(velocity);
				values.gDotG = g.dot(g);
				// fmaxf drops NaN, which would let a broken force look converged. Report it as infinite instead
				const float forceSq = force.lenSquared();
				values.maxForceSq = isfinite(forceSq) ? forceSq : __int_as_float(0x7f800000); // +inf
				values.nParticles = 1;
			}
		}

		Sums* const simulationBlockSums = blockSums + size_t(simulationId) * gridDim.x;
		const Sums blockTotal = BlockReduce(tempStorage).Reduce(values, CombineSums{});
		if (threadIdx.x == 0) {
			simulationBlockSums[blockIdx.x] = blockTotal;
			__threadfence();
			isLastBlock = atomicAdd(&nBlocksDone[simulationId], 1u) == gridDim.x - 1;
		}
		__syncthreads();
		if (!isLastBlock) return;

		// Every other block has fenced its sums before counting itself done. Read past L1, which may be stale
		__threadfence();
		Sums partial{};
		for (int block = threadIdx.x; block < gridDim.x; block += BlockSize) {
			const Sums* const source = &simulationBlockSums[block];
			Sums other;
			other.forceDotVelocity = __ldcg(&source->forceDotVelocity);
			other.velocityDotVelocity = __ldcg(&source->velocityDotVelocity);
			other.gDotVelocity = __ldcg(&source->gDotVelocity);
			other.gDotG = __ldcg(&source->gDotG);
			other.maxForceSq = __ldcg(&source->maxForceSq);
			other.nParticles = __ldcg(&source->nParticles);
			partial = CombineSums{}(partial, other);
		}
		__syncthreads(); // tempStorage is reused
		const Sums total = BlockReduce(tempStorage).Reduce(partial, CombineSums{});
		if (threadIdx.x == 0) {
			if (total.nParticles > 0)
				Decide(config, total, states[simulationId]);
			nBlocksDone[simulationId] = 0;
		}
	}

	// Same launch layout as SuperclusterIntegrateKernel: blockDim = (16, 4, 1)
	template <typename BoundaryCondition>
	__global__ void UpdateKernel(const Config config, const SimState* const states, SuperCluster* const superClusters, const SuperClusterMeta* const scMeta,
		PersistentCluster* const pclusters, const Float3* const preconditionedForce, ParticleState* const particles, int nSuperclusters, Float3 boxSize) {
		constexpr int nScsPerBlock = 4;
		const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? static_cast<int>(blockIdx.x * nScsPerBlock + threadIdx.y) : -1;

		__shared__ Float3 p0s[nScsPerBlock];
		if (threadIdx.x == 0) {
			p0s[threadIdx.y] = scIdGlobal == -1 ? Float3{} : superClusters[scIdGlobal].Position(0);
			BoundaryCondition::applyBCNM(p0s[threadIdx.y], boxSize);
		}
		__syncthreads();

		if (scIdGlobal == -1) return;
		const int pidInPcluster = scMeta[scIdGlobal].indexInPcluster[threadIdx.x];
		const int pcIdGlobal = scMeta[scIdGlobal]._pclusterIds[threadIdx.x];
		const int pidGlobal = scMeta[scIdGlobal].globalParticleIds[threadIdx.x];
		if (pidGlobal == -1) return;

		const SimState& state = states[scMeta[scIdGlobal].simulationId];
		const int slot = pcIdGlobal * PersistentCluster::maxParticles + pidInPcluster;
		const float dt = state.dt;

		Float3 pos = superClusters[scIdGlobal].Position(threadIdx.x);
		BoundaryCondition::applyHyperposNM(p0s[threadIdx.y], pos, boxSize);
		const Float3 g = preconditionedForce[slot];
		Float3 velocity = particles[slot].velocity;

		if (state.reset) {
			pos = pos - velocity * (0.5f * dt);
			velocity = Float3{};
		}
		velocity = velocity + g * dt;
		velocity = velocity * state.velocityScale + g * state.forceScale;
		const float maxSpeed = config.maxDisplacement / dt;
		if (velocity.lenSquared() > maxSpeed * maxSpeed)
			velocity = velocity * (maxSpeed / velocity.len());
		pos = pos + velocity * dt;

		particles[slot].velocity = velocity;
		superClusters[scIdGlobal].posX[threadIdx.x] = pos.x;
		superClusters[scIdGlobal].posY[threadIdx.x] = pos.y;
		superClusters[scIdGlobal].posZ[threadIdx.x] = pos.z;
		pclusters[pcIdGlobal].pqd[pidInPcluster].position = pos;
	}
}
