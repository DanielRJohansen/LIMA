#pragma once

#include "LimaTypes.cuh"

#include <vector>

struct PersistentCluster;
struct PersistentClusterMeta;
struct BondGroups;

// Energy minimization is preconditioned FIRE 2.0 (Guénolé et al. 2020, Comput. Mater. Sci. 175, 109584). Instead of F/m the
// dynamics are driven by G = P^-1 F, where P is
//  - a per-atom estimate of the diagonal of the Hessian: the bonded stiffness of the atom plus a nonbonded floor, and
//  - for small molecules contained in a single pcluster (water, ions), a rigid-body term, since translating or rotating
//    those molecules only stretches nonbonded interactions, which are much softer than their bonds.
// With this P every mode has roughly unit curvature, so the timestep is unitless and ~1.
// Host/device-shared data lives here, kernels in EnergyMinimization.cuh.
namespace EM {
	struct Config {
		// FIRE 2.0. dt is unitless because G is a displacement [nm]
		float dtStart = 0.3f;
		float dtMax = 1.f;
		float dtMin = 0.006f;
		int nDelay = 20;
		float dtGrow = 1.1f;
		float dtShrink = 0.5f;
		float alphaStart = 0.25f;
		float alphaShrink = 0.99f;
		float alphaMin = 0.03f;				// Keeps local oscillations damped while the rest of the system still descends

		float maxDisplacement = 0.02f;		// [nm] per atom per step

		// Preconditioner
		float nonbondedStiffness = 3e7f;	// [J/mol/nm^2] Added to the bonded stiffness of every atom, the only stiffness of ions
		float rigidStiffness = 3e8f;		// [J/mol/nm^2] Stiffness of rigid-body motion of single-pcluster molecules
	};

	// Diagonal preconditioner, 1/K per pcluster slot [mol nm^2/J]. Empty slots get 0
	std::vector<float> ComputeInverseStiffness(const std::vector<PersistentCluster>& pclusters, const std::vector<PersistentClusterMeta>& pclusterMeta,
		const BondGroups& bonds, Float3 boxSize, float nonbondedStiffness);

	// 1 for pclusters that contain an entire molecule, i.e. no bonds to other pclusters
	std::vector<uint8_t> FindWholeMoleculePclusters(size_t nPclusters, const BondGroups& bonds);

	// Sums over the particles of one simulation, first per block and then in a fixed order, so the result is deterministic
	// and independent of which other simulations share the batch
	struct Sums {
		double forceDotVelocity = 0.;		// Σ F·v, the power
		double velocityDotVelocity = 0.;	// Σ v·v
		double gDotVelocity = 0.;			// Σ G·v
		double gDotG = 0.;					// Σ G·G
		float maxForceSq = 0.f;				// max |F|^2 [(J/mol/nm)^2]
		int nParticles = 0;
	};

	struct SimState {
		float maxForce = 0.f;				// [J/mol/nm] of the latest step

		float dt = 0.f;
		float alpha = 0.f;
		int nPositive = 0;
		int iteration = 0;

		// Decided per step, consumed by the update kernel
		int reset = 0;						// Power was not positive: step half back and restart from rest
		float velocityScale = 0.f;			// Velocity mixing: v = velocityScale * v + forceScale * G
		float forceScale = 0.f;
	};

	struct ParticleState {
		Float3 velocity{};					// [nm]
	};
}
