#pragma once

#include "LimaTypes.cuh"
#include "Bodies.cuh"

// Deterministic accumulation of all forces on each particle in MD: the NB, bonded, PME and SNF kernels convert their
// partial forces to 64-bit fixed point and sum them with integer atomics. Unlike float atomics, integer addition is
// associative, so the sum is bitwise independent of the order the atomics arrive in, and the kernels may run concurrently.
// The integrate kernel takes each particle's sum, leaving the accumulator zeroed for the next step.
// Not used in EM, where forces can exceed the fixed point range.
struct ForceAccumulator {
	static constexpr float scale = 16777216.f;		// 2^24 -> resolution 6e-8 J/mol/nm, range +-5.5e11 J/mol/nm
	static constexpr float scaleInv = 1.f / scale;	// Power of 2, so ToFloat is exactly the rounded sum

	// SoA, indexed by the particle's slot scId * SuperCluster::maxParticles + indexInSupercluster.
	// potE is nullptr on steps that don't log data, and is then neither summed nor touched
	unsigned long long* fx = nullptr;
	unsigned long long* fy = nullptr;
	unsigned long long* fz = nullptr;
	unsigned long long* potE = nullptr;

	__device__ static unsigned long long ToFixed(float v) { return static_cast<unsigned long long>(llrintf(v * scale)); }
	__device__ static float ToFloat(unsigned long long v) { return static_cast<float>(static_cast<long long>(v)) * scaleInv; }

	// withPotE must be potE != nullptr. Compile time for the NB kernel, where the check costs
	template <bool withPotE>
	__device__ void Add(int slot, const ForceEnergy& fe) const {
		AddFixed<withPotE>(slot, ToFixed(fe.force.x), ToFixed(fe.force.y), ToFixed(fe.force.z), withPotE ? ToFixed(fe.potE) : 0);
	}
	__device__ void Add(int slot, const ForceEnergy& fe) const {
		if (potE) Add<true>(slot, fe);
		else Add<false>(slot, fe);
	}

	// For sources that already sum in the same fixed point, adding their sum exactly
	template <bool withPotE>
	__device__ void AddFixed(int slot, unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long e) const {
		atomicAdd(&fx[slot], x);
		atomicAdd(&fy[slot], y);
		atomicAdd(&fz[slot], z);
		if constexpr (withPotE)
			atomicAdd(&potE[slot], e);
	}
	__device__ void AddFixed(int slot, unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long e) const {
		if (potE) AddFixed<true>(slot, x, y, z, e);
		else AddFixed<false>(slot, x, y, z, e);
	}

	// Returns the particle's summed forces and zeroes its slot, so the accumulator is ready for the next step
	template <bool withPotE>
	__device__ ForceEnergy Take(int slot) const {
		const ForceEnergy fe{ Float3{ ToFloat(fx[slot]), ToFloat(fy[slot]), ToFloat(fz[slot]) }, withPotE ? ToFloat(potE[slot]) : 0.f };
		fx[slot] = 0;
		fy[slot] = 0;
		fz[slot] = 0;
		if constexpr (withPotE)
			potE[slot] = 0;
		return fe;
	}
};
