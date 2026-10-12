#pragma once

#include "LimaTypes.cuh"
#include "Bodies.cuh"

// Deterministic accumulation of forces on each particle in MD: the NB, PME and SNF kernels convert their
// partial forces to 64-bit fixed point and sum them with integer atomics. Unlike float atomics, integer addition is
// associative, so the sum is bitwise independent of the order the atomics arrive in, and the kernels may run concurrently.
// Each particle's first bondgroup stores its sum in a separate supercluster-slot plane. Integration reads it
// coalesced and gathers only the remaining appearances, all as integers before the single float conversion.
// The integrate kernel combines both sums before converting to float, leaving the accumulator zeroed for the next step.
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

	// Gather the separately stored bonded sums before the single float conversion. No rounding-order change.
	// Zero the global accumulator slot; bonded entries are overwritten by their owning groups next step.
	template <bool withPotE>
	__device__ ForceEnergy Take(int slot, const ForceAccumulator& primary, const ulonglong4* const bonded, const int* const refs, int nSlots) const {
		unsigned long long x = fx[slot] + primary.fx[slot], y = fy[slot] + primary.fy[slot], z = fz[slot] + primary.fz[slot];
		unsigned long long e = withPotE ? potE[slot] + primary.potE[slot] : 0;
#pragma unroll
		for (int i = 0; i < 3; i++) {
			const int index = refs[i * nSlots + slot];
			if (index < 0) break;
			const ulonglong2 xy = __ldg(reinterpret_cast<const ulonglong2*>(&bonded[index]));
			x += xy.x;
			y += xy.y;
			if constexpr (withPotE) {
				const ulonglong2 ze = __ldg(reinterpret_cast<const ulonglong2*>(&bonded[index].z));
				z += ze.x;
				e += ze.y;
			}
			else z += __ldg(&bonded[index].z);
		}
		const ForceEnergy fe{ Float3{ ToFloat(x), ToFloat(y), ToFloat(z) }, withPotE ? ToFloat(e) : 0.f };
		fx[slot] = 0;
		fy[slot] = 0;
		fz[slot] = 0;
		if constexpr (withPotE)
			potE[slot] = 0;
		return fe;
	}
};
