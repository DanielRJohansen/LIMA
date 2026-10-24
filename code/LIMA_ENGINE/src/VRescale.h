#pragma once

#include "GpuRuntime.h"

#include <cmath>
#include <cstdint>

// std::mt19937_64 as the standard defines it, so its output is the same, usable on the device. Trivially copyable, so its state can be
// copied to and from the device
struct Mt19937_64 {
	static constexpr int n = 312, m = 156;
	uint64_t state[n];
	int index = n;

	__host__ __device__ explicit Mt19937_64(uint64_t seed = 5489u) {
		state[0] = seed;
		for (int i = 1; i < n; i++)
			state[i] = 6364136223846793005ull * (state[i - 1] ^ (state[i - 1] >> 62)) + i;
	}
	__host__ __device__ uint64_t operator()() {
		if (index >= n) {
			constexpr uint64_t upper = 0xFFFFFFFF80000000ull, lower = 0x7FFFFFFFull;
			for (int i = 0; i < n; i++) {
				const uint64_t x = (state[i] & upper) | (state[(i + 1) % n] & lower);
				state[i] = state[(i + m) % n] ^ (x >> 1) ^ ((x & 1) ? 0xB5026F5AA96619E9ull : 0ull);
			}
			index = 0;
		}
		uint64_t y = state[index++];
		y ^= (y >> 29) & 0x5555555555555555ull;
		y ^= (y << 17) & 0x71D67FFFEDA60000ull;
		y ^= (y << 37) & 0xFFF7EEE000000000ull;
		y ^= y >> 43;
		return y;
	}
};

// The v-rescale thermostat (Bussi, Donadio & Parrinello 2007, GROMACS' tcoupl = v-rescale): each coupling rescales the kinetic
// energy to a random value that relaxes to the reference temperature with time constant tau_t, such that it samples the canonical
// distribution. Deterministic: the generator's output is fixed by the standard, and the normal and gamma deviates are computed here
// in a fixed order rather than by std's distributions, which differ between standard libraries. Runs on the host and the device,
// whose math functions may differ in the last bit
class VRescale {
	Mt19937_64 rng;

	__host__ __device__ double Uniform() { return (static_cast<double>(rng() >> 11) + 0.5) * 0x1.0p-53; }	// In (0, 1)
	__host__ __device__ double Normal() {	// Box-Muller
		const double u1 = Uniform();
		const double u2 = Uniform();
		return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
	}
	// Marsaglia & Tsang 2000, for shape >= 1
	__host__ __device__ double Gamma(double shape) {
		const double d = shape - 1.0 / 3.0;
		const double c = 1.0 / sqrt(9.0 * d);
		while (true) {
			const double x = Normal();
			const double v1 = 1.0 + c * x;
			if (v1 <= 0.0)
				continue;
			const double v = v1 * v1 * v1;
			const double u = Uniform();
			if (u < 1.0 - 0.0331 * x * x * x * x || log(u) < 0.5 * x * x + d * (1.0 - v + log(v)))
				return d * v;
		}
	}
	// The sum of n squared standard normal deviates, as GROMACS' vrescale_sumnoises
	__host__ __device__ double SumNoises(double n) {
		if (n < 2.0) {
			double sum = 0.0;
			for (long i = 0; i < lround(n); i++) {
				const double g = Normal();
				sum += g * g;
			}
			return sum;
		}
		return 2.0 * Gamma(0.5 * n);
	}

public:
	__host__ __device__ explicit VRescale(uint64_t seed = 1) : rng(seed) {}

	// The factor for the velocities, from the current temperature [K]. tauInCouplings is tau_t in units of the coupling interval.
	// As GROMACS' vrescale_resamplekin, with the kinetic energies in temperature units
	__host__ __device__ double Factor(double temperature, double referenceTemperature, int64_t degreesOfFreedom, double tauInCouplings) {
		if (!(temperature > 0.0) || degreesOfFreedom <= 1)
			return 1.0;
		const double decay = tauInCouplings > 0.1 ? exp(-1.0 / tauInCouplings) : 0.0;
		const double n = static_cast<double>(degreesOfFreedom);
		const double r = Normal();
		const double noises = SumNoises(n - 1.0);
		const double newTemperature = temperature + (1.0 - decay) * (referenceTemperature * (noises + r * r) / n - temperature)
			+ 2.0 * r * sqrt(temperature * referenceTemperature / n * (1.0 - decay) * decay);
		return sqrt((newTemperature > 0.0 ? newTemperature : 0.0) / temperature);
	}
};
