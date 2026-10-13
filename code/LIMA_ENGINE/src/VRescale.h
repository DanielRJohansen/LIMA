#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>

// The v-rescale thermostat (Bussi, Donadio & Parrinello 2007, GROMACS' tcoupl = v-rescale): each coupling rescales the kinetic
// energy to a random value that relaxes to the reference temperature with time constant tau_t, such that it samples the canonical
// distribution. Deterministic: std::mt19937_64's output is fixed by the standard, and the normal and gamma deviates are computed here
// in a fixed order rather than by std's distributions, which differ between standard libraries
class VRescale {
	std::mt19937_64 rng;

	double Uniform() { return (static_cast<double>(rng() >> 11) + 0.5) * 0x1.0p-53; }	// In (0, 1)
	double Normal() {	// Box-Muller
		const double u1 = Uniform();
		const double u2 = Uniform();
		return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * 3.14159265358979323846 * u2);
	}
	// Marsaglia & Tsang 2000, for shape >= 1
	double Gamma(double shape) {
		const double d = shape - 1.0 / 3.0;
		const double c = 1.0 / std::sqrt(9.0 * d);
		while (true) {
			const double x = Normal();
			const double v1 = 1.0 + c * x;
			if (v1 <= 0.0)
				continue;
			const double v = v1 * v1 * v1;
			const double u = Uniform();
			if (u < 1.0 - 0.0331 * x * x * x * x || std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v)))
				return d * v;
		}
	}
	// The sum of n squared standard normal deviates, as GROMACS' vrescale_sumnoises
	double SumNoises(double n) {
		if (n < 2.0) {
			double sum = 0.0;
			for (long i = 0; i < std::lround(n); i++) {
				const double g = Normal();
				sum += g * g;
			}
			return sum;
		}
		return 2.0 * Gamma(0.5 * n);
	}

public:
	explicit VRescale(uint64_t seed = 1) : rng(seed) {}

	// The factor for the velocities, from the current temperature [K]. tauInCouplings is tau_t in units of the coupling interval.
	// As GROMACS' vrescale_resamplekin, with the kinetic energies in temperature units
	double Factor(double temperature, double referenceTemperature, int64_t degreesOfFreedom, double tauInCouplings) {
		if (!(temperature > 0.0) || degreesOfFreedom <= 1)
			return 1.0;
		const double decay = tauInCouplings > 0.1 ? std::exp(-1.0 / tauInCouplings) : 0.0;
		const double n = static_cast<double>(degreesOfFreedom);
		const double r = Normal();
		const double noises = SumNoises(n - 1.0);
		const double newTemperature = temperature + (1.0 - decay) * (referenceTemperature * (noises + r * r) / n - temperature)
			+ 2.0 * r * std::sqrt(temperature * referenceTemperature / n * (1.0 - decay) * decay);
		return std::sqrt(std::max(newTemperature, 0.0) / temperature);
	}
};
