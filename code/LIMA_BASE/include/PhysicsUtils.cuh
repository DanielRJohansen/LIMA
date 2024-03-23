#pragma once

namespace PhysicsUtils {
	__device__ __host__ constexpr static float calcKineticEnergy(const float velocity, const float mass) {
		return 0.5f * mass * velocity * velocity;
	}

	// Calculate mean speed of particles. [K], [kg/mol]
	//sqrt((8RT)/(piM))
	static float tempToVelocity(double temperature /*[K]*/, double mass /*[kg/mol]*/) {
		return sqrtf(static_cast<float>(3.0 * BOLTZMANNCONSTANT * temperature / (mass / AVOGADROSNUMBER)));
	}
	// http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/kintem.html

	static float constexpr kineticEnergyToTemperature(long double kineticEnergy /*[J]*/, int numParticles) {
		const double temperature = kineticEnergy * (2.0 / 3.0) / (BOLTZMANNCONSTANT * numParticles);
		return static_cast<float>(temperature);
	}
}