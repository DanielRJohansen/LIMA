#pragma once

namespace PhysicsUtils {

	/// <summary></summary>
	/// <param name="speed">[m/s]</param>
	/// <param name="mass">[kg/mol]</param>
	/// <returns>[J/mol]</returns>
	__device__ __host__ constexpr static float calcKineticEnergy(const float speed, const float mass) {
		return 0.5f * mass * speed * speed;
	}

	// Calculate mean speed of particles. [K], [kg/mol]
	__device__ __host__ static float tempToVelocity(double temperature /*[K]*/, double mass /*[kg/mol]*/) {
		return sqrtf(static_cast<float>(3.0 * BOLTZMANNCONSTANT * temperature / (mass / AVOGADROSNUMBER)));
	}
	// http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/kintem.html

	// For multiple particles
	inline float constexpr kineticEnergyToTemperature(long double kineticEnergy /*[J/mol]*/, int numParticles) {
		const double temperature = kineticEnergy * (2.0 / 3.0) / (BOLTZMANNCONSTANT * numParticles);
		return static_cast<float>(temperature);
	}

	// For single particle
	__device__ __host__ inline constexpr float kineticEnergyToTemperature(float kineticEnergy /*[J/mol]*/) {
		return kineticEnergy * (2.0f / 3.0f) / (BOLTZMANNCONSTANT * AVOGADROSNUMBER);
	}

	/// <summary></summary>
	/// <param name="myCharge">[C/mol]</param>
	/// <param name="otherCharge">[C/mol]</param>
	/// <param name="diff">[nm]</param>
	/// <returns>[GN/mol]</returns>
	__device__ __host__ inline Float3 CalcCoulumbForce(const float myCharge, const float otherCharge, const Float3& diff /*self - other*/) 
	{
		const float modifiedCoulombConstant = COULOMBCONSTANT * 1e9 * 1e9 / AVOGADROSNUMBER * 1e-9;	// [GN/mol nm^2 / (C/mol)^2]

		return diff.norm() * modifiedCoulombConstant * (myCharge * otherCharge) / diff.lenSquared();
	}

	// <summary></summary>
	// <param name="myCharge">[C/mol]</param>
	// <param name="otherCharge">[C/mol]</param>
	// <param name="distance">[nm]</param>
	// <returns>[GJ/mol]</returns>
	__device__ __host__ inline constexpr float CalcCoulumbPotential(const float myCharge, const float otherCharge, const float distance) 
	{
		const float modifiedCoulombConstant = COULOMBCONSTANT * 1e9 / AVOGADROSNUMBER * 1e-9;	// [GN/mol * m * nm / (C/mol)^2]

		return modifiedCoulombConstant * (myCharge * otherCharge) / distance;
	}
}