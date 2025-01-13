#pragma once

namespace PhysicsUtils {

	constexpr float modifiedCoulombConstant = COULOMBCONSTANT / NANO / AVOGADROSNUMBER * KILO * KILO;	// [J/mol*nm / (kilo C/mol)^2]

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

	// For single particle - this is only correct for IDEAL GAS, make obsolete!
	__device__ __host__ inline constexpr float kineticEnergyToTemperature(float kineticEnergy /*[J/mol]*/) {
		return kineticEnergy * (2.0f / 3.0f) / (BOLTZMANNCONSTANT * AVOGADROSNUMBER);
	}
	
	__device__ __host__ inline constexpr float kineticEnergyToTemperature(double totalKineticEnergy /*[J/mol]*/, int64_t degreesOfFreedom) {
		return totalKineticEnergy * (2.0 / static_cast<double>(degreesOfFreedom)) / (BOLTZMANNCONSTANT * AVOGADROSNUMBER);
	}

	/// <summary>Slow host version, faster version in PhysicsUtilsDevice</summary>
	/// <param name="myCharge">[kilo C/mol]</param>	// TODO: These can probably be half for performance gains
	/// <param name="otherCharge">[kilo C/mol]</param>
	/// <param name="diff">self-other [nm]</param>
	/// <returns>[J/mol/nm]</returns>	
	__host__ inline Float3 CalcCoulumbForce(const float myCharge, const float otherCharge, const Float3& diff) 
	{
		return diff.norm() * modifiedCoulombConstant * (myCharge * otherCharge) / diff.lenSquared();
	}

	// <summary>Slow host version, faster version in PhysicsUtilsDevice</summary>
	// <param name="myCharge">[kilo C/mol]</param>
	// <param name="otherCharge">[kilo C/mol]</param>
	// <param name="distance">[nm]</param>
	// <returns>[J/mol]</returns>
	__host__ inline constexpr float CalcCoulumbPotential(const float myCharge, const float otherCharge, const float distance) 
	{		
		return modifiedCoulombConstant * (myCharge * otherCharge) / distance;
	}

	/// <summary>
	/// Computes the ewald splitting parameter Kappa
	/// </summary>
	/// <returns>[1/nm]</returns>
	constexpr float CalcEwaldkappa(float cutoffNM) {
		return 3.f / cutoffNM; // erfc(3) will yield a tolerance of 2e-5.. 
	}
}