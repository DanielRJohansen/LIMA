#pragma once

#include "PhysicsUtils.cuh"

// Functions optimized for CUDA
namespace PhysicsUtilsDevice {

	/// <summary></summary>
	/// <param name="myCharge">[kilo C/mol]</param>	// TODO: These can probably be half for performance gains
	/// <param name="otherCharge">[kilo C/mol]</param>
	/// <param name="diff">self-other [nm]</param>
	/// <returns>[1/l N/mol]</returns>
	__device__ __host__ inline Float3 CalcCoulumbForce(const float myCharge, const float otherCharge, const Float3& diff) 
	{
		const float modifiedCoulombConstant = COULOMBCONSTANT /NANO / NANO / AVOGADROSNUMBER * LIMA * KILO * KILO;	// [1/l N/mol nm^2 / (kilo C/mol)^2]

		const float invLen = rsqrtf(diff.lenSquared());                  // Computes 1 / sqrt(lenSquared)
		const float invLenCubed = invLen * invLen * invLen;       // Computes (1 / |diff|^3)

		return diff * (modifiedCoulombConstant * myCharge * otherCharge * invLenCubed);
	}

	// <summary></summary>
	// <param name="myCharge">[kilo C/mol]</param>
	// <param name="otherCharge">[kilo C/mol]</param>
	// <param name="diff">[nm]</param>
	// <returns>[J/mol]</returns>
	__device__ __host__ inline float CalcCoulumbPotential(const float myCharge, const float otherCharge, const Float3& diff)
	{
		// N * m = J
		const float modifiedCoulombConstant = COULOMBCONSTANT / NANO / AVOGADROSNUMBER * KILO * KILO;	// [J/mol * nm / (kilo C/mol)^2] 

		return modifiedCoulombConstant * (myCharge * otherCharge) * rsqrtf(diff.lenSquared());
	}
}