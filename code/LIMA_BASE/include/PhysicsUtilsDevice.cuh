#pragma once

#include "PhysicsUtils.cuh"

// Functions optimized for CUDA
namespace PhysicsUtilsDevice {
	using PhysicsUtils::modifiedCoulombConstant_Force;
	using PhysicsUtils::modifiedCoulombConstant_Potential;

	// <summary>Calculate the force without multiplying the coulumbConstant, so called must do that!!</summary>	
	/// <param name="myCharge">[kilo C/mol]</param>	// TODO: These can probably be half for performance gains
	/// <param name="otherCharge">[kilo C/mol]</param>
	/// <param name="diff">self-other [nm]</param>
	/// <returns>[1/l N/mol]</returns>
	__device__ inline Float3 CalcCoulumbForce_optim(const float myCharge, const float otherCharge, const Float3& diff) 
	{		
		const float invLen = rsqrtf(diff.lenSquared());                  // Computes 1 / sqrt(lenSquared)
		const float invLenCubed = invLen * invLen * invLen;       // Computes (1 / |diff|^3)

		const Float3 force = diff * myCharge * otherCharge * invLenCubed;
#ifdef FORCE_NAN_CHECK
		if (force.isNan())
			force.print('E');
#endif

		return force;
	}

	// <summary>Calculate the potential without multiplying the coulumbConstant, so called must do that!!</summary>
	// <param name="myCharge">[kilo C/mol]</param>
	// <param name="otherCharge">[kilo C/mol]</param>
	// <param name="diff">[nm]</param>
	// <returns>[J/mol   /   modifiedCoulombConstant_Force ]</returns>
	//constexpr float modifiedCoulombConstant_Force = 1.f;
	__device__ inline float CalcCoulumbPotential_optim(const float myCharge, const float otherCharge, const Float3& diff)
	{
		// N * m = J		
		return (myCharge * otherCharge) * rsqrtf(diff.lenSquared());
	}
}