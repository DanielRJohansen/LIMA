#pragma once

#include "PhysicsUtils.cuh"

// Functions optimized for CUDA
namespace PhysicsUtilsDevice {
	using PhysicsUtils::modifiedCoulombConstant;


	__device__ inline Float3 CalcCoulumbForce_optim(const float chargeProduct, const Float3& diff)
	{
		const float invLen = rsqrtf(diff.lenSquared());                  // Computes 1 / sqrt(lenSquared)
		const float invLenCubed = invLen * invLen * invLen;       // Computes (1 / |diff|^3)

		Float3 force = diff * chargeProduct * invLenCubed;
#ifdef FORCE_NAN_CHECK
		if (force.isNan())
			force.print('E');
#endif
		if constexpr (ENABLE_ERFC_FOR_EWALD) {
			const float erfcTerm = erfc(diff.len() * ewaldkappa_device);
			const float scalar = erfcTerm + 2.f * ewaldkappa_device / sqrt(PI) * diff.len() * exp(-ewaldkappa_device * ewaldkappa_device * diff.lenSquared());
			force *= scalar;
		}

		return force;
	}

	// <summary>Calculate the potential without multiplying the coulumbConstant, so called must do that!!</summary>
	// <param name="myCharge">[kilo C/mol]</param>
	// <param name="otherCharge">[kilo C/mol]</param>
	// <param name="diff">[nm]</param>
	// <returns>[J/mol   /   modifiedCoulombConstant ]</returns>
	//constexpr float modifiedCoulombConstant = 1.f;
	__device__ inline float CalcCoulumbPotential_optim(const float myCharge, const float otherCharge, const Float3& diff)
	{
		float potential = (myCharge * otherCharge) * rsqrtf(diff.lenSquared());
		if constexpr (ENABLE_ERFC_FOR_EWALD) {
			potential *= erfc(diff.len() * ewaldkappa_device);
		}
		
		return potential;
	}
}