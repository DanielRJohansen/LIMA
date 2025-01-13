#pragma once

#include "KernelConstants.cuh"
#include "DeviceAlgorithms.cuh"

namespace LAL {
	__device__ void CalcBspline(float f, float* w) {
		if constexpr (!USE_PRECOMPUTED_BSPLINES) {
			w[0] = (1.f - f) * (1.f - f) * (1.f - f) / 6.f;
			w[1] = (4.f - 6.f * f * f + 3.f * f * f * f) / 6.f;
			w[2] = (1.f + 3.f * f + 3.f * f * f - 3.f * f * f * f) / 6.f;
			w[3] = (f * f * f) / 6.f;
		}
		else {
			const float alpha = f * (BSPLINE_LUT_SIZE - 1);
			const int index = static_cast<int>(floor(alpha));
			const int idxUp = min(index + 1, BSPLINE_LUT_SIZE - 1);
			const float frac = alpha - index;   // in [0,1)
			w[0] = lerp(bsplineTable_device[index], bsplineTable_device[idxUp], frac);
			w[1] = lerp(bsplineTable_device[BSPLINE_LUT_SIZE + index], bsplineTable_device[BSPLINE_LUT_SIZE + idxUp], frac);

			float alphaInv = (1.f - f) * (BSPLINE_LUT_SIZE - 1);
			int idxInv = (int)floor(alphaInv);
			int idxInvUp = min(idxInv + 1, BSPLINE_LUT_SIZE - 1);
			float fracInv = alphaInv - idxInv;
			w[2] = lerp(bsplineTable_device[BSPLINE_LUT_SIZE + idxInv], bsplineTable_device[BSPLINE_LUT_SIZE + idxInvUp], fracInv);
			w[3] = lerp(bsplineTable_device[idxInv], bsplineTable_device[idxInvUp], fracInv);
		}		
	}
}


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
			if constexpr (!USE_PRECOMPUTED_ERFCSCALARS) {
				const float erfcTerm = erfc(diff.len() * ewaldkappa_device);
				const float scalar = erfcTerm + 2.f * ewaldkappa_device / sqrt(PI) * diff.len() * exp(-ewaldkappa_device * ewaldkappa_device * diff.lenSquared());
				force *= scalar;
			}
			else {
				const float distanceInArray = fminf(diff.len() * cutoffNmReciprocal_device * ERFC_LUT_SIZE - 1, ERFC_LUT_SIZE - 1);
				const int index = static_cast<int>(std::floor(distanceInArray));
				const int indexNext = std::min(index + 1, ERFC_LUT_SIZE - 1);
				const float frac = distanceInArray - static_cast<float>(index);
				const float scalar = LAL::lerp(erfcForcescalarTable_device[index], erfcForcescalarTable_device[indexNext], frac);// optim: look into using std::lerp

				force *= scalar;
			}
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
			if constexpr (!USE_PRECOMPUTED_ERFCSCALARS) {

				potential *= erfc(diff.len() * ewaldkappa_device);
			}
			else {
				const float distanceInArray = fminf(diff.len() * cutoffNmReciprocal_device * ERFC_LUT_SIZE - 1, ERFC_LUT_SIZE - 1);
				const int index = static_cast<int>(std::floor(distanceInArray));
				const int indexNext = std::min(index + 1, ERFC_LUT_SIZE - 1);
				const float frac = distanceInArray - static_cast<float>(index);
				const float scalar = LAL::lerp(erfcPotentialscalarTable_device[index], erfcForcescalarTable_device[indexNext], frac);// optim: look into using std::lerp
				potential *= scalar;
			}
		}

		return potential;
	}
}