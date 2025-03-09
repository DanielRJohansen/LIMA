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
			const int N = DeviceConstants::BSPLINE_LUT_SIZE;
			const float alpha = f * (N - 1);
			const int index = static_cast<int>(floor(alpha));
			const int idxUp = min(index + 1, N - 1);
			const float frac = alpha - index;   // in [0,1)
			w[0] = lerp(DeviceConstants::bsplineTable[index], DeviceConstants::bsplineTable[idxUp], frac);
			w[1] = lerp(DeviceConstants::bsplineTable[N + index], DeviceConstants::bsplineTable[N + idxUp], frac);

			float alphaInv = (1.f - f) * (N - 1);
			int idxInv = (int)floor(alphaInv);
			int idxInvUp = min(idxInv + 1, N - 1);
			float fracInv = alphaInv - idxInv;
			w[2] = lerp(DeviceConstants::bsplineTable[N + idxInv], DeviceConstants::bsplineTable[N + idxInvUp], fracInv);
			w[3] = lerp(DeviceConstants::bsplineTable[idxInv], DeviceConstants::bsplineTable[idxInvUp], fracInv);
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
				const float erfcTerm = erfc(diff.len() * DeviceConstants::ewaldKappa);
				const float scalar = erfcTerm + 2.f * DeviceConstants::ewaldKappa / sqrt(PI) * diff.len() * exp(-DeviceConstants::ewaldKappa * DeviceConstants::ewaldKappa * diff.lenSquared());
				force *= scalar;
			}
			else {
				const int N = DeviceConstants::ERFC_LUT_SIZE;
				const float distanceInArray = fminf(diff.len() * DeviceConstants::cutoffNmReciprocal * N - 1, N - 1);
				const int index = static_cast<int>(std::floor(distanceInArray));
				const int indexNext = std::min(index + 1, N - 1);
				const float frac = distanceInArray - static_cast<float>(index);
				const float scalar = LAL::lerp(DeviceConstants::erfcForcescalarTable[index], DeviceConstants::erfcForcescalarTable[indexNext], frac);// optim: look into using std::lerp

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
	__device__ inline float CalcCoulumbPotential_optim(const float chargeProduct, const Float3& diff)
	{
		float potential = (chargeProduct) * rsqrtf(diff.lenSquared());
		if constexpr (ENABLE_ERFC_FOR_EWALD) {
			if constexpr (!USE_PRECOMPUTED_ERFCSCALARS) {

				potential *= erfc(diff.len() * DeviceConstants::ewaldKappa);
			}
			else {
				const int N = DeviceConstants::ERFC_LUT_SIZE;
				const float distanceInArray = fminf(diff.len() * DeviceConstants::cutoffNmReciprocal * N - 1, N - 1);
				const int index = static_cast<int>(std::floor(distanceInArray));
				const int indexNext = std::min(index + 1, N - 1);
				const float frac = distanceInArray - static_cast<float>(index);
				const float scalar = LAL::lerp(DeviceConstants::erfcPotentialscalarTable[index], DeviceConstants::erfcForcescalarTable[indexNext], frac);// optim: look into using std::lerp
				potential *= scalar;
			}
		}

		return potential;
	}
}