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


// If we're gonna use this we need to deal with the license. But we dont need 0-4, so rather remake it for out purposes..
	__forceinline__ __device__ float raw_ex2(float a)
	{
		float r;
		asm("ex2.approx.ftz.f32 %0,%1;" : "=f"(r) : "f"(a));
		return r;
	}
	/* Approximate erfc(x) on [0, 4] with maximum absolute error of 1.15291e-7,
	   maximum relative error of 3.36619e-4, and maximum ulp error of 5621.588.
	*/
	__forceinline__ __device__ float fasterfc(float a)
	{
		float t;
		t = -1.64611265e-6f;   // -0x1.b9e000p-20
		t = fmaf(t, a, 2.95254722e-5f);  //  0x1.ef5af0p-16
		t = fmaf(t, a, -2.33422339e-4f);  // -0x1.e985aap-13
		t = fmaf(t, a, 1.04246172e-3f);  //  0x1.11466cp-10
		t = fmaf(t, a, -2.55015842e-3f);  // -0x1.4e411ep-9
		t = fmaf(t, a, 3.19798535e-4f);  //  0x1.4f5544p-12
		t = fmaf(t, a, 2.76054665e-2f);  //  0x1.c449b8p-6
		t = fmaf(t, a, -1.48274124e-1f);  // -0x1.2faa58p-3
		t = fmaf(t, a, -9.18447673e-1f);  // -0x1.d63ec6p-1
		t = fmaf(t, a, -1.62790680e+0f);  // -0x1.a0be80p+0
		t = t * a;
		return raw_ex2(t);
	}


	/*inline float CalcErfcScalar(float dist, float distSq) {
		const float erfcTerm = erfc(dist * DeviceConstants::ewaldKappaHardcoded);
		const float scalar = erfcTerm + 2.f * DeviceConstants::ewaldKappaHardcoded / PI_sqrt * dist * exp(-DeviceConstants::ewaldKappaHardcoded * DeviceConstants::ewaldKappaHardcoded * distSq);
		return scalar;
	}*/
	__device__ inline float CalcErfcScalar(float dist, float distSq) {
		//float kappa = 3.f / 1.2f;
		float kappa = DeviceConstants::ewaldKappa;
		float erfcTerm = fasterfc(dist * kappa);
		//const float erfcTerm = erfc(dist * kappa);
		float scalar = erfcTerm + 2.f * kappa / PI_sqrt * dist * exp(-kappa * kappa * distSq);

		return scalar;
	}

	/// <summary>
	/// Calculate the force without multiplying the coulumbConstant, so caller must do that!!
	/// </summary>
	/// <param name="chargeProduct"></param>
	/// <param name="diff"></param>
	/// <returns>[]</returns>
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
				/*const float erfcTerm = erfc(diff.len() * DeviceConstants::ewaldKappa);
				const float scalar = erfcTerm + 2.f * DeviceConstants::ewaldKappa / sqrt(PI) * diff.len() * exp(-DeviceConstants::ewaldKappa * DeviceConstants::ewaldKappa * diff.lenSquared());
				force *= scalar;*/

				float len = 1.f / invLen;
				force *= CalcErfcScalar(len, len*len);
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

	__device__ inline Float3 CalcCoulumbForce_optim(const float chargeProduct, const Float3& diff, const float distSq)
	{
		const float invLen = rsqrtf(distSq);                  // Computes 1 / sqrt(lenSquared)
		const float invLenCubed = invLen * invLen * invLen;       // Computes (1 / |diff|^3)

		Float3 force = diff * chargeProduct * invLenCubed;
#ifdef FORCE_NAN_CHECK
		if (force.isNan())
			force.print('E');
#endif
		if constexpr (ENABLE_ERFC_FOR_EWALD) {
			if constexpr (!USE_PRECOMPUTED_ERFCSCALARS) {
				force *= CalcErfcScalar(1.f / invLen, distSq);
			}
			else {
				const int N = DeviceConstants::ERFC_LUT_SIZE;
				const float distanceInArray = fminf(1.f/invLen * DeviceConstants::cutoffNmReciprocal * N - 1, N - 1);
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