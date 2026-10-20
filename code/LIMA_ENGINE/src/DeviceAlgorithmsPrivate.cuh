#pragma once

#include "DeviceAlgorithms.cuh"
#include "PhysicsUtils.cuh"
#include "platform/DeviceIntrinsics.cuh"

#include <algorithm>

namespace LAL {
	//static_assert(USE_PRECOMPUTED_BSPLINES == false, "Precomputed B-spline LUT is disabled. Please set USE_PRECOMPUTED_BSPLINES to false.");
	__device__ inline void CalcBspline(float f, float* w) {
		/*
		// Disabled precomputed B-spline LUT. Keep this implementation for potential future use.
		if constexpr (USE_PRECOMPUTED_BSPLINES) {
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
		*/		

		w[0] = (1.f - f) * (1.f - f) * (1.f - f) / 6.f;
		w[1] = (4.f - 6.f * f * f + 3.f * f * f * f) / 6.f;
		w[2] = (1.f + 3.f * f + 3.f * f * f - 3.f * f * f * f) / 6.f;
		w[3] = (f * f * f) / 6.f;
	}

	// Derivatives of CalcBspline's weights with respect to f
	__device__ inline void CalcBsplineDerivative(float f, float* dw) {
		dw[0] = -0.5f * (1.f - f) * (1.f - f);
		dw[1] = -2.f * f + 1.5f * f * f;
		dw[2] = 0.5f + f - 1.5f * f * f;
		dw[3] = 0.5f * f * f;
	}

	template<int order>
	__device__ __forceinline__
		float EvalPoly(const float x, const std::array<float, order>& coeffs)
	{
		static_assert(order > 0);
		float acc = coeffs[order - 1];
#pragma unroll
		for (int i = order - 2; i >= 0; --i)
			acc = acc * x + coeffs[i];
		return acc;
	}
}


// Functions optimized for CUDA
namespace PhysicsUtilsDevice {
	using PhysicsUtils::modifiedCoulombConstant;


// If we're gonna use this we need to deal with the license. But we dont need 0-4, so rather remake it for out purposes..
	__forceinline__ __device__ float raw_ex2(float a)
	{
		return DeviceIntrinsics::Exp2Approx(a);
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

	// Ewald real space force scalar with x = kappa*dist: erfc(x) + 2x/sqrt(pi) * exp(-x^2) = exp(-x^2) * R(x), R(x) = erfcx(x) + 2x/sqrt(pi),
	// with erfcx(x) = exp(x^2) * erfc(x). So a single ex2 serves both terms, and erfc for the potential is exp(-x^2) * (R(x) - 2x/sqrt(pi)).
	// erfcx is smooth, so a polynomial fits R on x in [0, 3.05] (kappa*cutoff = 3, see CalcEwaldkappa). Beyond that it is wrong,
	// so the cutoff must mask those pairs. R(0) = 1 exactly, so kappa = 0 (no PME) gives plain Coulomb.
	__device__ inline float EwaldGaussian(float distSq, float ewaldKappa) {
		constexpr float log2e = 1.4426950408889634f;
		return raw_ex2(distSq * (-ewaldKappa * ewaldKappa * log2e)); // The product of constants is hoisted out of the pair loops
	}
	// R(x). Minimax (Lawson-iterated Chebyshev) fit of (R(x)-1)/x, weighted for the relative error of R.
	// Max relative error of the force scalar 5.4e-6 (fasterfc + exp: 2.9e-6)
	__device__ inline float EwaldForcePoly(float x) {
		float r = -1.642438147e-4f;
		r = fmaf(r, x, 2.629942494e-3f);
		r = fmaf(r, x, -1.860453375e-2f);
		r = fmaf(r, x, 7.769382000e-2f);
		r = fmaf(r, x, -2.180595994e-1f);
		r = fmaf(r, x, 4.498292804e-1f);
		r = fmaf(r, x, -7.343890071e-1f);
		r = fmaf(r, x, 9.968287945e-1f);
		r = fmaf(r, x, 1.945109834e-4f);
		return fmaf(r, x, 1.f);
	}
	// erfc(x) = EwaldGaussian * ErfcxFromForcePoly. Max relative error of erfc 1.0e-4, near the cutoff where erfc ~ 2e-5 (fasterfc: 7.0e-5)
	__device__ inline float ErfcxFromForcePoly(float forcePoly, float x) {
		return fmaf(x, -2.f / PI_sqrt, forcePoly);
	}

	//static_assert(USE_PRECOMPUTED_ERFCSCALARS == false, "Precomputed ERFC scalar LUT is disabled. Please set USE_PRECOMPUTED_ERFCSCALARS to false.");
	__device__ inline float CalcErfcScalar(float dist, float distSq, float ewaldKappa) {
		/*
		// Disabled precomputed ERFC scalar LUT. Keep this implementation for potential future use.
		if constexpr (USE_PRECOMPUTED_ERFCSCALARS) {
			const int N = DeviceConstants::ERFC_LUT_SIZE;
			const float distanceInArray = fminf(dist * DeviceConstants::cutoffNmReciprocal * N - 1, N - 1);
			const int index = static_cast<int>(std::floor(distanceInArray));
			const int indexNext = std::min(index + 1, N - 1);
			const float frac = distanceInArray - static_cast<float>(index);
			const float scalar = LAL::lerp(DeviceConstants::erfcForcescalarTable[index], DeviceConstants::erfcForcescalarTable[indexNext], frac);// optim: look into using std::lerp
			
			return scalar;
			
		}
		*/
		return EwaldGaussian(distSq, ewaldKappa) * EwaldForcePoly(dist * ewaldKappa);
	}

	__device__ inline Float3 CalcCoulumbForceTrueImplementation(const float chargeProduct, const Float3& diff, const float distSq, const float ewaldKappa)
	{
		const float invLen = rsqrtf(distSq);                  // Computes 1 / sqrt(lenSquared)
		const float invLenCubed = invLen * invLen * invLen;       // Computes (1 / |diff|^3)

		Float3 force = diff * chargeProduct * invLenCubed * modifiedCoulombConstant;
		if constexpr (FORCE_CHECKS) {
			if (force.isNan())
				force.print('E');
		}
		if constexpr (ENABLE_ERFC_FOR_EWALD) {
			float len = 1.f / invLen;
			force *= CalcErfcScalar(len, distSq, ewaldKappa);
		}

		return force;
	}


	__device__ inline Float3 CalcCoulumbForce(const float chargeProduct, const Float3& diff, float distSq, const float ewaldKappa) {
		return CalcCoulumbForceTrueImplementation(chargeProduct, diff, diff.lenSquared(), ewaldKappa);
	}

	__device__ inline Float3 CalcCoulumbForce(const float chargeProduct, const Float3& diff, const float ewaldKappa) {
		return CalcCoulumbForce(chargeProduct, diff, diff.lenSquared(), ewaldKappa);
	}



	// <summary>Calculate the potential without multiplying the coulumbConstant, so called must do that!!</summary>
	// <param name="myCharge">[kilo C/mol]</param>
	// <param name="otherCharge">[kilo C/mol]</param>
	// <param name="diff">[nm]</param>
	// <returns>[J/mol   /   modifiedCoulombConstant ]</returns>
	//constexpr float modifiedCoulombConstant = 1.f;
	__device__ inline float CalcCoulumbPotentialTrueImplementation(const float chargeProduct, const float distSq, const float ewaldKappa)
	{
		float dist = sqrtf(distSq);
		float potential = (chargeProduct) * 1.f/dist;
		if constexpr (ENABLE_ERFC_FOR_EWALD) {
			/*
			// Disabled precomputed ERFC potential LUT. Keep this implementation for potential future use.
			if constexpr (USE_PRECOMPUTED_ERFCSCALARS) {
				const int N = DeviceConstants::ERFC_LUT_SIZE;
				const float distanceInArray = fminf(dist * DeviceConstants::cutoffNmReciprocal * N - 1, N - 1);
				const int index = static_cast<int>(std::floor(distanceInArray));
				const int indexNext = std::min(index + 1, N - 1);
				const float frac = distanceInArray - static_cast<float>(index);
				const float scalar = LAL::lerp(DeviceConstants::erfcPotentialscalarTable[index], DeviceConstants::erfcForcescalarTable[indexNext], frac);// optim: look into using std::lerp
				potential *= scalar;
			}
			else
			*/
			{
				potential *= erfc(dist * ewaldKappa);
			}
		}

		return potential;
	}

	// Returns [J/mol], with 0.5 to account for splitting the potential between the 2 particles
	__device__ inline float CalcCoulumbPotential(const float chargeProduct, const float distSq, const float ewaldKappa)
	{
		return CalcCoulumbPotentialTrueImplementation(chargeProduct, distSq, ewaldKappa) * modifiedCoulombConstant * 0.5f;
	}

}
