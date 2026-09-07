#pragma once

#include "DeviceAlgorithms.cuh"
#include "PhysicsUtils.cuh"

namespace LAL {
	//static_assert(USE_PRECOMPUTED_BSPLINES == false, "Precomputed B-spline LUT is disabled. Please set USE_PRECOMPUTED_BSPLINES to false.");
	__device__ void CalcBspline(float f, float* w) {
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

	template<int order>
	__device__ __forceinline__
		float EvalPoly(const float x, const std::array<float, order>& coeffs)
	{
		float acc = coeffs[order];
#pragma unroll
		for (int i = order - 1; i >= 0; --i)
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
		

		if constexpr (ERFC_USE_CHEBYSHEV_APPROXIMATION) {
			constexpr float a0 = 1.0021710689;
			constexpr float a1 = -0.1573428973;
			constexpr float a2 = 2.6808707411;
			constexpr float a3 = -29.9956062655;
			constexpr float a4 = 57.2335637829;
			constexpr float a5 = -28.4043932874;
			constexpr float a6 = -22.0595037569;
			constexpr float a7 = 27.4957110932;
			constexpr float a8 = -7.7900877300;

			float x = dist;

			float scalar =  (((((((a8 * x + a7) * x + a6) * x + a5) * x + a4) * x + a3) * x + a2) * x + a1) * x + a0;

			// TODO: WARNING: DANGER: This is being called with dist > 2, figure out where that comes from!!
			//if constexpr (FORCE_CHECKS) {
			//	if (scalar < -0.1f || scalar > 1.1f)
			//		printf("Scalar out of bounds: %f  dist: %f\n", scalar, dist);
			//}

			// handle small numeric errors
			scalar = std::clamp(scalar, 0.f, 1.f);

			return scalar;
		}
		else {
			float erfcTerm = fasterfc(dist * ewaldKappa);
			//const float erfcTerm = erfc(dist * kappa);
			float scalar = erfcTerm + 2.f * ewaldKappa / PI_sqrt * dist * exp(-ewaldKappa * ewaldKappa * distSq);

			return scalar;
		}
	}

	__device__ inline Float3 CalcCoulumbForceTrueImplementation(const float chargeProduct, const Float3& diff, const float distSq, const float ewaldKappa)
	{
		const float invLen = rsqrtf(distSq);                  // Computes 1 / sqrt(lenSquared)
		const float invLenCubed = invLen * invLen * invLen;       // Computes (1 / |diff|^3)

		Float3 force = diff * chargeProduct * invLenCubed * modifiedCoulombConstant;
#ifdef FORCE_NAN_CHECK
		if (force.isNan())
			force.print('E');
#endif
		if constexpr (ENABLE_ERFC_FOR_EWALD) {
			float len = 1.f / invLen;
			force *= CalcErfcScalar(len, distSq, ewaldKappa);
		}

		return force;
	}


	__device__ inline Float3 CalcCoulumbForceChebyshevPiecewise(
		const float chargeProduct, const Float3& diff, const float distSq, const float ewaldKappa)
	{
		if (distSq < 0.1f || distSq > (1.2f*1.2f)) {
			return CalcCoulumbForceTrueImplementation(chargeProduct, diff, distSq, ewaldKappa);
		}

		const float domainCutoff = 0.5f;
		static constexpr std::array<float, 9> coeffsNeardomain{
			2801.6250509480,
			-59908.3311363743,
			614658.4828306455,
			-3747455.1890168283,
			14453034.7755872533,
			-35615387.2337769344,
			54376097.6685821563,
			-46859432.2084176466,
			17421550.6577042267,
		}; // The magnitude of these is quite an issue...

		static constexpr std::array<float, 9> coeffsFardomain{
			229.0823566196,
			-1589.6481393717,
			4956.6701841007,
			-8980.2657392603,
			10263.5436292686,
			-7537.0226890896,
			3459.8949570903,
			-905.3826409340,
			103.2156813554,
		};

		const auto& a = distSq < domainCutoff ? coeffsNeardomain : coeffsFardomain;
		const float invLenCubedTimesErfcScalarApprox = LAL::EvalPoly<9>(distSq, a);

		return diff * chargeProduct * invLenCubedTimesErfcScalarApprox;
	}

	// TODO: Include modified coulumb constant here
	__device__ inline Float3 CalcCoulumbForce(const float chargeProduct, const Float3& diff, float distSq, const float ewaldKappa) {
		if constexpr (COULUMB_USE_CHEBYSHEV_APPROXIMATION)			
			return CalcCoulumbForceChebyshevPiecewise(chargeProduct, diff, distSq, ewaldKappa);
		else
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

	__device__ inline float CalcCoulumbPotentialChebyshevPiecewise(
		const float chargeProduct, const float distSq, const float ewaldKappa)
	{
		if (distSq < 0.1f || distSq >(1.2f * 1.2f)) {
			return CalcCoulumbPotentialTrueImplementation(chargeProduct, distSq, ewaldKappa) * modifiedCoulombConstant * 0.5f;
		}

		// TODO: It's kind of an issue that these scalars are so large. It will lead to some serious imprecision
		// We could, as we're computing these params scale them to be <1, and then rescale the result accordingly.
		const float domainCutoff = 0.5f;
		static constexpr std::array<float, 10> coeffsNeardomain{
			34.0379978922,
			-661.4914286200,
			6882.6123761172,
			-45966.2677339184,
			206170.4678880951,
			-626888.8898421292,
			1274173.0836184307,
			-1656858.2279262797,
			1245661.5552110448,
			-411659.2654445015,
		};

		static constexpr std::array<float, 10> coeffsFardomain{
			7.3918879318,
			-54.9064293050,
			188.6443996713,
			-388.6586874083,
			524.1132887011,
			-476.3149712651,
			290.2059096899,
			-113.8773273801,
			26.0453693742,
			-2.6403995315,
		};

		const auto& a = distSq < domainCutoff ? coeffsNeardomain : coeffsFardomain;
		const float potentialApproximation = LAL::EvalPoly<10>(distSq, a);

		return chargeProduct * potentialApproximation;
	}

	__device__ inline float CalcCoulumbPotential(const float chargeProduct, const float distSq, const float ewaldKappa)
	{
		if constexpr (COULUMB_USE_CHEBYSHEV_APPROXIMATION)
			return CalcCoulumbPotentialChebyshevPiecewise(chargeProduct, distSq, ewaldKappa);
		else
			return CalcCoulumbPotentialTrueImplementation(chargeProduct, distSq, ewaldKappa);
	}

}
