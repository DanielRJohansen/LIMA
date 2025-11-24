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


	__device__ inline float CalcErfcScalar(float dist, float distSq) {
		if constexpr (USE_PRECOMPUTED_ERFCSCALARS) {
			const int N = DeviceConstants::ERFC_LUT_SIZE;
			const float distanceInArray = fminf(dist * DeviceConstants::cutoffNmReciprocal * N - 1, N - 1);
			const int index = static_cast<int>(std::floor(distanceInArray));
			const int indexNext = std::min(index + 1, N - 1);
			const float frac = distanceInArray - static_cast<float>(index);
			const float scalar = LAL::lerp(DeviceConstants::erfcForcescalarTable[index], DeviceConstants::erfcForcescalarTable[indexNext], frac);// optim: look into using std::lerp
			
			return scalar;
			
		}
		else if constexpr (ERFC_USE_CHEBYSHEV_APPROXIMATION) {
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
			float kappa = DeviceConstants::ewaldKappa;
			float erfcTerm = fasterfc(dist * kappa);
			//const float erfcTerm = erfc(dist * kappa);
			float scalar = erfcTerm + 2.f * kappa / PI_sqrt * dist * exp(-kappa * kappa * distSq);

			return scalar;
		}
	}

	__device__ inline Float3 CalcCoulumbForceTrueImplementation(const float chargeProduct, const Float3& diff, const float distSq)
	{
		const float invLen = rsqrtf(distSq);                  // Computes 1 / sqrt(lenSquared)
		const float invLenCubed = invLen * invLen * invLen;       // Computes (1 / |diff|^3)

		Float3 force = diff * chargeProduct * invLenCubed;
#ifdef FORCE_NAN_CHECK
		if (force.isNan())
			force.print('E');
#endif
		if constexpr (ENABLE_ERFC_FOR_EWALD) {
			float len = 1.f / invLen;
			force *= CalcErfcScalar(len, distSq);
		}

		return force;
	}

	//__device__ inline Float3 CalcCoulumbForceChebyshev(const float chargeProduct, const Float3& diff, const float distSq)
	//{
	//	// ApproximationCutoff
	//	if (distSq < 0.1f || distSq > (1.2f*1.2f)) {
	//		return CalcCoulumbForceTrueImplementation(chargeProduct, diff, distSq);
	//	}


	//	constexpr float a0 = 81.8869829271;
	//	constexpr float a1 = -994.1834195934;
	//	constexpr float a2 = 5251.1528499177;
	//	constexpr float a3 = -15345.7257415252;
	//	constexpr float a4 = 26873.9011166293;
	//	constexpr float a5 = -28831.2495036422;
	//	constexpr float a6 = 18538.9981271726;
	//	constexpr float a7 = -6553.1557863429;
	//	constexpr float a8 = 978.3199500712;
	//	
	//	const float invLenCubedTimesErfcScalarApprox =
	//		fmaf(distSq, fmaf(distSq,
	//			fmaf(distSq, fmaf(distSq, fmaf(distSq, fmaf(distSq, fmaf(distSq, fmaf(distSq, a8, a7),
	//				a6), a5), a4), a3), a2), a1), a0);	
	//	
	//	return diff * chargeProduct * invLenCubedTimesErfcScalarApprox;
	//}


	__device__ inline Float3 CalcCoulumbForceChebyshevPiecewise(
		const float chargeProduct, const Float3& diff, const float distSq) 
	{
		if (distSq < 0.1f || distSq > (1.2f*1.2f)) {
			return CalcCoulumbForceTrueImplementation(chargeProduct, diff, distSq);
		}

		const float domainCutoff = 0.5f;
		static constexpr std::array<float, 9> coeffsNeardomain{
			187.7238724824,
			-4014.1788106412,
			41185.4079517427,
			-251099.5537387842,
			968430.6814281681,
			-2386421.5553119550,
			3643489.5602366733,
			-3139832.7458249126,
			1167337.1327849999,
		};

		static constexpr std::array<float, 9> coeffsFardomain{
			15.3497439236,
			-106.5149330052,
			332.1234300295,
			-601.7258661806,
			687.7123528098,
			-505.0208577000,
			231.8314792008,
			-60.6654824767,
			6.9160030527,
		};

		const auto& a = distSq < domainCutoff ? coeffsNeardomain : coeffsFardomain;
		const float invLenCubedTimesErfcScalarApprox = LAL::EvalPoly<9>(distSq, a);

		return diff * chargeProduct * invLenCubedTimesErfcScalarApprox;
	}


	__device__ inline Float3 CalcCoulumbForce(const float chargeProduct, const Float3& diff, float distSq) {
		if constexpr (COULUMB_USE_CHEBYSHEV_APPROXIMATION)			
			return CalcCoulumbForceChebyshevPiecewise(chargeProduct, diff, distSq);
		else
			return CalcCoulumbForceTrueImplementation(chargeProduct, diff, diff.lenSquared());
	}

	__device__ inline Float3 CalcCoulumbForce(const float chargeProduct, const Float3& diff) {
		return CalcCoulumbForce(chargeProduct, diff, diff.lenSquared());
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