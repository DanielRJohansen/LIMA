#pragma once

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "EngineUtils.cuh"
#include "PhysicsUtilsDevice.cuh"
#include "DeviceAlgorithmsPrivate.cuh"

#include <cfloat>

namespace LJ {
	enum CalcLJOrigin { ComComIntra, ComComInter, ComSol, SolCom, SolSolIntra, SolSolInter, Pairbond, PP };


	__device__ static const char* calcLJOriginString[] = {
		"ComComIntra", "ComComInter", "ComSol", "SolCom", "SolSolIntra", "SolSolInter", "Pairbond", "PP"
	};


	__device__ void calcLJForceOptimLogErrors(Float3 diff, float sigma, float epsilon, float s, int emvariant, float forceScalar, int originSelect, int pid0, int pid1, const char* const* originStrings) {
		printf(
			"LJ: "
			"diff: %8.3f %8.3f %8.3f  "
			"dist %8.6f  "
			"sigma: %4.2f  "
			"eps: %11.6f  "
			"s %9.6f  "
			"emvariant %2d  "
			"forceMagnitude %14.6f  "
			"origin %-12s  "
			"pIds: %2d %2d\n",
			diff.x, diff.y, diff.z,
			diff.len(),
			sigma,
			epsilon,
			s,
			emvariant,
			forceScalar * diff.len() * 24.f,
			originStrings[originSelect],
			pid0, pid1
		);
	}


	/// <summary></summary>
	/// <param name="diff">other minus self, attractive direction [nm]</param>
	/// <param name="dist_sq_reciprocal"></param>
	/// <param name="potE">[J/mol]</param>
	/// <param name="sigma">[nm]</param>
	/// <param name="epsilon">[J/mol]</param>
	/// <returns>Force [1/24 J/mol/nm] on p0. Caller must multiply with scalar 24. to get correct result</returns>
	template<bool computePotE, bool emvariant>
	__device__ inline Float3 calcLJForceOptim(const Float3& diff, const float dist_sq_reciprocal, float& potE, const float sigma, const float epsilon,
		CalcLJOrigin originSelect, /*For debug only*/
		int pid0 = -1, int pid1 = -1) {

		//return Float3{ diff.x > 0 ? 1.f/24.f : -1.f/24.f, 0.f, 0.f};


		//if (!(min(pid0, pid1) == 4 && max(pid0, pid1) == 83))
		//	return {};

		if constexpr (!ENABLE_LJ) {
			return {};
		}

		// Directly from book
		float s = (sigma * sigma) * dist_sq_reciprocal;								// [nm^2]/[nm^2] -> unitless	// OPTIM: Only calculate sigma_squared, since we never use just sigma
		s = s * s * s;
		float force_scalar = epsilon * s * dist_sq_reciprocal * (1.f - 2.f * s);	// Attractive when positive		[(kg*nm^2)/(nm^2*ns^2*mol)] ->----------------------	[(kg)/(ns^2*mol)]	

		if constexpr (emvariant)
			force_scalar = fmaxf(fminf(force_scalar, 1e+20), -1e+20); // Necessary to avoid inf * 0 = NaN

		const Float3 force = diff * force_scalar;

		if constexpr (FORCE_CHECKS) {
			/*int debugPid = 2251;
			bool eitherIsInvalid = pid0 == -1 || pid1 == -1;*/
			bool debugThis = false;// (pid0 == debugPid || pid1 == debugPid) && !eitherIsInvalid;


			if (force.isNan() || debugThis) {
					calcLJForceOptimLogErrors(diff, sigma, epsilon, s, emvariant, force_scalar, originSelect, pid0, pid1, calcLJOriginString);
				/*printf("LJ is nan. diff: %f %f %f dist %f sigma: %f eps: %f s %f emvariant %d forceScalar %f origin %s pIds: %d %d\n",
					diff.x, diff.y, diff.z, diff.len(), sigma, epsilon, s, emvariant, force_scalar, calcLJOriginString[(int)originSelect], pid0, pid1);*/
			}
		}

		if constexpr (computePotE && ENABLE_POTE) {
			potE += 4.f * epsilon * s * (s - 1.f) * 0.5f;	// 0.5 to account for splitting the potential between the 2 particles
		}

		if constexpr (emvariant)
			return EngineUtils::ForceActivationFunction(-1, force, 100.f);


		return force;	// [1/24 J/mol/nm]
	}

	// The EM pair interaction, on unscaled particle parameters. diff is from p0 to p1, the returned fe is on p0, invert to get fe on p1.
	// Separate from ComputePairNB, since the force activation function is applied to the LJ force alone.
	// The caller must skip masked pairs (padding particles may overlap other particles), rather than discard the result
	template<bool computePotE>
	__device__ inline ForceEnergy ComputePairNBEm(const Float3& diff, float sigma, float epsilon, float chargeProduct, float ewaldKappa)
	{
		ForceEnergy fe{}; // on p0

		fe.force = calcLJForceOptim<computePotE, true>(diff, 1.f / diff.lenSquared(), fe.potE, sigma, epsilon, CalcLJOrigin::PP) * 24.f;

		if constexpr (ENABLE_ES_SR) {
			if (chargeProduct != 0.f) {
				fe.force += PhysicsUtilsDevice::CalcCoulumbForce(chargeProduct, -diff, ewaldKappa);
				if constexpr (computePotE)
					fe.potE += PhysicsUtilsDevice::CalcCoulumbPotential(chargeProduct, diff.lenSquared(), ewaldKappa);
			}
		}

		if constexpr (FORCE_CHECKS) {
			if (fe.force.isNan() || isnan(fe.potE))
				printf("PP NB EM is nan. diff: %f %f %f  sigma: %f  eps: %f  charge product: %f  distance %f\n",
					diff.x, diff.y, diff.z, sigma, epsilon, chargeProduct, diff.len());
		}

		return fe;
	}

	// ComputePairNB expects both particles' epsilonSqrt and charge to be pre-scaled with these,
	// so the pair product directly yields 24*epsilon and modifiedCoulombConstant*chargeProduct
	constexpr float ljEpsilonSqrtScale = 4.898979485566356f; // sqrt(24)
	__device__ inline float CoulombChargeScale() { return sqrtf(PhysicsUtilsDevice::modifiedCoulombConstant); }
	__device__ inline void PrescaleNBParams(float& epsilonSqrt, float& charge) {
		epsilonSqrt *= ljEpsilonSqrtScale;
		charge *= CoulombChargeScale();
	}

	// The MD pair interaction. diff is from p0 to p1, the returned fe is on p0, invert to get fe on p1.
	// epsilonTimes24 and chargeProductTimesK are products of particle parameters pre-scaled with PrescaleNBParams
	// Branchless: LJ and coulomb share a single rsqrt, and are combined into a single scalar applied to diff.
	// Masked pairs are discarded with a select, since padding particles may overlap other particles and produce NaN.
	template<bool computePotE>
	__device__ inline ForceEnergy ComputePairNB(const Float3& diff, float sigma, float epsilonTimes24, float chargeProductTimesK, bool masked, float ewaldKappa)
	{
		ForceEnergy fe{}; // on p0

		const float distSq = diff.lenSquared();
		const float distInv = rsqrtf(distSq);
		const float distSqInv = 1.f / distSq; // Not distInv^2, the steep LJ terms amplify the rsqrt error into noticeable energy drift

		float forceScalar = 0.f; // force = diff * forceScalar, attractive when positive

		if constexpr (ENABLE_LJ) {
			float s = (sigma * sigma) * distSqInv;
			s = s * s * s;
			forceScalar += epsilonTimes24 * s * distSqInv * (1.f - 2.f * s);	// [J/mol/nm^2]

			if constexpr (computePotE && ENABLE_POTE)
				fe.potE += epsilonTimes24 * (1.f / 12.f) * s * (s - 1.f);	// 4*eps*s*(s-1) * 0.5 to account for splitting the potential between the 2 particles
		}

		if constexpr (ENABLE_ES_SR) {
			// Repulsive for equal charges, hence the negation
			float coulombScalar = -chargeProductTimesK * distInv * distSqInv;
			if constexpr (ENABLE_ERFC_FOR_EWALD && !ERFC_USE_CHEBYSHEV_APPROXIMATION) {
				// As CalcErfcScalar, but keeping erfc for the potential
				const float dist = distSq * distInv;
				const float erfcTerm = PhysicsUtilsDevice::fasterfc(dist * ewaldKappa);
				coulombScalar *= erfcTerm + 2.f * ewaldKappa / PI_sqrt * dist * exp(-ewaldKappa * ewaldKappa * distSq);
				if constexpr (computePotE)
					fe.potE += chargeProductTimesK * distInv * erfcTerm * 0.5f; // 0.5 to account for splitting the potential between the 2 particles
			}
			else {
				if constexpr (ENABLE_ERFC_FOR_EWALD)
					coulombScalar *= PhysicsUtilsDevice::CalcErfcScalar(distSq * distInv, distSq, ewaldKappa);
				if constexpr (computePotE)
					fe.potE += PhysicsUtilsDevice::CalcCoulumbPotentialTrueImplementation(chargeProductTimesK, distSq, ewaldKappa) * 0.5f; // 0.5 to account for splitting the potential between the 2 particles
			}
			forceScalar += coulombScalar;
		}

		fe.force = diff * (masked ? 0.f : forceScalar);
		fe.potE = masked ? 0.f : fe.potE;

		if constexpr (FORCE_CHECKS) {
			if (fe.force.isNan() || isnan(fe.potE))
				printf("PP NB is nan. diff: %f %f %f  sigma: %f  eps*24: %f  charge product: %f  distance %f\n",
					diff.x, diff.y, diff.z, sigma, epsilonTimes24, chargeProductTimesK, diff.len());
		}
		return fe;
	}
}
