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
			return EngineUtils::ForceActivationFunction(force, 100.f);


		return force;	// [1/24 J/mol/nm]
	}



		// Specific to solvent kernel	
	template<bool computePotE, bool emvariant>
	__device__ void ComputeSolventToSolventLJForcesIntrablock(ForceEnergy& fe, 
		const uint8_t& myAtomtype, const Float3& myPosition,
		const ParticleQuickData* const queryParticles,
		const NonbondedInteractionParams& precomputedOO,
		const float* const charges, /*{O, H}*/
		int nParticles, int nParticlesThisBatch, int batchOffset) 
	{	
		if (threadIdx.x >= nParticles)
			return;

		for (int queryIndexRel = 0; queryIndexRel < nParticlesThisBatch; queryIndexRel++) {
			//const bool sameAtom = index == queryIndex; 
			const int queryIndexAbs = batchOffset + queryIndexRel;
			const bool sameMolecule = threadIdx.x / 3 == queryIndexAbs / 3; //Not needed since H_epsilon is 0... // This is only valid for water-like molecules with 3 atoms each
			if (sameMolecule) // TODO: THis is worng. We should still do ES inside molecules, just not LJ
				continue;

			const Float3 diff = Float3(queryParticles[queryIndexRel].relPos) - myPosition;
			const float distSq = diff.lenSquared();

			if (EngineUtils::isOutsideCutoff(distSq)) { continue; }

			if (myAtomtype == 0 && queryParticles[queryIndexRel].atomType == 0) {
				fe.force += calcLJForceOptim<computePotE, emvariant>(diff, 1. / distSq, fe.potE,
					precomputedOO.sigma, precomputedOO.epsilon,
					CalcLJOrigin::SolSolIntra,
					threadIdx.x, queryIndexRel
				) * 24.f;
			}
			if constexpr (ENABLE_ES_SR) {
				const float chargeProduct = charges[myAtomtype] * charges[queryParticles[queryIndexRel].atomType];
				fe.force += PhysicsUtilsDevice::CalcCoulumbForce(chargeProduct, -diff, distSq);
				if constexpr (computePotE)
					fe.potE += PhysicsUtilsDevice::CalcCoulumbPotential(chargeProduct, distSq);
			}
		}
	}

	template<bool computePotE, bool emvariant>
	__device__ void ComputeSolventToSolventLJForcesInterblock(
		ForceEnergy& fe, 
		const uint8_t& myAtomtype, const Float3& myPosition,
		const ParticleQuickData* __restrict__ const queryParticles,
		const NonbondedInteractionParams precomputedOO, const float* const charges, /*{O, H}*/
		int nParticles, int nParticlesQueryThisBatch, float cutoffNmSq)
	{
		if (threadIdx.x >= nParticles)
			return;

		for (int queryIndex = 0; queryIndex < nParticlesQueryThisBatch; queryIndex++) {

			const Float3 diff = Float3(queryParticles[queryIndex].relPos) - myPosition;
			const float distSq = diff.lenSquared();
			if (distSq > cutoffNmSq)
				continue;


			if (myAtomtype == 0 && queryParticles[queryIndex].atomType == 0) {
				fe.force += calcLJForceOptim<computePotE, emvariant>(diff, 1. / distSq, fe.potE,
					precomputedOO.sigma, precomputedOO.epsilon,
					CalcLJOrigin::SolSolInter,
					threadIdx.x, queryIndex
				) * 24.f;
			}

			if constexpr (ENABLE_ES_SR) {
				const float chargeProduct = charges[myAtomtype] * charges[queryParticles[queryIndex].atomType];
				fe.force += PhysicsUtilsDevice::CalcCoulumbForce(chargeProduct, -diff, distSq);
				if constexpr (computePotE)
					fe.potE += PhysicsUtilsDevice::CalcCoulumbPotential(chargeProduct, distSq);
			}
		}
	}


	// Returns fe on p0, invert to get fe on p1
	template<bool computePotE, bool emvariant>
	__device__ ForceEnergy ComputeParticleParticleNB(const Float3& pos0, const Float3& pos1, const LJParameters& lj0, const LJParameters& lj1, const float& charge0, const float& charge1, int p0ParticleGlobalId, int p1ParticleGlobalId) 
	{
		ForceEnergy fe{}; // on p0
		
		const Float3 diff = pos1 - pos0;
		
		if (lj0.epsilonSqrt != -1.f && lj1.epsilonSqrt != -1.f) {
			//diff.print('d');
			fe.force = calcLJForceOptim<computePotE, emvariant>(diff, 1. / diff.lenSquared(), fe.potE,
				CalcSigma(lj0.sigmaHalf, lj1.sigmaHalf),
				CalcEpsilon(lj0.epsilonSqrt, lj1.epsilonSqrt),
				//precomputedOO.sigma, precomputedOO.epsilon,
				CalcLJOrigin::PP,
				p0ParticleGlobalId, p1ParticleGlobalId
			) * 24.f;
		}

		if constexpr (ENABLE_ES_SR) {
			if (!isnan(charge0) && !isnan(charge1)) {
				const float chargeProduct = charge0 * charge1;
				fe.force += PhysicsUtilsDevice::CalcCoulumbForce(chargeProduct, -diff);
				if constexpr (computePotE)
					fe.potE += PhysicsUtilsDevice::CalcCoulumbPotential(chargeProduct, diff.lenSquared());
			}
		}


		if constexpr (FORCE_CHECKS) {
			//if (fe.force.isNan() || isnan(fe.potE)) {
			//	printf("PP NB is nan. diff: %f %f %f  sigma: %f %f  eps: %f %f charge: %f %f distance %f\n",
			//		diff.x, diff.y, diff.z,
			//		p0.params.sigmaHalf, p1.params.sigmaHalf,
			//		p0.params.epsilonSqrt, p1.params.epsilonSqrt,
			//		p0.params.charge, p1.params.charge,
			//		diff.len());
			//}
		}

		return fe;
	}
}
