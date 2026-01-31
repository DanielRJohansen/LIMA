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
			"forceScalar %14.6f  "
			"origin %-12s  "
			"pIds: %2d %2d\n",
			diff.x, diff.y, diff.z,
			diff.len(),
			sigma,
			epsilon,
			s,
			emvariant,
			forceScalar,
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
			if (force.isNan() || false) {
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



	// For intraCompound or bonded-to compounds	
	template<bool computePotE, bool emvariant>
	__device__ Float3 computeCompoundCompoundLJForces(const Float3& self_pos, uint8_t atomtype_self, float& potE_sum,
		const Float3* const neighbor_positions, int neighbor_n_particles, const uint8_t* const atom_types,
		const BondedParticlesLUT* const bonded_particles_lut, CalcLJOrigin ljorigin, const ForceField_NB& forcefield,
		float chargeSelf, const float* const charges,
		const uint32_t* globalParticleIds
	)
	{
		Float3 force(0.f);
		Float3 electrostaticForce{};
		float electrostaticPotential{};

		for (int neighborparticle_id = 0; neighborparticle_id < neighbor_n_particles; neighborparticle_id++) {

			// If thread's assoc. particle is bonded to the particle in neighborcompound, continue
			if (bonded_particles_lut->get(threadIdx.x, neighborparticle_id)) { continue; }

			const int neighborparticle_atomtype = atom_types[neighborparticle_id];

			const Float3 diff = (neighbor_positions[neighborparticle_id] - self_pos);
			const float dist_sq_reciprocal = 1.f / diff.lenSquared();

			std::array<int, 2> globalPids = globalParticleIds == nullptr ? std::array<int, 2>{-1, -1} : std::array<int, 2>{(int)globalParticleIds[threadIdx.x], (int)globalParticleIds[neighborparticle_id]};

			force += calcLJForceOptim<computePotE, emvariant>(diff, dist_sq_reciprocal, potE_sum,
				calcSigma(atomtype_self, neighborparticle_atomtype, forcefield), calcEpsilon(atomtype_self, neighborparticle_atomtype, forcefield),
				ljorigin,
				globalPids[0], globalPids[1]
			);

			if constexpr (ENABLE_ES_SR) {
				electrostaticForce += PhysicsUtilsDevice::CalcCoulumbForce(chargeSelf * charges[neighborparticle_id], -diff);
				if constexpr (computePotE)
					electrostaticPotential += PhysicsUtilsDevice::CalcCoulumbPotential(chargeSelf * charges[neighborparticle_id], diff.lenSquared());
			}
		}

		potE_sum += electrostaticPotential;
		return force * 24.f + electrostaticForce;
	}

	// For non bonded-to compounds
	template<bool computePotE, bool emvariant>
    __device__ inline Float3 computeCompoundCompoundLJForces(const Float3& self_pos, float& potE_sum,
        const Float3* const neighbor_positions, const int neighbor_n_particles,
        const float chargeSelf, const float* const chargeNeighbors,
        const ForceField_NB::ParticleParameters& myParams, const ForceField_NB::ParticleParameters* const neighborParams)
	{
		Float3 force(0.f);
		Float3 electrostaticForce{};
		float electrostaticPotential{};
        const float cutoff_recip = DeviceConstants::cutoffNmSquaredReciprocal;

		for (int neighborparticle_id = 0; neighborparticle_id < neighbor_n_particles; neighborparticle_id++) {
			
            const Float3 diff = (neighbor_positions[neighborparticle_id] - self_pos);
            const float dist_sq_reciprocal = 1.f / diff.lenSquared();
            if (!EngineUtils::isOutsideCutoff_recip(dist_sq_reciprocal, cutoff_recip)) {
				force += calcLJForceOptim<computePotE, emvariant>(diff, dist_sq_reciprocal, potE_sum,
                    myParams.sigmaHalf + neighborParams[neighborparticle_id].sigmaHalf,
                    myParams.epsilonSqrt * neighborParams[neighborparticle_id].epsilonSqrt,
					CalcLJOrigin::ComComInter
				);

				//printf("OLD sigma %f %f eps %f %f charge %f %f dist %f\n", myParams.sigmaHalf, neighborParams[neighborparticle_id].sigmaHalf, myParams.epsilonSqrt, neighborParams[neighborparticle_id].epsilonSqrt, chargeSelf, chargeNeighbors[neighborparticle_id], diff.len());

				if constexpr (ENABLE_ES_SR) {
					electrostaticForce += PhysicsUtilsDevice::CalcCoulumbForce(chargeSelf * chargeNeighbors[neighborparticle_id], -diff);
					if constexpr (computePotE && ENABLE_POTE)
						electrostaticPotential += PhysicsUtilsDevice::CalcCoulumbPotential(chargeSelf * chargeNeighbors[neighborparticle_id], diff.lenSquared());
				}
			}
		}		

		potE_sum += electrostaticPotential;
		return force * 24.f + electrostaticForce;
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


	template<bool computePotE, bool emvariant>
	__device__ Float3 computeSolventToCompoundLJForces(const Float3& self_pos, float myCharge, const int n_particles, const Float3* const positions, float& potE_sum, const uint8_t atomtype_self,
		const ForceField_NB& forcefield, const ForcefieldTinymol& forcefieldTinymol_shared, const uint8_t* const tinymolTypeIds) {	// Specific to solvent kernel
		Float3 force{};
		Float3 electrostaticForce{};
		float electrostaticPotential{};

		for (int i = 0; i < n_particles; i++) {

			const Float3 diff = (positions[i] - self_pos);
			const float dist_sq_reciprocal = 1.f / diff.lenSquared();
			if (EngineUtils::isOutsideCutoff_recip(dist_sq_reciprocal)) { continue; }



			force += calcLJForceOptim<computePotE, emvariant>(diff, dist_sq_reciprocal, potE_sum,
				CalcSigma(forcefield.particle_parameters[atomtype_self].sigmaHalf, forcefieldTinymol_shared.types[tinymolTypeIds[i]].sigmaHalf),
				CalcEpsilon(forcefield.particle_parameters[atomtype_self].epsilonSqrt, forcefieldTinymol_shared.types[tinymolTypeIds[i]].epsilonSqrt),
				CalcLJOrigin::SolCom,
				atomtype_self, -1
			);

			if constexpr (ENABLE_ES_SR) {
				const float chargeProduct = myCharge * forcefieldTinymol_shared.types[tinymolTypeIds[i]].charge;
				electrostaticForce += PhysicsUtilsDevice::CalcCoulumbForce(chargeProduct, -diff);
				if constexpr (computePotE)
					electrostaticPotential += PhysicsUtilsDevice::CalcCoulumbPotential(chargeProduct, diff.lenSquared());
			}
		}

		potE_sum += electrostaticPotential;
		return force * 24.f + electrostaticForce;
	}
	
	template<bool computePotE, bool emvariant>
	__device__ Float3 computeCompoundToSolventLJForces(const Float3& self_pos, const int n_particles, const Float3* const positions,
		float& potE_sum, const uint8_t* atomtypes_others, const int sol_id, const ForcefieldTinymol& forcefieldTinymol_shared, const uint8_t tinymolTypeId,
		const float* const charges)
	{
		Float3 force(0.f);
		Float3 electrostaticForce{};
		float electrostaticPotential{};

		for (int i = 0; i < n_particles; i++) {
			 
			const Float3 diff = (positions[i] - self_pos);
			const float dist_sq_reciprocal = 1.f / diff.lenSquared();
			if (EngineUtils::isOutsideCutoff_recip(dist_sq_reciprocal)) { continue; }

			const auto& otherType = DeviceConstants::forcefield.particle_parameters[atomtypes_others[i]];

			force += calcLJForceOptim<computePotE, emvariant>(diff, dist_sq_reciprocal, potE_sum,
				CalcSigma(forcefieldTinymol_shared.types[tinymolTypeId].sigmaHalf, otherType.sigmaHalf),
				CalcEpsilon(forcefieldTinymol_shared.types[tinymolTypeId].epsilonSqrt, otherType.epsilonSqrt),
				CalcLJOrigin::ComSol,
				sol_id, -1
			);

			if constexpr (ENABLE_ES_SR) {
				const float chargeProduct = forcefieldTinymol_shared.types[tinymolTypeId].charge * charges[i];
				electrostaticForce += PhysicsUtilsDevice::CalcCoulumbForce(chargeProduct, -diff);
				if constexpr (computePotE)
					electrostaticPotential += PhysicsUtilsDevice::CalcCoulumbPotential(chargeProduct, diff.lenSquared());
			}
		}

		potE_sum += electrostaticPotential;
		return force * 24.f + electrostaticForce;
	}

	// Returns fe on p0, invert to get fe on p1
	template<bool computePotE, bool emvariant>
	__device__ ForceEnergy ComputeParticleParticleNB(const PData& p0, const PData& p1, int p0ParticleGlobalId, int p1ParticleGlobalId) 
	{
		ForceEnergy fe{}; // on p0
		
		//const Float3 diff = Float3(queryParticles[queryIndex].relPos) - myPosition;
		const Float3 diff = p1.position - p0.position;
		
		if (p0.params.epsilonSqrt != -1.f && p1.params.epsilonSqrt != -1.f) {
			//diff.print('d');
			fe.force = calcLJForceOptim<computePotE, emvariant>(diff, 1. / diff.lenSquared(), fe.potE,
				CalcSigma(p0.params.sigmaHalf, p1.params.sigmaHalf),
				CalcEpsilon(p0.params.epsilonSqrt, p1.params.epsilonSqrt),
				//precomputedOO.sigma, precomputedOO.epsilon,
				CalcLJOrigin::PP,
				p0ParticleGlobalId, p1ParticleGlobalId
			) * 24.f;


			//if (fe.force.len() > 10000.f) {
			//	printf("p0 %d p1 %d force %f %f %f p0 %f %f %f p1 %f %f %f dist %f sigma %f %f eps %f %f\n", p0ParticleGlobalId, p1ParticleGlobalId,
			//		fe.force.x, fe.force.y, fe.force.z, p0.position.x, p0.position.y, p0.position.z, p1.position.x, p1.position.y, p1.position.z, 
			//		diff.len(), p0.params.sigmaHalf, p1.params.sigmaHalf, p0.params.epsilonSqrt, p1.params.epsilonSqrt);
			//}

			//fe.force.print('f');	
			//printf("\nNEW sigma %f %f eps %f %f charge %f %f dist %f\n", p0.params.sigmaHalf, p1.params.sigmaHalf, p0.params.epsilonSqrt, p1.params.epsilonSqrt, p0.params.charge, p1.params.charge, diff.len());
		}

		if constexpr (ENABLE_ES_SR) {
			if (!isnan(p0.params.charge) && !isnan(p1.params.charge)) {
				const float chargeProduct = p0.params.charge * p1.params.charge;
				fe.force += PhysicsUtilsDevice::CalcCoulumbForce(chargeProduct, -diff);
				if constexpr (computePotE)
					fe.potE += PhysicsUtilsDevice::CalcCoulumbPotential(chargeProduct, diff.lenSquared());
			}
		}


		if constexpr (FORCE_CHECKS) {
			if (fe.force.isNan() || isnan(fe.potE)) {
				printf("PP NB is nan. diff: %f %f %f  sigma: %f %f  eps: %f %f charge: %f %f distance %f\n",
					diff.x, diff.y, diff.z,
					p0.params.sigmaHalf, p1.params.sigmaHalf,
					p0.params.epsilonSqrt, p1.params.epsilonSqrt,
					p0.params.charge, p1.params.charge,
					diff.len());
			}
		}

		return fe;
	}
}
