#pragma once

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "EngineUtils.cuh"
#include "PhysicsUtilsDevice.cuh"
#include <cfloat>

namespace LJ {
	// __constant__ mem version
	__device__ inline float calcSigma(uint8_t atomtype1, uint8_t atomtype2) {
		return (forcefield_device.particle_parameters[atomtype1].sigma + forcefield_device.particle_parameters[atomtype2].sigma) * 0.5f;
	}
	// __shared__ mem version
	__device__ inline float calcSigma(uint8_t atomtype1, uint8_t atomtype2, const ForceField_NB& forcefield) {
		return (forcefield.particle_parameters[atomtype1].sigma + forcefield.particle_parameters[atomtype2].sigma) * 0.5f;
	}

	__device__ inline float calcEpsilon(uint8_t atomtype1, uint8_t atomtype2) {
		return sqrtf(forcefield_device.particle_parameters[atomtype1].epsilon * forcefield_device.particle_parameters[atomtype2].epsilon);
	}
	__device__ inline float calcEpsilon(uint8_t atomtype1, uint8_t atomtype2, const ForceField_NB& forcefield) {
		//return 1.4f;
		return __fsqrt_rn(forcefield.particle_parameters[atomtype1].epsilon * forcefield.particle_parameters[atomtype2].epsilon);
	}

	__device__ inline float CalcSigmaTinymol(uint8_t tinymolType1, uint8_t tinymolType2, const ForcefieldTinymol& forcefield) {
		return (forcefield.types[tinymolType1].sigma + forcefield.types[tinymolType2].sigma) * 0.5f;
	}
	__device__ inline float CalcEpsilonTinymol(uint8_t tinymolType1, uint8_t tinymolType2, const ForcefieldTinymol& forcefield) {
		return sqrtf(forcefield.types[tinymolType1].epsilon * forcefield.types[tinymolType2].epsilon);
	}

	__device__ inline float CalcSigma(float s1, float s2) {
		return (s1 + s2) * 0.5f;
	}
	__device__ inline float CalcEpsilon(float e1, float e2) {
		return sqrtf(e1 * e2);
	}

	enum CalcLJOrigin { ComComIntra, ComComInter, ComSol, SolCom, SolSolIntra, SolSolInter, Pairbond };


	__device__ static const char* calcLJOriginString[] = {
		"ComComIntra", "ComComInter", "ComSol", "SolCom", "SolSolIntra", "SolSolInter"
	};


	__device__ void calcLJForceOptimLogErrors(float s, float epsilon, Float3 force, CalcLJOrigin originSelect, float distNM, Float3 diff, float force_scalar, float sigma, int type1, int type2) {
		//auto pot = 4. * epsilon * s * (s - 1.f) * 0.5;
		//if (force.len() > 1.f || pot > 1e+8) {
		if (distNM < 0.05f) {
			//printf("\nBlock %d thread %d\n", blockIdx.x, threadIdx.x);
			////((*pos1 - *pos0) * force_scalar).print('f');
			//pos0.print('0');
			//pos1.print('1');
			printf("\nLJ Force %s: dist nm %f force %f sigma %f epsilon %f t1 %d t2 %d\n",
				calcLJOriginString[(int)originSelect], distNM, (diff * force_scalar).len(), sigma, epsilon, type1, type2);
		}
	}

#ifdef ENABLE_LJ

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
		int type1 = -1, int type2 = -1) {


		// Directly from book
		float s = (sigma * sigma) * dist_sq_reciprocal;								// [nm^2]/[nm^2] -> unitless	// OPTIM: Only calculate sigma_squared, since we never use just sigma
		s = s * s * s;
		float force_scalar = epsilon * s * dist_sq_reciprocal * (1.f - 2.f * s);	// Attractive when positive		[(kg*nm^2)/(nm^2*ns^2*mol)] ->----------------------	[(kg)/(ns^2*mol)]	

		if constexpr (emvariant)
			force_scalar = fmaxf(fminf(force_scalar, 1e+20), -1e+20); // Necessary to avoid inf * 0 = NaN

		const Float3 force = diff * force_scalar;
#ifdef FORCE_NAN_CHECK
		if (force.isNan()) {
			printf("LJ is nan. diff: %f %f %f  sigma: %f  eps: %f s %f distSqRecip %f emvariant %d forceScalar %f firstPart %f\n",
				   diff.x, diff.y, diff.z, sigma, epsilon, s, dist_sq_reciprocal, emvariant, force_scalar, epsilon * s * dist_sq_reciprocal);
		}
#endif

		if constexpr (computePotE && ENABLE_POTE) {
			potE += 4.f * epsilon * s * (s - 1.f) * 0.5f;	// 0.5 to account for splitting the potential between the 2 particles
		}

		if constexpr (emvariant)
			return EngineUtils::ForceActivationFunction(force, 100.f);
		
#if defined LIMASAFEMODE
		calcLJForceOptimLogErrors(s, epsilon, force, originSelect, diff.len(), diff, force_scalar, sigma, type1, type2);
#endif

		return force;	// [1/24 J/mol/nm]
	}
#else 
	__device__ static Float3 calcLJForceOptim(const Float3&, float, float& , float , float , CalcLJOrigin, int, int) {
		return Float3{};
	}
#endif






	// For intraCompound or bonded-to compounds	
	template<bool computePotE, bool emvariant>
	__device__ Float3 computeCompoundCompoundLJForces(const Float3& self_pos, uint8_t atomtype_self, float& potE_sum,
		const Float3* const neighbor_positions, int neighbor_n_particles, const uint8_t* const atom_types,
		const BondedParticlesLUT* const bonded_particles_lut, CalcLJOrigin ljorigin, const ForceField_NB& forcefield, 
        float chargeSelf, const float* const charges)
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

			//float a = calcEpsilon(atomtype_self, neighborparticle_atomtype, forcefield);
			//float b = __half2float(nonbondedInteractionParams_device[atomtype_self * ForceField_NB::MAX_TYPES + neighborparticle_atomtype].epsilon);
			//if (std::abs(a - b) / a > 1e-5) {
			//	printf("Epsilon mismatch %f %f\n", a, b);
			//}

			force += calcLJForceOptim<computePotE, emvariant>(diff, dist_sq_reciprocal, potE_sum,
				calcSigma(atomtype_self, neighborparticle_atomtype, forcefield), calcEpsilon(atomtype_self, neighborparticle_atomtype, forcefield),
				ljorigin,
				threadIdx.x, neighborparticle_id
			);

			if constexpr (ENABLE_ES_SR) {
				electrostaticForce += PhysicsUtilsDevice::CalcCoulumbForce_optim(chargeSelf* charges[neighborparticle_id], -diff);
				electrostaticPotential += PhysicsUtilsDevice::CalcCoulumbPotential_optim(chargeSelf, charges[neighborparticle_id], diff);
			}
		}

		potE_sum += electrostaticPotential * PhysicsUtilsDevice::modifiedCoulombConstant_Potential * 0.5f;
		return force * 24.f + electrostaticForce * PhysicsUtilsDevice::modifiedCoulombConstant_Force;
	}

	// For non bonded-to compounds
	template<bool computePotE, bool emvariant>
	__device__ inline Float3 computeCompoundCompoundLJForces(const Float3& self_pos, const uint8_t atomtype_self, float& potE_sum,
		const Float3* const neighbor_positions, const int neighbor_n_particles, const uint8_t* const atom_types, const ForceField_NB& forcefield, 
        const float chargeSelf, const float* const chargeNeighbors,
        const ForceField_NB::ParticleParameters& myParams, const ForceField_NB::ParticleParameters* const neighborParams, int& util)
	{
		Float3 force(0.f);
		Float3 electrostaticForce{};
		float electrostaticPotential{};
		//int hits = 0;

		// TODO: i dont have any unittests that test whether this works as i expect. Would require a compound with 2 atoms not in the same gridnode
		const Float3 selfRelOffset{
			fabsf(self_pos.x) > static_cast<float>(BoxGrid::blocksizeNM) / 2.f ? copysignf(static_cast<float>(BoxGrid::blocksizeNM), self_pos.x) : 0.f,
			fabsf(self_pos.y) > static_cast<float>(BoxGrid::blocksizeNM) / 2.f ? copysignf(static_cast<float>(BoxGrid::blocksizeNM), self_pos.y) : 0.f,
			fabsf(self_pos.z) > static_cast<float>(BoxGrid::blocksizeNM) / 2.f ? copysignf(static_cast<float>(BoxGrid::blocksizeNM), self_pos.z) : 0.f
		};

		for (int neighborparticle_id = 0; neighborparticle_id < neighbor_n_particles; neighborparticle_id++) {
			
            const Float3 diff = (neighbor_positions[neighborparticle_id] - self_pos);
            const float dist_sq_reciprocal = 1.f / diff.lenSquared();
			if (!EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) {
				//hits++;
                //const NonbondedInteractionParams params = nonbondedInteractionParams_device[static_cast<int>(atomtype_self) * ForceField_NB::MAX_TYPES + neighborparticle_atomtype];
				force += calcLJForceOptim<computePotE, emvariant>(diff, dist_sq_reciprocal, potE_sum,
                    (myParams.sigma + neighborParams[neighborparticle_id].sigma) * 0.5f,
                    __fsqrt_rn(myParams.epsilon * neighborParams[neighborparticle_id].epsilon),
                    /*calcSigma(atomtype_self, neighborparticle_atomtype, forcefield),
                    calcEpsilon(atomtype_self, neighborparticle_atomtype, forcefield),*/
                    //nonbondedInteractionParams_device[static_cast<int>(atomtype_self) * ForceField_NB::MAX_TYPES + neighborparticle_atomtype].sigma,
                    //nonbondedInteractionParams_device[static_cast<int>(atomtype_self) * ForceField_NB::MAX_TYPES + neighborparticle_atomtype].epsilon,
                    //params.sigma, params.epsilon,
					//calcSigma(atomtype_self, neighborparticle_atomtype), calcEpsilon(atomtype_self, neighborparticle_atomtype),
					CalcLJOrigin::ComComInter
					//global_id_self, neighbor_compound->particle_global_ids[neighborparticle_id]
				);

				electrostaticForce += PhysicsUtilsDevice::CalcCoulumbForce_optim(chargeSelf * chargeNeighbors[neighborparticle_id], -diff);
				if constexpr (computePotE && ENABLE_POTE)
					electrostaticPotential += PhysicsUtilsDevice::CalcCoulumbPotential_optim(chargeSelf, chargeNeighbors[neighborparticle_id], diff);
			}

			if constexpr (ENABLE_ES_SR) {
				if ((neighbor_positions[neighborparticle_id] - selfRelOffset).LargestMagnitudeElement() < static_cast<float>(BoxGrid::blocksizeNM) * 1.5f) // TODO: magic nr
				{
                    //electrostaticForce += PhysicsUtilsDevice::CalcCoulumbForce_optim(chargeSelf, chargeNeighbors[neighborparticle_id], -diff * LIMA_TO_NANO);
     //               electrostaticForce += PhysicsUtilsDevice::CalcCoulumbForce_optim(chargeSelf * chargeNeighbors[neighborparticle_id], -diff);
					//if constexpr (computePotE && ENABLE_POTE)
					//	electrostaticPotential += PhysicsUtilsDevice::CalcCoulumbPotential_optim(chargeSelf, chargeNeighbors[neighborparticle_id], diff);
				}
			}
		}
		//atomicAdd(&util, hits);
		


		potE_sum += electrostaticPotential * PhysicsUtilsDevice::modifiedCoulombConstant_Potential * 0.5f;
		return force * 24.f + electrostaticForce * PhysicsUtilsDevice::modifiedCoulombConstant_Force;
	}

	// Specific to solvent kernel	
	template<bool computePotE, bool emvariant>
	__device__ Float3 computeSolventToSolventLJForces(const Float3& relpos_self, const uint8_t tinymolTypeIdSelf, const Float3* const relpos_others, int n_elements, const bool exclude_own_index, float& potE_sum,
		const ForcefieldTinymol& forcefieldTinymol_shared, const uint8_t* const tinymolTypeIds) {
		Float3 force{};

		for (int i = 0; i < n_elements; i++) {
			// If computing within block, dont compute force against thread's solvent
			if (exclude_own_index && threadIdx.x == i) { continue; }

			const Float3 diff = (relpos_others[i] - relpos_self);
			const float dist_sq_reciprocal = 1.f / diff.lenSquared();
			if (EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) { continue; }

			force += calcLJForceOptim<computePotE, emvariant>(diff, dist_sq_reciprocal, potE_sum,				
				CalcSigmaTinymol(tinymolTypeIdSelf, tinymolTypeIds[i], forcefieldTinymol_shared),
				CalcEpsilonTinymol(tinymolTypeIdSelf, tinymolTypeIds[i], forcefieldTinymol_shared),
				exclude_own_index ? CalcLJOrigin::SolSolIntra : CalcLJOrigin::SolSolInter,
				threadIdx.x, i
			);
		}
		return force * 24.f;
	}

	// False fix the hardcoded template params here
	template<bool computePotE, bool emvariant>
	__device__ Float3 computeSolventToCompoundLJForces(const Float3& self_pos, const int n_particles, const Float3* const positions, float& potE_sum, const uint8_t atomtype_self,
		const ForceField_NB& forcefield, const ForcefieldTinymol& forcefieldTinymol_shared, const uint8_t* const tinymolTypeIds) {	// Specific to solvent kernel
		Float3 force{};

		const float mySigma = forcefield.particle_parameters[atomtype_self].sigma;
		const float myEps = forcefield.particle_parameters[atomtype_self].epsilon;		

		for (int i = 0; i < n_particles; i++) {

			const Float3 diff = (positions[i] - self_pos);
			const float dist_sq_reciprocal = 1.f / diff.lenSquared();
			if (EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) { continue; }



			force += calcLJForceOptim<computePotE, emvariant>(diff, dist_sq_reciprocal, potE_sum,
				CalcSigma(mySigma, forcefieldTinymol_shared.types[tinymolTypeIds[i]].sigma),
				CalcEpsilon(myEps, forcefieldTinymol_shared.types[tinymolTypeIds[i]].epsilon),
				CalcLJOrigin::SolCom,
				atomtype_self, -1
			);
		}
		return force * 24.f;
	}
	
	template<bool computePotE, bool emvariant>
	__device__ Float3 computeCompoundToSolventLJForces(const Float3& self_pos, const int n_particles, const Float3* const positions,
		float& potE_sum, const uint8_t* atomtypes_others, const int sol_id, const ForcefieldTinymol& forcefieldTinymol_shared, const uint8_t tinymolTypeId)
	{
		Float3 force(0.f);

		for (int i = 0; i < n_particles; i++) {
			 
			const Float3 diff = (positions[i] - self_pos);
			const float dist_sq_reciprocal = 1.f / diff.lenSquared();
			if (EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) { continue; }

			force += calcLJForceOptim<computePotE, emvariant>(diff, dist_sq_reciprocal, potE_sum,
				CalcSigma(forcefieldTinymol_shared.types[tinymolTypeId].sigma, forcefield_device.particle_parameters[atomtypes_others[i]].sigma),
				CalcEpsilon(forcefieldTinymol_shared.types[tinymolTypeId].epsilon, forcefield_device.particle_parameters[atomtypes_others[i]].epsilon),
				CalcLJOrigin::ComSol,
				sol_id, -1
			);
		}
		return force * 24.f;
	}
}
