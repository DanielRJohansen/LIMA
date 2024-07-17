#pragma once

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "EngineUtils.cuh"






namespace ShortRangeElectrostatics {

	/// <summary>
	/// 
	/// </summary>
	/// <param name="posA">Relative to it's node</param>
	/// <param name="posB">Relative to A's node</param>
	/// <returns></returns>
	__device__ inline bool withinRange(const Float3& posA, const Float3& posB) {
		const Int3 nodeDiff = LIMAPOSITIONSYSTEM::PositionToNodeIndex(posA) - LIMAPOSITIONSYSTEM::PositionToNodeIndex(posB);

		return nodeDiff.MaxAbsElement() <= 2;
	}

} // namespace ShortRangeElectrostatics

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

	enum CalcLJOrigin { ComComIntra, ComComInter, ComSol, SolCom, SolSolIntra, SolSolInter };


	__device__ static const char* calcLJOriginString[] = {
		"ComComIntra", "ComComInter", "ComSol", "SolCom", "SolSolIntra", "SolSolInter"
	};


	__device__ void calcLJForceOptimLogErrors(float s, float epsilon, Float3 force, CalcLJOrigin originSelect, float dist, Float3 diff, float force_scalar, float sigma, float type1, float type2) {
		auto pot = 4. * epsilon * s * (s - 1.f) * 0.5;
		//if (force.len() > 1.f || pot > 1e+8) {
		const float distNM = dist / NANO_TO_LIMA;
		if (distNM < 0.05f) {
			//printf("\nBlock %d thread %d\n", blockIdx.x, threadIdx.x);
			////((*pos1 - *pos0) * force_scalar).print('f');
			//pos0.print('0');
			//pos1.print('1');
			printf("\nLJ Force %s: dist nm %f force %f sigma %f epsilon %f t1 %d t2 %d\n",
				calcLJOriginString[(int)originSelect], distNM, (diff * force_scalar).len(), sigma / NANO_TO_LIMA, epsilon, type1, type2);
		}
	}

#ifdef ENABLE_LJ

	// This function does not add the 24 scalar, the caller function must do so!
	// Calculates LJ force on p0	(attractive to p1. Negative values = repulsion )//
	// Returns force in J/mol*M??
	__device__ static Float3 calcLJForceOptim(const Float3& diff, const float dist_sq_reciprocal, float& potE, const float sigma /*[nm]*/, const float epsilon /*[(kg*nm^2)/(ns^2*mol)]*/,
		CalcLJOrigin originSelect, /*For debug only*/
		int type1 = -1, int type2 = -1) {


		// Directly from book
		float s = (sigma * sigma) * dist_sq_reciprocal;								// [nm^2]/[nm^2] -> unitless	// OPTIM: Only calculate sigma_squared, since we never use just sigma
		s = s * s * s;
		const float force_scalar = epsilon * s * dist_sq_reciprocal * (1.f - 2.f * s);	// Attractive. Negative, when repulsive		[(kg*nm^2)/(nm^2*ns^2*mol)] ->----------------------	[(kg)/(ns^2*mol)]	

		const Float3 force = diff * force_scalar;

		if constexpr (CALC_POTE) {
			potE += 2.f * epsilon * s * (s - 1.f);	// 2.f instead of 2.f to account for 2 particles doing the same calculation
		}
#if defined LIMASAFEMODE
		calcLJForceOptimLogErrors(s, epsilon, force, originSelect, diff.len(), diff, force_scalar, sigma, type1, type2);
#endif

		return force;	// GN/mol [(kg*nm)/(ns^2*mol)]
	}
#else 
	__device__ static Float3 calcLJForceOptim(const Float3&, float, float& , float , float , CalcLJOrigin, int, int) {
		return Float3{};
	}
#endif






	// For intraCompound or bonded-to compounds
	__device__ Float3 computeCompoundCompoundLJForces(const Float3& self_pos, uint8_t atomtype_self, float& potE_sum,
		const Float3* const neighbor_positions, int neighbor_n_particles, const uint8_t* const atom_types,
		const BondedParticlesLUT* const bonded_particles_lut, CalcLJOrigin ljorigin, const ForceField_NB& forcefield)
	{
		Float3 force(0.f);
		Float3 electrostaticForce{};
		for (int neighborparticle_id = 0; neighborparticle_id < neighbor_n_particles; neighborparticle_id++) {

			// If thread's assoc. particle is bonded to the particle in neighborcompound, continue
			if (bonded_particles_lut->get(threadIdx.x, neighborparticle_id)) { continue; }

			const int neighborparticle_atomtype = atom_types[neighborparticle_id];

			const Float3 diff = (neighbor_positions[neighborparticle_id] - self_pos);
			const float dist_sq_reciprocal = 1.f / diff.lenSquared();
			if (!EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) {	// TODO: no need to do check for intracompound
				force += calcLJForceOptim(diff, dist_sq_reciprocal, potE_sum,
					calcSigma(atomtype_self, neighborparticle_atomtype, forcefield), calcEpsilon(atomtype_self, neighborparticle_atomtype, forcefield),
					ljorigin,
					threadIdx.x, neighborparticle_id
				);
			}

			//electrostaticForce += PhysicsUtils::CalcCoulumbForce()

		}
		return force * 24.f;
	}

	// For non bonded-to compounds
	__device__ Float3 computeCompoundCompoundLJForces(const Float3& self_pos, uint8_t atomtype_self, float& potE_sum,
		const Float3* const neighbor_positions, int neighbor_n_particles, const uint8_t* const atom_types, const ForceField_NB& forcefield, float chargeSelf, const half* const chargeNeighbors)
	{
		Float3 force(0.f);
		Float3 electrostaticForce{};

		for (int neighborparticle_id = 0; neighborparticle_id < neighbor_n_particles; neighborparticle_id++) {
			const int neighborparticle_atomtype = atom_types[neighborparticle_id];

			const Float3 diff = (neighbor_positions[neighborparticle_id] - self_pos);
			const float dist_sq_reciprocal = 1.f / diff.lenSquared();
			if (!EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) {
				force += calcLJForceOptim(diff, dist_sq_reciprocal, potE_sum,
					calcSigma(atomtype_self, neighborparticle_atomtype, forcefield), calcEpsilon(atomtype_self, neighborparticle_atomtype, forcefield),
					CalcLJOrigin::ComComInter
					//global_id_self, neighbor_compound->particle_global_ids[neighborparticle_id]
				);
			}

			if constexpr (ENABLE_ELECTROSTATICS) {
				electrostaticForce += PhysicsUtils::CalcCoulumbForce(chargeSelf, chargeNeighbors[neighborparticle_id], -diff * LIMA_TO_NANO);
				potE_sum += PhysicsUtils::CalcCoulumbPotential(chargeSelf, chargeNeighbors[neighborparticle_id], diff.len() * LIMA_TO_NANO) * 0.5f;
			}
		}
		return force * 24.f + electrostaticForce;
	}


	__device__ Float3 computeSolventToSolventLJForces(const Float3& relpos_self, const Float3* const relpos_others, int n_elements, const bool exclude_own_index, float& potE_sum) {	// Specific to solvent kernel
		Float3 force{};

		for (int i = 0; i < n_elements; i++) {

			// If computing within block, dont compute force against thread's solvent
			if (exclude_own_index && threadIdx.x == i) { continue; }

			const Float3 diff = (relpos_others[i] - relpos_self);
			const float dist_sq_reciprocal = 1.f / diff.lenSquared();
			if (EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) { continue; }

			force += calcLJForceOptim(diff, dist_sq_reciprocal, potE_sum,
				forcefield_device.particle_parameters[ATOMTYPE_SOLVENT].sigma,
				forcefield_device.particle_parameters[ATOMTYPE_SOLVENT].epsilon,
				exclude_own_index ? CalcLJOrigin::SolSolIntra : CalcLJOrigin::SolSolInter,
				threadIdx.x, i
			);
		}
		return force * 24.f;
	}

	__device__ Float3 computeSolventToCompoundLJForces(const Float3& self_pos, const int n_particles, const Float3* const positions, float& potE_sum, const uint8_t atomtype_self,
		const ForceField_NB& forcefield) {	// Specific to solvent kernel
		Float3 force{};
		for (int i = 0; i < n_particles; i++) {

			const Float3 diff = (positions[i] - self_pos);
			const float dist_sq_reciprocal = 1.f / diff.lenSquared();
			if (EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) { continue; }

			force += calcLJForceOptim(diff, dist_sq_reciprocal, potE_sum,
				calcSigma(atomtype_self, ATOMTYPE_SOLVENT, forcefield),
				calcEpsilon(atomtype_self, ATOMTYPE_SOLVENT, forcefield),
				CalcLJOrigin::SolCom,
				atomtype_self, ATOMTYPE_SOLVENT
			);
		}
		return force * 24.f;
	}

	__device__ Float3 computeCompoundToSolventLJForces(const Float3& self_pos, const int n_particles, const Float3* const positions,
		float& potE_sum, const uint8_t* atomtypes_others, const int sol_id)
	{
		Float3 force(0.f);
		for (int i = 0; i < n_particles; i++) {

			const Float3 diff = (positions[i] - self_pos);
			const float dist_sq_reciprocal = 1.f / diff.lenSquared();
			if (EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) { continue; }

			force += calcLJForceOptim(diff, dist_sq_reciprocal, potE_sum,
				calcSigma(ATOMTYPE_SOLVENT, atomtypes_others[i]),
				calcEpsilon(ATOMTYPE_SOLVENT, atomtypes_others[i]),
				CalcLJOrigin::ComSol,
				sol_id, -1
			);
		}
		return force * 24.f;
	}
}