#pragma once

#include <cuda_runtime.h>

#include "LimaTypes.cuh"
#include "Constants.cuh"
#include "Bodies.cuh"



namespace LimaForcecalc {

// ------------------------------------------------------------------------------------------- BONDED FORCES -------------------------------------------------------------------------------------------//

__device__ void calcPairbondForces(Float3* pos_a, Float3* pos_b, PairBond* bondtype, Float3* results, float* potE) {
	// Calculates bond force on both particles					//
	// Calculates forces as J/mol*M								//
	// kb [J/(mol*lm^2)]
	Float3 difference = *pos_a - *pos_b;						//	[lm]
	float error = difference.len() - bondtype->b0;				//	[lm]

	*potE += 0.5f * bondtype->kb * (error * error);				// [J/mol]

	float force_scalar = -bondtype->kb * error;				//	[J/(mol*lm)] = [kg/(mol*s^2)]

	difference = difference.norm();								// dif_unit_vec, but shares variable with dif
	results[0] = difference * force_scalar;						// [kg * lm / (mol*ls^2)] = [lN]
	results[1] = difference * force_scalar * -1;				// [kg * lm / (mol*ls^2)] = [lN]

	if (0.5f * bondtype->kb * (error * error) > 100000 || true) {
		//printf("error %f    pote %f    force %f\n", error, 0.5f * bondtype->kb * (error * error), force_scalar);
		//printf("Dif %f\n",difference.len());
		//printf("len %f\n", (*pos_a - *pos_b).len());
		//results[0].print('r');
		//results[1].print('R');
		//pos_a->print('a');
		//pos_b->print('b');
	}
}


__device__ void calcAnglebondForces(Float3* pos_left, Float3* pos_middle, Float3* pos_right, AngleBond* angletype, Float3* results, float* potE) {
	const Float3 v1 = (*pos_left - *pos_middle).norm();
	const Float3 v2 = (*pos_right - *pos_middle).norm();
	const Float3 normal = v1.cross(v2).norm();	// Poiting towards y, when right is pointing toward x

	const Float3 inward_force_direction1 = (v1.cross(normal * -1.f)).norm();
	const Float3 inward_force_direction2 = (v2.cross(normal)).norm();

	const float angle = Float3::getAngle(v1, v2);					// [rad]
	const float error = angle - angletype->theta_0;				// [rad]

	// Simple implementation
	//*potE += angletype->k_theta * error * error * 0.5f;		// Energy [J/mol]0
	//float torque = angletype->k_theta * (error);				// Torque [J/(mol*rad)]

	// Correct implementation
	*potE += -angletype->k_theta * (cosf(error) - 1.f);		// Energy [J/mol]
	const float torque = angletype->k_theta * sinf(error);	// Torque [J/(mol*rad)]

	results[0] = inward_force_direction1 * (torque / (*pos_left - *pos_middle).len());
	results[2] = inward_force_direction2 * (torque / (*pos_right - *pos_middle).len());
	results[1] = (results[0] + results[2]) * -1.f;


	if (threadIdx.x == 2 && blockIdx.x == 0 || 1) {
		//printf("\nangle %f error %f force %f t0 %f kt %f\n", angle, error, force_scalar, angletype->theta_0, angletype->k_theta);
	}
}
__device__ void calcDihedralbondForces(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	Float3 normal1 = (*pos_left - *pos_lm).cross((*pos_rm - *pos_lm)).norm();
	Float3 normal2 = (*pos_lm - *pos_rm).cross((*pos_right - *pos_rm)).norm();
	// Both vectors point "left" (looking from lm to rm). 

	float torsion = Float3::getAngle(normal2, normal1);
	//const float torsion = normal1.getAngleSigned(normal2);

	const bool angle_is_negative = (normal2.dot(*pos_left - *pos_lm)) > 0.f;
	torsion = angle_is_negative ? torsion * -1.f : torsion;

	normal2 *= -1;
	// Now  normal2 is flipped meaning both vectors point inward when 0 < torsion < 3.14, and outwards otherwise

	*potE += dihedral->k_phi * (1.f + cosf(dihedral->n * torsion - dihedral->phi_0));

	float torque = -dihedral->k_phi * (dihedral->n * sinf(dihedral->n * torsion - dihedral->phi_0));

	if (0) {
		normal1.print('1');
		normal2.print('2');
		pos_left->print('L');
		pos_lm->print('l');
		pos_rm->print('r');
		pos_right->print('R');
		printf("angle neg %d\n", angle_is_negative);
		//printf("torsion %f      ref %f     error %f     force: %f\n", torsion, dihedral->phi_0, error, force_scalar);
		float pot = dihedral->k_phi * (1 + cosf(dihedral->n * torsion - dihedral->phi_0));
		printf("torsion %f     torque: %f    pot %f     phi_0 %f k_phi %f\n", torsion, torque, pot, dihedral->phi_0, dihedral->k_phi);
	}



	results[0] = normal1 * (torque / (*pos_left - *pos_lm).len());
	results[3] = normal2 * (torque / (*pos_right - *pos_rm).len());
	// Not sure about the final two forces, for now we'll jsut split the sum of opposite forces between them.
	//results[1] = (results[0] + results[3]) * -1.f * 0.5;
	//results[2] = (results[0] + results[3]) * -1.f * 0.5;
	results[1] = (results[0]) * -1.f;
	results[2] = (results[3]) * -1.f;
}


// ------------------------------------------------------------------------------------------- LJ Forces -------------------------------------------------------------------------------------------//

__device__ static Float3 calcLJForce(const Float3* pos0, const Float3* pos1, float* data_ptr, float* potE, const float sigma, const float epsilon, int type1 = -1, int type2 = -1) {
	// Calculates LJ force on p0	(attractive to p1. Negative values = repulsion )//
	// input positions in cartesian coordinates [nm]
	// sigma [nm]
	// epsilon [J/mol]->[(kg*nm^2)/(ns^2*mol)]
	// Returns force in J/mol*M		?????????????!?!?//

	

	// Directly from book
	float dist_sq = (*pos1 - *pos0).lenSquared();
	float s = (sigma * sigma) / dist_sq;								// [nm^2]/[nm^2] -> unitless	// OPTIM: Only calculate sigma_squared, since we never use just sigma
	s = s * s * s;
	float force_scalar = 24.f * epsilon * s / dist_sq * (1.f - 2.f * s);// *FEMTO_TO_LIMA* FEMTO_TO_LIMA;	// Attractive. Negative, when repulsive		[(kg*nm^2)/(nm^2*ns^2*mol)] ->----------------------	[(kg)/(ns^2*mol)]	

	*potE += 4. * epsilon * s * (s - 1.f) * 0.5;


	if (((*pos1 - *pos0) * force_scalar).len() > 100 || true) {
		//printf("\nBlock %d thread %d type1 %d\n", blockIdx.x, threadIdx.x, type1);
		//((*pos1 - *pos0) * force_scalar).print('f');
		//pos0->print('0');
		//pos1->print('1');
		//printf("dist nm %f force %f\n", sqrt(dist_sq) / NANO_TO_LIMA, ((*pos1 - *pos0) * force_scalar).len());
	}
	//printf("\ndist: %d %f force %f\n", threadIdx.x, (*pos0 - *pos1).len() / NANO_TO_LIMA, ((*pos1 - *pos0) * force_scalar).len());

#ifdef LIMA_VERBOSE
	//if (threadIdx.x == 0 && blockIdx.x == 0) {
	//	float dist = (*pos0 - *pos1).len() * NORMALIZER;
	//	printf("\ndist %f force %f pot %f, sigma %f, s %f distsq %f eps %.10f\n", dist, force_scalar * (*pos1 - *pos0).len(), *potE, sigma * NORMALIZER, s, dist_sq, epsilon);
	//}
	pos0->print('0');
	pos0->print('1');
}
#endif
return (*pos1 - *pos0) * force_scalar;										// GN/mol [(kg*nm)/(ns^2*mol)]
	}


}	// End of namespace LimaForcecalc




