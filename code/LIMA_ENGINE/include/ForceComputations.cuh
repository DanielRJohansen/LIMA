#pragma once

//#include <cuda_runtime.h>

#include "LIMA_BASE/include/LimaTypes.cuh"
#include "LIMA_BASE/include/Constants.cuh"
#include "LIMA_BASE/include/Bodies.cuh"



namespace LimaForcecalc {

// ------------------------------------------------------------------------------------------- BONDED FORCES -------------------------------------------------------------------------------------------//

__device__ void calcPairbondForces(Float3* pos_a, Float3* pos_b, SingleBond* bondtype, Float3* results, float* potE) {
	// Calculates bond force on both particles					//
	// Calculates forces as J/mol*M								//
	// kb [J/(mol*lm^2)]
	Float3 difference = *pos_a - *pos_b;						//	[lm]
	const float error = difference.len() - bondtype->b0;				//	[lm]

	*potE += 0.5f * bondtype->kb * (error * error);				// [J/mol]

	const float force_scalar = -bondtype->kb * error;				//	[J/(mol*lm)] = [kg/(mol*s^2)]

	Float3 dir = difference.norm();								// dif_unit_vec, but shares variable with dif
	results[0] = dir * force_scalar;							// [kg * lm / (mol*ls^2)] = [lN]
	results[1] = dir * force_scalar * -1.f;						// [kg * lm / (mol*ls^2)] = [lN]

#ifdef LIMASAFEMODE
	if (results[0].len() > 0.5f) {
		printf("\nSingleBond: dist %f error: %f [nm] b0 %f [nm] force %f\n", difference.len()/NANO_TO_LIMA, error / NANO_TO_LIMA, bondtype->b0 / NANO_TO_LIMA, force_scalar);
	}
#endif
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

#ifdef LIMASAFEMODE
	if (results[0].len() > 0.1f) {
		printf("\nAngleBond: angle %f [rad] error %f [rad] force %f t0 [rad] %f kt %f\n", angle, error, results[0], angletype->theta_0, angletype->k_theta);
	}
#endif
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
		printf("\ntorsion %f [rad]    torque: %f    pot %f     phi_0 %f [rad] k_phi %f\n", 
			torsion, torque, pot, dihedral->phi_0, dihedral->k_phi);
	}



	results[0] = normal1 * (torque / (*pos_left - *pos_lm).len());
	results[3] = normal2 * (torque / (*pos_right - *pos_rm).len());
	// Not sure about the final two forces, for now we'll jsut split the sum of opposite forces between them.
	//results[1] = (results[0] + results[3]) * -1.f * 0.5;
	//results[2] = (results[0] + results[3]) * -1.f * 0.5;
	results[1] = (results[0]) * -1.f;
	results[2] = (results[3]) * -1.f;

#ifdef LIMASAFEMODE
	if (results[0].len() > 0.5f) {
		printf("\nDihedralBond: torsion %f [rad] torque: %f phi_0 [rad] %f k_phi %f\n", 
			torsion, torque, dihedral->phi_0, dihedral->k_phi);
	}
#endif
}


// ------------------------------------------------------------------------------------------- LJ Forces -------------------------------------------------------------------------------------------//
enum CalcLJOrigin { ComComIntra, ComComInter, ComSol, SolCom, SolSolIntra, SolSolInter };


__device__ static const char* calcLJOriginString[] = {
	"ComComIntra", "ComComInter", "ComSol", "SolCom", "SolSolIntra", "SolSolInter"
};

__device__ static Float3 calcLJForce(const Float3* pos0, const Float3* pos1, float* data_ptr, float* potE, const float sigma, const float epsilon, 
	CalcLJOrigin originSelect, /*For debug only*/
	int type1 = -1, int type2 = -1) {
	// Calculates LJ force on p0	(attractive to p1. Negative values = repulsion )//
	// input positions in cartesian coordinates [nm]
	// sigma [nm]
	// epsilon [J/mol]->[(kg*nm^2)/(ns^2*mol)]
	// Returns force in J/mol*M		?????????????!?!?//


	// Directly from book
	const float dist_sq = (*pos1 - *pos0).lenSquared();
	float s = (sigma * sigma) / dist_sq;								// [nm^2]/[nm^2] -> unitless	// OPTIM: Only calculate sigma_squared, since we never use just sigma
	s = s * s * s;
	const float force_scalar = 24.f * epsilon * s / dist_sq * (1.f - 2.f * s);	// Attractive. Negative, when repulsive		[(kg*nm^2)/(nm^2*ns^2*mol)] ->----------------------	[(kg)/(ns^2*mol)]	

	*potE += 4. * epsilon * s * (s - 1.f) * 0.5;

	const Float3 force = (*pos1 - *pos0) * force_scalar;

#ifdef LIMASAFEMODE
	if (force > 0.5f) {
		//printf("\nBlock %d thread %d\n", blockIdx.x, threadIdx.x);
		////((*pos1 - *pos0) * force_scalar).print('f');
		//pos0->print('0');
		//pos1->print('1');
		printf("\nLJ Force %s: dist nm %f force %f sigma %f t1 %d t2 %d\n",calcLJOriginString[(int)originSelect], sqrt(dist_sq) / NANO_TO_LIMA, ((*pos1 - *pos0) * force_scalar).len(), sigma/NANO_TO_LIMA, type1, type2);
	}
#endif

	return force;	// GN/mol [(kg*nm)/(ns^2*mol)]
	}


}	// End of namespace LimaForcecalc




