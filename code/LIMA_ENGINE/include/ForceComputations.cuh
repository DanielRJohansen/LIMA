#pragma once

//#include <cuda_runtime.h>

#include "LIMA_BASE/include/LimaTypes.cuh"
#include "LIMA_BASE/include/Constants.cuh"
#include "LIMA_BASE/include/Bodies.cuh"



namespace LimaForcecalc {

// ------------------------------------------------------------------------------------------- BONDED FORCES -------------------------------------------------------------------------------------------//

__device__ void calcSinglebondForces(Float3* pos_a, Float3* pos_b, SingleBond* bondtype, Float3* results, float* potE) {
	// Calculates bond force on both particles					//
	// Calculates forces as J/mol*M								//
	// kb [J/(mol*lm^2)]
	Float3 difference = *pos_a - *pos_b;						//	[lm]
	const double error = difference.len() - (double)bondtype->b0;				//	[lm]

	*potE = 0.5f * bondtype->kb * (error * error);				// [J/mol]
	const float force_scalar = -bondtype->kb * error;				//	[J/(mol*lm)] = [kg/(mol*s^2)]

	const double error_fm = error / (double)PICO_TO_LIMA;
	//*potE += 0.5 * (double)bondtype->kb * (error_fm * error);				// [J/mol]
	//const double force_scalar = (double) -bondtype->kb * error_fm;				//	[J/(mol*lm)] = [kg/(mol*s^2)]

	Float3 dir = difference.norm();								// dif_unit_vec, but shares variable with dif
	results[0] = dir * force_scalar;							// [kg * lm / (mol*ls^2)] = [lN]
	results[1] = -dir * force_scalar;						// [kg * lm / (mol*ls^2)] = [lN]

#ifdef LIMASAFEMODE
	if (abs(error) > bondtype->b0/2.f || false) {
		printf("\nSingleBond: dist %f error: %f [nm] b0 %f [nm] kb %.10f [J/mol] force %f\n", difference.len()/NANO_TO_LIMA, error / NANO_TO_LIMA, bondtype->b0 / NANO_TO_LIMA, bondtype->kb, force_scalar);
		printf("errfm %f\n", error_fm);
		printf("pot %f\n", *potE);
	}
#endif
}

__device__ void calcAnglebondForces(Float3* pos_left, Float3* pos_middle, Float3* pos_right, const AngleBond const* angletype, Float3* results, float* potE) {
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
		printf("\nAngleBond: angle %f [rad] error %f [rad] force %f t0 %f [rad] kt %f\n", angle, error, results[0].len(), angletype->theta_0, angletype->k_theta);
	}
#endif
}

// From resource: https://nosarthur.github.io/free%20energy%20perturbation/2017/02/01/dihedral-force.html
// Greatly inspired by OpenMD's CharmmDihedral algorithm
__device__ void calcDihedralbondForces(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	Float3 r12 = (*pos_lm - *pos_left) / NANO_TO_LIMA;
	Float3 r23 = (*pos_rm - *pos_lm) / NANO_TO_LIMA;
	Float3 r34 = (*pos_right - *pos_rm) / NANO_TO_LIMA;

	Float3 A = r12.cross(r23);
	const float rAinv = 1.f / A.len();
	Float3 B = r23.cross(r34);
	const float rBinv = 1.f / B.len();
	Float3 C = r23.cross(A);
	const float rCinv = 1.f / C.len();

	const float cos_phi = A.dot(B) * (rAinv * rBinv);
	const float sin_phi = C.dot(B) * (rCinv * rBinv);
	const float torsion = -atan2(sin_phi, cos_phi);


	*potE = (double)dihedral->k_phi * (1. + cos((double)dihedral->n * torsion - (double)dihedral->phi_0));
	const float torque = dihedral->k_phi * (dihedral->n * sin(dihedral->n * torsion - dihedral->phi_0)) / NANO_TO_LIMA;



	B = B * rBinv;
	Float3 f1, f2, f3;
	if (fabs(sin_phi) > 0.1f) {
		A = A * rAinv;

		const Float3 dcosdA = (A * cos_phi - B) * rAinv;
		const Float3 dcosdB = (B * cos_phi - A) * rBinv;

		const float k = torque / sin_phi;	// Wtf is k????

		f1 = r23.cross(dcosdA) * k;
		f3 = -r23.cross(dcosdB) * k;
		f2 = (r34.cross(dcosdB) - r12.cross(dcosdA)) * k;
	}
	else {
		C = C * rCinv;

		const Float3 dsindC = (C * sin_phi - B) * rCinv;
		const Float3 dsindB = (B * sin_phi - C) * rBinv;

		const float k = -torque / cos_phi;

		// TODO: This is ugly, fix it
		f1 = Float3{
			((r23.y * r23.y + r23.z * r23.z) * dsindC.x - r23.x * r23.y * dsindC.y - r23.x * r23.z * dsindC.z),
			((r23.z * r23.z + r23.x * r23.x) * dsindC.y - r23.y * r23.z * dsindC.z - r23.y * r23.x * dsindC.x),
			((r23.x * r23.x + r23.y * r23.y) * dsindC.z - r23.z * r23.x * dsindC.x - r23.z * r23.y * dsindC.y)
		} * k;
		

		f3 = dsindB.cross(r23) * k;

		f2 = Float3{
			(-(r23.y * r12.y + r23.z * r12.z) * dsindC.x + (2.f * r23.x * r12.y - r12.x * r23.y) * dsindC.y + (2.f * r23.x * r12.z - r12.x * r23.z) * dsindC.z + dsindB.z * r34.y - dsindB.y * r34.z),
			(-(r23.z * r12.z + r23.x * r12.x) * dsindC.y + (2.f * r23.y * r12.z - r12.y * r23.z) * dsindC.z + (2.f * r23.y * r12.x - r12.y * r23.x) * dsindC.x + dsindB.x * r34.z - dsindB.z * r34.x),
			(-(r23.x * r12.x + r23.y * r12.y) * dsindC.z + (2.f * r23.z * r12.x - r12.z * r23.x) * dsindC.x + (2.f * r23.z * r12.y - r12.z * r23.y) * dsindC.y + dsindB.y * r34.x - dsindB.x * r34.y)
		} * k;
	}

	results[0] = f1;
	results[1] = f2-f1;
	results[2] = f3-f2;
	results[3] = -f3;

#ifdef LIMASAFEMODE
	Float3 force_spillover = Float3{};
	for (int i = 0; i < 4; i++) {
		force_spillover += results[i];
	}
	if (force_spillover.len() > 0.00000001) {
		force_spillover.print('s');
	}
#endif
}

__device__ void calcImproperdihedralbondForces(Float3* i, Float3* j, Float3* k, Float3* l, ImproperDihedralBond* improper, Float3* results, float* potE) {
	Float3 normal1 = (*j - *i).cross((*k - *i)).norm();	// i is usually the center on
	Float3 normal2 = (*j - *l).cross((*k - *l)).norm();	// l is the outsider
	return;	// DISABLED untill made stable

	float angle = Float3::getAngle(normal2, normal1); 
	if (angle < 0) {
		printf("angle:: %f", angle);
	}

	const bool angle_is_negative = (normal1.dot(*l- *i)) > 0.f;
	if (angle_is_negative) {
		angle = -angle;
	}


	if (angle > PI || angle < -PI) {
		printf("Anlg too large!! %f\n\n\n\n", angle);
	}



	const float error = angle - improper->psi_0;

	*potE += 0.5f * improper->k_psi * (error*error);
	const float torque = improper->k_psi * (angle - improper->psi_0);

	if (0) {
		//normal1.print('1');
		//normal2.print('2');
		//i->print('i');
		//j->print('j');
		//k->print('k');
		//l->print('l');
		float pot = 0.5f * improper->k_psi * (error * error);
		printf("\tangle %f [rad]    torque: %f    pot %f     phi_0 %f [rad] k_phi %f\n\n",
			angle, torque, pot, improper->psi_0, improper->k_psi);
	}

	// Tongue in cheek here. Each particle follows the vector that will minize the potential
	// at this current step. I am fairly sure this is correct. TODO: Clear with Ali
	//results[3] = normal1 * (torque / (*l - *i).len());	// l
	//
	//results[1] = normal2 * (torque / (*j - *i).len()); // j
	//results[2] = normal2 * (torque / (*k - *i).len()); // k
	//
	//results[0] = -(results[3] + results[1] + results[2]);	// i restores quilibrium of force

	// This is the simple way, always right-ish
	results[3] = normal2 * (torque / (*l - *i).len());	// l

	results[1] = normal1 * (torque / (*j - *i).len()); // j
	results[2] = normal1 * (torque / (*k - *i).len()); // j
	
	results[0] = -(results[3] + results[1] + results[2]);	// i restores quilibrium of force


#ifdef LIMASAFEMODE
	if (results[0].len() > 0.5f) {
		printf("\nImproperdihedralBond: angle %f [rad] torque: %f psi_0 [rad] %f k_psi %f\n",
			angle, torque, improper->psi_0, improper->k_psi);
	}
#endif
}

// ------------------------------------------------------------------------------------------- LJ Forces -------------------------------------------------------------------------------------------//
enum CalcLJOrigin { ComComIntra, ComComInter, ComSol, SolCom, SolSolIntra, SolSolInter };


__device__ static const char* calcLJOriginString[] = {
	"ComComIntra", "ComComInter", "ComSol", "SolCom", "SolSolIntra", "SolSolInter"
};

template<bool em_variant=0>
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

	*potE += 4. * epsilon * s * (s - 1.f) * 0.5;	// 0.5 to account for 2 particles doing the same calculation
	const Float3 force = (*pos1 - *pos0) * force_scalar;

#ifdef LIMASAFEMODE
	if constexpr (!em_variant) {	// During EM dt is lower, so large forces are not a problem
		auto pot = 4. * epsilon * s * (s - 1.f) * 0.5;
		if (force.len() > 1.f || pot > 1e+8) {
			//printf("\nBlock %d thread %d\n", blockIdx.x, threadIdx.x);
			////((*pos1 - *pos0) * force_scalar).print('f');
			//pos0->print('0');
			//pos1->print('1');
			printf("\nLJ Force %s: dist nm %f force %f sigma %f t1 %d t2 %d\n", 
				calcLJOriginString[(int)originSelect], sqrt(dist_sq) / NANO_TO_LIMA, ((*pos1 - *pos0) * force_scalar).len(), sigma / NANO_TO_LIMA, type1, type2);
		}
	}	
#endif

	return force;	// GN/mol [(kg*nm)/(ns^2*mol)]
	}
}	// End of namespace LimaForcecalc