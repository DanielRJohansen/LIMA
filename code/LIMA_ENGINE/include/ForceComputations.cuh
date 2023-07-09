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
__device__ void calcDihedralbondForces1(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	Float3 grad_l, grad_lm, grad_rm, grad_r;

	float cosphi_1, cosphi_2;

	{
		const Float3 r1 = (*pos_lm - *pos_left).norm();
		const Float3 r2 = (*pos_rm - *pos_lm).norm();
		const Float3 r3 = (*pos_right - *pos_rm).norm();

		const Float3 m_hat = (r1.cross(r2)).norm();
		const Float3 n_hat = (r2.cross(r3)).norm();

		grad_l = (m_hat.cross(n_hat).cross(m_hat).cross(r2)) / (r1.cross(r2)).len();

		grad_lm = ((n_hat.cross(m_hat).cross(n_hat).cross(r3)) / (r2.cross(r3).len())) 
			- ((m_hat.cross(n_hat).cross(m_hat).cross(r1 + r2)) / (r1.cross(r2).len()));

		cosphi_1 = n_hat.dot(m_hat);

		//r1.print('m');
	}
	


	// Mirror it for right side
	{
		const Float3 r1 = (*pos_rm - *pos_right).norm();
		const Float3 r2 = (*pos_lm - *pos_rm).norm();
		const Float3 r3 = (*pos_left - *pos_lm).norm();

		const Float3 m_hat = (r1.cross(r2)).norm();
		const Float3 n_hat = (r2.cross(r3)).norm();

		grad_r = (m_hat.cross(n_hat).cross(m_hat).cross(r2)) / (r1.cross(r2)).len();

		grad_rm = ((n_hat.cross(m_hat).cross(n_hat).cross(r3)) / (r2.cross(r3).len()))
			- ((m_hat.cross(n_hat).cross(m_hat).cross(r1 + r2)) / (r1.cross(r2).len()));

		cosphi_2 = n_hat.dot(m_hat);
	}
	
	float torsion = std::acosf(cosphi_1);
	
	grad_l.print('l');
	printf("Lens %f %f %f %f\n", grad_l.len(), grad_lm.len(), grad_rm.len(), grad_r.len());


	*potE = (double)dihedral->k_phi * (1. + cos((double)dihedral->n * torsion - (double)dihedral->phi_0));
	const double torque = -dihedral->n * dihedral->k_phi * (dihedral->n * sin(dihedral->n * torsion - dihedral->phi_0));

	results[0] = grad_l * torque;
	results[1] = grad_lm * torque;
	results[2] = grad_rm * torque;
	results[3] = grad_r * torque;

}

// From resource: https://nosarthur.github.io/free%20energy%20perturbation/2017/02/01/dihedral-force.html
__device__ void calcDihedralbondForces(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	Float3 grad_l, grad_lm, grad_rm, grad_r;

	float cosphi_1, cosphi_2;
	float normaldot;

	{
		const Float3 r1 = (*pos_lm - *pos_left);
		const Float3 r2 = (*pos_rm - *pos_lm);
		const Float3 r3 = (*pos_right - *pos_rm);

		const Float3 m_hat = (r1.cross(r2)).norm();
		const Float3 n_hat = (r2.cross(r3)).norm();

		grad_l = (n_hat.cross(m_hat).cross(n_hat).cross(r2)) / (r1.cross(r2)).len();

		grad_lm = ((n_hat.cross(m_hat).cross(n_hat).cross(r3)) / (r2.cross(r3).len()))
			- ((m_hat.cross(n_hat).cross(m_hat).cross(r1 + r2)) / (r1.cross(r2).len()));

		cosphi_1 = n_hat.dot(m_hat);

		normaldot = n_hat.dot(r1);
		//r1.print('m');
	}



	// Mirror it for right side
	{
		const Float3 r1 = (*pos_rm - *pos_right);
		const Float3 r2 = (*pos_lm - *pos_rm);
		const Float3 r3 = (*pos_left - *pos_lm);

		const Float3 m_hat = (r1.cross(r2)).norm();
		const Float3 n_hat = (r2.cross(r3)).norm();

		grad_r = (m_hat.cross(n_hat).cross(m_hat).cross(r2)) / (r1.cross(r2)).len();

		((n_hat.cross(m_hat).cross(n_hat).cross(r2)) / (r1.cross(r2)).len()).print('a');
		grad_r.print('b');
		grad_rm = ((n_hat.cross(m_hat).cross(n_hat).cross(r3)) / (r2.cross(r3).len()))
			- ((m_hat.cross(n_hat).cross(m_hat).cross(r1 + r2)) / (r1.cross(r2).len()));

		cosphi_2 = n_hat.dot(m_hat);
	}


	printf("len ratio = %f mag ratio %f\n", (abs(pos_left->distToLine(*pos_lm, *pos_rm)) / abs(pos_right->distToLine(*pos_lm, *pos_rm))),	grad_r.len() / grad_l.len());

	//grad_l = grad_l.norm();
	//grad_r = grad_r.norm();
	//grad_lm = grad_lm.norm();
	//grad_rm = grad_rm.norm();

	float torsion = std::acosf(CPPD::min(cosphi_1, 1.f));
	//if (normaldot > 0.f) { torsion = -torsion; }


	//grad_l.print('l');
	//printf("Lens %f %f %f %f\n", grad_l.len(), grad_lm.len(), grad_rm.len(), grad_r.len());
	printf("Torsion %f cosphi %f is_1 %d\n", torsion, cosphi_1, cosphi_1 = 1.f);



	

	*potE = (double)dihedral->k_phi*1000.f * (1. + cos((double)dihedral->n * torsion - (double)dihedral->phi_0));
	double torque = -dihedral->n * dihedral->k_phi * 1000.f * (dihedral->n * sin(dihedral->n * torsion - dihedral->phi_0));

	if ((double)dihedral->k_phi * (1. + cos((double)dihedral->n * torsion - (double)dihedral->phi_0)) < 0.f) {
		printf("This is really not supposed to happen!\n");
		torque *= 100000.f;	// Force an error
	}

	//results[0] = grad_l * (torque / (*pos_lm - *pos_left).len());
	//results[3] = grad_r * (torque / (*pos_rm - *pos_right).len());
	//results[0] = grad_l * (torque / abs(pos_left->distToLine(*pos_lm, *pos_rm)));
	//results[3] = grad_r * (torque / abs(pos_right->distToLine(*pos_lm, *pos_rm)));
	results[0] = grad_l * torque;
	results[3] = grad_r * torque;


	float force_scalar = -(results[0] + results[3]).len();
	//results[1] = grad_lm * force_scalar;
	//results[2] = grad_rm * force_scalar;
	results[1] = grad_lm * torque;
	results[2] = grad_rm * torque;


	Float3 force_spillover = Float3{};
	for (int i = 0; i < 4; i++) {
		force_spillover += results[i];
	}
	force_spillover.print('s');
	results[0].print('l');
	for (int i = 0; i < 4; i++) {
		//results[i] -= force_spillover * 0.25f;
	}

	//force_spillover.print('f');
	//grad_lm.print('l');
	//grad_rm.print('r');

	{
		const Float3 normal1 = (*pos_left - *pos_lm).norm_d().cross((*pos_rm - *pos_lm).norm_d()).norm_d();
		Float3 normal2 = (*pos_lm - *pos_rm).norm_d().cross((*pos_right - *pos_rm).norm_d()).norm_d();
		// Both vectors point "left" (looking from lm to rm). 



		double torsion2 = Float3::getAngle(normal2, normal1);

		const bool angle_is_negative = (normal2.dot(*pos_left - *pos_lm)) > 0.f;
		if (angle_is_negative) {
			torsion2 = -torsion2;
		}

		//printf("my %f his %f normalsdot %f graddot %f gradr_dot_m\n", torsion2, acosf(cosphi_1), normaldot, grad_l.dot(grad_r), grad_r.dot(m_hat);
	}
	

}

__device__ void calcDihedralbondForces2(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	const Float3 normal1 = (*pos_left - *pos_lm).norm_d().cross((*pos_rm - *pos_lm).norm_d()).norm_d();
	Float3 normal2 = (*pos_lm - *pos_rm).norm_d().cross((*pos_right - *pos_rm).norm_d()).norm_d();
	// Both vectors point "left" (looking from lm to rm). 

	double torsion = Float3::getAngle(normal2, normal1);
	//const float torsion = normal1.getAngleSigned(normal2);

	const bool angle_is_negative = (normal2.dot(*pos_left - *pos_lm)) > 0.f;
	if (angle_is_negative) {
		torsion = -torsion;
	}

	//torsion = angle_is_negative ? torsion * -1.f : torsion;

	normal2 *= -1;
	// Now  normal2 is flipped meaning both vectors point inward when 0 < torsion < 3.14, and outwards otherwise
	

	// This is according to chatgpt
	*potE = (double)dihedral->k_phi * (1. + cos((double) dihedral->n * torsion - (double)dihedral->phi_0));
	const double torque =  -dihedral->n * dihedral->k_phi * (dihedral->n * sin(dihedral->n * torsion - dihedral->phi_0));




	if (0) {
		//normal1.print('1');
		//normal2.print('2');
		//pos_left->print('L');
		//pos_lm->print('l');
		//pos_rm->print('r');
		//pos_right->print('R');
		//printf("angle neg %d\n", angle_is_negative);
		//printf("torsion %f      ref %f     error %f     force: %f\n", torsion, dihedral->phi_0, error, force_scalar);
		float pot = dihedral->k_phi * (1 + cosf(dihedral->n * torsion - dihedral->phi_0));
		printf("\ntorsion %f [rad]    torque: %f    pot %f     phi_0 %f [rad] k_phi %f\n", 
			torsion, torque, pot, dihedral->phi_0, dihedral->k_phi);
	}



	//results[0] = normal1 * (torque / (*pos_left - *pos_lm).len());
	//results[3] = normal2 * (torque / (*pos_right - *pos_rm).len());
	results[0] = normal1 * (torque / pos_left->distToLine(*pos_lm, *pos_rm));
	results[3] = normal2 * (torque / pos_right->distToLine(*pos_lm, *pos_rm));

	// Not sure about the final two forces, for now we'll jsut split the sum of opposite forces between them.
	results[1] = (results[0] + results[3]) * -1.f * 0.5;
	results[2] = (results[0] + results[3]) * -1.f * 0.5;
	//results[1] = -results[0];
	//results[2] = -results[3];

#ifdef LIMASAFEMODE
	if (results[0].len() > 0.5f) {
		printf("\nDihedralBond: torsion %f [rad] torque: %f phi_0 [rad] %f k_phi %f\n", 
			torsion, torque, dihedral->phi_0, dihedral->k_phi);
	}
#endif
}


__device__ void calcImproperdihedralbondForces(Float3* i, Float3* j, Float3* k, Float3* l, ImproperDihedralBond* improper, Float3* results, float* potE) {
	Float3 normal1 = (*j - *i).cross((*k - *i)).norm();	// i is usually the center on
	Float3 normal2 = (*j - *l).cross((*k - *l)).norm();	// l is the outsider


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




