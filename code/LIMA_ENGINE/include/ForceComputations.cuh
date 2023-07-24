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
__device__ void calcDihedralbondForces5(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	return;
	Float3 grad_l, grad_lm, grad_rm, grad_r;
	float torsion;
	
	{
		const Float3 r1 = (*pos_lm - *pos_left);
		const Float3 r2 = (*pos_rm - *pos_lm);
		const Float3 r3 = (*pos_right - *pos_rm);

		const Float3 m_hat = (r1.cross(r2)).norm();
		const Float3 n_hat = (r2.cross(r3)).norm();

		//grad_l = (n_hat.cross(m_hat).cross(n_hat).cross(r2)) / (r1.cross(r2)).len();
		grad_r = (m_hat.cross(n_hat).cross(m_hat).cross(r2)) / (r1.cross(r2)).len();
		grad_lm = ((n_hat.cross(m_hat).cross(n_hat).cross(r3)) / (r2.cross(r3).len()))
			- ((m_hat.cross(n_hat).cross(m_hat).cross(r1 + r2)) / (r1.cross(r2).len()));

		torsion = abs(Float3::getAngle(m_hat, n_hat));
	}



	// Mirror it for right side
	{
		const Float3 r1 = (*pos_rm - *pos_right);
		const Float3 r2 = (*pos_lm - *pos_rm);
		const Float3 r3 = (*pos_left - *pos_lm);

		const Float3 m_hat = (r1.cross(r2)).norm();
		const Float3 n_hat = (r2.cross(r3)).norm();

		//grad_r = (n_hat.cross(m_hat).cross(n_hat).cross(r2)) / (r1.cross(r2)).len();
		grad_l = (m_hat.cross(n_hat).cross(m_hat).cross(r2)) / (r1.cross(r2)).len();

		//((n_hat.cross(m_hat).cross(n_hat).cross(r2)) / (r1.cross(r2)).len()).print('a');
		//grad_r.print('b');
		grad_rm = ((n_hat.cross(m_hat).cross(n_hat).cross(r3)) / (r2.cross(r3).len()))
			- ((m_hat.cross(n_hat).cross(m_hat).cross(r1 + r2)) / (r1.cross(r2).len()));
	}


	//	printf("len ratio = %f mag ratio %f\n", (abs(pos_left->distToLine(*pos_lm, *pos_rm)) / abs(pos_right->distToLine(*pos_lm, *pos_rm))),	grad_r.len() / grad_l.len());


	//grad_l.print('1');
	//grad_lm.print('2');
	//grad_rm.print('3');
	//grad_r.print('4');
	//(grad_l + grad_lm + grad_rm + grad_r).print('S');

	//grad_l = grad_l.norm();
	//grad_r = grad_r.norm();
	//grad_lm = grad_lm.norm();
	//grad_rm = grad_rm.norm();

	float dtl = abs(pos_left->distToLine(*pos_lm, *pos_rm));
	float dtl2 = abs(pos_right->distToLine(*pos_lm, *pos_rm));
	printf("grad_l %.12f myscalar %.12f ratio %f\n", grad_l.len(), 1.f / dtl, grad_l.len() / (1.f / dtl));
	printf("grad_r %.12f myscalar %.12f ratio %f\n", grad_r.len(), 1.f / dtl2, grad_r.len() / (1.f / dtl2));

	*potE = (double)dihedral->k_phi * (1. + cos((double)dihedral->n * torsion - (double)dihedral->phi_0));
	double torque = -dihedral->n * dihedral->k_phi * (dihedral->n * sin(dihedral->n * torsion - dihedral->phi_0));

	printf("torsion %f torque %f\n", torsion, torque);

	if ((double)dihedral->k_phi * (1. + cos((double)dihedral->n * torsion - (double)dihedral->phi_0)) < 0.f) {
		printf("This is really not supposed to happen!\n");
		torque *= 100000.f;	// Force an error
	}



	results[0] = grad_l * (torque);
	results[1] = grad_lm * (torque);
	results[2] = grad_rm * (torque);
	results[3] = grad_r * (torque);

	//results[0] = grad_l * (torque / (*pos_lm - *pos_left).len());
	//results[3] = grad_r * (torque / (*pos_rm - *pos_right).len());
	//results[0] = grad_l * (torque / abs(pos_left->distToLine(*pos_lm, *pos_rm)));
	//results[3] = grad_r * (torque / abs(pos_right->distToLine(*pos_lm, *pos_rm)));
	//results[0] = grad_l * torque;
	//results[3] = grad_r * torque;

	// Caclculate the scaler needed for the interior gradients to neutralize the exterior gradients
	//const Float3 ext = -(results[0] + results[3]);	// negative sum of exterior atoms
	//
	////const float scalar_lm = (ext.y * grad_rm.x - grad_rm.y * ext.x) / (grad_rm.x * grad_lm.y + grad_lm.x); wrong
	//const float scalar_lm = (ext.y * grad_rm.x - grad_rm.y * ext.x) / (grad_rm.x - grad_lm.y - grad_rm.y * grad_lm.x);
	//const float scalar_rm = (ext.x - grad_lm.x * scalar_lm) / grad_rm.x;
	//
	//results[1] = grad_lm * scalar_lm;
	//results[2] = grad_rm * scalar_rm;



	Float3 force_spillover = Float3{};
	for (int i = 0; i < 4; i++) {
		force_spillover += results[i];
	}
	//ext.print('e');
	force_spillover.print('s');
	//results[0].print('L');
	//results[1].print('l');
	//results[2].print('r');
	//results[3].print('R');
	printf("\n");
	//for (int i = 0; i < 4; i++) {
	//	//results[i] -= force_spillover * 0.25f;
	//}

	//force_spillover.print('f');
	//grad_lm.print('l');
	//grad_rm.print('r');

	{
		//const Float3 normal1 = (*pos_left - *pos_lm).norm_d().cross((*pos_rm - *pos_lm).norm_d()).norm_d();
		//Float3 normal2 = (*pos_lm - *pos_rm).norm_d().cross((*pos_right - *pos_rm).norm_d()).norm_d();
		//// Both vectors point "left" (looking from lm to rm). 



		//double torsion2 = Float3::getAngle(normal2, normal1);

		//const bool angle_is_negative = (normal2.dot(*pos_left - *pos_lm)) > 0.f;
		//if (angle_is_negative) {
		//	torsion2 = -torsion2;
		//}

		//printf("my %f his %f normalsdot %f graddot %f gradr_dot_m\n", torsion2, acosf(cosphi_1), normaldot, grad_l.dot(grad_r), grad_r.dot(m_hat);
	}
}


__device__ void calcDihedralbondForces8(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	Float3 r12 = *pos_lm - *pos_left;
	Float3 r23 = *pos_rm - *pos_lm;
	Float3 r34 = *pos_right - *pos_rm;

	r12 = r12 / NANO_TO_LIMA;
	r23 = r23 / NANO_TO_LIMA;
	r34 = r34 / NANO_TO_LIMA;



	// First get the correct direction the exterior atoms should move in
	Float3 dir1 = r23.cross(r12).norm();
	Float3 dir4 = r23.cross(r34).norm();
	if (dir1.dot(r34) > 0.f) {
		dir4 = -dir4;
	}
	else {
		dir1 = -dir1;
	}
	const Float3 delta1 = dir1 / pos_left->distToLine(*pos_lm, *pos_rm);
	const Float3 delta4 = dir4 / pos_right->distToLine(*pos_lm, *pos_rm);

	r34.print('r');
	delta1.print('d');

	const Float3 fuckisthis = r34.cross(delta4) - r12.cross(delta1);

	Float3 delta2 = (fuckisthis - delta1).norm();
	Float3 delta3 = (-delta4 - fuckisthis).norm();



	Float3 ext = -(delta1 + delta4);
	const float scalar_lm = (ext.y * delta3.x - delta3.y * ext.x) / (delta3.x - delta2.y - delta3.y * delta2.x);
	const float scalar_rm = (ext.x - delta2.x * scalar_lm) / delta3.x;

	delta2 *= scalar_lm;
	delta3 *= scalar_rm;

	//delta2 = delta2 - r23.norm() * r23.dot(delta2) / r23.len();
	//delta3 = delta3 - r23.norm() * r23.dot(delta3) / r23.len();

	const float torsion = Float3::getAngle(delta1, delta4);
	*potE = (double)dihedral->k_phi * (1. + cos((double)dihedral->n * torsion - (double)dihedral->phi_0));
	const double torque = -dihedral->k_phi * (dihedral->n * sin(dihedral->n * torsion - dihedral->phi_0));
	//float torque = 10000000.f;

	results[0] = delta1 * (torque);
	results[1] = delta2 * (torque);
	results[2] = delta3 * (torque);
	results[3] = delta4 * (torque);

	results[0].print('1');
	results[1].print('2');
	results[2].print('3');
	results[3].print('4');

	Float3 force_spillover = Float3{};
	for (int i = 0; i < 4; i++) {
		force_spillover += results[i];
	}
	if (force_spillover.len() > 0.0001) {
		
	}
	force_spillover.print('s');
	printf("\n");

}

__device__ void calcDihedralbondForces7(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	const Float3 r21 = *pos_left - *pos_lm;
	const Float3 r32 = *pos_lm - *pos_rm;
	const Float3 r43 = *pos_rm - *pos_right;
	//return;
	Float3 A = r21.cross(r32);
	const float rA = A.len();
	Float3 B = r32.cross(r43);
	const float rB = B.len();

	// If this is very small, the torsion angle is undefined
	if (rA * rB < 0.0000001f) { 
		printf("returning\n");
		return; }

	A = A.norm();
	B = B.norm();

	float cos_phi = A.dot(B);
	if (cos_phi > 1.0) cos_phi = 1.0;
	if (cos_phi < -1.0) cos_phi = -1.0;

	const float torsion = acosf(cos_phi);


	
	// Simple implementation
	//const float error = torsion - dihedral->phi_0;				// [rad]
	//*potE += dihedral->k_phi * error * error * 0.5f;		// Energy [J/mol]0
	//float torque = dihedral->k_phi * (error);				// Torque [J/(mol*rad)]

	*potE = (double)dihedral->k_phi * (1. + cos((double)dihedral->n * torsion - (double)dihedral->phi_0));
	const double torque = -dihedral->k_phi * (dihedral->n * sin(dihedral->n * torsion - dihedral->phi_0));

	Float3 dcosdA = (A - B * cos_phi) / rA;
	Float3 dcosdB = (B - A * cos_phi) / rB;

	const Float3 d1 = r32.cross(dcosdA);
	const Float3 d2 = (r43.cross(dcosdB) - r21.cross(dcosdA));
	const Float3 d3 = dcosdB.cross(r32);

	float scalar = 1.f / d1.len();
	//d1*= scalar;
	//d2 *= scalar;
	//d3 *= scalar;
	//Float3 f1 = r32.cross(dcosdA) * torque;
	//Float3 f2 = (r43.cross(dcosdB) - r21.cross(dcosdA)) * torque;


	//printf("cosphi %f torsion %f\n", cos_phi, torsion);
	//f1.print('1');
	//f2.print('2');
	//f3.print('3');
	
	// Remove the component along rotational axis.
	Float3 g1 = d3 - d2;
	Float3 g2 = d2 - d1;

	g1 = g1 - r32.norm() * r32.dot(g1)/r32.len();
	g2 = g2 - r32.norm() * r32.dot(g2)/r32.len();

	if (r32.dot(g1) + r32.dot(g2) > 0.0001) {
		printf("dots %f %f\n", r32.dot(g1), r32.dot(g2));
	}
	
	results[3] = d1 * torque;
	results[2] = g2 * torque;
	results[1] = g1 * torque;
	results[0] = -d3 * torque;

	

	float error = torsion - dihedral->phi_0;
	if (dihedral->n == 3) {
		error = CPPD::min(abs(error), abs(torsion - (dihedral->phi_0 + PI * 2.f / 3.f)));
		error = CPPD::min(abs(error), abs(torsion - (dihedral->phi_0 - PI * 2.f / 3.f)));
	}

	//if (threadIdx.x == 1) {
		//printf("cosphi %f torsion %f torque %f error %f\n", cos_phi, torsion, torque, error);
		//pos_left->print('L');
		//pos_lm->print('l');
		//pos_rm->print('r');
		//pos_right->print('R');

		//printf("ra %f rb %f\n", rA, rB);
	//}
	

	
	Float3 force_spillover = Float3{};
	for (int i = 0; i < 4; i++) {
		force_spillover += results[i];
	}
	if (force_spillover.len() > 0.0001) {
		force_spillover.print('s');
	}
	

	float dtl = abs(pos_left->distToLine(*pos_lm, *pos_rm));
	float dtl2 = abs(pos_right->distToLine(*pos_lm, *pos_rm));
	//printf("grad_l %.12f myscalar %.12f ratio %f\n", r32.cross(dcosdA).len(), 1.f / dtl, r32.cross(dcosdA).len() / (1.f / dtl));
	//printf("grad_r %.12f myscalar %.12f ratio %f\n", dcosdB.cross(r32).len(), 1.f / dtl2, dcosdB.cross(r32).len() / (1.f / dtl2));

	float leverarm_correction = 1.f / abs(pos_left->distToLine(*pos_lm, *pos_rm))
		/ CPPD::max(r32.cross(dcosdA).len(),0.00001f);	// This is the wrong one, they are still mixed up

	float ratio = r32.cross(dcosdA).len() / (1.f / dtl);
	for (int i = 0; i < 4; i++) {
		//results[i] *= leverarm_correction;
	}
	
	if (threadIdx.x == 1) {
		dcosdA.print('d');
		A.print('A');
		printf("lever %f rA %f\n", leverarm_correction, rA);
		results[3].print('3');
		results[2].print('2');
		results[1].print('1');
		results[0].print('0');
	}
	
}

__device__ void calcDihedralbondForces(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	Float3 r12 = *pos_lm - *pos_left;
	Float3 r23 = *pos_rm - *pos_lm;
	Float3 r34 = *pos_right - *pos_rm;
	r12 = r12 / NANO_TO_LIMA;
	r23 = r23 / NANO_TO_LIMA;
	r34 = r34 / NANO_TO_LIMA;

	Float3 A = r12.cross(r23);
	float rAinv = 1.f / A.len();
	Float3 B = r23.cross(r34);
	float rBinv = 1.f / B.len();
	Float3 C = r23.cross(A);
	float rCinv = 1.f / C.len();

	float cos_phi = A.dot(B) * (rAinv * rBinv);
	float sin_phi = C.dot(B) * (rCinv * rBinv);

	float phi = -atan2(sin_phi, cos_phi);


	//printf("sin %f cos %f phi %f\n", sin_phi, cos_phi, phi);

	*potE = (double)dihedral->k_phi * (1. + cos((double)dihedral->n * phi - (double)dihedral->phi_0));
	float torque = dihedral->k_phi * (dihedral->n * sin(dihedral->n * phi - dihedral->phi_0)) / NANO_TO_LIMA;
	//torque = 1.f;




	//r12.print('a');
	//r23.print('b');
	//r34.print('c');

	//A.print('A');
	//B.print('B');
	//C.print('C');


	B = B * rBinv;
	Float3 f1, f2, f3;
	if (fabs(sin_phi) > 0.1f) {
		A = A * rAinv;

		Float3 dcosdA = (A * cos_phi - B) * rAinv;
		Float3 dcosdB = (B * cos_phi - A) * rBinv;

		torque = torque / sin_phi;

		f1 = Float3{
			torque * (r23.y * dcosdA.z - r23.z * dcosdA.y),
			torque * (r23.z * dcosdA.x - r23.x * dcosdA.z),
			torque * (r23.x * dcosdA.y - r23.y * dcosdA.x)
		};

		f3 = Float3{
			torque * (r23.z * dcosdB.y - r23.y * dcosdB.z),
			torque * (r23.x * dcosdB.z - r23.z * dcosdB.x),
			torque * (r23.y * dcosdB.x - r23.x * dcosdB.y)
		};

		f2 = Float3{
			torque * (r12.z * dcosdA.y - r12.y * dcosdA.z + r34.y * dcosdB.z - r34.z * dcosdB.y),
			torque * (r12.x * dcosdA.z - r12.z * dcosdA.x + r34.z * dcosdB.x - r34.x * dcosdB.z),
			torque * (r12.y * dcosdA.x - r12.x * dcosdA.y + r34.x * dcosdB.y - r34.y * dcosdB.x)
		};

	}
	else {
		C = C * rCinv;

		Float3 dsindC = (C * sin_phi - B) * rCinv;
		Float3 dsindB = (B * sin_phi - C) * rBinv;

		torque = -torque / cos_phi;

		f1 = Float3{
			torque * ((r23.y * r23.y + r23.z * r23.z) * dsindC.x - r23.x * r23.y * dsindC.y - r23.x * r23.z * dsindC.z),
			torque * ((r23.z * r23.z + r23.x * r23.x) * dsindC.y - r23.y * r23.z * dsindC.z - r23.y * r23.x * dsindC.x),
			torque * ((r23.x * r23.x + r23.y * r23.y) * dsindC.z - r23.z * r23.x * dsindC.x - r23.z * r23.y * dsindC.y)
		};
		

		f3 = dsindB.cross(r23) * torque;

		f2 = Float3{
			torque * (-(r23.y * r12.y + r23.z * r12.z) * dsindC.x + (2.f * r23.x * r12.y - r12.x * r23.y) * dsindC.y + (2.f * r23.x * r12.z - r12.x * r23.z) * dsindC.z + dsindB.z * r34.y - dsindB.y * r34.z),
			torque * (-(r23.z * r12.z + r23.x * r12.x) * dsindC.y + (2.f * r23.y * r12.z - r12.y * r23.z) * dsindC.z + (2.f * r23.y * r12.x - r12.y * r23.x) * dsindC.x + dsindB.x * r34.z - dsindB.z * r34.x),
			torque * (-(r23.x * r12.x + r23.y * r12.y) * dsindC.z + (2.f * r23.z * r12.x - r12.z * r23.x) * dsindC.x + (2.f * r23.z * r12.y - r12.z * r23.y) * dsindC.y + dsindB.y * r34.x - dsindB.x * r34.y)
		};
	}

	//A.print('A');
	//B.print('B');
	//C.print('C');

	//f1.print('q');
	//f2.print('w');
	//f3.print('e');

	//printf("\ntorsion %f [rad]    torque: %f    pot %f     phi_0 %f [rad] k_phi %f\n", phi, torque, *potE, dihedral->phi_0, dihedral->k_phi);

	results[0] = f1;
	results[1] = f2-f1;
	results[2] = f3-f2;
	results[3] = -f3;

	/*results[0].print('1');
	results[1].print('2');
	results[2].print('3');
	results[3].print('4');*/

	Float3 force_spillover = Float3{};
	for (int i = 0; i < 4; i++) {
		force_spillover += results[i];
	}
	if (force_spillover.len() > 0.00000001) {
		force_spillover.print('s');
	}
}


__device__ void calcDihedralbondForces6(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
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
	*potE = (double)dihedral->k_phi * (1. + cos((double)dihedral->n * torsion - (double)dihedral->phi_0));
	const double torque = -dihedral->n * dihedral->k_phi * (dihedral->n * sin(dihedral->n * torsion - dihedral->phi_0));




	if (0) {
		//normal1.print('1');
		//normal2.print('2');
		//
		//
		//
		//;
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




