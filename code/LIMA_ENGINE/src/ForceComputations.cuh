#pragma once

//#include <cuda_runtime.h>

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "EngineUtils.cuh"


namespace LimaForcecalc {

// ------------------------------------------------------------------------------------------- BONDED FORCES -------------------------------------------------------------------------------------------//

__device__ void calcSinglebondForces(const Float3& pos_a, const Float3& pos_b, const SingleBond& bondtype, Float3* results, float& potE, bool bridgekernel) {
	// Calculates bond force on both particles					//
	// Calculates forces as J/mol*M								//
	// kb [J/(mol*lm^2)]
	const Float3 difference = pos_a - pos_b;						//	[lm]
	const float error = difference.len() - bondtype.b0;				//	[lm]

	if constexpr (CALC_POTE) {
		potE = 0.5f * bondtype.kb * (error * error);				// [J/mol]
	}
	const float force_scalar = -bondtype.kb * error;				//	[J/(mol*lm)] = [kg/(mol*s^2)]

	//const double error_fm = error / PICO_TO_LIMA;
	//*potE += 0.5 * (double)bondtype->kb * (error_fm * error);				// [J/mol]
	//const double force_scalar = (double) -bondtype->kb * error_fm;				//	[J/(mol*lm)] = [kg/(mol*s^2)]

	const Float3 dir = difference.norm();								// dif_unit_vec, but shares variable with dif
	results[0] = dir * force_scalar;							// [kg * lm / (mol*ls^2)] = [lN]
	results[1] = -dir * force_scalar;						// [kg * lm / (mol*ls^2)] = [lN]

#if defined LIMASAFEMODE
	if (abs(error) > bondtype.b0/2.f || 0) {
		//std::cout << "SingleBond : " << kernelname << " dist " << difference.len() / NANO_TO_LIMA;
		printf("\nSingleBond: bridge %d dist %f error: %f [nm] b0 %f [nm] kb %.10f [J/mol] force %f\n", bridgekernel, difference.len() / NANO_TO_LIMA, error / NANO_TO_LIMA, bondtype.b0 / NANO_TO_LIMA, bondtype.kb, force_scalar);
		pos_a.print('a');
		pos_b.print('b');
		//printf("errfm %f\n", error_fm);
		//printf("pot %f\n", *potE);
	}
#endif
}

__device__ void calcAnglebondForces(const Float3& pos_left, const Float3& pos_middle, const Float3& pos_right, const AngleBond& angletype, Float3* results, float& potE) {
	const Float3 v1 = (pos_left - pos_middle).norm();
	const Float3 v2 = (pos_right - pos_middle).norm();
	const Float3 normal = v1.cross(v2).norm();	// Poiting towards y, when right is pointing toward x

	const Float3 inward_force_direction1 = (v1.cross(normal * -1.f)).norm();
	const Float3 inward_force_direction2 = (v2.cross(normal)).norm();

	const float angle = Float3::getAngleOfNormVectors(v1, v2);
	const float error = angle - angletype.theta_0;				// [rad]

	// Simple implementation
	if constexpr (CALC_POTE) {
		potE = angletype.k_theta * error * error * 0.5f;		// Energy [J/mol]0
	}
	const float torque = angletype.k_theta * (error);				// Torque [J/(mol*rad)]

	// Correct implementation
	//potE = -angletype.k_theta * (cosf(error) - 1.f);		// Energy [J/mol]
	//const float torque = angletype.k_theta * sinf(error);	// Torque [J/(mol*rad)]

	results[0] = inward_force_direction1 * (torque / (pos_left - pos_middle).len());
	results[2] = inward_force_direction2 * (torque / (pos_right - pos_middle).len());
	results[1] = (results[0] + results[2]) * -1.f;



#if defined LIMASAFEMODE
	if (results[0].len() > 0.1f) {
		printf("\nAngleBond: angle %f [rad] error %f [rad] force %f t0 %f [rad] kt %f\n", angle, error, results[0].len(), angletype.theta_0, angletype.k_theta);
	}
#endif
}

// From resource: https://nosarthur.github.io/free%20energy%20perturbation/2017/02/01/dihedral-force.html
// Greatly inspired by OpenMD's CharmmDihedral algorithm
__device__ void calcDihedralbondForces(const Float3& pos_left, const Float3& pos_lm, const Float3& pos_rm, const Float3& pos_right, const DihedralBond& dihedral, Float3* results, float& potE) {
	const Float3 r12 = (pos_lm - pos_left);
	const Float3 r23 = (pos_rm - pos_lm);
	const Float3 r34 = (pos_right - pos_rm);

	Float3 A = r12.cross(r23);
	const float rAinv = 1.f / A.len();
	Float3 B = r23.cross(r34);
	const float rBinv = 1.f / B.len();
	Float3 C = r23.cross(A);
	const float rCinv = 1.f / C.len();

	const float cos_phi = A.dot(B) * (rAinv * rBinv);
	const float sin_phi = C.dot(B) * (rCinv * rBinv);
	const float torsion = -atan2(sin_phi, cos_phi);

	if constexpr (CALC_POTE) {
		potE = __half2float(dihedral.k_phi) * (1. + cos(__half2float(dihedral.n) * torsion - __half2float(dihedral.phi_0)));
	}
	const float torque = __half2float(dihedral.k_phi) * (__half2float(dihedral.n) * sin(__half2float(dihedral.n) * torsion - __half2float(dihedral.phi_0))) / NANO_TO_LIMA;



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

#if defined LIMASAFEMODE
	Float3 force_spillover = Float3{};
	for (int i = 0; i < 4; i++) {
		force_spillover += results[i];
	}
	if (force_spillover.len()*10000.f > results[0].len()) {
		force_spillover.print('s');
		results[0].print('0');
	}

	if (isnan(potE) && r12.len() == 0) {
		printf("Bad torsion: Block %d t %d\n", blockIdx.x, threadIdx.x);
		//printf("torsion %f torque %f\n", torsion, torque);
		//r12.print('1');
		//pos_left.print('L');
		//pos_lm.print('l');
		potE = 6969696969.f;
	}
#endif
}

// https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html
// Plane described by i,j,k, and l is out of plane, connected to i
__device__ void calcImproperdihedralbondForces(const Float3& i, const Float3& j, const Float3& k, const Float3& l, const ImproperDihedralBond& improper, Float3* results, float& potE) {
	const Float3 ij_norm = (j - i).norm();
	const Float3 ik_norm = (k - i).norm();
	const Float3 il_norm = (l - i).norm();
	const Float3 lj_norm = (j - l).norm();
	const Float3 lk_norm = (k - l).norm();

	const Float3 plane_normal = (ij_norm).cross((ik_norm)).norm();	// i is usually the center on
	const Float3 plane2_normal = (lj_norm.cross(lk_norm)).norm();


	float angle = Float3::getAngleOfNormVectors(plane_normal, plane2_normal);
	const bool angle_is_negative = (plane_normal.dot(il_norm)) > 0.f;
	if (angle_is_negative) {
		angle = -angle;
	}

	const float error = angle - improper.psi_0;

	if constexpr (CALC_POTE) {
		potE = 0.5f * improper.k_psi * (error * error);
	}
	const float torque = improper.k_psi * (angle - improper.psi_0) * LIMA_TO_NANO;

	// This is the simple way, always right-ish
	results[3] = plane2_normal * (torque / l.distToLine(j, k));
	results[0] = -plane_normal * (torque / i.distToLine(j, k));

	const Float3 residual = -(results[0] + results[3]);
	const float ratio = (j - i).len() / ((j - i).len() + (k - i).len());
	results[1] = residual * (1.f-ratio);
	results[2] = residual * (ratio);


#if defined LIMASAFEMODE
	Float3 force_spillover = Float3{};
	for (int i = 0; i < 4; i++) {
		force_spillover += results[i];
	}
	if (force_spillover.len()*10000.f > results[0].len()) {
		force_spillover.print('s');
	}
	if (angle > PI || angle < -PI) {
		printf("Anlg too large!! %f\n\n\n\n", angle);
	}
	if (results[0].len() > 0.5f) {
		printf("\nImproperdihedralBond: angle %f [rad] torque: %f psi_0 [rad] %f k_psi %f\n",
			angle, torque, improper.psi_0, improper.k_psi);
	}
#endif
}


// ------------------------------------------------------------------------------------------- LJ Forces -------------------------------------------------------------------------------------------//
enum CalcLJOrigin { ComComIntra, ComComInter, ComSol, SolCom, SolSolIntra, SolSolInter };


__device__ static const char* calcLJOriginString[] = {
	"ComComIntra", "ComComInter", "ComSol", "SolCom", "SolSolIntra", "SolSolInter"
};

// This function does not add the 24 scalar, the caller fucntion must do so!
__device__ static Float3 calcLJForceOptim(const Float3& diff, const float dist_sq_reciprocal, float& potE, const float sigma, const float epsilon,
	CalcLJOrigin originSelect, /*For debug only*/
	int type1 = -1, int type2 = -1) {
	// Calculates LJ force on p0	(attractive to p1. Negative values = repulsion )//
	// input positions in cartesian coordinates [nm]
	// sigma [nm]
	// epsilon [J/mol]->[(kg*nm^2)/(ns^2*mol)]
	// Returns force in J/mol*M		?????????????!?!?//

	// Directly from book
	float s = (sigma * sigma) * dist_sq_reciprocal;								// [nm^2]/[nm^2] -> unitless	// OPTIM: Only calculate sigma_squared, since we never use just sigma
	s = s * s * s;
	const float force_scalar = epsilon * s * dist_sq_reciprocal * (1.f - 2.f * s);	// Attractive. Negative, when repulsive		[(kg*nm^2)/(nm^2*ns^2*mol)] ->----------------------	[(kg)/(ns^2*mol)]	

	const Float3 force = diff * force_scalar;

	if constexpr (CALC_POTE) {
		//potE += 4.f * epsilon * s * (s - 1.f) * 0.5f;	// 0.5 to account for 2 particles doing the same calculation
		potE += 2.f * epsilon * s * (s - 1.f);
	}
#if defined LIMASAFEMODE
	auto pot = 4. * epsilon * s * (s - 1.f) * 0.5;
	if (force.len() > 1.f || pot > 1e+8) {
		//printf("\nBlock %d thread %d\n", blockIdx.x, threadIdx.x);
		////((*pos1 - *pos0) * force_scalar).print('f');
		//pos0.print('0');
		//pos1.print('1');
		/*printf("\nLJ Force %s: dist nm %f force %f sigma %f epsilon %f t1 %d t2 %d\n",
			calcLJOriginString[(int)originSelect], sqrt(dist_sq) / NANO_TO_LIMA, ((pos1 - pos0) * force_scalar).len(), sigma / NANO_TO_LIMA, epsilon, type1, type2);*/
	}
#endif

	return force;	// GN/mol [(kg*nm)/(ns^2*mol)]
}


















__device__ void cudaAtomicAdd(Float3& target, const Float3& add) {
	atomicAdd(&target.x, add.x);
	atomicAdd(&target.y, add.y);
	atomicAdd(&target.z, add.z);
}

// ------------------------------------------------------------ Forcecalc handlers ------------------------------------------------------------ //

// only works if n threads >= n bonds
__device__ Float3 computeSinglebondForces(const SingleBond* const singlebonds, const int n_singlebonds, const Float3* const positions,
	Float3* const forces_interim, float* const potentials_interim, float* const potE, int bridgekernel)
{
	// First clear the buffer which will store the forces.
	forces_interim[threadIdx.x] = Float3(0.f);
	potentials_interim[threadIdx.x] = 0.f;
	__syncthreads();

	for (int bond_offset = 0; (bond_offset * blockDim.x) < n_singlebonds; bond_offset++) {
		const SingleBond* pb = nullptr;
		Float3 forces[2] = { Float3{}, Float3{} };
		float potential = 0.f;
		const int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < n_singlebonds) {
			pb = &singlebonds[bond_index];

			LimaForcecalc::calcSinglebondForces(
				positions[pb->atom_indexes[0]],
				positions[pb->atom_indexes[1]],
				*pb,
				forces,
				potential,
				bridgekernel
			);
		}

		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && pb != nullptr) {
				for (int i = 0; i < 2; i++) {
					forces_interim[pb->atom_indexes[i]] += forces[i];
					potentials_interim[pb->atom_indexes[i]] += potential * 0.5f;
				}
			}
			__syncthreads();
		}
	}

	*potE += potentials_interim[threadIdx.x];
	const Float3 force = forces_interim[threadIdx.x];
	__syncthreads();

	return force;
}


__device__ Float3 computeAnglebondForces(const AngleBond* const anglebonds, const int n_anglebonds, const Float3* const positions,
	Float3* const forces_interim, float* const potentials_interim, float* const potE)
{
	// First clear the buffer which will store the forces.
	forces_interim[threadIdx.x] = Float3(0.f);
	potentials_interim[threadIdx.x] = 0.f;
	__syncthreads();

	for (int bond_offset = 0; (bond_offset * blockDim.x) < n_anglebonds; bond_offset++) {
		const AngleBond* ab = nullptr;
		Float3 forces[3] = { Float3{}, Float3{}, Float3{} };
		float potential = 0.f;
		const int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < n_anglebonds) {
			ab = &anglebonds[bond_index];

			LimaForcecalc::calcAnglebondForces(
				positions[ab->atom_indexes[0]],
				positions[ab->atom_indexes[1]],
				positions[ab->atom_indexes[2]],
				*ab,
				forces,
				potential
			);
		}


		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && ab != nullptr) {
				for (int i = 0; i < ab->n_atoms; i++) {
					forces_interim[ab->atom_indexes[i]] += forces[i];
					potentials_interim[ab->atom_indexes[i]] += potential / 3.f;
				}
			}
			__syncthreads();
		}
	}

	*potE += potentials_interim[threadIdx.x];
	const Float3 force = forces_interim[threadIdx.x];
	__syncthreads();

	return force;
}


__device__ Float3 computeDihedralForces(const DihedralBond* const dihedrals, const int n_dihedrals, const Float3* const positions,
	Float3* const forces_interim, float* const potentials_interim, float* const potE)
{
	// First clear the buffer which will store the forces.
	forces_interim[threadIdx.x] = Float3(0.f);
	potentials_interim[threadIdx.x] = 0.f;
	__syncthreads();

	for (int bond_offset = 0; (bond_offset * blockDim.x) < n_dihedrals; bond_offset++) {
		const DihedralBond* db = nullptr;
		Float3 forces[4] = { Float3{}, Float3{}, Float3{}, Float3{} };
		float potential = 0.f;
		const int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < n_dihedrals) {
			db = &dihedrals[bond_index];
			LimaForcecalc::calcDihedralbondForces(
				positions[db->atom_indexes[0]] / NANO_TO_LIMA,
				positions[db->atom_indexes[1]] / NANO_TO_LIMA,
				positions[db->atom_indexes[2]] / NANO_TO_LIMA,
				positions[db->atom_indexes[3]] / NANO_TO_LIMA,
				*db,
				forces,
				potential
			);


			if constexpr (USE_ATOMICS_FOR_BONDS_RESULTS) {
				for (int i = 0; i < db->n_atoms; i++) {
					cudaAtomicAdd(forces_interim[db->atom_indexes[i]], forces[i]);
					atomicAdd(&potentials_interim[db->atom_indexes[i]], potential * 0.25f);
				}
			}
		}

		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && db != nullptr) {
				for (int i = 0; i < 4; i++) {
					forces_interim[db->atom_indexes[i]] += forces[i];
					potentials_interim[db->atom_indexes[i]] += potential * 0.25f;
				}
			}
			__syncthreads();
		}
	}

	*potE += potentials_interim[threadIdx.x];
	const Float3 force = forces_interim[threadIdx.x];
	__syncthreads();

	return force;
}

__device__ Float3 computeImproperdihedralForces(const ImproperDihedralBond* const impropers, const int n_impropers, const Float3* const positions,
	Float3* const forces_interim, float* const potentials_interim, float* const potE)
{
	__syncthreads();

	// First clear the buffer which will store the forces.
	forces_interim[threadIdx.x] = Float3(0.f);
	potentials_interim[threadIdx.x] = 0.f;
	__syncthreads();

	for (int bond_offset = 0; (bond_offset * blockDim.x) < n_impropers; bond_offset++) {
		const ImproperDihedralBond* db = nullptr;
		Float3 forces[4] = { Float3{}, Float3{}, Float3{}, Float3{} };
		float potential = 0.f;
		const int bond_index = threadIdx.x + bond_offset * blockDim.x;


		if (bond_index < n_impropers) {
			db = &impropers[bond_index];

			LimaForcecalc::calcImproperdihedralbondForces(
				positions[db->atom_indexes[0]] / NANO_TO_LIMA,
				positions[db->atom_indexes[1]] / NANO_TO_LIMA,
				positions[db->atom_indexes[2]] / NANO_TO_LIMA,
				positions[db->atom_indexes[3]] / NANO_TO_LIMA,
				*db,
				forces,
				potential
			);

			if constexpr (USE_ATOMICS_FOR_BONDS_RESULTS) {
				for (int i = 0; i < db->n_atoms; i++) {
					cudaAtomicAdd(forces_interim[db->atom_indexes[i]], forces[i]);
					atomicAdd(&potentials_interim[db->atom_indexes[i]], potential * 0.25f);
				}
			}
		}

		if constexpr (!USE_ATOMICS_FOR_BONDS_RESULTS) {
			for (int i = 0; i < blockDim.x; i++) {
				if (threadIdx.x == i && db != nullptr) {
					for (int i = 0; i < db->n_atoms; i++) {
						forces_interim[db->atom_indexes[i]] += forces[i];
						potentials_interim[db->atom_indexes[i]] += potential * 0.25f;
					}
				}
				__syncthreads();
			}
		}
	}

	*potE += potentials_interim[threadIdx.x];
	const Float3 force = forces_interim[threadIdx.x];
	__syncthreads();

	return force;
}



}	// End of namespace LimaForcecalc
