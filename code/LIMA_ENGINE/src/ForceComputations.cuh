#pragma once

//#include <cuda_runtime.h>

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "EngineUtils.cuh"


namespace LimaForcecalc 
{


template <bool energyMinimize>
__device__ inline void calcSinglebondForces(const Float3& pos_a, const Float3& pos_b, const SingleBond::Parameters& bondParams, Float3* results, float& potE, bool bridgekernel) {
	// Calculates bond force on both particles					
	// Calculates forces as J/mol*M								
	const Float3 difference = pos_a - pos_b;						// [nm]
	const float error = difference.len() - bondParams.b0;				// [nm]

	if constexpr (ENABLE_POTE) {
		potE = 0.5f * bondParams.kb * (error * error);				// [J/mol]
	}
	float force_scalar = -bondParams.kb * error;				// [J/mol/nm]

	// In EM mode we might have some VERY long bonds, to avoid explosions, we cap the error used to calculate force to 2*b0
	// Note that we still get the correct value for potE
	if constexpr (energyMinimize) {
		if (error > bondParams.b0 * 2.f) {
			force_scalar = -bondParams.kb * bondParams.b0 * 2.f;
		}
	}

	const Float3 dir = difference.norm();							// dif_unit_vec, but shares variable with dif
	results[0] = dir * force_scalar;								// [kg * nm / (mol*ls^2)] = [1/n N]
	results[1] = -dir * force_scalar;								// [kg * nm / (mol*ls^2)] = [1/n N]

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

__device__ inline void calcAnglebondForces(const Float3& pos_left, const Float3& pos_middle, const Float3& pos_right, const AngleUreyBradleyBond& angletype, Float3* results, float& potE) {
	const Float3 v1 = (pos_left - pos_middle).norm();
	const Float3 v2 = (pos_right - pos_middle).norm();
	const Float3 normal = v1.cross(v2).norm();	// Poiting towards y, when right is pointing toward x

	const Float3 inward_force_direction1 = (v1.cross(normal * -1.f)).norm();
	const Float3 inward_force_direction2 = (v2.cross(normal)).norm();

	const float angle = Float3::getAngleOfNormVectors(v1, v2);
	const float error = angle - angletype.params.theta0;				// [rad]

	// Simple implementation
	if constexpr (ENABLE_POTE) {
		potE = angletype.params.kTheta * error * error * 0.5f;		// Energy [J/mol]
	}
	const float torque = angletype.params.kTheta * (error);				// Torque [J/(mol*rad)]

	// Correct implementation
	//potE = -angletype.k_theta * (cosf(error) - 1.f);		// Energy [J/mol]
	//const float torque = angletype.k_theta * sinf(error);	// Torque [J/(mol*rad)]

	results[0] = inward_force_direction1 * (torque / (pos_left - pos_middle).len());
	results[2] = inward_force_direction2 * (torque / (pos_right - pos_middle).len());
	results[1] = (results[0] + results[2]) * -1.f;

	// UreyBradley potential
	if constexpr (ENABLE_UREYBRADLEY) {
		const Float3 difference = pos_left - pos_right;						// [nm]
		const float error = difference.len() - angletype.params.ub0;		// [nm]

		if constexpr (ENABLE_POTE) {
			potE += 0.5f * angletype.params.kUB * (error * error);			// [J/mol]
		}
		float force_scalar = -angletype.params.kUB * error;					// [J/mol/nm] 

		// In EM mode we might have some VERY long bonds, to avoid explosions, we cap the error used to calculate force to 2*b0
		// Note that we still get the correct value for potE
	/*	if constexpr (energyMinimize) {
			if (error > bondParams.b0 * 2.f) {
				force_scalar = -bondParams.kb * bondParams.b0 * 2.f;
			}
		}*/

		//printf("UB: %f Angleforce %f\n", force_scalar, results[0].len());

		const Float3 dir = difference.norm();							// dif_unit_vec, but shares variable with dif

		//printf("%f %f %f %f %f\n", force_scalar, error, angletype.params.kUB, angletype.params.ub0, results[0].len());


		results[0] += dir * force_scalar;								// [J/mol/nm] = [1/lima N]
		results[2] += -dir * force_scalar;


		//results[0] += -inward_force_direction1 * force_scalar;
		//results[2] += -inward_force_direction2 * force_scalar;
	}
#if defined LIMASAFEMODE
	if (results[0].len() > 0.1f) {
		printf("\nAngleBond: angle %f [rad] error %f [rad] force %f t0 %f [rad] kt %f\n", angle, error, results[0].len(), angletype.theta_0, angletype.k_theta);
	}
#endif
}

// From resource: https://nosarthur.github.io/free%20energy%20perturbation/2017/02/01/dihedral-force.html
// Greatly inspired by OpenMD's CharmmDihedral algorithm
__device__ inline void calcDihedralbondForces(const Float3& pos_left, const Float3& pos_lm, const Float3& pos_rm, const Float3& pos_right, 
	const DihedralBond& dihedral, Float3* results, float& potE) {
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

	//if constexpr (ENABLE_POTE) {
	//	potE = __half2float(dihedral.params.k_phi) * (1. + cos(__half2float(dihedral.params.n) * torsion - __half2float(dihedral.params.phi_0)));
	//}
	//const float torque = __half2float(dihedral.params.k_phi) * (__half2float(dihedral.params.n) * sin(__half2float(dihedral.params.n) * torsion 
	//	- __half2float(dihedral.params.phi_0))) / NANO_TO_LIMA;

	if constexpr (ENABLE_POTE) {
		potE = dihedral.params.k_phi * (1. + cos(dihedral.params.n * torsion - dihedral.params.phi_0));
	}
	const float torque = dihedral.params.k_phi * (dihedral.params.n * sin(dihedral.params.n * torsion - dihedral.params.phi_0));

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
__device__ inline void calcImproperdihedralbondForces(const Float3& i, const Float3& j, const Float3& k, const Float3& l, const ImproperDihedralBond& improper, Float3* results, float& potE) {
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

	const float error = angle - improper.params.psi_0;

	if constexpr (ENABLE_POTE) {
		potE = 0.5f * improper.params.k_psi * (error * error);
	}
	const float torque = improper.params.k_psi * (angle - improper.params.psi_0);

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












__device__ inline void cudaAtomicAdd(Float3& target, const Float3& add) {
	atomicAdd(&target.x, add.x);
	atomicAdd(&target.y, add.y);
	atomicAdd(&target.z, add.z);
}

// ------------------------------------------------------------ Forcecalc handlers ------------------------------------------------------------ //

// only works if n threads >= n bonds
template<bool energyMinimization>
__device__ inline Float3 computeSinglebondForces(const SingleBond* const singlebonds, const int n_singlebonds, const Float3* const positions,
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

			LimaForcecalc::calcSinglebondForces<energyMinimization>(
				positions[pb->atom_indexes[0]],
				positions[pb->atom_indexes[1]],
				pb->params,
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


__device__ inline Float3 computeAnglebondForces(const AngleUreyBradleyBond* const anglebonds, const int n_anglebonds, const Float3* const positions,
	Float3* const forces_interim, float* const potentials_interim, float* const potE)
{
	// First clear the buffer which will store the forces.
	forces_interim[threadIdx.x] = Float3(0.f);
	potentials_interim[threadIdx.x] = 0.f;
	__syncthreads();

	for (int bond_offset = 0; (bond_offset * blockDim.x) < n_anglebonds; bond_offset++) {
		const AngleUreyBradleyBond* ab = nullptr;
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
				for (int i = 0; i < ab->nAtoms; i++) {
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


__device__ inline Float3 computeDihedralForces(const DihedralBond* const dihedrals, const int n_dihedrals, const Float3* const positions,
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
				positions[db->atom_indexes[0]],
				positions[db->atom_indexes[1]],
				positions[db->atom_indexes[2]],
				positions[db->atom_indexes[3]],
				*db,
				forces,
				potential
			);


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

	//if (threadIdx.x == 0 )
	//	for (int i = 0; i < blockDim.x; i++) 
	//		printf("%f\n", potentials_interim[i]);

	return force;
}

__device__ inline Float3 computeImproperdihedralForces(const ImproperDihedralBond* const impropers, const int n_impropers, const Float3* const positions,
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
				positions[db->atom_indexes[0]],
				positions[db->atom_indexes[1]],
				positions[db->atom_indexes[2]],
				positions[db->atom_indexes[3]],
				*db,
				forces,
				potential
			);

			if constexpr (USE_ATOMICS_FOR_BONDS_RESULTS) {
				for (int i = 0; i < db->nAtoms; i++) {
					cudaAtomicAdd(forces_interim[db->atom_indexes[i]], forces[i]);
					atomicAdd(&potentials_interim[db->atom_indexes[i]], potential * 0.25f);
				}
			}
		}

		if constexpr (!USE_ATOMICS_FOR_BONDS_RESULTS) {
			for (int i = 0; i < blockDim.x; i++) {
				if (threadIdx.x == i && db != nullptr) {
					for (int i = 0; i < db->nAtoms; i++) {
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
