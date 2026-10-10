#pragma once

//#include <cuda_runtime.h>

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "EngineUtils.cuh"

#include "LennardJonesInteractions.cuh"


namespace LimaForcecalc 
{


template <bool energyMinimize>
__device__ inline void calcSinglebondForces(const Float3& p0, const Float3& p1, const SingleBond::Parameters& bondParams, Float3* results, float& potE, bool bridgekernel, int id0=-1, int id1=-1) {
	// Calculates bond force on both particles					
	// Calculates forces as J/mol*M								
	const Float3 difference = p0 - p1;						// [nm]
	const float error = difference.len() - bondParams.b0;				// [nm]

	if constexpr (ENABLE_POTE) {
		potE = 0.5f * bondParams.kb * (error * error);				// [J/mol]
	}
	float force_scalar = -bondParams.kb * error;				// [J/mol/nm]
	//printf("ForceScalar %f error %f\n", force_scalar, error);
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

	if constexpr (FORCE_CHECKS) {
		//printf("p0 %d %f %f %f p1 %d %f %f %f dist %f\n", id0, p0.x, p0.y, p0.z, id1, p1.x, p1.y, p1.z, difference.len());
		if (results[0].isNan()) {
			printf("Singlebond produces NAN force: kb %f b0 %f dist %f\n", bondParams.kb, bondParams.b0, difference.len());
		}
	}

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

template <bool energyMinimize>
__device__ inline void calcAnglebondForces(const Float3& pos_left, const Float3& pos_middle, const Float3& pos_right, const AngleUreyBradleyBond& angletype, Float3* results, float& potE) {
	const Float3 v1 = (pos_left - pos_middle).norm();
	const Float3 v2 = (pos_right - pos_middle).norm();
	Float3 normal = v1.cross(v2).norm();	// Poiting towards y, when right is pointing toward x
	if (energyMinimize && normal.lenSquared() == 0.f) {
		// Collinear: the bending plane is undefined, so the bending force is zero, matching GROMACS.
		// EM must still escape a straight angle away from theta0, so pick any plane containing v1.
		const bool alongX = fabsf(v1.x) >= 0.9f;
		normal = v1.cross(Float3{ alongX ? 0.f : 1.f, alongX ? 1.f : 0.f, 0.f }).norm();
	}

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


		if constexpr (FORCE_CHECKS) {
			if (results[0].isNan()) {
				printf("Anglebond produces NAN force: kTheta %f theta0 %f distance %f\n", angletype.params.kTheta, angletype.params.theta0, difference.len());
			}
		}
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








// ------------------------------------------------------------ Forcecalc handlers ------------------------------------------------------------ //

// Accumulates the per-particle bond results of a bondgroup in shared memory. Both are deterministic:
// Serial: float sums, added one thread at a time in a fixed order. Used in EM, where forces may exceed the fixed point range
struct SerialBondAccumulator {
	static constexpr bool parallel = false;
	static constexpr int nThreads = 32;
	__device__ static int ThreadIndex() { return threadIdx.x & 31; }
	__device__ static void Sync() { __syncwarp(); }
	float4* feInterrims;

	__device__ void Init(int index) const { feInterrims[index] = float4{ 0, 0, 0, 0 }; }
	__device__ void Add(int index, const Float3& force, float potE) const {
		feInterrims[index] = ::Add(feInterrims[index], make_float4(force.x, force.y, force.z, potE));
	}
	__device__ ForceEnergy Get(int index) const {
		const float4 v = feInterrims[index];
		return ForceEnergy{ Float3{ v.x, v.y, v.z }, v.w };
	}
};
// Parallel: 64-bit fixed point summed with shared memory integer atomics, which are order independent, so all threads can add at once.
// Same scale as ForceAccumulator: resolution 6e-8, range +-5.5e11
template <bool withPotE>
struct FixedPointBondAccumulatorT {
	static constexpr bool parallel = true;
	static constexpr int nThreads = 32;
	__device__ static int ThreadIndex() { return threadIdx.x & 31; }
	__device__ static void Sync() { __syncwarp(); }
	static constexpr float scale = 16777216.f;		// 2^24
	static constexpr float scaleInv = 1.f / scale;	// Power of 2, so Get is exactly the rounded sum
	unsigned int* words; // [low components][high components], each component has stride particles
	int stride;

	__device__ static unsigned long long ToFixed(float v) { return static_cast<unsigned long long>(llrintf(v * scale)); }
	__device__ static float ToFloat(unsigned long long v) { return static_cast<float>(static_cast<long long>(v)) * scaleInv; }

	__device__ void Init(int index) const {
		for (int i = 0; i < 3 + withPotE; i++) {
			words[i * stride + index] = 0;
			words[(3 + withPotE + i) * stride + index] = 0;
		}
	}
	// Shared 64-bit atomicAdd uses a CAS loop. Two native 32-bit atomics sum the same bits:
	// each low-word wrap contributes one carry to the high word, regardless of arrival order.
	// Readers must wait for all additions (ScatterBondResults supplies the warp barrier).
	__device__ void AddFixed(int index, unsigned long long value) const {
		const unsigned int low = static_cast<unsigned int>(value);
		const unsigned int previous = atomicAdd(&words[index], low);
		atomicAdd(&words[(3 + withPotE) * stride + index], static_cast<unsigned int>(value >> 32) + (previous > UINT_MAX - low));
	}
	__device__ unsigned long long GetFixed(int index) const {
		return static_cast<unsigned long long>(words[index]) | static_cast<unsigned long long>(words[(3 + withPotE) * stride + index]) << 32;
	}
	__device__ void Add(int index, const Float3& force, float potE) const {
		AddFixed(index, ToFixed(force.x));
		AddFixed(stride + index, ToFixed(force.y));
		AddFixed(2 * stride + index, ToFixed(force.z));
		if constexpr (withPotE) AddFixed(3 * stride + index, ToFixed(potE));
	}
	__device__ ForceEnergy Get(int index) const {
		return ForceEnergy{ Float3{ ToFloat(GetFixed(index)), ToFloat(GetFixed(stride + index)), ToFloat(GetFixed(2 * stride + index)) }, withPotE ? ToFloat(GetFixed(3 * stride + index)) : 0.f };
	}
	template <typename ForceAccumulator>
	__device__ void StorePrimary(const ForceAccumulator& target, int index, int slot) const {
		target.fx[slot] = GetFixed(index);
		target.fy[slot] = GetFixed(stride + index);
		target.fz[slot] = GetFixed(2 * stride + index);
		if constexpr (withPotE) target.potE[slot] = GetFixed(3 * stride + index);
	}
	// One group owns each output entry, so global stores need no atomic read-modify-write.
	__device__ void StoreTo(ulonglong4* const target, int index, int outputIndex) const {
		target[outputIndex] = make_ulonglong4(GetFixed(index), GetFixed(stride + index), GetFixed(2 * stride + index), withPotE ? GetFixed(3 * stride + index) : 0);
	}
};
using FixedPointBondAccumulator = FixedPointBondAccumulatorT<true>;

// Adds each thread's bond results to the accumulator, then syncs so the bonds buffer may be reused
template <typename Accumulator, typename AddResults>
__device__ inline void ScatterBondResults(const Accumulator& acc, bool hasBond, AddResults addResults) {
	if constexpr (Accumulator::parallel) {
		if (hasBond)
			addResults();
		acc.Sync();
	}
	else {
		for (int tid = 0; tid < Accumulator::nThreads; tid++) {
			if (acc.ThreadIndex() == tid && hasBond)
				addResults();
			acc.Sync();
		}
	}
}

// A bondgroup may contain several molecules (small molecules are packed together), which can be far apart.
// So periodic boundaries are applied per bond, placing each atom at the image nearest the bond's first atom
template <typename BoundaryCondition, int n>
__device__ inline void LoadBondPositions(const Float3* const positions, const uint8_t(&ids)[n], Float3(&out)[n], const Float3& boxSize, const Float3& boxSizeInv) {
	out[0] = positions[ids[0]];
	for (int i = 1; i < n; i++) {
		out[i] = positions[ids[i]];
		BoundaryCondition::ApplyHyperpos(out[0], out[i], boxSize, boxSizeInv);
	}
}

// only works if n threads >= n bonds
template<typename BoundaryCondition, typename Accumulator, bool energyMinimization>
__device__ inline void computeSinglebondForces(const SingleBond* const singlebonds, const int n_singlebonds, const Float3* const positions,	const Accumulator& acc, int bridgekernel,
	const Float3& boxSize, const Float3& boxSizeInv)
{
	for (int bond_offset = 0; (bond_offset * Accumulator::nThreads) < n_singlebonds; bond_offset++) {
		const SingleBond* pb = nullptr;
		Float3 forces[2] = { Float3{}, Float3{} };
		float potential = 0.f;
		const int bond_index = acc.ThreadIndex() + bond_offset * Accumulator::nThreads;

		if (bond_index < n_singlebonds) {
			pb = &singlebonds[bond_index];

			Float3 pos[SingleBond::nAtoms];
			LoadBondPositions<BoundaryCondition>(positions, pb->idInBondgroup, pos, boxSize, boxSizeInv);
			LimaForcecalc::calcSinglebondForces<energyMinimization>(
				pos[0],
				pos[1],
				pb->params,
				forces,
				potential,
				bridgekernel,
				(int)pb->idInBondgroup[0],
				(int)pb->idInBondgroup[1]
			);
		}

		ScatterBondResults(acc, pb != nullptr, [&] {
				acc.Add(pb->idInBondgroup[0], forces[0], potential * 0.5f);
				acc.Add(pb->idInBondgroup[1], forces[1], potential * 0.5f);
			});
	}
}

template<typename BoundaryCondition, typename Accumulator>
__device__ inline void computePairbondForces(const PairBond* const pairbonds, const int n_pairbonds, const Float3* const positions,	const Accumulator& acc,
	const Float3& boxSize, const Float3& boxSizeInv)
{
	for (int bond_offset = 0; (bond_offset * Accumulator::nThreads) < n_pairbonds; bond_offset++) {
		const PairBond* pb = nullptr;
		Float3 forces[2] = { Float3{}, Float3{} };
		float potential = 0.f;
		const int bond_index = acc.ThreadIndex() + bond_offset * Accumulator::nThreads;

		if (bond_index < n_pairbonds) {
			pb = &pairbonds[bond_index];

			Float3 pos[PairBond::nAtoms];
			LoadBondPositions<BoundaryCondition>(positions, pb->atom_indexes, pos, boxSize, boxSizeInv);
			const Float3 diff = pos[1] - pos[0];
			const float distSqReciprocal = 1.f / diff.lenSquared();

			const Float3 forceOnLeft = LJ::calcLJForceOptim<true, false>(diff, distSqReciprocal, potential, pb->params.sigma, pb->params.epsilon, LJ::CalcLJOrigin::Pairbond) * 24.f;
			forces[0] = forceOnLeft;
			forces[1] = -forceOnLeft;
		}

		ScatterBondResults(acc, pb != nullptr, [&] {
				acc.Add(pb->atom_indexes[0], forces[0], potential); // No *0.5f here, since LJ computes the pot per atom already;
				acc.Add(pb->atom_indexes[1], forces[1], potential);
			});
	}
}

template<typename BoundaryCondition, typename Accumulator, bool energyMinimization>
__device__ inline void computeAnglebondForces(const AngleUreyBradleyBond* const anglebonds, const int n_anglebonds, const Float3* const positions, const Accumulator& acc,
	const Float3& boxSize, const Float3& boxSizeInv)
{
	for (int bond_offset = 0; (bond_offset * Accumulator::nThreads) < n_anglebonds; bond_offset++) {
		const AngleUreyBradleyBond* ab = nullptr;
		Float3 forces[3] = { Float3{}, Float3{}, Float3{} };
		float potential = 0.f;
		const int bond_index = acc.ThreadIndex() + bond_offset * Accumulator::nThreads;

		if (bond_index < n_anglebonds) {
			ab = &anglebonds[bond_index];

			Float3 pos[AngleUreyBradleyBond::nAtoms];
			LoadBondPositions<BoundaryCondition>(positions, ab->atom_indexes, pos, boxSize, boxSizeInv);
			LimaForcecalc::calcAnglebondForces<energyMinimization>(
				pos[0],
				pos[1],
				pos[2],
				*ab,
				forces,
				potential
			);
		}


		ScatterBondResults(acc, ab != nullptr, [&] {
				acc.Add(ab->atom_indexes[0], forces[0], potential / 3.f);
				acc.Add(ab->atom_indexes[1], forces[1], potential / 3.f);
				acc.Add(ab->atom_indexes[2], forces[2], potential / 3.f);
			});
	}
}


template<typename BoundaryCondition, typename Accumulator>
__device__ inline void computeDihedralForces(const DihedralBond* const dihedrals, const int n_dihedrals, const Float3* const positions,	const Accumulator& acc,
	const Float3& boxSize, const Float3& boxSizeInv)
{
	for (int bond_offset = 0; (bond_offset * Accumulator::nThreads) < n_dihedrals; bond_offset++) {
		const DihedralBond* db = nullptr;
		Float3 forces[4] = { Float3{}, Float3{}, Float3{}, Float3{} };
		float potential = 0.f;
		const int bond_index = acc.ThreadIndex() + bond_offset * Accumulator::nThreads;

		if (bond_index < n_dihedrals) {
			db = &dihedrals[bond_index];
			Float3 pos[DihedralBond::nAtoms];
			LoadBondPositions<BoundaryCondition>(positions, db->atom_indexes, pos, boxSize, boxSizeInv);
			LimaForcecalc::calcDihedralbondForces(
				pos[0],
				pos[1],
				pos[2],
				pos[3],
				*db,
				forces,
				potential
			);
		}

		ScatterBondResults(acc, db != nullptr, [&] {
				acc.Add(db->atom_indexes[0], forces[0], potential * 0.25f);
				acc.Add(db->atom_indexes[1], forces[1], potential * 0.25f);
				acc.Add(db->atom_indexes[2], forces[2], potential * 0.25f);
				acc.Add(db->atom_indexes[3], forces[3], potential * 0.25f);
			});
	}
}

template<typename BoundaryCondition, typename Accumulator>
__device__ inline void computeImproperdihedralForces(const ImproperDihedralBond* const impropers, const int n_impropers, const Float3* const positions,	const Accumulator& acc,
	const Float3& boxSize, const Float3& boxSizeInv)
{
	for (int bond_offset = 0; (bond_offset * Accumulator::nThreads) < n_impropers; bond_offset++) {
		const ImproperDihedralBond* db = nullptr;
		Float3 forces[4] = { Float3{}, Float3{}, Float3{}, Float3{} };
		float potential = 0.f;
		const int bond_index = acc.ThreadIndex() + bond_offset * Accumulator::nThreads;


		if (bond_index < n_impropers) {
			db = &impropers[bond_index];

			Float3 pos[ImproperDihedralBond::nAtoms];
			LoadBondPositions<BoundaryCondition>(positions, db->atom_indexes, pos, boxSize, boxSizeInv);
			LimaForcecalc::calcImproperdihedralbondForces(
				pos[0],
				pos[1],
				pos[2],
				pos[3],
				*db,
				forces,
				potential
			);
		}

		ScatterBondResults(acc, db != nullptr, [&] {
				acc.Add(db->atom_indexes[0], forces[0], potential * 0.25f);
				acc.Add(db->atom_indexes[1], forces[1], potential * 0.25f);
				acc.Add(db->atom_indexes[2], forces[2], potential * 0.25f);
				acc.Add(db->atom_indexes[3], forces[3], potential * 0.25f);
			});
	}
}



}	// End of namespace LimaForcecalc
