//#pragma once - this file must NOT be included multiple times

#include "Engine.cuh"
#include "ForceComputations.cuh"
#include "KernelWarnings.cuh"
#include "EngineUtils.cuh"

#include "BoundaryCondition.cuh"
#include "DeviceAlgorithms.cuh"
#include "Utilities.h"

//#include <cuda/pipeline>
#include "LennardJonesInteractions.cuh"



















// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//


template <typename BoundaryCondition, bool energyMinimize>
__global__ void PclusterSnfKernel(const PersistentCluster* const pc, const PersistentClusterMeta* const pcMeta, UniformElectricField electricField, ForceEnergy* const forceEnergy, int pclusterOffset, int nPclusters) {

	const int workId = blockIdx.x * blockDim.x + threadIdx.x;
	if (workId >= nPclusters) return;
	const int pcId = pclusterOffset + workId;

	
	for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
		const int pidGlobal = pcMeta[pcId].particleIdsGlobal[pid];
		if (pidGlobal == -1)
			continue;

		float charge = pc[pcId].pqd[pid].params.charge;
		Float3 force = electricField.GetForce(charge);

		forceEnergy[pcId * PersistentCluster::maxParticles + pid] = ForceEnergy{ force, 0.f };
	}
}

__global__ void ElasticPositionsForceKernel(const PersistentCluster* const pc, const PersistentClusterMeta* const pcMeta, const Float3* const elasticPositions, ForceEnergy* const forceEnergy, int pclusterOffset, int nPclusters, Float3 boxSize) {
	const int workId = blockIdx.x * blockDim.x + threadIdx.x;
	if (workId >= nPclusters) return;
	const int pcId = pclusterOffset + workId;


	for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
		const int pidGlobal = pcMeta[pcId].particleIdsGlobal[pid];
		if (pidGlobal == -1)
			continue;

		const float mass = pcMeta[pcId].mass[pid];
		const Float3 position = pc[pcId].pqd[pid].position;
		const Float3 ep = elasticPositions[pidGlobal];

		Float3 elasticPosition{
			isnan(ep.x) ? position.x : ep.x,
			isnan(ep.y) ? position.y : ep.y,
			isnan(ep.z) ? position.z : ep.z
		};
		PeriodicBoundaryCondition::applyHyperposNM(position, elasticPosition, boxSize);

		const Float3 difference = elasticPosition - position;
		const float distSq = difference.lenSquared();
		const float dist = difference.len();
		const Float3 forceDirection = distSq < 0.00001f ? Float3{ 0.f } : difference.norm();

		
		// Magnitude =  Coeff * (e^(wx^2) - 1) / (e^(wx^2) + 1)  // Coeff controls magnitude, w controls gradient
		const float coefficient = 1000000.f; // [kJ/mol/nm]
		const float w = 6;
		float eTerm = expf(w * dist);
		const float forceMagnitude = mass * coefficient * (eTerm - 1.0f) / (eTerm + 1.0f); 

		const float potentialEnergy = logf(eTerm + 1.0f) * coefficient;

		forceEnergy[pcId * PersistentCluster::maxParticles + pid] = ForceEnergy{ forceDirection * forceMagnitude, potentialEnergy };
	}
}





static const int THREADS_PER_BONDSGROUPSKERNEL = 64;
template <typename BoundaryCondition, bool emVariant>
__global__ void BondgroupsKernel(const BondGroupsDevice bondGroups, const BoxState boxState, ForceEnergy* const forceEnergiesOut, const PersistentCluster* const pclusters, Float3 boxSize, Float3 boxSizeInv) {
	__shared__ Float3 positions[THREADS_PER_BONDSGROUPSKERNEL];

	// MD accumulates in parallel with deterministic fixed point atomics, EM serially in float, since its forces may exceed the fixed point range
	using Accumulator = std::conditional_t<emVariant, LimaForcecalc::SerialBondAccumulator, LimaForcecalc::FixedPointBondAccumulator>;
	__shared__ unsigned long long accumulatorBuffer[4 * THREADS_PER_BONDSGROUPSKERNEL]; // Large enough for either
	Accumulator acc;
	if constexpr (emVariant) acc = Accumulator{ reinterpret_cast<float4*>(accumulatorBuffer) };
	else acc = Accumulator{ accumulatorBuffer, THREADS_PER_BONDSGROUPSKERNEL };

	static const int batchSize = THREADS_PER_BONDSGROUPSKERNEL;
	static const int largestBondBytesize = std::max(sizeof(AngleUreyBradleyBond), sizeof(DihedralBond));
	__shared__ char _bondsBuffer[largestBondBytesize * batchSize];	

	const BondGroup* const bondGroup = &bondGroups.groups[blockIdx.x];

	acc.Init(threadIdx.x);

	// Fetch positions. Periodic boundaries are applied per bond, since a group may contain several molecules far apart
	if (threadIdx.x < bondGroup->nParticles) {
		const BondGroup::ParticleRef pRef = bondGroups.particles[bondGroup->indexOfFirstParticle + threadIdx.x];
		positions[threadIdx.x] = pclusters[pRef.pcid].pqd[pRef.pid].position;  //boxState.compoundsRelposNm[pRef.compoundId * MAX_COMPOUND_PARTICLES + pRef.localIdInCompound] + relShift;
	}
	__syncthreads();

	
	{
		__syncthreads();
		SingleBond* bondsBuffer = reinterpret_cast<SingleBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nSinglebonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nSinglebonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroups.singlebonds[bondGroup->indexOfFirstSinglebond + bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeSinglebondForces<BoundaryCondition, Accumulator, emVariant>(bondsBuffer, std::min(batchSize, bondGroup->nSinglebonds - batchStart), positions, acc, 0, boxSize, boxSizeInv);
		}
	}

	{
		__syncthreads();
		AngleUreyBradleyBond* bondsBuffer = reinterpret_cast<AngleUreyBradleyBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nAnglebonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nAnglebonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroups.anglebonds[bondGroup->indexOfFirstAnglebond + bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeAnglebondForces<BoundaryCondition, Accumulator, emVariant>(bondsBuffer, std::min(batchSize, bondGroup->nAnglebonds - batchStart), positions, acc, boxSize, boxSizeInv);
		}
	}

	{
		__syncthreads();
		DihedralBond* bondsBuffer = reinterpret_cast<DihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nDihedralbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nDihedralbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroups.dihedralbonds[bondGroup->indexOfFirstDihedralbond + bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeDihedralForces<BoundaryCondition, Accumulator>(bondsBuffer, std::min(batchSize, bondGroup->nDihedralbonds - batchStart), positions, acc, boxSize, boxSizeInv);
		}
	}

	{
		__syncthreads();
		ImproperDihedralBond* bondsBuffer = reinterpret_cast<ImproperDihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nImproperdihedralbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nImproperdihedralbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroups.improperdihedralbonds[bondGroup->indexOfFirstImproperdihedralbond + bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeImproperdihedralForces<BoundaryCondition, Accumulator>(bondsBuffer, std::min(batchSize, bondGroup->nImproperdihedralbonds - batchStart), positions, acc, boxSize, boxSizeInv);
		}
	}


	{
		__syncthreads();
		// TODO: i have no clue if pairbonds should also compute SR electrostatics?
		PairBond* bondsBuffer = reinterpret_cast<PairBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nPairbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nPairbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroups.pairbonds[bondGroup->indexOfFirstPairbond + bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computePairbondForces<BoundaryCondition, Accumulator>(bondsBuffer, std::min(batchSize, bondGroup->nPairbonds - batchStart), positions, acc, boxSize, boxSizeInv);
		}
	}

	if (threadIdx.x < bondGroup->nParticles)
		forceEnergiesOut[bondGroup->indexOfFirstParticle + threadIdx.x] = acc.Get(threadIdx.x);
}



// Deterministic accumulation of NB forces in MD: each partial force is converted to 64-bit fixed point
// and summed with integer atomics. Unlike float atomics, integer addition is associative, so the sum is bitwise
// independent of the order the atomics arrive in. The accumulator is small enough to stay L2 resident.
// Not used in EM, where forces can exceed the fixed point range.
struct NbForceAccumulator {
	static constexpr float scale = 16777216.f;		// 2^24 -> resolution 6e-8 J/mol/nm, range +-5.5e11 J/mol/nm
	static constexpr float scaleInv = 1.f / scale;	// Power of 2, so ToFloat is exactly the rounded sum

	// SoA, indexed by scId * SuperCluster::maxParticles + particleIndex. potE is only zeroed/used on logging steps
	unsigned long long* fx = nullptr;
	unsigned long long* fy = nullptr;
	unsigned long long* fz = nullptr;
	unsigned long long* potE = nullptr;

	__device__ static unsigned long long ToFixed(float v) { return static_cast<unsigned long long>(llrintf(v * scale)); }
	__device__ static float ToFloat(unsigned long long v) { return static_cast<float>(static_cast<long long>(v)) * scaleInv; }

	template <bool withPotE>
	__device__ void Add(int index, const ForceEnergy& fe) const {
		atomicAdd(&fx[index], ToFixed(fe.force.x));
		atomicAdd(&fy[index], ToFixed(fe.force.y));
		atomicAdd(&fz[index], ToFixed(fe.force.z));
		if constexpr (withPotE)
			atomicAdd(&potE[index], ToFixed(fe.potE));
	}

	template <bool withPotE>
	__device__ ForceEnergy Load(int index) const {
		return ForceEnergy{ Float3{ ToFloat(fx[index]), ToFloat(fy[index]), ToFloat(fz[index]) }, withPotE ? ToFloat(potE[index]) : 0.f };
	}
};

namespace NbNonlocal {
	// This lane's row of an entry's noInteractions with an own quarter, rotated so bit k is j particle (lane+k)&3. All set if
	// the entry does not interact with the quarter
	__device__ inline uint32_t RotatedNoInteractions(uint32_t noInteractions, uint32_t ownQuarterMask, int ownQuarter, int lane) {
		const uint32_t row = (noInteractions >> (lane * 4)) & 0xFu;
		return ((ownQuarterMask >> ownQuarter) & 1) ? ((row | (row << 4)) >> lane) & 0xFu : 0xFu;
	}

	// entry.noInteractions[ownQuarter], for a runtime ownQuarter. Indexing the array would put entry in local memory
	__device__ inline uint32_t NoInteractionsOfQuarter(const QuarterEntry& entry, int ownQuarter) {
		const uint64_t all = uint64_t(entry.noInteractions[0]) | uint64_t(entry.noInteractions[1]) << 16
			| uint64_t(entry.noInteractions[2]) << 32 | uint64_t(entry.noInteractions[3]) << 48;
		return static_cast<uint32_t>(all >> (ownQuarter * 16)) & 0xFFFFu;
	}

	template <bool withPotE>
	__device__ inline ForceEnergy ShuffleForceEnergy(const ForceEnergy& fe, int srcLane) {
		ForceEnergy out{};
		out.force.x = __shfl_sync(0xFFFFFFFFu, fe.force.x, srcLane, 4);
		out.force.y = __shfl_sync(0xFFFFFFFFu, fe.force.y, srcLane, 4);
		out.force.z = __shfl_sync(0xFFFFFFFFu, fe.force.z, srcLane, 4);
		if constexpr (withPotE)
			out.potE = __shfl_sync(0xFFFFFFFFu, fe.potE, srcLane, 4);
		return out;
	}

	// accOwn[quarter] += fe, with selects, as a runtime index would move accOwn to local memory
	template <bool withPotE>
	__device__ inline void AddToQuarter(ForceEnergy (&accOwn)[4], int quarter, const ForceEnergy& fe) {
#pragma unroll
		for (int q = 0; q < 4; q++) {
			const bool hit = q == quarter;
			accOwn[q].force.x += hit ? fe.force.x : 0.f;
			accOwn[q].force.y += hit ? fe.force.y : 0.f;
			accOwn[q].force.z += hit ? fe.force.z : 0.f;
			if constexpr (withPotE)
				accOwn[q].potE += hit ? fe.potE : 0.f;
		}
	}
}

// Nonbonded forces. Each warp computes the QuarterEntries of one supercluster, see EmitQuarterEntriesKernel.
// blockDim = 64, 2 superclusters per block. computePotE must match logData of the following SuperclusterIntegrateKernel.
//
// The warp's 8 groups of 4 lanes each take one entry: lane l stages particle l of the entry's j quarter in shared memory,
// then for every own quarter in the entry's mask, holds own particle l of that quarter, and in iteration k pairs it with
// j particle (l+k)&3. Own forces thus accumulate in the lane that owns the particle for the whole task, and in MD the j forces
// in accJ[k] only need routing to their lanes once per entry. Loading the j quarter and flushing its forces once per entry
// rather than once per 4x4 block matters: the kernel is close to L2 bound.
// Entries are sorted by ownQuarterMask, so the groups of a warp mostly share a mask, and the warp loops over their union.
// Pairs beyond the cutoff are masked, so the forces do not depend on how the task builder grouped the pairs.
//
// MD sums the forces with NbForceAccumulator. EM forces can exceed its range, so EM stores each entry's j forces in the entry's
// quarter of the SCResult of its (own, j) supercluster pair, and the own forces in the own supercluster's first SCResult.
// SuperclusterIntegrateKernel sums them in a fixed order. In EM every quarter of such a pair has an entry, also quarters without
// pairs in range, so every result is written each step. Entries of the own supercluster add their j forces to the own forces.
//
// __launch_bounds__: 64 registers avoid spills (16 blocks per SM), the potE variants need 80 (12 blocks per SM).
// After larger changes, check ncu for local memory traffic and revisit the bounds.
template <typename BoundaryCondition, bool energyMinimize, bool computePotE>
__global__ void __launch_bounds__(64, computePotE ? 12 : 16) NbNonlocalKernel(const SuperCluster* const superClusters,
	const QuarterEntryTask* const tasks, const QuarterEntry* const entries,
	const NbForceAccumulator nbForceAcc /*MD only*/, SCResult* const results /*EM only*/, const int* const entryResultIndices /*EM only*/,
	const SuperClusterMeta* const scMeta, Float3 boxSize, Float3 boxSizeInv, float ewaldKappa, float cutoffSq, int nSuperclusters)
{
	using namespace NbNonlocal;
	constexpr int nGroups = 8;
	const int warpInBlock = threadIdx.x >> 5;
	const int scId = blockIdx.x * 2 + warpInBlock;
	if (scId >= nSuperclusters)
		return; // Only warp level syncs below
	const int laneInWarp = threadIdx.x & 31;
	const int group = laneInWarp >> 2;
	const int lane = laneInWarp & 3;

	// In MD, charge and epsilon are pre-scaled, see LJ::PrescaleNBParams
	__shared__ float4 ownPosCharge[2][SuperCluster::maxParticles];
	__shared__ float2 ownSigmaEpsilon[2][SuperCluster::maxParticles];
	__shared__ float4 jPosCharge[2][nGroups][4];
	__shared__ float2 jSigmaEpsilon[2][nGroups][4];

	if (laneInWarp < SuperCluster::maxParticles) {
		float4 pq = superClusters[scId].posCharge[laneInWarp];
		float2 se = superClusters[scId].sigmaEpsilon[laneInWarp];
		if constexpr (!energyMinimize)
			LJ::PrescaleNBParams(se.y, pq.w);
		ownPosCharge[warpInBlock][laneInWarp] = pq;
		ownSigmaEpsilon[warpInBlock][laneInWarp] = se;
	}
	const QuarterEntryTask task = tasks[scId];
	__syncwarp();

	const Float3 reference{ ownPosCharge[warpInBlock][0].x, ownPosCharge[warpInBlock][0].y, ownPosCharge[warpInBlock][0].z };
	ForceEnergy accOwn[4]{};

	const int end = task.firstEntry + task.nEntries;
	for (int base = task.firstEntry; base < end; base += nGroups) {
		const int e = base + group;
		const bool activeGroup = e < end;
		const QuarterEntry entry = entries[activeGroup ? e : base];
		const uint32_t ownQuarterMask = activeGroup ? entry.ownQuarterMask : 0u;

		const int jIndex = entry.jQuarter * 4 + lane;
		float4 pqJ = superClusters[entry.jScId].posCharge[jIndex];
		float2 seJ = superClusters[entry.jScId].sigmaEpsilon[jIndex];
		const bool validJ = seJ.y != -1.f;
		if constexpr (!energyMinimize)
			LJ::PrescaleNBParams(seJ.y, pqJ.w);
		BoundaryCondition::ApplyHyperpos(reference, pqJ.x, pqJ.y, pqJ.z, boxSize, boxSizeInv);
		jPosCharge[warpInBlock][group][lane] = pqJ;
		jSigmaEpsilon[warpInBlock][group][lane] = seJ;
		__syncwarp();

		const uint32_t unionMask = __reduce_or_sync(0xFFFFFFFFu, ownQuarterMask);
		ForceEnergy jForce{}; // On j particle jIndex

		if constexpr (!energyMinimize) {
			ForceEnergy accJ[4]{};
#pragma unroll
			for (int ownQuarter = 0; ownQuarter < 4; ownQuarter++) {
				if (!((unionMask >> ownQuarter) & 1))
					continue;
				const uint32_t noInteractions = RotatedNoInteractions(entry.noInteractions[ownQuarter], ownQuarterMask, ownQuarter, lane);
				const float4 pqI = ownPosCharge[warpInBlock][ownQuarter * 4 + lane];
				const float2 seI = ownSigmaEpsilon[warpInBlock][ownQuarter * 4 + lane];
#pragma unroll
				for (int k = 0; k < 4; k++) {
					const int jLocal = (lane + k) & 3;
					const float4 pqJ = jPosCharge[warpInBlock][group][jLocal];
					const float2 seJ = jSigmaEpsilon[warpInBlock][group][jLocal];
					const Float3 diff{ pqI.x - pqJ.x, pqI.y - pqJ.y, pqI.z - pqJ.z };
					const bool masked = ((noInteractions >> k) & 1) | (diff.lenSquared() >= cutoffSq);
					ForceEnergy fe = LJ::ComputePairNB<computePotE>(diff, LJ::CalcSigma(seI.x, seJ.x), LJ::CalcEpsilon(seI.y, seJ.y), pqI.w * pqJ.w, masked, ewaldKappa);
					accJ[k] += fe;
					accOwn[ownQuarter] += fe.InvertForce();
				}
			}

			// Lane m collects accJ[k] from lane (m-k)&3, which computed j particle m in iteration k
			jForce = accJ[0];
#pragma unroll
			for (int k = 1; k < 4; k++)
				jForce += ShuffleForceEnergy<computePotE>(accJ[k], (lane - k) & 3);
		}
		else {
			// The EM pair interaction is much larger, so the loops are not unrolled, keeping one copy of it and the registers it needs.
			// The j forces are therefore routed every iteration, and the own forces added to their quarter with selects
#pragma unroll 1
			for (int ownQuarter = 0; ownQuarter < 4; ownQuarter++) {
				if (!((unionMask >> ownQuarter) & 1))
					continue;
				const uint32_t noInteractions = RotatedNoInteractions(NoInteractionsOfQuarter(entry, ownQuarter), ownQuarterMask, ownQuarter, lane);
				const float4 pqI = ownPosCharge[warpInBlock][ownQuarter * 4 + lane];
				const float2 seI = ownSigmaEpsilon[warpInBlock][ownQuarter * 4 + lane];
				ForceEnergy ownForce{};
#pragma unroll 1
				for (int k = 0; k < 4; k++) {
					const int jLocal = (lane + k) & 3;
					const float4 pqJ = jPosCharge[warpInBlock][group][jLocal];
					const float2 seJ = jSigmaEpsilon[warpInBlock][group][jLocal];
					const Float3 diff{ pqI.x - pqJ.x, pqI.y - pqJ.y, pqI.z - pqJ.z };
					const bool masked = ((noInteractions >> k) & 1) | (diff.lenSquared() >= cutoffSq);
					ForceEnergy fe = masked ? ForceEnergy{}
						: LJ::ComputePairNBEm<computePotE>(diff, LJ::CalcSigma(seI.x, seJ.x), LJ::CalcEpsilon(seI.y, seJ.y), pqI.w * pqJ.w, ewaldKappa);
					ownForce += fe.InvertForce();
					// Lane m collects the force lane (m-k)&3 computed on j particle m
					jForce += ShuffleForceEnergy<computePotE>(fe, (lane - k) & 3);
				}
				AddToQuarter<computePotE>(accOwn, ownQuarter, ownForce);
			}
		}

		if constexpr (energyMinimize) {
			if (activeGroup) {
				if (entry.jScId == scId) {
					AddToQuarter<computePotE>(accOwn, entry.jQuarter, jForce);
				}
				else {
					results[entryResultIndices[e]].Store<computePotE>(jIndex, jForce);
				}
			}
		}
		else {
			if (activeGroup && validJ)
				nbForceAcc.Add<computePotE>(entry.jScId * SuperCluster::maxParticles + jIndex, jForce);
		}
		__syncwarp(); // Before the j quarter is restaged
	}

	// Sum the own forces over the groups. A butterfly gives all groups the bitwise identical sum
#pragma unroll
	for (int q = 0; q < 4; q++) {
		for (int offset = 4; offset < 32; offset <<= 1) {
			accOwn[q].force.x += __shfl_xor_sync(0xFFFFFFFFu, accOwn[q].force.x, offset);
			accOwn[q].force.y += __shfl_xor_sync(0xFFFFFFFFu, accOwn[q].force.y, offset);
			accOwn[q].force.z += __shfl_xor_sync(0xFFFFFFFFu, accOwn[q].force.z, offset);
			if constexpr (computePotE)
				accOwn[q].potE += __shfl_xor_sync(0xFFFFFFFFu, accOwn[q].potE, offset);
		}
	}
	// Group q writes own quarter q
#pragma unroll
	for (int q = 0; q < 4; q++) {
		if (group != q)
			continue;
		if constexpr (energyMinimize)
			results[scMeta[scId].resultsStartIndex].Store<computePotE>(q * 4 + lane, accOwn[q]);
		else if (superClusters[scId].Valid(q * 4 + lane))
			nbForceAcc.Add<computePotE>(scId * SuperCluster::maxParticles + q * 4 + lane, accOwn[q]);
	}
}

// blockDim=(16, 4, 1)
template<typename BoundaryCondition, bool emvariant, bool logData>
__global__ void SuperclusterIntegrateKernel(const ForceEnergyInterims forceEnergies, Float3* const emForces /*Only available in EM*/, int data_logging_interval, 
	const SCResult* const scResults /*Only used in EM*/, const NbForceAccumulator nbForceAcc /*Only used in MD*/,
	SuperCluster* superClusters, const SuperClusterMeta* const scMeta, PersistentCluster* const pclusters, const PersistentClusterMeta* const pcMeta, PersistentclusterInterimState* const pcStates, 
	int64_t step, const IntegrationSimulationData* simulationData, int nSuperclusters, float* forcesMagnitudeSquaredBuffer, /*Only available in EM*/
	Float3 boxSize, Float3* fixedParticleMovementBuffer, Float3* forceMaskBuffer, const Rotation* fixedParticleRotationBuffer,
	Float3* trajBuffer, float* potEBuffer, float* velocityBuffer, Float3* forceBuffer /*Only available in LIVEEDIT*/  /*,
const ForceEnergy* const nbForceenergy*/) {

	const int nScsPerBlock = 4;

	const int scIdLocal = threadIdx.y;
	const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? static_cast<int>(blockIdx.x * nScsPerBlock + threadIdx.y) : -1;


	//__shared__ Float3 positions[SuperCluster::nParticles * nScsPerBlock];
	__shared__ Float3 p0s[nScsPerBlock];
	__shared__ int resultsStartIndex[nScsPerBlock];
	__shared__ int nResults[nScsPerBlock];
	__shared__ IntegrationSimulationData simulations[nScsPerBlock];


	if (threadIdx.x == 0) {
		resultsStartIndex[threadIdx.y] = scIdGlobal == -1 ? -1 : scMeta[scIdGlobal].resultsStartIndex;
		nResults[threadIdx.y] = scIdGlobal == -1 ? 0 : scMeta[scIdGlobal].nResults;
		simulations[threadIdx.y] = scIdGlobal == -1 ? IntegrationSimulationData{} : simulationData[scMeta[scIdGlobal].simulationId];
		p0s[threadIdx.y] = scIdGlobal == -1 ? Float3{} : superClusters[scIdGlobal].Position(0);

		// By applying BC here, we dont need to wait for thread0 later in the kernel
		BoundaryCondition::applyBCNM(p0s[threadIdx.y], boxSize);// TODO: We should use either SC CoM, or a particle close to the middle..
	}
	__syncthreads();

	//const int pidInPcluster = threadIdx.x % 4;									// Always safe
	const int pidInPcluster = scIdGlobal == -1 ? -1 : scMeta[scIdGlobal].indexInPcluster[threadIdx.x];	// May be -1
	const int pcIdGlobal = scIdGlobal == -1 ? -1 :  scMeta[scIdGlobal]._pclusterIds[threadIdx.x];		// May be -1
	const int pidGlobal = scIdGlobal == -1 ? -1 : scMeta[scIdGlobal].globalParticleIds[threadIdx.x];	// May be -1
	//const int pidGlobal = pcIdGlobal == -1 ? -1 : pcMeta[pcIdGlobal].particleIdsGlobal[pidInPcluster];

	if (pidGlobal == -1)
		return;// NO SYNCS AFTER THIS!

	if constexpr (INDEXING_CHECKS) {
		if (pidInPcluster == -1 || pcIdGlobal == -1 || pidGlobal == -1) {
			printf("Indexing check failed in SuperclusterIntegrateKernel! scIdGlobal %d, pidInPcluster %d, pcIdGlobal %d, pidGlobal %d\n", scIdGlobal, pidInPcluster, pcIdGlobal, pidGlobal);
		}
	}

	Float3 pos = superClusters[scIdGlobal].Position(threadIdx.x);
	BoundaryCondition::applyHyperposNM(p0s[threadIdx.y], pos, boxSize);
	const auto& simulation = simulations[scIdLocal];

	// Collect ForceEnergy from all sources
	ForceEnergy fe{};
	// Gather from NB kernels
	if constexpr (emvariant) {
		for (int i = resultsStartIndex[scIdLocal]; i < resultsStartIndex[scIdLocal] + nResults[scIdLocal]; i++) {
			const ForceEnergy result = scResults[i].Load<logData>(threadIdx.x);
			KernelHelpersWarnings::ForceCheck(result.force);
			fe += result;
		}
	}
	else {
		const ForceEnergy result = nbForceAcc.Load<logData>(scIdGlobal * SuperCluster::maxParticles + threadIdx.x);
		KernelHelpersWarnings::ForceCheck(result.force);
		fe += result;
	}
//	fe += nbForceenergy[scIdGlobal * SuperCluster::nParticles + threadIdx.x];
	// Gather from the bondgroups this particle appears in
	{
		const BondgroupRefManager& beRefs = pcMeta[pcIdGlobal].bondgroupReferences[pidInPcluster];
		for (int i = 0; i < beRefs.nBondgroupApperances; i++)
			fe += forceEnergies.forceEnergiesBondgroups[beRefs.bondgroupApperances[i].indexInForceEnergiesBondgroups];
	}
	fe += forceEnergies.snf[pcIdGlobal * PersistentCluster::maxParticles + pidInPcluster]; // TODO: Should this be compiletime, or maybe launch param decided to ignore? Since most simulations would ignore...
	fe += forceEnergies.pme[pcIdGlobal * PersistentCluster::maxParticles + pidInPcluster]; // TODO: OPTIM: These should follow SC layout, not PC


	// Write the force to global buffer, only needed to monitor simulations
	forcesMagnitudeSquaredBuffer[pidGlobal] = fe.force.lenSquared();

	// ------------------------------------------------------------ Integration --------------------------------------------------------------- //	
	float speed = 0.f;

	const float mass = pcMeta[pcIdGlobal].mass[pidInPcluster];

	// Energy minimize. The particles are moved by EM::UpdateKernel once the step is decided for the whole simulation
	if constexpr (emvariant) {
		// Overlapping particles in unminimized structures can produce infinite forces. Push them apart in a pseudorandom
		// direction, with a force large enough to still dominate the max force
		const bool finite = isfinite(fe.force.lenSquared());
		emForces[pcIdGlobal * PersistentCluster::maxParticles + pidInPcluster] = finite ? fe.force : EngineUtils::GenerateRandomForce(pidGlobal) * 1e9f;
	}
	else {

		if (forceMaskBuffer) {
			fe.force = fe.force * forceMaskBuffer[pidGlobal];
		}

		const Float3 forcePrev = pcStates[pcIdGlobal].forces_prev[pidInPcluster];
		const Float3 velPrev = pcStates[pcIdGlobal].vels_prev[pidInPcluster];
		Float3 vel_now = EngineUtils::integrateVelocityVVS(velPrev, forcePrev, fe.force, simulation.dt, mass);
		Float3 pos_now = EngineUtils::IntegratePositionVVS(pos, vel_now, fe.force, mass, simulation.dt);

		if constexpr (FORCE_CHECKS) {
			if ((pos_now - pos).len() > 0.5f) {
				printf("Warning: Particle %d in PC %d moved %f nm in one step.\n", pidInPcluster, pcIdGlobal, (pos_now - pos).len());
			}
		}

		// TODO: This should be happening in the EM variant, but i cant get that working properly
		if (fixedParticleMovementBuffer != nullptr) {
			Float3 fixedMovement = fixedParticleMovementBuffer[pidGlobal];
			if (fixedMovement.lenSquared() > 0) {
				fe.force = Float3{};
				vel_now = Float3{};
				pos_now = pos + fixedMovement;
			}				
		}

		if (fixedParticleRotationBuffer != nullptr) {
			fe.force = Float3{};
			vel_now = Float3{}; 			
			const Rotation& rotation = fixedParticleRotationBuffer[pidGlobal];
			BoundaryCondition::applyHyperposNM(rotation.center, pos_now, boxSize);
			LAL::RotatePoint(pos_now, rotation.center, rotation.rotation);
			BoundaryCondition::applyHyperposNM(p0s[threadIdx.y], pos_now, boxSize);
			//LAL::RotatePoint(pos_now, Float3{}, Float3{ 0.001f, 0.f, 0.f });
		}
	
		pos = pos_now;// Save pos locally, but only push to box as this kernel ends

		Float3 velScaled;
		velScaled = vel_now * simulation.thermostatScalar;

		pcStates[pcIdGlobal].forces_prev[pidInPcluster] = fe.force;
		pcStates[pcIdGlobal].vels_prev[pidInPcluster] = velScaled;

		speed = velScaled.len();
	}

	// ------------------------------------------------------------ Boundary Condition --------------------------------------------------------------- //	

	//BoundaryCondition::applyHyperposNM(p0s[threadIdx.y], pos);

	EngineUtils::LogPclusterData(pcIdGlobal - simulation.pclusters.offset, pidInPcluster, step, data_logging_interval, pos, fe.potE, fe.force, speed,
		simulation.pclusters.count * PersistentCluster::maxParticles, trajBuffer + simulation.logOffset, potEBuffer + simulation.logOffset,
		velocityBuffer + simulation.logOffset, forceBuffer + simulation.logOffset);

	superClusters[scIdGlobal].SetPosition(threadIdx.x, pos);
	pclusters[pcIdGlobal].pqd[pidInPcluster].position = pos;
}


