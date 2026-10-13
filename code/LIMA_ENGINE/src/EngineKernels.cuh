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
#include "ForceAccumulator.cuh"



















// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//


// Finds pcluster particles in the superclusters, which hold the current positions
struct ParticleSlots {
	const SuperCluster* superClusters = nullptr;
	const int* pclusterParticleSlots = nullptr;

	__device__ int Slot(int pcId, int pid) const { return pclusterParticleSlots[pcId * PersistentCluster::maxParticles + pid]; }
	__device__ Float3 Position(int slot) const { return superClusters[slot / SuperCluster::maxParticles].Position(slot % SuperCluster::maxParticles); }
};

// Where the SNF kernels put a particle's force: added to forceAcc in MD, stored in the pcluster layout snf buffer in EM
struct SnfForceOutput {
	ForceEnergy* snf = nullptr;					// EM only
	ForceAccumulator forceAcc{};				// MD only

	__device__ void Put(const ParticleSlots& slots, int pcId, int pid, const ForceEnergy& fe) const {
		if (snf) snf[pcId * PersistentCluster::maxParticles + pid] = fe;
		else forceAcc.Add(slots.Slot(pcId, pid), fe);
	}
};

template <typename BoundaryCondition, bool energyMinimize>
__global__ void PclusterSnfKernel(const PersistentCluster* const pc, const PersistentClusterMeta* const pcMeta, UniformElectricField electricField, const SnfForceOutput output,
	const ParticleSlots particleSlots, int pclusterOffset, int nPclusters) {

	const int workId = blockIdx.x * blockDim.x + threadIdx.x;
	if (workId >= nPclusters) return;
	const int pcId = pclusterOffset + workId;

	
	for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
		const int pidGlobal = pcMeta[pcId].particleIdsGlobal[pid];
		if (pidGlobal == -1)
			continue;

		float charge = pc[pcId].pqd[pid].params.charge;
		Float3 force = electricField.GetForce(charge);

		output.Put(particleSlots, pcId, pid, ForceEnergy{ force, 0.f });
	}
}

__global__ void ElasticPositionsForceKernel(const PersistentClusterMeta* const pcMeta, const Float3* const elasticPositions, const SnfForceOutput output,
	const ParticleSlots particleSlots, int pclusterOffset, int nPclusters, Float3 boxSize) {
	const int workId = blockIdx.x * blockDim.x + threadIdx.x;
	if (workId >= nPclusters) return;
	const int pcId = pclusterOffset + workId;


	for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
		const int pidGlobal = pcMeta[pcId].particleIdsGlobal[pid];
		if (pidGlobal == -1)
			continue;

		const float mass = pcMeta[pcId].mass[pid];
		const Float3 position = particleSlots.Position(particleSlots.Slot(pcId, pid));
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

		output.Put(particleSlots, pcId, pid, ForceEnergy{ forceDirection * forceMagnitude, potentialEnergy });
	}
}





static const int THREADS_PER_BONDSGROUPSKERNEL = 128;
static constexpr int BONDGROUPS_PER_BLOCK = THREADS_PER_BONDSGROUPSKERNEL / 32;
static constexpr int BONDGROUP_MAX_PARTICLES = 64;
static_assert(BONDGROUP_MAX_PARTICLES <= 64, "Compact bond records reserve 6 bits per atom id");
static_assert(LimaForcecalc::FixedPointBondAccumulator::scale == ForceAccumulator::scale, "Integration combines bonded and nonbonded sums in the same fixed point");
// withCmap: the CMAP code raises the kernel from 62 to 84 registers, which cuts occupancy, so batches without CMAP use a variant without it
template <typename BoundaryCondition, bool emVariant, bool computePotE, bool withCmap>
__global__ void BondgroupsKernel(const CompactBondGroupsDevice bondGroups, int nBondgroups, const BoxState boxState, ForceEnergy* const forceEnergiesOut /*EM only*/,
	const ForceAccumulator primaryBondForces, ulonglong4* const extraBondForces, const int* const extraResultIndices /*MD only*/, const SuperCluster* const superClusters, const int* const particleSlots, Float3 boxSize, Float3 boxSizeInv) {
	const int groupInBlock = threadIdx.x / 32;
	const int groupId = blockIdx.x * BONDGROUPS_PER_BLOCK + groupInBlock;
	if (groupId >= nBondgroups) return; // Uniform per warp; no block barriers below
	const int lane = threadIdx.x & 31;
	__shared__ Float3 positionsBuffer[BONDGROUPS_PER_BLOCK][BONDGROUP_MAX_PARTICLES];
	Float3* const positions = positionsBuffer[groupInBlock];

	// MD accumulates in parallel with deterministic fixed point atomics, EM serially in float, since its forces may exceed the fixed point range
	using Accumulator = std::conditional_t<emVariant, LimaForcecalc::SerialBondAccumulator, LimaForcecalc::FixedPointBondAccumulatorT<computePotE>>;
	__shared__ __align__(16) unsigned long long accumulatorBuffers[BONDGROUPS_PER_BLOCK][(emVariant ? 2 : computePotE ? 4 : 3) * BONDGROUP_MAX_PARTICLES];
	auto* accumulatorBuffer = accumulatorBuffers[groupInBlock];
	Accumulator acc;
	if constexpr (emVariant) acc = Accumulator{ reinterpret_cast<float4*>(accumulatorBuffer) };
	else acc = Accumulator{ reinterpret_cast<unsigned int*>(accumulatorBuffer), BONDGROUP_MAX_PARTICLES };

	static const int batchSize = 32;
	static const int largestBondBytesize = std::max(sizeof(AngleUreyBradleyBond), sizeof(DihedralBond));
	static_assert(sizeof(CmapBond) <= largestBondBytesize);
	__shared__ __align__(16) char bondsBuffers[BONDGROUPS_PER_BLOCK][largestBondBytesize * batchSize];
	char* const _bondsBuffer = bondsBuffers[groupInBlock];

	const BondGroup* const bondGroup = &bondGroups.groups[groupId];

	int slots[2]{ -1, -1 };
#pragma unroll
	for (int i = 0; i < 2; i++) {
		const int pid = lane + i * 32;
		acc.Init(pid);
		if (pid < bondGroup->nParticles) {
			const int slot = particleSlots[bondGroup->indexOfFirstParticle + pid];
			slots[i] = slot;
			const float4 pq = __ldg(&superClusters[slot / SuperCluster::maxParticles].posCharge[slot % SuperCluster::maxParticles]);
			positions[pid] = Float3{ pq.x, pq.y, pq.z };
		}
	}
	acc.Sync();

	
	{
		acc.Sync();
		SingleBond* bondsBuffer = reinterpret_cast<SingleBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nSinglebonds; batchStart += batchSize) {
			if (batchStart + lane < bondGroup->nSinglebonds) {
				const int bondIndex = batchStart + lane;
				bondGroups.singlebonds.Load(bondGroup->indexOfFirstSinglebond + bondIndex, bondsBuffer[lane]);
			}
			acc.Sync();

			LimaForcecalc::computeSinglebondForces<BoundaryCondition, Accumulator, emVariant>(bondsBuffer, std::min(batchSize, bondGroup->nSinglebonds - batchStart), positions, acc, 0, boxSize, boxSizeInv);
		}
	}

	{
		acc.Sync();
		AngleUreyBradleyBond* bondsBuffer = reinterpret_cast<AngleUreyBradleyBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nAnglebonds; batchStart += batchSize) {
			if (batchStart + lane < bondGroup->nAnglebonds) {
				const int bondIndex = batchStart + lane;
				bondGroups.anglebonds.Load(bondGroup->indexOfFirstAnglebond + bondIndex, bondsBuffer[lane]);
			}
			acc.Sync();

			LimaForcecalc::computeAnglebondForces<BoundaryCondition, Accumulator, emVariant>(bondsBuffer, std::min(batchSize, bondGroup->nAnglebonds - batchStart), positions, acc, boxSize, boxSizeInv);
		}
	}

	{
		acc.Sync();
		DihedralBond* bondsBuffer = reinterpret_cast<DihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nDihedralbonds; batchStart += batchSize) {
			if (batchStart + lane < bondGroup->nDihedralbonds) {
				const int bondIndex = batchStart + lane;
				bondGroups.dihedralbonds.Load(bondGroup->indexOfFirstDihedralbond + bondIndex, bondsBuffer[lane]);
			}
			acc.Sync();

			LimaForcecalc::computeDihedralForces<BoundaryCondition, Accumulator>(bondsBuffer, std::min(batchSize, bondGroup->nDihedralbonds - batchStart), positions, acc, boxSize, boxSizeInv);
		}
	}

	{
		acc.Sync();
		ImproperDihedralBond* bondsBuffer = reinterpret_cast<ImproperDihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nImproperdihedralbonds; batchStart += batchSize) {
			if (batchStart + lane < bondGroup->nImproperdihedralbonds) {
				const int bondIndex = batchStart + lane;
				bondGroups.improperdihedralbonds.Load(bondGroup->indexOfFirstImproperdihedralbond + bondIndex, bondsBuffer[lane]);
			}
			acc.Sync();

			LimaForcecalc::computeImproperdihedralForces<BoundaryCondition, Accumulator>(bondsBuffer, std::min(batchSize, bondGroup->nImproperdihedralbonds - batchStart), positions, acc, boxSize, boxSizeInv);
		}
	}

	if constexpr (withCmap) {
		acc.Sync();
		CmapBond* bondsBuffer = reinterpret_cast<CmapBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nCmapbonds; batchStart += batchSize) {
			if (batchStart + lane < bondGroup->nCmapbonds)
				bondsBuffer[lane] = bondGroups.cmapbonds[bondGroup->indexOfFirstCmapbond + batchStart + lane];
			acc.Sync();

			LimaForcecalc::computeCmapForces<BoundaryCondition, Accumulator>(bondsBuffer, std::min(batchSize, bondGroup->nCmapbonds - batchStart), positions, acc,
				bondGroups.cmapGridPoints, bondGroups.cmapGridSpacing, boxSize, boxSizeInv);
		}
	}


	{
		acc.Sync();
		// TODO: i have no clue if pairbonds should also compute SR electrostatics?
		PairBond* bondsBuffer = reinterpret_cast<PairBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nPairbonds; batchStart += batchSize) {
			if (batchStart + lane < bondGroup->nPairbonds) {
				const int bondIndex = batchStart + lane;
				bondGroups.pairbonds.Load(bondGroup->indexOfFirstPairbond + bondIndex, bondsBuffer[lane]);
			}
			acc.Sync();

			LimaForcecalc::computePairbondForces<BoundaryCondition, Accumulator>(bondsBuffer, std::min(batchSize, bondGroup->nPairbonds - batchStart), positions, acc, boxSize, boxSizeInv);
		}
	}

#pragma unroll
	for (int i = 0; i < 2; i++) {
		const int pid = lane + i * 32;
		if (pid < bondGroup->nParticles) {
			if constexpr (emVariant)
				forceEnergiesOut[bondGroup->indexOfFirstParticle + pid] = acc.Get(pid);
			else {
				const int extra = extraResultIndices[bondGroup->indexOfFirstParticle + pid];
				if (extra < 0) acc.StorePrimary(primaryBondForces, pid, slots[i]);
				else acc.StoreTo(extraBondForces, pid, extra);
			}
		}
	}
}



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
// MD sums the forces with ForceAccumulator. EM forces can exceed its range, so EM stores each entry's j forces in the entry's
// quarter of the SCResult of its (own, j) supercluster pair, and the own forces in the own supercluster's first SCResult.
// SuperclusterIntegrateKernel sums them in a fixed order. In EM every quarter of such a pair has an entry, also quarters without
// pairs in range, so every result is written each step. Entries of the own supercluster add their j forces to the own forces.
//
// __launch_bounds__: 64 registers avoid spills (16 blocks per SM), the potE variants need 80 (12 blocks per SM).
// After larger changes, check ncu for local memory traffic and revisit the bounds.
template <typename BoundaryCondition, bool energyMinimize, bool computePotE>
__global__ void __launch_bounds__(64, computePotE ? 12 : 16) NbNonlocalKernel(const SuperCluster* const superClusters,
	const QuarterEntryTask* const tasks, const QuarterEntry* const entries,
	const ForceAccumulator forceAcc /*MD only*/, SCResult* const results /*EM only*/, const int* const entryResultIndices /*EM only*/,
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
				forceAcc.Add<computePotE>(entry.jScId * SuperCluster::maxParticles + jIndex, jForce);
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
			forceAcc.Add<computePotE>(scId * SuperCluster::maxParticles + q * 4 + lane, accOwn[q]);
	}
}

// Particles edited live, indexed by global particle id. Each is nullptr unless in use
struct LiveEditBuffers {
	const Float3* fixedParticleMovement = nullptr;
	const Float3* forceMask = nullptr;
	const Rotation* fixedParticleRotation = nullptr;

	__host__ __device__ bool Any() const { return fixedParticleMovement || forceMask || fixedParticleRotation; }
};

struct LogBuffers {
	Float3* traj = nullptr;
	float* potE = nullptr;
	float* velocity = nullptr;
	Float3* force = nullptr;
	int loggingInterval = 0;
	int64_t step = 0;

	__device__ void Log(const IntegrationSimulationData& simulation, int pcId, int indexInPc, Float3 position, const ForceEnergy& fe, float speed) const {
		EngineUtils::LogPclusterData(pcId - simulation.pclusters.offset, indexInPc, step, loggingInterval, position, fe.potE, fe.force, speed,
			simulation.pclusters.count * PersistentCluster::maxParticles, traj + simulation.logOffset, potE + simulation.logOffset,
			velocity + simulation.logOffset, force + simulation.logOffset);
	}
};

// Velocity Verlet step of each particle, with the forces all force kernels summed in forceAcc. A particle's position and integration
// state are in its supercluster slot, see ParticleIntegrationState. One supercluster per 16 threads, blockDim = (16, 4, 1)
template<typename BoundaryCondition, bool logData>
__global__ void SuperclusterIntegrateKernel(const ForceAccumulator forceAcc, const ForceAccumulator primaryBondForces, const ulonglong4* const extraBondForces, const int* const bondReferences, SuperCluster* const superClusters, ParticleIntegrationState* const states,
	const SuperClusterMeta* const scMeta, const IntegrationSimulationData* const simulationData, int nSuperclusters,
	Float3 boxSize, float* const forcesMagnitudeSquared, const LiveEditBuffers liveEdit, const LogBuffers log)
{
	const int scId = blockIdx.x * blockDim.y + threadIdx.y;
	if (scId >= nSuperclusters)
		return;
	const int slot = scId * SuperCluster::maxParticles + threadIdx.x;
	ParticleIntegrationState state = states[slot];
	if (state.mass == 0.f) // Unused slot
		return;

	const IntegrationSimulationData& simulation = simulationData[scMeta[scId].simulationId];
	const int pcId = state.pclusterParticle / PersistentCluster::maxParticles;
	const int indexInPc = state.pclusterParticle % PersistentCluster::maxParticles;
	const int pidGlobal = liveEdit.Any() ? scMeta[scId].globalParticleIds[threadIdx.x] : -1; // Only needed for live edits

	// The supercluster is kept in the periodic image of its first particle
	Float3 p0 = superClusters[scId].Position(0);
	BoundaryCondition::applyBCNM(p0, boxSize);
	Float3 pos = superClusters[scId].Position(threadIdx.x);
	BoundaryCondition::applyHyperposNM(p0, pos, boxSize);

	ForceEnergy fe = forceAcc.Take<logData>(slot, primaryBondForces, extraBondForces, bondReferences, nSuperclusters * SuperCluster::maxParticles);
	KernelHelpersWarnings::ForceCheck(fe.force);
	// Monitoring reads the force magnitudes. Engine::StoreIntegrationStates computes them from forcePrev, unless live edits may change it
	if (liveEdit.Any())
		forcesMagnitudeSquared[pidGlobal] = fe.force.lenSquared();

	if (liveEdit.forceMask)
		fe.force = fe.force * liveEdit.forceMask[pidGlobal];

	Float3 velocity = EngineUtils::integrateVelocityVVS(state.velocity, state.forcePrev, fe.force, simulation.dt, state.mass);
	Float3 posNext = EngineUtils::IntegratePositionVVS(pos, velocity, fe.force, state.mass, simulation.dt);

	if constexpr (FORCE_CHECKS) {
		if ((posNext - pos).len() > 0.5f)
			printf("Warning: Particle %d in PC %d moved %f nm in one step.\n", indexInPc, pcId, (posNext - pos).len());
	}

	if (liveEdit.fixedParticleMovement) {
		const Float3 fixedMovement = liveEdit.fixedParticleMovement[pidGlobal];
		if (fixedMovement.lenSquared() > 0) {
			fe.force = Float3{};
			velocity = Float3{};
			posNext = pos + fixedMovement;
		}
	}
	if (liveEdit.fixedParticleRotation) {
		fe.force = Float3{};
		velocity = Float3{};
		const Rotation& rotation = liveEdit.fixedParticleRotation[pidGlobal];
		BoundaryCondition::applyHyperposNM(rotation.center, posNext, boxSize);
		LAL::RotatePoint(posNext, rotation.center, rotation.rotation);
		BoundaryCondition::applyHyperposNM(p0, posNext, boxSize);
	}

	state.velocity = velocity * simulation.thermostatScalar;
	state.forcePrev = fe.force;
	states[slot] = state;

	if constexpr (logData)
		log.Log(simulation, pcId, indexInPc, posNext, fe, state.velocity.len());

	superClusters[scId].SetPosition(threadIdx.x, posNext);
}

namespace Constraints {
	constexpr int maxIterations = 100;
	// Relative error of the squared constraint length. Tight, as the energy drift of rigid water depends on it
	constexpr float tolerance = 1e-6f;
	// Relative error of the velocity along a constraint. Looser, as its float rounding floor in coupled constraints is close to 1e-6
	constexpr float velocityTolerance = 1e-5f;

	// Set when a solve fails, so the warning is printed once per process rather than once per group and step
	static __device__ int failureReported = 0;
	__device__ inline void ReportFailure(const char* what, int iterations) {
		if (atomicExch(&failureReported, 1) == 0)
			printf("Warning: %s did not converge in %d iterations, or a constraint degenerated. The system is likely unstable, eg. dt too large\n", what, iterations);
	}

	// The particles of a pcluster in its supercluster, -1 for the unused indices
	struct Group {
		int slots[PersistentCluster::maxParticles];
	};

	// Removes the velocity components along the constraints (RATTLE), at the constrained positions. Sequential, so deterministic.
	// A single constraint is exact after one pass. Coupled ones iterate until the velocity along each constraint is small
	// relative to the speeds of its particles, as an absolute measure would stall at the rounding of fast particles
	__device__ inline void ConstrainVelocities(const PclusterConstraints& constraints, const Group& group, const Float3* const positions,
		Float3* const velocities, const float* const invMasses) {
		for (int iteration = 0; iteration < maxIterations; iteration++) {
			float maxError = 0.f;
			for (int c = 0; c < constraints.nConstraints; c++) {
				const int a = group.slots[constraints.particles[c] & 3], b = group.slots[constraints.particles[c] >> 2];
				const Float3 r = positions[a] - positions[b];
				const Float3 dv = velocities[a] - velocities[b];
				const float rr = r.dot(r);
				const float rv = r.dot(dv);
				const float invMassSum = invMasses[a] + invMasses[b];
				if (invMassSum == 0.f)	// Both particles fixed by live edits
					continue;
				const float k = rv / (rr * invMassSum);
				velocities[a] -= r * (k * invMasses[a]);
				velocities[b] += r * (k * invMasses[b]);
				// Relative to the particles' speeds: the relative velocity itself can be cancelled to rounding noise
				maxError = fmaxf(maxError, fabsf(rv) / (sqrtf(rr) * (velocities[a].len() + velocities[b].len()) + 1e-30f));
			}
			if (constraints.nConstraints == 1 || maxError < velocityTolerance)
				return;
		}
		ReportFailure("RATTLE", maxIterations);
	}

	// SHAKE: moves the particles' steps so the new positions satisfy the constraints, along the constraint directions at the old
	// positions. Works on the steps rather than the positions, as those are small, so the solution is precise. Sequential
	__device__ inline void ConstrainSteps(const PclusterConstraints& constraints, const Group& group, const Float3* const positions,
		Float3* const steps, const float* const invMasses) {
		for (int iteration = 0; iteration < maxIterations; iteration++) {
			float maxError = 0.f;
			for (int c = 0; c < constraints.nConstraints; c++) {
				const int a = group.slots[constraints.particles[c] & 3], b = group.slots[constraints.particles[c] >> 2];
				const Float3 rOld = positions[a] - positions[b];
				const Float3 rNew = rOld + steps[a] - steps[b];
				const float dd = constraints.distance[c] * constraints.distance[c];
				const float diff = dd - rNew.dot(rNew);
				const float invMassSum = invMasses[a] + invMasses[b];
				const float alignment = rNew.dot(rOld);
				if (invMassSum == 0.f)	// Both particles fixed by live edits
					continue;
				if (alignment <= 0.f) {	// The step turned the constraint by 90 degrees or more, so there is no solution along rOld
					ReportFailure("SHAKE", iteration);
					return;
				}
				const float g = diff / (2.f * alignment * invMassSum);
				steps[a] += rOld * (g * invMasses[a]);
				steps[b] -= rOld * (g * invMasses[b]);
				maxError = fmaxf(maxError, fabsf(diff) / dd);
			}
			if (maxError < tolerance)
				return;
		}
		ReportFailure("SHAKE", maxIterations);
	}
}

// Puts the positions on the constraints before the first MD step, like GROMACS constrains the starting coordinates. Otherwise
// the first step would, and its jump in energy would show in the energy of the run. blockDim = (16, 4, 1)
template<typename BoundaryCondition>
__global__ void ConstrainPositionsKernel(SuperCluster* const superClusters, const ParticleIntegrationState* const states,
	const PclusterConstraints* const pclusterConstraints, int nSuperclusters, Float3 boxSize)
{
	constexpr int n = SuperCluster::maxParticles;
	__shared__ Float3 positions[4][n];
	__shared__ Float3 steps[4][n];
	__shared__ float invMasses[4][n];
	__shared__ int pclusterParticles[4][n];
	const int scInBlock = threadIdx.y;
	const int scId = blockIdx.x * blockDim.y + threadIdx.y;
	const bool inRange = scId < nSuperclusters;
	ParticleIntegrationState state{};
	state.mass = 0.f;
	if (inRange)
		state = states[scId * n + threadIdx.x];
	const bool active = inRange && state.mass != 0.f;

	Float3 pos{}, p0{};
	if (active) {
		p0 = superClusters[scId].Position(0);
		BoundaryCondition::applyBCNM(p0, boxSize);
		pos = superClusters[scId].Position(threadIdx.x);
		BoundaryCondition::applyHyperposNM(p0, pos, boxSize);
	}
	positions[scInBlock][threadIdx.x] = pos;
	steps[scInBlock][threadIdx.x] = Float3{};
	invMasses[scInBlock][threadIdx.x] = active ? 1.f / state.mass : 0.f;
	pclusterParticles[scInBlock][threadIdx.x] = active ? state.pclusterParticle : -1;
	__syncthreads();

	const int pcId = state.pclusterParticle / PersistentCluster::maxParticles;
	if (active && state.pclusterParticle % PersistentCluster::maxParticles == 0) {
		PclusterConstraints constraints = pclusterConstraints[pcId];
		Constraints::Group group{ { -1, -1, -1, -1 } };
		for (int i = 0; i < n; i++) {
			const int particle = pclusterParticles[scInBlock][i];
			if (particle >= 0 && particle / PersistentCluster::maxParticles == pcId)
				group.slots[particle % PersistentCluster::maxParticles] = i;
		}
		for (int c = 0; c < constraints.nConstraints; c++)
			if (group.slots[constraints.particles[c] & 3] < 0 || group.slots[constraints.particles[c] >> 2] < 0)
				constraints.nConstraints = 0;
		if (constraints.nConstraints > 0)
			Constraints::ConstrainSteps(constraints, group, positions[scInBlock], steps[scInBlock], invMasses[scInBlock]);
	}
	__syncthreads();

	if (active)
		superClusters[scId].SetPosition(threadIdx.x, pos + steps[scInBlock][threadIdx.x]);
}

// SuperclusterIntegrateKernel with distance constraints (constraints = h-bonds), as RATTLE: the velocities are made orthogonal to
// the constraints, and SHAKE puts the positions after the step back on them. The position correction is added to forcePrev as
// the constraint force, so the next step's half kick includes it, while the stored velocity stays the constrained full-step
// velocity. A constrained group never spans pclusters, and each pcluster's constraints are solved by the thread of its first
// particle in a fixed order, so the result is deterministic. One supercluster per 16 threads, blockDim = (16, 4, 1)
template<typename BoundaryCondition, bool logData>
__global__ void SuperclusterIntegrateConstrainedKernel(const ForceAccumulator forceAcc, const ForceAccumulator primaryBondForces, const ulonglong4* const extraBondForces, const int* const bondReferences, SuperCluster* const superClusters, ParticleIntegrationState* const states,
	const SuperClusterMeta* const scMeta, const IntegrationSimulationData* const simulationData, int nSuperclusters,
	Float3 boxSize, float* const forcesMagnitudeSquared, const LiveEditBuffers liveEdit, const LogBuffers log, const PclusterConstraints* const pclusterConstraints)
{
	constexpr int n = SuperCluster::maxParticles;
	__shared__ Float3 positions[4][n];	// Before the step
	__shared__ Float3 velocities[4][n];	// Then the steps
	__shared__ float invMasses[4][n];
	__shared__ int pclusterParticles[4][n];
	const int scInBlock = threadIdx.y;
	const int scId = blockIdx.x * blockDim.y + threadIdx.y;
	const int slot = scId * n + threadIdx.x;

	// No early returns, as the threads synchronize below
	const bool inRange = scId < nSuperclusters;
	ParticleIntegrationState state{};
	state.mass = 0.f;
	if (inRange)
		state = states[slot];
	const bool active = inRange && state.mass != 0.f;	// Unused slots have no mass

	const int pcId = state.pclusterParticle / PersistentCluster::maxParticles;
	const int indexInPc = state.pclusterParticle % PersistentCluster::maxParticles;
	const IntegrationSimulationData* const simulation = active ? &simulationData[scMeta[scId].simulationId] : nullptr;
	const int pidGlobal = active ? scMeta[scId].globalParticleIds[threadIdx.x] : -1;

	Float3 pos{}, velocity{};
	ForceEnergy fe{};
	Float3 p0{};
	if (active) {
		// The supercluster is kept in the periodic image of its first particle
		p0 = superClusters[scId].Position(0);
		BoundaryCondition::applyBCNM(p0, boxSize);
		pos = superClusters[scId].Position(threadIdx.x);
		BoundaryCondition::applyHyperposNM(p0, pos, boxSize);

		fe = forceAcc.Take<logData>(slot, primaryBondForces, extraBondForces, bondReferences, nSuperclusters * n);
		KernelHelpersWarnings::ForceCheck(fe.force);
		// The constraint force goes into forcePrev, so monitoring gets the measured force here
		forcesMagnitudeSquared[pidGlobal] = fe.force.lenSquared();
		if (liveEdit.forceMask)
			fe.force = fe.force * liveEdit.forceMask[pidGlobal];
		velocity = EngineUtils::integrateVelocityVVS(state.velocity, state.forcePrev, fe.force, simulation->dt, state.mass);
	}
	// Particles moved by live edits take their given step. The constraints treat them as infinitely heavy, so only their
	// partners are corrected
	const bool fixed = active && (liveEdit.fixedParticleRotation
		|| (liveEdit.fixedParticleMovement && liveEdit.fixedParticleMovement[pidGlobal].lenSquared() > 0));
	if (fixed) {
		fe.force = Float3{};
		velocity = Float3{};
	}
	positions[scInBlock][threadIdx.x] = pos;
	velocities[scInBlock][threadIdx.x] = velocity;
	invMasses[scInBlock][threadIdx.x] = active && !fixed ? 1.f / state.mass : 0.f;
	pclusterParticles[scInBlock][threadIdx.x] = active ? state.pclusterParticle : -1;
	__syncthreads();

	PclusterConstraints constraints{};
	Constraints::Group group{ { -1, -1, -1, -1 } };
	if (active && indexInPc == 0) {
		constraints = pclusterConstraints[pcId];
		for (int i = 0; i < n; i++) {
			const int particle = pclusterParticles[scInBlock][i];
			if (particle >= 0 && particle / PersistentCluster::maxParticles == pcId)
				group.slots[particle % PersistentCluster::maxParticles] = i;
		}
		// Pclusters are placed in superclusters whole, so this only guards against a broken layout
		for (int c = 0; c < constraints.nConstraints; c++)
			if (group.slots[constraints.particles[c] & 3] < 0 || group.slots[constraints.particles[c] >> 2] < 0) {
				printf("Constrained particles of pcluster %d are not in the same supercluster\n", pcId);
				constraints.nConstraints = 0;
			}
		if (constraints.nConstraints > 0)
			Constraints::ConstrainVelocities(constraints, group, positions[scInBlock], velocities[scInBlock], invMasses[scInBlock]);
	}
	__syncthreads();

	velocity = velocities[scInBlock][threadIdx.x];
	Float3 step{};
	if (active)	// As IntegratePositionVVS, but without adding pos, which would cost the step most of its precision
		step = velocity * simulation->dt + fe.force * (0.5f / state.mass * simulation->dt * simulation->dt);
	if (fixed) {
		Float3 posNext = pos + step;
		if (liveEdit.fixedParticleMovement && liveEdit.fixedParticleMovement[pidGlobal].lenSquared() > 0)
			posNext = pos + liveEdit.fixedParticleMovement[pidGlobal];
		if (liveEdit.fixedParticleRotation) {
			const Rotation& rotation = liveEdit.fixedParticleRotation[pidGlobal];
			BoundaryCondition::applyHyperposNM(rotation.center, posNext, boxSize);
			LAL::RotatePoint(posNext, rotation.center, rotation.rotation);
			BoundaryCondition::applyHyperposNM(p0, posNext, boxSize);
		}
		step = posNext - pos;
	}
	__syncthreads();	// Everyone has read their velocity before the steps replace them
	velocities[scInBlock][threadIdx.x] = step;
	__syncthreads();

	if (constraints.nConstraints > 0)
		Constraints::ConstrainSteps(constraints, group, positions[scInBlock], velocities[scInBlock], invMasses[scInBlock]);
	__syncthreads();

	if (!active)
		return;
	const Float3 constrainedStep = velocities[scInBlock][threadIdx.x];
	const Float3 posNext = pos + constrainedStep;
	// The correction to the step is dt^2/2m of the constraint force. Fixed particles are not corrected, so they get none
	const Float3 constraintForce = (constrainedStep - step) * (2.f * state.mass / (simulation->dt * simulation->dt));

	state.velocity = velocity * simulation->thermostatScalar;
	state.forcePrev = fe.force + constraintForce;
	states[slot] = state;

	if constexpr (logData)
		log.Log(*simulation, pcId, indexInPc, posNext, fe, state.velocity.len());

	superClusters[scId].SetPosition(threadIdx.x, posNext);
}

// EM: sums each particle's forces from the force kernels into emForces, for EM::UpdateKernel to move the particles once the
// step is decided for the whole simulation. One supercluster per 16 threads, blockDim = (16, 4, 1)
template<typename BoundaryCondition, bool logData>
__global__ void EmCollectForcesKernel(const ForceEnergyInterims forceEnergies, const SCResult* const scResults, Float3* const emForces,
	SuperCluster* const superClusters, const SuperClusterMeta* const scMeta, PersistentCluster* const pclusters, const PersistentClusterMeta* const pcMeta,
	const IntegrationSimulationData* const simulationData, int nSuperclusters, Float3 boxSize, float* const forcesMagnitudeSquared, const LogBuffers log)
{
	const int scId = blockIdx.x * blockDim.y + threadIdx.y;
	if (scId >= nSuperclusters)
		return;
	const int pcId = scMeta[scId]._pclusterIds[threadIdx.x];
	const int indexInPc = scMeta[scId].indexInPcluster[threadIdx.x];
	const int pidGlobal = scMeta[scId].globalParticleIds[threadIdx.x];
	if (pidGlobal == -1)
		return;

	// The supercluster is kept in the periodic image of its first particle
	Float3 p0 = superClusters[scId].Position(0);
	BoundaryCondition::applyBCNM(p0, boxSize);
	Float3 pos = superClusters[scId].Position(threadIdx.x);
	BoundaryCondition::applyHyperposNM(p0, pos, boxSize);

	ForceEnergy fe{};
	// Gather from NB kernels
	const int resultsStart = scMeta[scId].resultsStartIndex;
	for (int i = resultsStart; i < resultsStart + scMeta[scId].nResults; i++) {
		const ForceEnergy result = scResults[i].Load<logData>(threadIdx.x);
		KernelHelpersWarnings::ForceCheck(result.force);
		fe += result;
	}
	// Gather from the bondgroups this particle appears in
	const BondgroupRefManager& beRefs = pcMeta[pcId].bondgroupReferences[indexInPc];
	for (int i = 0; i < beRefs.nBondgroupApperances; i++)
		fe += forceEnergies.forceEnergiesBondgroups[beRefs.bondgroupApperances[i].indexInForceEnergiesBondgroups];
	fe += forceEnergies.snf[pcId * PersistentCluster::maxParticles + indexInPc];
	fe += forceEnergies.pme[pcId * PersistentCluster::maxParticles + indexInPc];

	forcesMagnitudeSquared[pidGlobal] = fe.force.lenSquared();

	// Overlapping particles in unminimized structures can produce infinite forces. Push them apart in a pseudorandom
	// direction, with a force large enough to still dominate the max force
	const bool finite = isfinite(fe.force.lenSquared());
	emForces[pcId * PersistentCluster::maxParticles + indexInPc] = finite ? fe.force : EngineUtils::GenerateRandomForce(pidGlobal) * 1e9f;

	if constexpr (logData)
		log.Log(simulationData[scMeta[scId].simulationId], pcId, indexInPc, pos, fe, 0.f);

	superClusters[scId].SetPosition(threadIdx.x, pos);
	pclusters[pcId].pqd[indexInPc].position = pos;
}


