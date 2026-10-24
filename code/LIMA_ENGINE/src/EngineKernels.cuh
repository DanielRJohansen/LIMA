//#pragma once - this file must NOT be included multiple times

#include "Engine.cuh"
#include "ForceComputations.cuh"
#include "KernelWarnings.cuh"
#include "EngineUtils.cuh"

#include "BoundaryCondition.cuh"
#include "DeviceAlgorithms.cuh"
#include "Utilities.h"

//#include <cuda/pipeline>
#include "platform/DeviceIntrinsics.cuh"
#include "LennardJonesInteractions.cuh"
#include "ForceAccumulator.cuh"
#include "Virial.cuh"



















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
//
// warpsPerGroup: in MD, few bondgroups (small systems) leave most of the GPU idle while each warp walks its group's bonds serially,
// so then up to 4 warps of a block share a group: each takes every warpsPerGroup'th batch of each bond type, adding to the group's
// shared accumulator. Its fixed point sums are exact, so the forces do not depend on the split. Part 0 loads and stores the group.
// The warps of an inactive group still reach the block barriers. EM accumulates serially in float, so it has 1
template <typename BoundaryCondition, bool emVariant, bool computePotE, bool computeVirial, bool withCmap>
__global__ void BondgroupsKernel(const CompactBondGroupsDevice bondGroups, int nBondgroups, int warpsPerGroup, const BoxState boxState, ForceEnergy* const forceEnergiesOut /*EM only*/,
	const ForceAccumulator primaryBondForces, ulonglong4* const extraBondForces, const int* const extraResultIndices /*MD only*/, const SuperCluster* const superClusters, const int* const particleSlots, Float3 boxSize, Float3 boxSizeInv,
	float ewaldKappa /*0 without PME*/, long long* const virialPartials /*MD pressure steps only*/) {
	const int warpInBlock = threadIdx.x / 32;
	const int groupInBlock = warpInBlock / warpsPerGroup;
	const int part = warpInBlock - groupInBlock * warpsPerGroup;
	const int groupId = blockIdx.x * (BONDGROUPS_PER_BLOCK / warpsPerGroup) + groupInBlock;
	const bool split = warpsPerGroup > 1;
	if (!split && groupId >= nBondgroups) return; // Uniform per warp; no block barriers without a split
	const bool active = groupId < nBondgroups;
	const auto GroupSync = [split]() { if (split) __syncthreads(); else __syncwarp(); };
	const int lane = threadIdx.x & 31;
	LIMA_SHARED_UNINIT(Float3, positionsBuffer, [BONDGROUPS_PER_BLOCK][BONDGROUP_MAX_PARTICLES]);
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
	static_assert(sizeof(PairBond) <= largestBondBytesize);
	__shared__ __align__(16) char bondsBuffers[BONDGROUPS_PER_BLOCK][largestBondBytesize * batchSize];
	char* const _bondsBuffer = bondsBuffers[warpInBlock];

	const BondGroup* const bondGroup = &bondGroups.groups[active ? groupId : 0];
	const bool loadsAndStores = active && part == 0;
	// The batches of each bond type this warp takes
	const int firstBatch = active ? part * batchSize : 1 << 30;
	const int batchStride = warpsPerGroup * batchSize;

	int slots[2]{ -1, -1 };
	if (loadsAndStores) {
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
		if constexpr (computeVirial) {
#pragma unroll
			for (int i = 0; i < 2; i++) {
				const int pid = lane + i * 32;
				if (pid > 0 && pid < bondGroup->nParticles)
					BoundaryCondition::applyHyperposNM(positions[0], positions[pid], boxSize);
			}
		}
	}
	GroupSync();

	
	{
		acc.Sync();
		SingleBond* bondsBuffer = reinterpret_cast<SingleBond*>(_bondsBuffer);
		for (int batchStart = firstBatch; batchStart < bondGroup->nSinglebonds; batchStart += batchStride) {
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
		for (int batchStart = firstBatch; batchStart < bondGroup->nAnglebonds; batchStart += batchStride) {
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
		for (int batchStart = firstBatch; batchStart < bondGroup->nDihedralbonds; batchStart += batchStride) {
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
		for (int batchStart = firstBatch; batchStart < bondGroup->nImproperdihedralbonds; batchStart += batchStride) {
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
		for (int batchStart = firstBatch; batchStart < bondGroup->nCmapbonds; batchStart += batchStride) {
			if (batchStart + lane < bondGroup->nCmapbonds)
				bondsBuffer[lane] = bondGroups.cmapbonds[bondGroup->indexOfFirstCmapbond + batchStart + lane];
			acc.Sync();

			LimaForcecalc::computeCmapForces<BoundaryCondition, Accumulator>(bondsBuffer, std::min(batchSize, bondGroup->nCmapbonds - batchStart), positions, acc,
				bondGroups.cmapGridPoints, bondGroups.cmapGridSpacing, boxSize, boxSizeInv);
		}
	}


	{
		acc.Sync();
		PairBond* bondsBuffer = reinterpret_cast<PairBond*>(_bondsBuffer);
		for (int batchStart = firstBatch; batchStart < bondGroup->nPairbonds; batchStart += batchStride) {
			if (batchStart + lane < bondGroup->nPairbonds) {
				const int bondIndex = batchStart + lane;
				bondGroups.pairbonds.Load(bondGroup->indexOfFirstPairbond + bondIndex, bondsBuffer[lane]);
			}
			acc.Sync();

			LimaForcecalc::computePairbondForces<BoundaryCondition, Accumulator>(bondsBuffer, std::min(batchSize, bondGroup->nPairbonds - batchStart), positions, acc, boxSize, boxSizeInv, ewaldKappa);
		}
	}

	GroupSync();
	if (!loadsAndStores) return; // No barriers below
	float virialLocal = 0.f;
#pragma unroll
	for (int i = 0; i < 2; i++) {
		const int pid = lane + i * 32;
		if (pid < bondGroup->nParticles) {
			const ForceEnergy result = acc.Get(pid);
			if constexpr (computeVirial)
				virialLocal += (positions[pid] - positions[0]).dot(result.force);
			if constexpr (emVariant)
				forceEnergiesOut[bondGroup->indexOfFirstParticle + pid] = result;
			else {
				const int extra = extraResultIndices[bondGroup->indexOfFirstParticle + pid];
				if (extra < 0) acc.StorePrimary(primaryBondForces, pid, slots[i]);
				else acc.StoreTo(extraBondForces, pid, extra);
			}
		}
	}
	if constexpr (computeVirial) {
		for (int offset = 16; offset > 0; offset >>= 1)
			virialLocal += __shfl_down_sync(LIMA_FULL_WARP_MASK, virialLocal, offset);
		if (lane == 0)
			virialPartials[groupId] = Virial::ToFixed(virialLocal);
	}
}



namespace NbNonlocal {
	// This lane's row of an entry's noInteractions with an own quarter: bit k is j particle (lane+k)&3. The task builder stores
	// the rows rotated, and all set for quarters the entry does not interact with
	__device__ inline uint32_t LaneNoInteractions(uint32_t noInteractions, int lane) {
		return (noInteractions >> (lane * 4)) & 0xFu;
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
		out.force.x = __shfl_sync(LIMA_FULL_WARP_MASK, fe.force.x, srcLane, 4);
		out.force.y = __shfl_sync(LIMA_FULL_WARP_MASK, fe.force.y, srcLane, 4);
		out.force.z = __shfl_sync(LIMA_FULL_WARP_MASK, fe.force.z, srcLane, 4);
		if constexpr (withPotE)
			out.potE = __shfl_sync(LIMA_FULL_WARP_MASK, fe.potE, srcLane, 4);
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

// Nonbonded forces. Each warp computes the QuarterEntries of one supercluster, see EmitQuarterEntriesKernel. In MD the entries of
// a small simulation's superclusters are split into contiguous chunks over several warps, see Engine::UploadNbChunks: such a
// simulation has too few superclusters to fill the GPU, and a warp's serial walk over its entries bounds the step. The grid has
// maxChunks warps per supercluster, and those beyond its simulation's chunk count exit.
// blockDim = 64, 2 warps per block. computePotE must match logData of the following SuperclusterIntegrateKernel.
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
template <typename BoundaryCondition, bool energyMinimize, bool computePotE, bool computeVirial>
__global__ void __launch_bounds__(64, computePotE ? 12 : 16) NbNonlocalKernel(const SuperCluster* const superClusters,
	const QuarterEntryTask* const tasks, const QuarterEntry* const entries,
	const ForceAccumulator forceAcc /*MD only*/, SCResult* const results /*EM only*/, const int* const entryResultIndices /*EM only*/,
	const SuperClusterMeta* const scMeta, Float3 boxSize, Float3 boxSizeInv, float ewaldKappa /*EM*/, const PhysicsUtilsDevice::EwaldRealSpace ewald /*MD*/,
	float cutoffSq, int nSuperclusters, int maxChunks /*1 in EM*/, const int* const simulationChunks /*MD only*/,
	long long* const virialPartials /*MD pressure steps only*/)
{
	using namespace NbNonlocal;
	constexpr int nGroups = 8;
	const int warpInBlock = threadIdx.x >> 5;
	const int warpId = blockIdx.x * 2 + warpInBlock;
	const int scId = warpId / maxChunks;
	const int chunk = warpId - scId * maxChunks;
	if (scId >= nSuperclusters)
		return; // Only warp level syncs below
	const int laneInWarp = threadIdx.x & 31;
	const int nChunks = energyMinimize || maxChunks == 1 ? 1 : simulationChunks[scMeta[scId].simulationId];
	if (chunk >= nChunks) {
		if constexpr (computeVirial)
			if (laneInWarp == 0) virialPartials[warpId] = 0;
		return;
	}
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
	float virialLocal = 0.f;

	// The warp's chunk of the task's entries, in whole groups
	const int groupsPerChunk = ((task.nEntries + nGroups - 1) / nGroups + nChunks - 1) / nChunks;
	const int begin = task.firstEntry + chunk * groupsPerChunk * nGroups;
	const int end = min(task.firstEntry + task.nEntries, begin + groupsPerChunk * nGroups);
	for (int base = begin; base < end; base += nGroups) {
		const int e = base + group;
		const bool activeGroup = e < end;
		QuarterEntry entry = entries[activeGroup ? e : base];
		const uint32_t ownQuarterMask = activeGroup ? entry.ownQuarterMask : 0u;
		if (!activeGroup)
			for (int q = 0; q < 4; q++)
				entry.noInteractions[q] = 0xFFFFu; // The group still runs the warp's quarters, with all pairs masked

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

		const uint32_t unionMask = __reduce_or_sync(LIMA_FULL_WARP_MASK, ownQuarterMask);
		ForceEnergy jForce{}; // On j particle jIndex

		if constexpr (!energyMinimize) {
			ForceEnergy accJ[4]{};
#pragma unroll
			for (int ownQuarter = 0; ownQuarter < 4; ownQuarter++) {
				if (!((unionMask >> ownQuarter) & 1))
					continue;
				const uint32_t noInteractions = LaneNoInteractions(entry.noInteractions[ownQuarter], lane);
				const float4 pqI = ownPosCharge[warpInBlock][ownQuarter * 4 + lane];
				const float2 seI = ownSigmaEpsilon[warpInBlock][ownQuarter * 4 + lane];
#pragma unroll
				for (int k = 0; k < 4; k++) {
					const int jLocal = (lane + k) & 3;
					const float4 pqJ = jPosCharge[warpInBlock][group][jLocal];
					const float2 seJ = jSigmaEpsilon[warpInBlock][group][jLocal];
					const Float3 diff{ pqI.x - pqJ.x, pqI.y - pqJ.y, pqI.z - pqJ.z };
					const bool masked = ((noInteractions >> k) & 1) | (diff.lenSquared() >= cutoffSq);
					ForceEnergy fe = LJ::ComputePairNB<computePotE>(diff, LJ::CalcSigma(seI.x, seJ.x), LJ::CalcEpsilon(seI.y, seJ.y), pqI.w * pqJ.w, masked, ewald);
					if constexpr (computeVirial)
					virialLocal -= diff.dot(fe.force);
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
				const uint32_t noInteractions = LaneNoInteractions(NoInteractionsOfQuarter(entry, ownQuarter), lane);
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
			accOwn[q].force.x += __shfl_xor_sync(LIMA_FULL_WARP_MASK, accOwn[q].force.x, offset);
			accOwn[q].force.y += __shfl_xor_sync(LIMA_FULL_WARP_MASK, accOwn[q].force.y, offset);
			accOwn[q].force.z += __shfl_xor_sync(LIMA_FULL_WARP_MASK, accOwn[q].force.z, offset);
			if constexpr (computePotE)
				accOwn[q].potE += __shfl_xor_sync(LIMA_FULL_WARP_MASK, accOwn[q].potE, offset);
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

	if constexpr (computeVirial) {
		for (int offset = 16; offset > 0; offset >>= 1)
			virialLocal += __shfl_down_sync(LIMA_FULL_WARP_MASK, virialLocal, offset);
		if (laneInWarp == 0)
			virialPartials[warpId] = Virial::ToFixed(virialLocal);
	}
}

// Particles edited live, indexed by global particle id. Each is nullptr unless in use
struct LiveEditBuffers {
	const Float3* fixedParticleMovement = nullptr;
	const Float3* forceMask = nullptr;
	const Rotation* fixedParticleRotation = nullptr;

	__host__ __device__ bool Any() const { return fixedParticleMovement || forceMask || fixedParticleRotation; }
};

// The step's index for the integration kernels. A replayed step graph passes the step it was captured at, so the counter is advanced
// and only raised to the passed step: that is the step itself when captured, and the previous step plus one when replayed
__global__ void AdvanceStepCounterKernel(int64_t* const counter, int64_t step) {
	*counter = max(*counter + 1, step);
}

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
	Float3 boxSize, float* const forcesMagnitudeSquared, const LiveEditBuffers liveEdit, const LogBuffers log, const CapacityOverflow overflow, const int64_t* const stepIndex,
	const ThermostatImpulse* const thermostatImpulses)
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

	// Negated, so non-finite steps are caught too
	if (!((posNext - pos).lenSquared() < CapacityOverflow::maxStep * CapacityOverflow::maxStep) && !liveEdit.Any())
		overflow.Report(CapacityOverflow::Unstable, scMeta[scId].globalParticleIds[threadIdx.x], static_cast<int>(*stepIndex), 0);

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

	const ThermostatImpulse impulse = thermostatImpulses[scMeta[scId].simulationId];
	state.velocity = velocity * (impulse.step == *stepIndex ? impulse.scalar : simulation.thermostatScalar);
	state.forcePrev = fe.force;
	states[slot] = state;

	if constexpr (logData)
		log.Log(simulation, pcId, indexInPc, posNext, fe, state.velocity.len());

	superClusters[scId].SetPosition(threadIdx.x, posNext);
}

namespace Constraints {
	constexpr int maxIterations = 100;
	// SHAKE's relative error of the squared constraint length
	constexpr float tolerance = 1e-6f;

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

	// Removes the velocity components along the constraints (RATTLE), at the constrained positions. The corrections are linear in
	// the constraints' multipliers, so with at most 3 constraints they are solved directly, as a 3x3 system, rather than iterated
	__device__ inline void ConstrainVelocities(const PclusterConstraints& constraints, const Group& group, const Float3* const positions,
		Float3* const velocities, const float* const invMasses) {
		constexpr int maxC = PclusterConstraints::maxConstraints;
		static_assert(maxC == 3);
		int a[maxC], b[maxC];
		Float3 r[maxC];
		for (int c = 0; c < constraints.nConstraints; c++) {
			a[c] = group.slots[constraints.particles[c] & 3];
			b[c] = group.slots[constraints.particles[c] >> 2];
			r[c] = positions[a[c]] - positions[b[c]];
		}
		// Row c: the change of the velocity along constraint c per unit multiplier of constraint k, which moves a[k] by +r[k]/m
		// and b[k] by -r[k]/m. Unused rows, and those of constraints between two particles fixed by live edits, are identity
		float m[maxC][maxC] = { { 1.f, 0.f, 0.f }, { 0.f, 1.f, 0.f }, { 0.f, 0.f, 1.f } };
		float rhs[maxC] = {};
		for (int c = 0; c < constraints.nConstraints; c++) {
			if (invMasses[a[c]] + invMasses[b[c]] == 0.f)
				continue;
			for (int k = 0; k < constraints.nConstraints; k++) {
				const float side = (a[c] == a[k] ? invMasses[a[c]] : a[c] == b[k] ? -invMasses[a[c]] : 0.f)
					- (b[c] == a[k] ? invMasses[b[c]] : b[c] == b[k] ? -invMasses[b[c]] : 0.f);
				m[c][k] = side * r[c].dot(r[k]);
			}
			rhs[c] = -r[c].dot(velocities[a[c]] - velocities[b[c]]);
		}
		// Cramer's rule
		const float det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
			+ m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
		if (!(fabsf(det) > 0.f)) {
			ReportFailure("RATTLE", 0);
			return;
		}
		float multipliers[maxC];
		for (int k = 0; k < maxC; k++) {
			float mk[maxC][maxC];
			for (int i = 0; i < maxC; i++)
				for (int j = 0; j < maxC; j++)
					mk[i][j] = j == k ? rhs[i] : m[i][j];
			multipliers[k] = (mk[0][0] * (mk[1][1] * mk[2][2] - mk[1][2] * mk[2][1]) - mk[0][1] * (mk[1][0] * mk[2][2] - mk[1][2] * mk[2][0])
				+ mk[0][2] * (mk[1][0] * mk[2][1] - mk[1][1] * mk[2][0])) / det;
		}
		for (int k = 0; k < constraints.nConstraints; k++) {
			velocities[a[k]] += r[k] * (multipliers[k] * invMasses[a[k]]);
			velocities[b[k]] -= r[k] * (multipliers[k] * invMasses[b[k]]);
		}
	}

	// A water made rigid by its 3 constraints: an oxygen and two hydrogens of equal mass, at equal distance from it. Slots in the group's
	// supercluster
	struct Water {
		int o, h1, h2;
		float distanceOH, distanceHH;	// [nm]
	};
	// Whether the constraints are such a water, with no particle fixed by live edits
	__device__ inline bool FindWater(const PclusterConstraints& constraints, const Group& group, const float* const invMasses, Water& water) {
		if (constraints.nConstraints != 3)
			return false;
		for (int hh = 0; hh < 3; hh++) {
			const int c1 = (hh + 1) % 3, c2 = (hh + 2) % 3;
			if (constraints.distance[c1] != constraints.distance[c2])
				continue;
			const int h1 = constraints.particles[hh] & 3, h2 = constraints.particles[hh] >> 2;
			const int a1 = constraints.particles[c1] & 3, b1 = constraints.particles[c1] >> 2;
			const int a2 = constraints.particles[c2] & 3, b2 = constraints.particles[c2] >> 2;
			const int o = (a1 == a2 || a1 == b2) ? a1 : (b1 == a2 || b1 == b2) ? b1 : -1;
			if (o < 0 || o == h1 || o == h2 || (a1 == o ? b1 : a1) == (a2 == o ? b2 : a2))
				continue;
			water = { group.slots[o], group.slots[h1], group.slots[h2], constraints.distance[c1], constraints.distance[hh] };
			return invMasses[water.o] > 0.f && invMasses[water.h1] > 0.f && invMasses[water.h1] == invMasses[water.h2];
		}
		return false;
	}

	// SETTLE (Miyamoto & Kollman 1992, as GROMACS' settle.cpp): moves a rigid water's steps so its positions after the step satisfy
	// its constraints, analytically. Expects the positions before the step to satisfy them, which the previous step ensured. Works
	// relative to the oxygen before the step, so the values stay small and precise. Returns false if the step is too large to solve
	__device__ inline bool Settle(const Water& water, const Float3* const positions, Float3* const steps, const float* const invMasses) {
		const float massO = 1.f / invMasses[water.o], massH = 1.f / invMasses[water.h1];
		const float totalMass = massO + 2.f * massH;
		// The canonical water around its center of mass: O at (0, ra, 0), the hydrogens at (-+rc, -rb, 0)
		const float rc = 0.5f * water.distanceHH;
		const float height = sqrtf(water.distanceOH * water.distanceOH - rc * rc);
		const float ra = 2.f * massH * height / totalMass;
		const float rb = height - ra;

		const Float3 b0 = positions[water.h1] - positions[water.o];
		const Float3 c0 = positions[water.h2] - positions[water.o];
		// The unconstrained positions after the step, around their center of mass
		const Float3 aNew = steps[water.o], bNew = b0 + steps[water.h1], cNew = c0 + steps[water.h2];
		const Float3 com = (aNew * massO + (bNew + cNew) * massH) * (1.f / totalMass);
		const Float3 a1 = aNew - com, b1 = bNew - com, c1 = cNew - com;

		// A frame with z normal to the water's plane before the step, and x perpendicular to that and the new oxygen
		Float3 ez = b0.cross(c0);
		Float3 ex = a1.cross(ez);
		Float3 ey = ez.cross(ex);
		ex = ex * rsqrtf(ex.lenSquared());
		ey = ey * rsqrtf(ey.lenSquared());
		ez = ez * rsqrtf(ez.lenSquared());
		const float xb0 = ex.dot(b0), yb0 = ey.dot(b0), xc0 = ex.dot(c0), yc0 = ey.dot(c0);
		const float za1 = ez.dot(a1);
		const float xb1 = ex.dot(b1), yb1 = ey.dot(b1), zb1 = ez.dot(b1);
		const float xc1 = ex.dot(c1), yc1 = ey.dot(c1), zc1 = ez.dot(c1);

		// The tilts out of the plane: phi of the oxygen, psi of the H-H axis
		const float sinPhi = za1 / ra;
		const float cosPhiSq = 1.f - sinPhi * sinPhi;
		if (!(cosPhiSq > 0.f))
			return false;
		const float cosPhi = sqrtf(cosPhiSq);
		const float sinPsi = (zb1 - zc1) / (2.f * rc * cosPhi);
		const float cosPsi = sqrtf(fmaxf(1.f - sinPsi * sinPsi, 0.f));
		const float ya2 = ra * cosPhi;
		const float xb2 = -rc * cosPsi;
		const float t1 = -rb * cosPhi, t2 = rc * sinPsi * sinPhi;
		const float yb2 = t1 - t2, yc2 = t1 + t2;

		// The rotation theta in the plane, from the conservation of angular momentum
		const float alpha = xb2 * (xb0 - xc0) + yb0 * yb2 + yc0 * yc2;
		const float beta = xb2 * (yc0 - yb0) + xb0 * yb2 + xc0 * yc2;
		const float gamma = xb0 * yb1 - xb1 * yb0 + xc0 * yc1 - xc1 * yc0;
		const float alphaBetaSq = alpha * alpha + beta * beta;
		const float rootSq = alphaBetaSq - gamma * gamma;
		if (!(rootSq > 0.f))
			return false;
		const float sinTheta = (alpha * gamma - beta * sqrtf(rootSq)) / alphaBetaSq;
		const float cosTheta = sqrtf(fmaxf(1.f - sinTheta * sinTheta, 0.f));

		const Float3 a3{ -ya2 * sinTheta, ya2 * cosTheta, za1 };
		const Float3 b3{ xb2 * cosTheta - yb2 * sinTheta, xb2 * sinTheta + yb2 * cosTheta, zb1 };
		const Float3 c3{ -xb2 * cosTheta - yc2 * sinTheta, -xb2 * sinTheta + yc2 * cosTheta, zc1 };
		const auto toBox = [&](const Float3& v) { return ex * v.x + ey * v.y + ez * v.z; };
		steps[water.o] = com + toBox(a3);
		steps[water.h1] = com + toBox(b3) - b0;
		steps[water.h2] = com + toBox(c3) - c0;
		return true;
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

namespace ConstrainedCell {
	constexpr int n = SuperCluster::maxParticles;
	constexpr int maxScs = SuperClustersControl::maxClustersPerBlock;

	// The particles of one pcluster are all in the superclusters of one cell, see ClusteringKernel, so the constrained kernels take a
	// cell per block, with thread (i, k) on particle i of the cell's k'th supercluster. The cell's superclusters have consecutive ids
	struct Cell {
		int firstScId = 0;
		int nScs = 0;
		__device__ Cell(const SuperClustersControl& scControl, int cell)
			: nScs(scControl.nSuperclustersInBlocks[cell]) {
			firstScId = nScs > 0 ? scControl.scIdsInBlocks[cell * maxScs] : 0;
		}
	};

	// The pcluster's particles as indices into the cell's arrays, ordered as in the pcluster. Returns false, after reporting it,
	// if a constrained particle is not in the cell, which only a broken layout would do
	__device__ inline bool FindGroup(const PclusterConstraints& constraints, int pcId, const int* const pclusterParticleSlots, const Cell& cell,
		Constraints::Group& group) {
		for (int k = 0; k < PersistentCluster::maxParticles; k++) {
			const int local = pclusterParticleSlots[pcId * PersistentCluster::maxParticles + k] - cell.firstScId * n;
			group.slots[k] = local >= 0 && local < cell.nScs * n ? local : -1;
		}
		for (int c = 0; c < constraints.nConstraints; c++)
			if (group.slots[constraints.particles[c] & 3] < 0 || group.slots[constraints.particles[c] >> 2] < 0) {
				printf("Constrained particles of pcluster %d are not in the same cell\n", pcId);
				return false;
			}
		return true;
	}
}

// Puts the positions on the constraints before the first MD step, like GROMACS constrains the starting coordinates. Otherwise
// the first step would, and its jump in energy would show in the energy of the run. gridDim = nGridnodes, blockDim = (16, 16, 1),
// see ConstrainedCell
template<typename BoundaryCondition>
__global__ void ConstrainPositionsKernel(const SuperClustersControl scControl, const ParticleIntegrationState* const states,
	const PclusterConstraints* const pclusterConstraints, const int* const pclusterParticleSlots, Float3 boxSize)
{
	using namespace ConstrainedCell;
	LIMA_SHARED_UNINIT(Float3, positions, [maxScs * n]);
	LIMA_SHARED_UNINIT(Float3, steps, [maxScs * n]);
	__shared__ float invMasses[maxScs * n];
	const Cell cell(scControl, blockIdx.x);
	if (cell.nScs == 0)
		return;
	const int local = threadIdx.y * n + threadIdx.x;
	const int scId = cell.firstScId + threadIdx.y;
	const bool inRange = threadIdx.y < cell.nScs;
	ParticleIntegrationState state{};
	state.mass = 0.f;
	if (inRange)
		state = states[scId * n + threadIdx.x];
	const bool active = inRange && state.mass != 0.f;

	// The cell is kept in the periodic image of its first particle
	Float3 pos{};
	if (active) {
		Float3 p0 = scControl.scData[cell.firstScId].Position(0);
		BoundaryCondition::applyBCNM(p0, boxSize);
		pos = scControl.scData[scId].Position(threadIdx.x);
		BoundaryCondition::applyHyperposNM(p0, pos, boxSize);
	}
	positions[local] = pos;
	steps[local] = Float3{};
	invMasses[local] = active ? 1.f / state.mass : 0.f;
	__syncthreads();

	const int pcId = state.pclusterParticle / PersistentCluster::maxParticles;
	if (active && state.pclusterParticle % PersistentCluster::maxParticles == 0) {
		const PclusterConstraints constraints = pclusterConstraints[pcId];
		Constraints::Group group{};
		if (constraints.nConstraints > 0 && FindGroup(constraints, pcId, pclusterParticleSlots, cell, group))
			Constraints::ConstrainSteps(constraints, group, positions, steps, invMasses);
	}
	__syncthreads();

	if (active)
		scControl.scData[scId].SetPosition(threadIdx.x, pos + steps[local]);
}

// SuperclusterIntegrateKernel with distance constraints (settled water, constraints = h-bonds), as RATTLE: the velocities are made
// orthogonal to the constraints, and SETTLE (water) or SHAKE puts the positions after the step back on them. The position correction is added to forcePrev as
// the constraint force, so the next step's half kick includes it, while the stored velocity stays the constrained full-step
// velocity. A constrained group never spans pclusters, and each pcluster's constraints are solved by the thread of its first
// particle in a fixed order, so the result is deterministic. gridDim = nGridnodes, blockDim = (16, 16, 1), see ConstrainedCell
template<typename BoundaryCondition, bool logData, bool computeVirial>
__global__ void SuperclusterIntegrateConstrainedKernel(const ForceAccumulator forceAcc, const ForceAccumulator primaryBondForces, const ulonglong4* const extraBondForces, const int* const bondReferences,
	const SuperClustersControl scControl, ParticleIntegrationState* const states, const IntegrationSimulationData* const simulationData, int nSuperclusters,
	Float3 boxSize, float* const forcesMagnitudeSquared, const LiveEditBuffers liveEdit, const LogBuffers log, const PclusterConstraints* const pclusterConstraints,
	const int* const pclusterParticleSlots, const CapacityOverflow overflow, const int64_t* const stepIndex, long long* const virialPartials,
	const ThermostatImpulse* const thermostatImpulses)
{
	using namespace ConstrainedCell;
	LIMA_SHARED_UNINIT(Float3, positions, [maxScs * n]);	// Before the step
	LIMA_SHARED_UNINIT(Float3, velocities, [maxScs * n]);	// Then the steps
	__shared__ float invMasses[maxScs * n];
	const Cell cell(scControl, blockIdx.x);
	if (cell.nScs == 0)
		return;
	const int local = threadIdx.y * n + threadIdx.x;
	const int scId = cell.firstScId + threadIdx.y;
	const int slot = scId * n + threadIdx.x;
	const SuperClusterMeta* const scMeta = scControl.scMeta;

	// No early returns, as the threads synchronize below
	const bool inRange = threadIdx.y < cell.nScs;
	ParticleIntegrationState state{};
	state.mass = 0.f;
	if (inRange)
		state = states[slot];
	const bool active = inRange && state.mass != 0.f;	// Unused slots have no mass

	const int pcId = state.pclusterParticle / PersistentCluster::maxParticles;
	const int indexInPc = state.pclusterParticle % PersistentCluster::maxParticles;
	const IntegrationSimulationData* const simulation = active ? &simulationData[scMeta[scId].simulationId] : nullptr;
	const int pidGlobal = active ? scMeta[scId].globalParticleIds[threadIdx.x] : -1;
	// The velocity scale of this step, see ThermostatImpulse
	float thermostatScalar = 1.f;
	if (active) {
		const ThermostatImpulse impulse = thermostatImpulses[scMeta[scId].simulationId];
		thermostatScalar = impulse.step == *stepIndex ? impulse.scalar : simulation->thermostatScalar;
	}

	Float3 pos{}, velocity{};
	ForceEnergy fe{};
	Float3 p0{};
	if (active) {
		// The cell is kept in the periodic image of its first particle. Its thread writes that position only after the syncs below
		p0 = scControl.scData[cell.firstScId].Position(0);
		BoundaryCondition::applyBCNM(p0, boxSize);
		pos = scControl.scData[scId].Position(threadIdx.x);
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
	positions[local] = pos;
	velocities[local] = velocity;
	invMasses[local] = active && !fixed ? 1.f / state.mass : 0.f;
	__syncthreads();

	PclusterConstraints constraints{};
	Constraints::Group group{};
	Constraints::Water water{};
	bool isWater = false;
	if (active && indexInPc == 0) {
		constraints = pclusterConstraints[pcId];
		if (constraints.nConstraints > 0 && !FindGroup(constraints, pcId, pclusterParticleSlots, cell, group))
			constraints.nConstraints = 0;
		if (constraints.nConstraints > 0) {
			Constraints::ConstrainVelocities(constraints, group, positions, velocities, invMasses);
			isWater = Constraints::FindWater(constraints, group, invMasses, water);
		}
	}
	__syncthreads();

	velocity = velocities[local];
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
	velocities[local] = step;
	__syncthreads();

	// SETTLE for rigid water, SHAKE for the rest, and for a water whose step is too large for SETTLE, which then reports it
	if (constraints.nConstraints > 0 && !(isWater && Constraints::Settle(water, positions, velocities, invMasses)))
		Constraints::ConstrainSteps(constraints, group, positions, velocities, invMasses);
	__syncthreads();

	const Float3 constrainedStep = velocities[local];
	const Float3 posNext = pos + constrainedStep;
	if (active && !(constrainedStep.lenSquared() < CapacityOverflow::maxStep * CapacityOverflow::maxStep) && !liveEdit.Any())
		overflow.Report(CapacityOverflow::Unstable, pidGlobal, static_cast<int>(*stepIndex), 0);
	// The correction to the step is dt^2/2m of the constraint force. Fixed particles are not corrected, so they get none
	const Float3 constraintForce = active
		? (constrainedStep - step) * (2.f * state.mass / (simulation->dt * simulation->dt)) : Float3{};
	float constraintVirial = active ? (pos - p0).dot(constraintForce) : 0.f;
	if constexpr (computeVirial) {
		for (int offset = 8; offset > 0; offset >>= 1)
			constraintVirial += __shfl_down_sync(LIMA_FULL_WARP_MASK, constraintVirial, offset, 16);
		if (inRange && threadIdx.x == 0)
			virialPartials[scId] = Virial::ToFixed(constraintVirial);
	}
	if (!active)
		return;

	state.velocity = velocity * thermostatScalar;
	state.forcePrev = fe.force + constraintForce;
	states[slot] = state;

	if constexpr (logData)
		log.Log(*simulation, pcId, indexInPc, posNext, fe, state.velocity.len());

	scControl.scData[scId].SetPosition(threadIdx.x, posNext);
}

// Isotropic pressure coupling preserves fractional coordinates. Integration states stay in the same slots and velocities are
// unchanged; the following neighbor-list rebuild moves the scaled particles to their new grid cells.
__global__ void ScaleBoxKernel(SuperCluster* superClusters, int nSuperclusters, float coordinateScale) {
	const int slot = blockIdx.x * blockDim.x + threadIdx.x;
	if (slot >= nSuperclusters * SuperCluster::maxParticles)
		return;
	const int scId = slot / SuperCluster::maxParticles;
	const int lane = slot % SuperCluster::maxParticles;
	Float3 position = superClusters[scId].Position(lane);
	superClusters[scId].SetPosition(lane, position * coordinateScale);
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


