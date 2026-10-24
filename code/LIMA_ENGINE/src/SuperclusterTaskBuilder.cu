// Builds the supercluster interaction tasks, Engine::MakeSuperClusterTasksGPU. A translation unit of its own, compiled
// in parallel with Engine.cu
//
// The neighbor search, FindNeighborsKernel, lists the 4x4 blocks of particle quarters with pairs in range, and
// EmitQuarterEntriesKernel turns them into the QuarterEntries of NbNonlocalKernel. EM additionally needs a result slot for every
// supercluster pair, see MakeNbTasksEM.

#include "Engine.cuh"
#include "TaskBuilderControl.cuh"
#include "SimulationData.h"
#include "BoundaryCondition.cuh"
#include "DeviceAlgorithms.cuh"
#include "platform/DeviceIntrinsics.cuh"
#include "LimaPositionSystem.cuh"
#include "Utilities.h"

#include <cfloat>

namespace {
	constexpr int maxNeighborsPerSc = TaskBuilderControlContents::maxTasksPerSc;
	constexpr int neighborSearchWarpsPerBlock = 4;

	// EM's lists include 4x4 blocks that are within this distance beyond the cutoff. MD's buffer is set from the verlet buffer
	// tolerance, see EngineBatchData::listBuffer [nm]
	constexpr float emListBuffer = 0.05f;
}

__device__ inline bool Warp_ScAreBonded(const SuperClusterMeta& sc0, const SuperClusterMeta& sc1, const PclustersBondedToPcluster* const pclustersBondedToPcluster) {
	constexpr WarpMask mask = LIMA_FULL_WARP_MASK;

	const int lane = threadIdx.x & 31;
	const int nPairs = sc0.nUniquePcIds * sc1.nUniquePcIds;

	bool bonded = false;

	for (int pairId = lane; pairId < nPairs; pairId += 32) {
		const int i = pairId / sc1.nUniquePcIds;
		const int j = pairId % sc1.nUniquePcIds;

		bonded |= pclustersBondedToPcluster[sc0.uniquePclusterIds[i]].Contains(sc1.uniquePclusterIds[j]);
	}

	return __any_sync(mask, bonded);
}

// Whether a pcluster is split over both superclusters, see ClusteringKernel
__device__ inline bool Warp_ScsSharePcluster(const SuperClusterMeta& sc0, const SuperClusterMeta& sc1) {
	const int lane = threadIdx.x & 31;
	const int own = lane < sc0.nUniquePcIds ? sc0.uniquePclusterIds[lane] : -1;
	bool shared = false;
	for (int j = 0; j < sc1.nUniquePcIds; j++)
		shared |= own == sc1.uniquePclusterIds[j];
	return __any_sync(LIMA_FULL_WARP_MASK, shared);
}

// Bits of TaskBuilderControl::scExternalBonds
constexpr uint8_t bondsToOtherScs = 1;		// A pcluster of the supercluster is bonded to a pcluster outside it
constexpr uint8_t splitPclusters = 2;		// A pcluster of the supercluster is split over it and other superclusters of its cell



// ------------------------------------------------------- Neighbor search ------------------------------------------------------- //

// gridDim = nGridnodes, blockDim = 32. Bounding sphere, cell, valid particles and external bonds of every supercluster, and the AABB of every cell
__global__ void SuperclusterBoundsKernel(const SuperClustersControl scControl, const PclustersBondedToPcluster* const pclustersBondedToPcluster,
	const PersistentClusterMeta* const pclusterMeta, float4* const scSpheres, int* const scCells, uint16_t* const scValidMasks, uint8_t* const scExternalBonds, float4* const cellMin, float4* const cellMax)
{
	static_assert(sizeof(PclustersBondedToPcluster) == 32 * sizeof(int), "One lane per value of a bonded set");
	const int cell = blockIdx.x;
	const int lane = threadIdx.x;
	const int n = scControl.nSuperclustersInBlocks[cell];
	Float3 lo{ FLT_MAX }, hi{ -FLT_MAX };
	for (int k = 0; k < n; k++) {
		const int scId = scControl.scIdsInBlocks[cell * SuperClustersControl::maxClustersPerBlock + k];
		const int i = lane & 15;
		const bool valid = scControl.scData[scId].Valid(i);
		const Float3 p = scControl.scData[scId].Position(i);
		Float3 sum = valid ? p : Float3{};
		float count = valid ? 1.f : 0.f;
		for (int offset = 8; offset > 0; offset >>= 1) {
			sum.x += __shfl_xor_sync(LIMA_FULL_WARP_MASK, sum.x, offset);
			sum.y += __shfl_xor_sync(LIMA_FULL_WARP_MASK, sum.y, offset);
			sum.z += __shfl_xor_sync(LIMA_FULL_WARP_MASK, sum.z, offset);
			count += __shfl_xor_sync(LIMA_FULL_WARP_MASK, count, offset);
		}
		const Float3 center = sum * (1.f / fmaxf(count, 1.f));
		float radius = valid ? (p - center).len() : 0.f;
		for (int offset = 8; offset > 0; offset >>= 1)
			radius = fmaxf(radius, __shfl_xor_sync(LIMA_FULL_WARP_MASK, radius, offset));
		const uint32_t validMask = __ballot_sync(LIMA_FULL_WARP_MASK, valid) & 0xFFFF;
		if (lane == 0) {
			scSpheres[scId] = float4{ center.x, center.y, center.z, radius };
			scCells[scId] = cell;
			scValidMasks[scId] = static_cast<uint16_t>(validMask);
		}
		// Whether any value in the bonded sets of the supercluster's pclusters is a pcluster outside it, or any of its pclusters is split
		// over superclusters. Superclusters with neither are bonded to no other supercluster, so FindNeighborsKernel can skip checking their pairs
		const int nUnique = scControl.scMeta[scId].nUniquePcIds;
		const int ownPcluster = lane < nUnique ? scControl.scMeta[scId].uniquePclusterIds[lane] : -1;
		bool externalBond = false, split = false;
		for (int u = 0; u < nUnique; u++) {
			const int bondedPcluster = pclustersBondedToPcluster[__shfl_sync(LIMA_FULL_WARP_MASK, ownPcluster, u)].Get(lane);
			bool internal = !PclustersBondedToPcluster::IsValue(bondedPcluster);
			for (int v = 0; v < nUnique; v++)
				internal |= __shfl_sync(LIMA_FULL_WARP_MASK, ownPcluster, v) == bondedPcluster;
			externalBond |= !internal;
		}
		const int pclusterOfSlot = valid && lane < SuperCluster::maxParticles ? scControl.scMeta[scId]._pclusterIds[i] : -1;
		const uint32_t samePcluster = __match_any_sync(LIMA_FULL_WARP_MASK, pclusterOfSlot);
		split = pclusterOfSlot >= 0 && __popc(samePcluster) < pclusterMeta[pclusterOfSlot].nParticles;
		externalBond = __any_sync(LIMA_FULL_WARP_MASK, externalBond);
		split = __any_sync(LIMA_FULL_WARP_MASK, split);
		if (lane == 0)
			scExternalBonds[scId] = (externalBond ? bondsToOtherScs : 0) | (split ? splitPclusters : 0);

		if (valid) {
			lo = Float3{ fminf(lo.x, p.x), fminf(lo.y, p.y), fminf(lo.z, p.z) };
			hi = Float3{ fmaxf(hi.x, p.x), fmaxf(hi.y, p.y), fmaxf(hi.z, p.z) };
		}
	}
	for (int offset = 16; offset > 0; offset >>= 1) {
		lo.x = fminf(lo.x, __shfl_xor_sync(LIMA_FULL_WARP_MASK, lo.x, offset));
		lo.y = fminf(lo.y, __shfl_xor_sync(LIMA_FULL_WARP_MASK, lo.y, offset));
		lo.z = fminf(lo.z, __shfl_xor_sync(LIMA_FULL_WARP_MASK, lo.z, offset));
		hi.x = fmaxf(hi.x, __shfl_xor_sync(LIMA_FULL_WARP_MASK, hi.x, offset));
		hi.y = fmaxf(hi.y, __shfl_xor_sync(LIMA_FULL_WARP_MASK, hi.y, offset));
		hi.z = fmaxf(hi.z, __shfl_xor_sync(LIMA_FULL_WARP_MASK, hi.z, offset));
	}
	if (lane == 0) {
		cellMin[cell] = float4{ lo.x, lo.y, lo.z, 0.f };
		cellMax[cell] = float4{ hi.x, hi.y, hi.z, 0.f };
	}
}

// blockDim = 32 * neighborSearchWarpsPerBlock, one warp per own supercluster.
// Lists the superclusters each supercluster owns an interaction with (queryId >= own id): those with a pair of particles within
// the list radius, and which of their 4x4 blocks of quarters have such a pair. Pairs within a self interaction count in one order only.
// entryCounts is the number of QuarterEntries EmitQuarterEntriesKernel will write, with allQuarters as passed to it.
// The cells from -cellRangeLo to cellRangeHi around the own cell are searched, but only those whose AABB is within reach, and of
// those only the superclusters whose bounding sphere is. Candidates are enumerated in a fixed order, so the lists are deterministic.
__global__ void __launch_bounds__(32 * neighborSearchWarpsPerBlock, 8) FindNeighborsKernel(const SuperClustersControl scControl, const float4* const scSpheres, const uint8_t* const scExternalBonds, const int* const scCells,
	const float4* const cellMin, const float4* const cellMax, const PclustersBondedToPcluster* const pclustersBondedToPcluster,
	ScNeighbor* const neighbors, int* const nNeighbors, int* const entryCounts, int* const overflow,
	int nSuperclusters, Int3 boxSize, Int3 cellRangeLo, Int3 cellRangeHi, float listRadius, bool allQuarters, Float3 boxSizeF, Float3 boxSizeInv)
{
	__shared__ float4 stage[neighborSearchWarpsPerBlock][SuperCluster::maxParticles];
	const int warpInBlock = threadIdx.x >> 5;
	const int lane = threadIdx.x & 31;
	const int scId = blockIdx.x * neighborSearchWarpsPerBlock + warpInBlock;
	if (scId >= nSuperclusters)
		return;

	const int ownIndex = lane & 15;
	const Float3 ownPos = scControl.scData[scId].Position(ownIndex);
	const bool ownValid = scControl.scData[scId].Valid(ownIndex);
	const float4 ownSphere = scSpheres[scId];
	const Float3 ownCenter{ ownSphere.x, ownSphere.y, ownSphere.z };
	const uint8_t ownExternalBonds = scExternalBonds[scId];
	const int ownCell = scCells[scId];
	const float listRadiusSq = listRadius * listRadius;
	// Padding is placed far away, own at +1e18 and staged candidates at -1e18, so the pair test needs no validity checks
	const Float3 ownTestPos = ownValid ? ownPos : Float3{ 1e18f };
	// Bit m: whether staged particle (lane >> 4) + 2m may pair with this own particle within a self interaction, which counts pairs in one order
	uint32_t selfMask = 0;
#pragma unroll
	for (int m = 0; m < 8; m++)
		if ((lane >> 4) + 2 * m > ownIndex) selfMask |= 1u << m;
	const uint32_t ownQuarterBit = 1u << ((ownIndex >> 2) * 4);

	const int cell = scCells[scId];
	const int nCellsPerSimulation = boxSize.x * boxSize.y * boxSize.z;
	const int gridOffset = (cell / nCellsPerSimulation) * nCellsPerSimulation;
	const NodeIndex cell3d = BoxGrid::Get3dIndex(cell - gridOffset, boxSize);
	const Int3 side{ cellRangeLo.x + cellRangeHi.x + 1, cellRangeLo.y + cellRangeHi.y + 1, cellRangeLo.z + cellRangeHi.z + 1 };
	const int nOffsets = side.x * side.y * side.z;

	int count = 0;
	int entryCount = 0;
	for (int base = 0; base < nOffsets; base += 32) {
		const int o = base + lane;
		int targetCell = 0;
		bool cellHit = false;
		if (o < nOffsets) {
			const NodeIndex relative{ o % side.x - cellRangeLo.x, (o / side.x) % side.y - cellRangeLo.y, o / (side.x * side.y) - cellRangeLo.z };
			NodeIndex target = cell3d + relative;
			PeriodicBoundaryCondition::applyBC(target, boxSize);
			targetCell = gridOffset + BoxGrid::Get1dIndex(target, boxSize);
			// Distance from the own sphere's center to the cell's AABB, in the image closest to the own center. Empty cells have lo > hi
			const float4 lo4 = cellMin[targetCell];
			const float4 hi4 = cellMax[targetCell];
			if (lo4.x <= hi4.x) {
				Float3 lo{ lo4.x, lo4.y, lo4.z }, hi{ hi4.x, hi4.y, hi4.z };
				const Float3 center = (lo + hi) * 0.5f;
				Float3 imageCenter = center;
				PeriodicBoundaryCondition::ApplyHyperpos(ownCenter, imageCenter, boxSizeF, boxSizeInv);
				lo += imageCenter - center;
				hi += imageCenter - center;
				const Float3 gap{ fmaxf(0.f, fmaxf(lo.x - ownCenter.x, ownCenter.x - hi.x)), fmaxf(0.f, fmaxf(lo.y - ownCenter.y, ownCenter.y - hi.y)),
					fmaxf(0.f, fmaxf(lo.z - ownCenter.z, ownCenter.z - hi.z)) };
				const float reach = listRadius + ownSphere.w;
				cellHit = gap.lenSquared() < reach * reach;
			}
		}

		uint32_t cellHits = __ballot_sync(LIMA_FULL_WARP_MASK, cellHit);
		while (cellHits) {
			const int src = __ffs(cellHits) - 1;
			cellHits &= cellHits - 1;
			const int candidateCell = __shfl_sync(LIMA_FULL_WARP_MASK, targetCell, src);
			const int nCandidates = scControl.nSuperclustersInBlocks[candidateCell];

			int candidateId = -1;
			bool coarseHit = false;
			float4 candidateSphere{};	// Center in the image nearest the own center
			if (lane < nCandidates) {
				candidateId = scControl.scIdsInBlocks[candidateCell * SuperClustersControl::maxClustersPerBlock + lane];
				if (candidateId >= scId) {
					const float4 s = scSpheres[candidateId];
					Float3 c{ s.x, s.y, s.z };
					PeriodicBoundaryCondition::ApplyHyperpos(ownCenter, c, boxSizeF, boxSizeInv);
					const float reach = listRadius + ownSphere.w + s.w;
					coarseHit = (c - ownCenter).lenSquared() < reach * reach;
					candidateSphere = float4{ c.x, c.y, c.z, s.w };
				}
			}

			uint32_t candidates = __ballot_sync(LIMA_FULL_WARP_MASK, coarseHit);
			while (candidates) {
				const int k = __ffs(candidates) - 1;
				candidates &= candidates - 1;
				const int queryId = __shfl_sync(LIMA_FULL_WARP_MASK, candidateId, k);
				const bool selfTask = queryId == scId;

				// Most candidates passing the sphere test have no pair in range. No pair is possible unless an own particle is within
				// reach of the candidate's sphere, which costs one test per own particle rather than loading the candidate's particles.
				// The slack keeps it conservative under rounding, so the lists are unchanged
				const Float3 candidateCenter{ __shfl_sync(LIMA_FULL_WARP_MASK, candidateSphere.x, k), __shfl_sync(LIMA_FULL_WARP_MASK, candidateSphere.y, k),
					__shfl_sync(LIMA_FULL_WARP_MASK, candidateSphere.z, k) };
				const float candidateReach = listRadius + __shfl_sync(LIMA_FULL_WARP_MASK, candidateSphere.w, k) + 1e-4f;
				if (!__any_sync(LIMA_FULL_WARP_MASK, ownValid && (ownPos - candidateCenter).lenSquared() < candidateReach * candidateReach))
					continue;

				if (lane < SuperCluster::maxParticles) {
					Float3 p = scControl.scData[queryId].Position(lane);
					PeriodicBoundaryCondition::ApplyHyperpos(ownCenter, p, boxSizeF, boxSizeInv);
					if (!scControl.scData[queryId].Valid(lane)) p = Float3{ -1e18f };
					stage[warpInBlock][lane] = float4{ p.x, p.y, p.z, 0.f };
				}
				__syncwarp();
				uint32_t inRange = 0;
#pragma unroll
				for (int m = 0; m < 8; m++) {
					const float4 q = stage[warpInBlock][(lane >> 4) + 2 * m];
					if ((ownTestPos - Float3{ q.x, q.y, q.z }).lenSquared() < listRadiusSq)
						inRange |= 1u << m;
				}
				inRange &= selfTask ? selfMask : 0xFFu;
				// Staged particles 2m and 2m+1 relative to the half warp are in query quarter m / 2
				uint32_t bits = 0;
#pragma unroll
				for (int jq = 0; jq < 4; jq++)
					if (inRange & (3u << (2 * jq))) bits |= ownQuarterBit << jq;
				bits = __reduce_or_sync(LIMA_FULL_WARP_MASK, bits);
				__syncwarp(); // Before the stage is reused

				if (bits) {
					// The pclusters split over superclusters are all within a cell
					const bool bonded = selfTask
						|| ((ownExternalBonds & bondsToOtherScs) && Warp_ScAreBonded(scControl.scMeta[scId], scControl.scMeta[queryId], pclustersBondedToPcluster))
						|| ((ownExternalBonds & splitPclusters) && scCells[queryId] == ownCell && Warp_ScsSharePcluster(scControl.scMeta[scId], scControl.scMeta[queryId]));
					if (lane == 0 && count < maxNeighborsPerSc)
						neighbors[scId * maxNeighborsPerSc + count] = ScNeighbor{ queryId, static_cast<uint16_t>(bits), static_cast<uint16_t>(bonded) };
					count++;
					entryCount += allQuarters && !selfTask ? 4 : __popc((bits | bits >> 4 | bits >> 8 | bits >> 12) & 0xF);
				}
			}
		}
	}

	if (lane == 0) {
		if (count > maxNeighborsPerSc) {
			atomicMax(overflow, count);
			count = maxNeighborsPerSc;
		}
		nNeighbors[scId] = count;
		entryCounts[scId] = entryCount;
	}
}



// ---------------------------------------------------------- Entries ----------------------------------------------------------- //

// Bit j: whether queryIds[j] is in the sorted set. One pass over the set for all four, ending at the first value beyond them
template <int size>
__device__ inline uint32_t BondedQueryMask(const StaticSet<size>& set, const int (&queryIds)[4]) {
	const int maxId = max(max(queryIds[0], queryIds[1]), max(queryIds[2], queryIds[3]));
	uint32_t mask = 0;
	for (int i = 0; i < size; i++) {
		const int value = set.Get(i);
		if (!StaticSet<size>::IsValue(value) || value > maxId)
			break;
#pragma unroll
		for (int j = 0; j < 4; j++)
			mask |= uint32_t(value == queryIds[j]) << j;
	}
	return mask;
}

// Order in which ownQuarterMask classes are emitted. Neighbouring classes differ in one bit (Gray code), so when a warp's
// 8 groups straddle two classes, the union of their masks is small. EM's entries without pairs in range (class 0) go last
__constant__ uint8_t quarterMaskEmitOrder[16] = { 1, 3, 2, 6, 7, 5, 4, 12, 13, 15, 14, 10, 11, 9, 8, 0 };

// gridDim = nSuperclusters, blockDim = 32. Writes the QuarterEntries of each supercluster: one per quarter of a neighbor with a pair
// within the list radius. Sorted by ownQuarterMask in quarterMaskEmitOrder, and within a class by (chunk of 32 neighbors, neighbor,
// query quarter), so the result is deterministic.
// EM (allQuarters) also writes entries for the quarters of other superclusters without pairs in range, so every quarter of their
// result gets written, the index of each entry's result, and the result range of each supercluster
template <bool allQuarters>
__global__ void EmitQuarterEntriesKernel(const SuperClustersControl scControl, const ScNeighbor* const neighbors,
	const int* const nNeighbors, const uint16_t* const scValidMasks, const ParticlesBondedToParticle* const particlesBondedToParticle,
	const int* const entryStarts, QuarterEntryTask* const entryTasks, QuarterEntry* const entries,
	const TaskBuilderControlContents tbContents /*EM only*/, int* const entryResultIndices /*EM only*/)
{
	constexpr uint32_t noEntry = 16;
	const int scId = blockIdx.x;
	const int lane = threadIdx.x;
	const int n = nNeighbors[scId];
	__shared__ int classOffsets[16];
	const uint32_t ownValid = scValidMasks[scId];

	// 4 bits per query quarter: the own quarters it has a pair within the list radius with
	auto OwnMasks = [](uint32_t quarterMask) {
		uint32_t ownMasks = 0;
#pragma unroll
		for (int jq = 0; jq < 4; jq++)
#pragma unroll
			for (int iq = 0; iq < 4; iq++)
				ownMasks |= ((quarterMask >> (iq * 4 + jq)) & 1u) << (jq * 4 + iq);
		return ownMasks;
	};
	// The class of an entry for a query quarter, or noEntry
	auto EntryClass = [scId](uint32_t ownMasks, int jq, int queryScId) {
		const uint32_t maskClass = (ownMasks >> (jq * 4)) & 0xF;
		return maskClass || (allQuarters && queryScId != scId) ? maskClass : noEntry;
	};

	// Count the entries of each class, then turn the counts into each class's first index
	if (lane < 16)
		classOffsets[lane] = 0;
	__syncwarp();
	for (int k = lane; k < n; k += 32) {
		const ScNeighbor neighbor = neighbors[scId * maxNeighborsPerSc + k];
		const uint32_t ownMasks = OwnMasks(neighbor.quarterMask);
		for (int jq = 0; jq < 4; jq++) {
			const uint32_t c = EntryClass(ownMasks, jq, neighbor.queryScId);
			if (c != noEntry) atomicAdd(&classOffsets[c], 1);
		}
	}
	__syncwarp();
	if (lane == 0) {
		int sum = entryStarts[scId];
		for (int r = 0; r < 16; r++) {
			const int c = quarterMaskEmitOrder[r];
			const int count = classOffsets[c];
			classOffsets[c] = sum;
			sum += count;
		}
		entryTasks[scId] = QuarterEntryTask{ entryStarts[scId], sum - entryStarts[scId] };
		if constexpr (allQuarters) {
			scControl.scMeta[scId].resultsStartIndex = tbContents.nResultsPrefixsum[scId];
			scControl.scMeta[scId].nResults = tbContents.nResults[scId];
		}
	}
	// The own particles' ids, for the bonded exclusions
	__shared__ int ownParticleIds[SuperCluster::maxParticles];
	if (lane < SuperCluster::maxParticles)
		ownParticleIds[lane] = scControl.scMeta[scId].globalParticleIds[lane];
	__syncwarp();

	for (int base = 0; base < n; base += 32) {
		const int k = base + lane;
		const bool validNeighbor = k < n;
		const ScNeighbor neighbor = validNeighbor ? neighbors[scId * maxNeighborsPerSc + k] : ScNeighbor{ 0, 0, 0 };
		const uint32_t ownMasks = OwnMasks(neighbor.quarterMask);
		const int queryScId = neighbor.queryScId;
		const bool selfTask = queryScId == scId;
		const uint32_t queryValid = validNeighbor ? scValidMasks[queryScId] : 0;

		// The result this supercluster writes the query's forces to: after the query's own result, the results of the
		// superclusters owning an interaction with it are ordered by id, so find this one's position with a binary search
		int resultIndex = -1;
		if constexpr (allQuarters) {
			if (validNeighbor && !selfTask) {
				const int* const owners = &tbContents.scIdsQueryNonowned[queryScId * TaskBuilderControlContents::maxTasksPerSc];
				int lo = 0, hi = tbContents.nInteractionsNonowned[queryScId];
				while (lo < hi) {
					const int mid = (lo + hi) / 2;
					if (owners[mid] < scId) lo = mid + 1;
					else hi = mid;
				}
				resultIndex = tbContents.nResultsPrefixsum[queryScId] + 1 + lo;
			}
		}

		// Within a class, order by (neighbor, query quarter), so a warp's groups get the quarters of the same neighbor and their
		// force atomics share a cache line (the accumulator is SoA, a supercluster's 16 x components are one line)
		// classLanes[c][jq]: the lanes whose entry for query quarter jq has class c
		__shared__ uint32_t classLanes[noEntry + 1][4];
		for (int i = lane; i < (noEntry + 1) * 4; i += 32)
			classLanes[i / 4][i % 4] = 0;
		__syncwarp();
		uint32_t classes[4];
		for (int jq = 0; jq < 4; jq++) {
			classes[jq] = validNeighbor ? EntryClass(ownMasks, jq, queryScId) : noEntry;
			const uint32_t sameClass = __match_any_sync(LIMA_FULL_WARP_MASK, classes[jq]);
			if (lane == __ffs(sameClass) - 1)
				classLanes[classes[jq]][jq] = sameClass;
		}
		__syncwarp();
		const uint32_t lowerLanes = (1u << lane) - 1;
		int dsts[4];
		for (int jq = 0; jq < 4; jq++) {
			int rank = 0;
			for (int jq2 = 0; jq2 < 4; jq2++)
				rank += __popc(classLanes[classes[jq]][jq2] & lowerLanes) + (jq2 < jq && classes[jq2] == classes[jq]);
			dsts[jq] = classes[jq] != noEntry ? classOffsets[classes[jq]] + rank : 0;
		}
		__syncwarp();
		if (lane < 16)
			classOffsets[lane] += __popc(classLanes[lane][0]) + __popc(classLanes[lane][1]) + __popc(classLanes[lane][2]) + __popc(classLanes[lane][3]);
		__syncwarp();

		for (int jq = 0; jq < 4; jq++) {
			const uint32_t maskClass = classes[jq];
			const int dst = dsts[jq];
			if (maskClass == noEntry)
				continue;

			const uint32_t queryValid4 = (queryValid >> (jq * 4)) & 0xF;
			int queryParticleIds[4]{ -1, -1, -1, -1 };
			if (neighbor.bonded)
				for (int jLocal = 0; jLocal < 4; jLocal++)
					queryParticleIds[jLocal] = scControl.scMeta[queryScId].globalParticleIds[jq * 4 + jLocal];
			QuarterEntry out{};
			out.jScId = queryScId;
			out.jQuarter = static_cast<uint8_t>(jq);
			out.ownQuarterMask = static_cast<uint8_t>(maskClass);
#pragma unroll
			for (int iq = 0; iq < 4; iq++) {
				uint32_t noInteractions = 0;
#pragma unroll
				for (int iLocal = 0; iLocal < 4; iLocal++) {
					const bool valid = (ownValid >> (iq * 4 + iLocal)) & 1;
					const uint32_t row = valid ? (~queryValid4 & 0xFu) : 0xFu;
					noInteractions |= row << (iLocal * 4);
				}
				if (selfTask && iq == jq)
					noInteractions |= 0xF731; // Pairs with jLocal <= iLocal are computed in the other order
				if (neighbor.bonded && ((maskClass >> iq) & 1)) {
					for (int iLocal = 0; iLocal < 4; iLocal++) {
						const int pidOwn = ownParticleIds[iq * 4 + iLocal];
						if (pidOwn != -1)
							noInteractions |= BondedQueryMask(particlesBondedToParticle[pidOwn], queryParticleIds) << (iLocal * 4);
					}
				}
				// Rotate row iLocal so bit k is query particle (iLocal+k)&3, the one lane iLocal pairs with in iteration k of NbNonlocalKernel
				uint32_t rotated = 0;
#pragma unroll
				for (int iLocal = 0; iLocal < 4; iLocal++) {
					const uint32_t row = (noInteractions >> (iLocal * 4)) & 0xFu;
					rotated |= (((row | (row << 4)) >> iLocal) & 0xFu) << (iLocal * 4);
				}
				out.noInteractions[iq] = static_cast<uint16_t>((maskClass >> iq) & 1 ? rotated : 0xFFFFu);
			}
			entries[dst] = out;
			if constexpr (allQuarters)
				entryResultIndices[dst] = resultIndex;
		}
	}
}



// ------------------------------------------------------------- EM -------------------------------------------------------------- //

// gridDim = nSuperclusters, blockDim = 32. Lists each supercluster as an owner of an interaction with each of its neighbors.
// TaskBuilderControl::Reset must have zeroed nInteractionsNonowned
__global__ void DistributeEmInteractionsKernel(const ScNeighbor* const neighbors, const int* const nNeighbors, TaskBuilderControlContents tbContents, int* const overflow) {
	const int scId = blockIdx.x;
	const int n = nNeighbors[scId];
	for (int k = threadIdx.x; k < n; k += blockDim.x) {
		const int queryScId = neighbors[scId * maxNeighborsPerSc + k].queryScId;
		if (queryScId == scId)
			continue;
		const int index = atomicAdd(&tbContents.nInteractionsNonowned[queryScId], 1);
		if (index < TaskBuilderControlContents::maxTasksPerSc)
			tbContents.scIdsQueryNonowned[queryScId * TaskBuilderControlContents::maxTasksPerSc + index] = scId;
		else
			atomicMax(overflow, index + 1);
	}
}

// blockDim = maxTasksPerSc. The owners are placed by atomics, so they must be sorted to be deterministic
__global__ void SortEmInteractionsKernel(TaskBuilderControlContents tbContents) {
	const int scId = blockIdx.x;

	LAL::Sort(&tbContents.scIdsQueryNonowned[scId * TaskBuilderControlContents::maxTasksPerSc], TaskBuilderControlContents::maxTasksPerSc, [](const int& id) {
		return id;
		});

	if (threadIdx.x == 0)
		tbContents.nResults[scId] = 1 + tbContents.nInteractionsNonowned[scId]; // 1 for its own forces, 1 per supercluster owning an interaction with it
}



// ------------------------------------------------------------ Host ------------------------------------------------------------- //

void Engine::FindSuperclusterNeighbors(cudaStream_t stream, float listRadius, bool allQuarters) {
	const int n = batch->nSuperclusters;
	const int nCells = batch->nGridnodes;
	const Int3 boxSize = batch->boxSize;
	const Float3 boxSizeF = batch->boxLength;
	auto& tb = *batch->taskbuilderControl;

	// Particles within the list radius are at most ceil(listRadius / cell length) cells away, plus 1 for particles protruding from
	// their cells. Cells are 1 nm, or slightly longer in a non-integer box. FindNeighborsKernel skips the cells out of reach. In
	// small boxes, every cell is visited once
	const Float3 cellLength = boxSizeF / NodeIndex(boxSize).toFloat3();
	const auto Range = [&](float length) { return static_cast<int>(std::ceil(listRadius / length)) + 1; };
	const Int3 range{ Range(cellLength.x), Range(cellLength.y), Range(cellLength.z) };
	const Int3 rangeLo{ std::min(range.x, (boxSize.x - 1) / 2), std::min(range.y, (boxSize.y - 1) / 2), std::min(range.z, (boxSize.z - 1) / 2) };
	const Int3 rangeHi{ std::min(range.x, boxSize.x / 2), std::min(range.y, boxSize.y / 2), std::min(range.z, boxSize.z / 2) };

	tb.scSpheres.Expand(n, 1.2);
	tb.scCells.Expand(n, 1.2);
	tb.scValidMasks.Expand(n, 1.2);
	tb.scExternalBonds.Expand(n, 1.2);
	tb.cellMin.Expand(nCells);
	tb.cellMax.Expand(nCells);
	tb.neighbors.Expand(size_t(n) * maxNeighborsPerSc, 1.2);
	tb.nNeighbors.Expand(n, 1.2);
	tb.entryCounts.Expand(n + 1, 1.2);
	tb.overflow.Expand(1);

	SuperclusterBoundsKernel<<<nCells, 32, 0, stream>>>(*batch->superClustersControl, tb.contents.pclustersBondedToPcluster, batch->pClusterMetaDevice.Get(), tb.scSpheres.Get(),
		tb.scCells.Get(), tb.scValidMasks.Get(), tb.scExternalBonds.Get(), tb.cellMin.Get(), tb.cellMax.Get());
	cudaMemsetAsync(tb.overflow.Get(), 0, sizeof(int), stream);
	FindNeighborsKernel<<<(n + neighborSearchWarpsPerBlock - 1) / neighborSearchWarpsPerBlock, 32 * neighborSearchWarpsPerBlock, 0, stream>>>(
		*batch->superClustersControl, tb.scSpheres.Get(), tb.scExternalBonds.Get(), tb.scCells.Get(), tb.cellMin.Get(), tb.cellMax.Get(),
		tb.contents.pclustersBondedToPcluster,
		tb.neighbors.Get(), tb.nNeighbors.Get(), tb.entryCounts.Get(), tb.overflow.Get(),
		n, boxSize, rangeLo, rangeHi, listRadius, allQuarters, boxSizeF, boxSizeF.Inv());
	LIMA_UTILS::genericErrorCheckNoSync("FindNeighborsKernel");
}

namespace {
	// overflow as read back after the task build: a count that exceeded maxNeighborsPerSc, or 0
	void ThrowIfNeighborOverflow(int overflow) {
		if (overflow > 0)
			throw std::runtime_error("A supercluster has " + std::to_string(overflow) + " neighbors, more than the capacity of " + std::to_string(maxNeighborsPerSc));
	}
}

void Engine::MakeNbTasksMD(cudaStream_t stream) {
	const int n = batch->nSuperclusters;
	auto& tb = *batch->taskbuilderControl;
	FindSuperclusterNeighbors(stream, batch->params.cutoff_nm + batch->listBuffer, false);

	tb.entryStarts.Expand(n + 1, 1.2);
	batch->quarterEntryTasksDevice.Expand(n, 1.2);

	cudaMemsetAsync(tb.entryCounts.Get() + n, 0, sizeof(int), stream);
	CubWrappers::ExclusiveScan(tb.entryCounts.Get(), tb.entryCounts.Get() + n + 1, tb.entryStarts.Get(), batch->scanTempStorage, stream);
	auto& readback = batch->rebuildReadback;
	readback.Reserve(2);
	CubWrappers::GatherInts(tb.overflow.Get(), nullptr, 1, readback.Device(), stream);
	CubWrappers::GatherInts(tb.entryStarts.Get() + n, nullptr, 1, readback.Device() + 1, stream);
	readback.Written(stream);

	// Queued before waiting for the counts, see RunClustering.
	// Integration clears the atomic planes. Primary bonded planes are overwritten by their owning groups.
	// Clear both on rebuild, since supercluster slots change and particles without bonds have no writer.
	batch->forceAccumulatorDevice.Expand(size_t(n) * SuperCluster::maxParticles * 8, 1.2);
	cudaMemsetAsync(batch->forceAccumulatorDevice.Get(), 0, sizeof(unsigned long long) * n * SuperCluster::maxParticles * 8, stream);
	readback.Wait(stream);
	ThrowIfNeighborOverflow(readback.Host()[0]);
	batch->nQuarterEntries = readback.Host()[1];
	batch->quarterEntriesDevice.Expand(std::max(batch->nQuarterEntries, 1), 1.2);

	EmitQuarterEntriesKernel<false><<<n, 32, 0, stream>>>(*batch->superClustersControl, tb.neighbors.Get(), tb.nNeighbors.Get(), tb.scValidMasks.Get(),
		tb.contents.particlesBondedToParticle, tb.entryStarts.Get(), batch->quarterEntryTasksDevice.Get(), batch->quarterEntriesDevice.Get(),
		tb.contents, nullptr);
	LIMA_UTILS::genericErrorCheckNoSync("EmitQuarterEntriesKernel");
	cudaStreamQuery(stream); // Submits the work queued after the sync at once, see RunClustering
}

// EM uses the same entries as MD, but stores its forces in SCResults rather than summing them with fixed point atomics, see
// NbNonlocalKernel. Each supercluster has a result for its own forces, followed by one per supercluster owning an interaction
// with it, ordered by owner id.
void Engine::MakeNbTasksEM(cudaStream_t stream) {
	const int n = batch->nSuperclusters;
	auto& tb = *batch->taskbuilderControl;
	FindSuperclusterNeighbors(stream, batch->params.cutoff_nm + emListBuffer, true);

	int neighborOverflow = 0;
	cudaMemcpyAsync(&neighborOverflow, tb.overflow.Get(), sizeof(int), cudaMemcpyDeviceToHost, stream);
	tb.AllocateEmBuffers();
	tb.Reset(stream); // Synchronizes the stream
	ThrowIfNeighborOverflow(neighborOverflow);
	DistributeEmInteractionsKernel<<<n, 32, 0, stream>>>(tb.neighbors.Get(), tb.nNeighbors.Get(), tb.contents, tb.overflow.Get());
	LIMA_UTILS::genericErrorCheckNoSync("DistributeEmInteractionsKernel");
	SortEmInteractionsKernel<<<n, TaskBuilderControlContents::maxTasksPerSc, 0, stream>>>(tb.contents);
	LIMA_UTILS::genericErrorCheckNoSync("SortEmInteractionsKernel");

	tb.entryStarts.Expand(n + 1, 1.2);
	batch->quarterEntryTasksDevice.Expand(n, 1.2);
	cudaMemsetAsync(tb.contents.nResults + n, 0, sizeof(int), stream);
	cudaMemsetAsync(tb.entryCounts.Get() + n, 0, sizeof(int), stream);
	CubWrappers::ExclusiveScan(tb.contents.nResults, tb.contents.nResults + n + 1, tb.contents.nResultsPrefixsum, batch->scanTempStorage, stream);
	CubWrappers::ExclusiveScan(tb.entryCounts.Get(), tb.entryCounts.Get() + n + 1, tb.entryStarts.Get(), batch->scanTempStorage, stream);

	auto& readback = batch->rebuildReadback;
	readback.Reserve(3);
	CubWrappers::GatherInts(tb.overflow.Get(), nullptr, 1, readback.Device(), stream);
	CubWrappers::GatherInts(tb.contents.nResultsPrefixsum + n, nullptr, 1, readback.Device() + 1, stream);
	CubWrappers::GatherInts(tb.entryStarts.Get() + n, nullptr, 1, readback.Device() + 2, stream);
	cudaStreamSynchronize(stream);
	const int overflow = readback.Host()[0];
	batch->nResults = readback.Host()[1];
	batch->nQuarterEntries = readback.Host()[2];
	if (overflow > 0)
		throw std::runtime_error("A supercluster is the neighbor of " + std::to_string(overflow) + " superclusters, more than the capacity of "
			+ std::to_string(TaskBuilderControlContents::maxTasksPerSc));
	batch->quarterEntriesDevice.Expand(std::max(batch->nQuarterEntries, 1), 1.2);
	batch->quarterEntryResultIndicesDevice.Expand(std::max(batch->nQuarterEntries, 1), 1.2);
	// scResultsDevice is expanded lazily in _deviceMaster

	EmitQuarterEntriesKernel<true><<<n, 32, 0, stream>>>(*batch->superClustersControl, tb.neighbors.Get(), tb.nNeighbors.Get(), tb.scValidMasks.Get(),
		tb.contents.particlesBondedToParticle, tb.entryStarts.Get(), batch->quarterEntryTasksDevice.Get(), batch->quarterEntriesDevice.Get(),
		tb.contents, batch->quarterEntryResultIndicesDevice.Get());
	LIMA_UTILS::genericErrorCheckNoSync("EmitQuarterEntriesKernel");
}

bool Engine::MakeSuperClusterTasksGPU(cudaStream_t stream) {
	InvalidateStepGraph();
	if (batch->nSuperclusters == 0)
		return true;

	if (!batch->taskbuilderControl || batch->taskbuilderControl->nSuperclustersUpperbound < batch->nSuperclusters)
		batch->taskbuilderControl = std::make_unique<TaskBuilderControl>(
			batch->nSuperclusters * 2, batch->particlesBondedToParticle, batch->pclustersBondedToPcluster);

	tasksBuiltForEm = batch->params.em_variant;
	if (batch->params.em_variant)
		MakeNbTasksEM(stream);
	else
		MakeNbTasksMD(stream);
	return true;
}
