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
#include "LimaPositionSystem.cuh"
#include "Utilities.h"

#include <cfloat>

namespace {
	constexpr int maxNeighborsPerSc = TaskBuilderControlContents::maxTasksPerSc;
	constexpr int neighborSearchWarpsPerBlock = 4;

	// The lists also include 4x4 blocks that are within this distance beyond the cutoff, so pairs moving within the cutoff
	// before the next task build are still computed. 0.05 nm misses as few pairs as the supercluster granular lists of
	// earlier versions did [nm]
	constexpr float listBuffer = 0.05f;
}

__device__ inline bool Warp_ScAreBonded(const SuperClusterMeta& sc0, const SuperClusterMeta& sc1, const PclustersBondedToPcluster* const pclustersBondedToPcluster) {
	constexpr unsigned int mask = 0xFFFFFFFFu;

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



// ------------------------------------------------------- Neighbor search ------------------------------------------------------- //

// gridDim = nGridnodes, blockDim = 32. Bounding sphere, cell, valid particles and external bonds of every supercluster, and the AABB of every cell
__global__ void SuperclusterBoundsKernel(const SuperClustersControl scControl, const PclustersBondedToPcluster* const pclustersBondedToPcluster,
	float4* const scSpheres, int* const scCells, uint16_t* const scValidMasks, uint8_t* const scExternalBonds, float4* const cellMin, float4* const cellMax)
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
			sum.x += __shfl_xor_sync(0xFFFFFFFFu, sum.x, offset);
			sum.y += __shfl_xor_sync(0xFFFFFFFFu, sum.y, offset);
			sum.z += __shfl_xor_sync(0xFFFFFFFFu, sum.z, offset);
			count += __shfl_xor_sync(0xFFFFFFFFu, count, offset);
		}
		const Float3 center = sum * (1.f / fmaxf(count, 1.f));
		float radius = valid ? (p - center).len() : 0.f;
		for (int offset = 8; offset > 0; offset >>= 1)
			radius = fmaxf(radius, __shfl_xor_sync(0xFFFFFFFFu, radius, offset));
		const uint32_t validMask = __ballot_sync(0xFFFFFFFFu, valid) & 0xFFFF;
		if (lane == 0) {
			scSpheres[scId] = float4{ center.x, center.y, center.z, radius };
			scCells[scId] = cell;
			scValidMasks[scId] = static_cast<uint16_t>(validMask);
		}
		// Whether any value in the bonded sets of the supercluster's pclusters is a pcluster outside it. Superclusters without are bonded
		// to no other supercluster, so FindNeighborsKernel can skip checking their pairs
		const int nUnique = scControl.scMeta[scId].nUniquePcIds;
		const int ownPcluster = lane < nUnique ? scControl.scMeta[scId].uniquePclusterIds[lane] : -1;
		bool externalBond = false;
		for (int u = 0; u < nUnique; u++) {
			const int bondedPcluster = pclustersBondedToPcluster[__shfl_sync(0xFFFFFFFFu, ownPcluster, u)].Get(lane);
			bool internal = !PclustersBondedToPcluster::IsValue(bondedPcluster);
			for (int v = 0; v < nUnique; v++)
				internal |= __shfl_sync(0xFFFFFFFFu, ownPcluster, v) == bondedPcluster;
			externalBond |= !internal;
		}
		externalBond = __any_sync(0xFFFFFFFFu, externalBond);
		if (lane == 0)
			scExternalBonds[scId] = externalBond;

		if (valid) {
			lo = Float3{ fminf(lo.x, p.x), fminf(lo.y, p.y), fminf(lo.z, p.z) };
			hi = Float3{ fmaxf(hi.x, p.x), fmaxf(hi.y, p.y), fmaxf(hi.z, p.z) };
		}
	}
	for (int offset = 16; offset > 0; offset >>= 1) {
		lo.x = fminf(lo.x, __shfl_xor_sync(0xFFFFFFFFu, lo.x, offset));
		lo.y = fminf(lo.y, __shfl_xor_sync(0xFFFFFFFFu, lo.y, offset));
		lo.z = fminf(lo.z, __shfl_xor_sync(0xFFFFFFFFu, lo.z, offset));
		hi.x = fmaxf(hi.x, __shfl_xor_sync(0xFFFFFFFFu, hi.x, offset));
		hi.y = fmaxf(hi.y, __shfl_xor_sync(0xFFFFFFFFu, hi.y, offset));
		hi.z = fmaxf(hi.z, __shfl_xor_sync(0xFFFFFFFFu, hi.z, offset));
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
__global__ void FindNeighborsKernel(const SuperClustersControl scControl, const float4* const scSpheres, const uint8_t* const scExternalBonds, const int* const scCells,
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
	const bool ownExternalBonds = scExternalBonds[scId];
	const float listRadiusSq = listRadius * listRadius;

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

		uint32_t cellHits = __ballot_sync(0xFFFFFFFFu, cellHit);
		while (cellHits) {
			const int src = __ffs(cellHits) - 1;
			cellHits &= cellHits - 1;
			const int candidateCell = __shfl_sync(0xFFFFFFFFu, targetCell, src);
			const int nCandidates = scControl.nSuperclustersInBlocks[candidateCell];

			int candidateId = -1;
			bool coarseHit = false;
			if (lane < nCandidates) {
				candidateId = scControl.scIdsInBlocks[candidateCell * SuperClustersControl::maxClustersPerBlock + lane];
				if (candidateId >= scId) {
					const float4 s = scSpheres[candidateId];
					Float3 c{ s.x, s.y, s.z };
					PeriodicBoundaryCondition::ApplyHyperpos(ownCenter, c, boxSizeF, boxSizeInv);
					const float reach = listRadius + ownSphere.w + s.w;
					coarseHit = (c - ownCenter).lenSquared() < reach * reach;
				}
			}

			uint32_t candidates = __ballot_sync(0xFFFFFFFFu, coarseHit);
			while (candidates) {
				const int k = __ffs(candidates) - 1;
				candidates &= candidates - 1;
				const int queryId = __shfl_sync(0xFFFFFFFFu, candidateId, k);
				const bool selfTask = queryId == scId;

				if (lane < SuperCluster::maxParticles) {
					Float3 p = scControl.scData[queryId].Position(lane);
					PeriodicBoundaryCondition::ApplyHyperpos(ownCenter, p, boxSizeF, boxSizeInv);
					stage[warpInBlock][lane] = float4{ p.x, p.y, p.z, scControl.scData[queryId].Valid(lane) ? 1.f : 0.f };
				}
				__syncwarp();
				uint32_t bits = 0;
#pragma unroll
				for (int m = 0; m < 8; m++) {
					const int queryIndex = (lane >> 4) + 2 * m;
					const float4 q = stage[warpInBlock][queryIndex];
					const bool inRange = (ownPos - Float3{ q.x, q.y, q.z }).lenSquared() < listRadiusSq;
					if (ownValid && q.w != 0.f && inRange && (!selfTask || queryIndex > ownIndex))
						bits |= 1u << ((ownIndex >> 2) * 4 + (queryIndex >> 2));
				}
				bits = __reduce_or_sync(0xFFFFFFFFu, bits);
				__syncwarp(); // Before the stage is reused

				if (bits) {
					const bool bonded = selfTask || (ownExternalBonds && Warp_ScAreBonded(scControl.scMeta[scId], scControl.scMeta[queryId], pclustersBondedToPcluster));
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

// Order in which ownQuarterMask classes are emitted. Neighbouring classes differ in one bit (Gray code), so when a warp's
// 8 groups straddle two classes, the union of their masks is small. EM's entries without pairs in range (class 0) go last
__constant__ uint8_t quarterMaskEmitOrder[16] = { 1, 3, 2, 6, 7, 5, 4, 12, 13, 15, 14, 10, 11, 9, 8, 0 };

// gridDim = nSuperclusters, blockDim = 32. Writes the QuarterEntries of each supercluster: one per quarter of a neighbor with a pair
// within the list radius. Sorted by ownQuarterMask in quarterMaskEmitOrder, and within a class by (chunk of 32 neighbors, query
// quarter, neighbor), so the result is deterministic.
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

		for (int jq = 0; jq < 4; jq++) {
			const uint32_t maskClass = validNeighbor ? EntryClass(ownMasks, jq, queryScId) : noEntry;
			const uint32_t sameClass = __match_any_sync(0xFFFFFFFFu, maskClass);
			const int rank = __popc(sameClass & ((1u << lane) - 1));
			const int dst = maskClass != noEntry ? classOffsets[maskClass] + rank : 0;
			__syncwarp();
			if (maskClass != noEntry && rank == 0)
				classOffsets[maskClass] += __popc(sameClass);
			__syncwarp();
			if (maskClass == noEntry)
				continue;

			const uint32_t queryValid4 = (queryValid >> (jq * 4)) & 0xF;
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
						const int pidOwn = scControl.scMeta[scId].globalParticleIds[iq * 4 + iLocal];
						if (pidOwn == -1) continue;
						for (int jLocal = 0; jLocal < 4; jLocal++) {
							const int pidQuery = scControl.scMeta[queryScId].globalParticleIds[jq * 4 + jLocal];
							if (pidQuery != -1 && particlesBondedToParticle[pidOwn].Contains(pidQuery))
								noInteractions |= 1u << (iLocal * 4 + jLocal);
						}
					}
				}
				out.noInteractions[iq] = static_cast<uint16_t>(noInteractions);
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
	const Float3 boxSizeF = NodeIndex(boxSize).toFloat3();
	auto& tb = *batch->taskbuilderControl;

	// Cells are 1 nm, so particles within the list radius are at most ceil(listRadius) cells away, plus 1 for particles
	// protruding from their cells. FindNeighborsKernel skips the cells out of reach. In small boxes, every cell is visited once
	const int range = static_cast<int>(std::ceil(listRadius)) + 1;
	const Int3 rangeLo{ std::min(range, (boxSize.x - 1) / 2), std::min(range, (boxSize.y - 1) / 2), std::min(range, (boxSize.z - 1) / 2) };
	const Int3 rangeHi{ std::min(range, boxSize.x / 2), std::min(range, boxSize.y / 2), std::min(range, boxSize.z / 2) };

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

	SuperclusterBoundsKernel<<<nCells, 32, 0, stream>>>(*batch->superClustersControl, tb.contents.pclustersBondedToPcluster, tb.scSpheres.Get(),
		tb.scCells.Get(), tb.scValidMasks.Get(), tb.scExternalBonds.Get(), tb.cellMin.Get(), tb.cellMax.Get());
	cudaMemsetAsync(tb.overflow.Get(), 0, sizeof(int), stream);
	FindNeighborsKernel<<<(n + neighborSearchWarpsPerBlock - 1) / neighborSearchWarpsPerBlock, 32 * neighborSearchWarpsPerBlock, 0, stream>>>(
		*batch->superClustersControl, tb.scSpheres.Get(), tb.scExternalBonds.Get(), tb.scCells.Get(), tb.cellMin.Get(), tb.cellMax.Get(),
		tb.contents.pclustersBondedToPcluster,
		tb.neighbors.Get(), tb.nNeighbors.Get(), tb.entryCounts.Get(), tb.overflow.Get(),
		n, boxSize, rangeLo, rangeHi, listRadius, allQuarters, boxSizeF, boxSizeF.Inv());
	LIMA_UTILS::genericErrorCheck(stream, "FindNeighborsKernel");
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
	FindSuperclusterNeighbors(stream, batch->params.cutoff_nm + listBuffer, false);

	tb.entryStarts.Expand(n + 1, 1.2);
	batch->quarterEntryTasksDevice.Expand(n, 1.2);
	// Integration clears the atomic planes. Primary bonded planes are overwritten by their owning groups.
	// Clear both on rebuild, since supercluster slots change and particles without bonds have no writer.
	batch->forceAccumulatorDevice.Expand(size_t(n) * SuperCluster::maxParticles * 8, 1.2);
	cudaMemsetAsync(batch->forceAccumulatorDevice.Get(), 0, sizeof(unsigned long long) * n * SuperCluster::maxParticles * 8, stream);

	cudaMemsetAsync(tb.entryCounts.Get() + n, 0, sizeof(int), stream);
	CubWrappers::ExclusiveScan(tb.entryCounts.Get(), tb.entryCounts.Get() + n + 1, tb.entryStarts.Get(), stream);
	int overflow = 0;
	cudaMemcpyAsync(&overflow, tb.overflow.Get(), sizeof(int), cudaMemcpyDeviceToHost, stream);
	cudaMemcpyAsync(&batch->nQuarterEntries, tb.entryStarts.Get() + n, sizeof(int), cudaMemcpyDeviceToHost, stream);
	cudaStreamSynchronize(stream);
	ThrowIfNeighborOverflow(overflow);
	batch->quarterEntriesDevice.Expand(std::max(batch->nQuarterEntries, 1), 1.2);

	EmitQuarterEntriesKernel<false><<<n, 32, 0, stream>>>(*batch->superClustersControl, tb.neighbors.Get(), tb.nNeighbors.Get(), tb.scValidMasks.Get(),
		tb.contents.particlesBondedToParticle, tb.entryStarts.Get(), batch->quarterEntryTasksDevice.Get(), batch->quarterEntriesDevice.Get(),
		tb.contents, nullptr);
	LIMA_UTILS::genericErrorCheck(stream, "EmitQuarterEntriesKernel");
}

// EM uses the same entries as MD, but stores its forces in SCResults rather than summing them with fixed point atomics, see
// NbNonlocalKernel. Each supercluster has a result for its own forces, followed by one per supercluster owning an interaction
// with it, ordered by owner id.
void Engine::MakeNbTasksEM(cudaStream_t stream) {
	const int n = batch->nSuperclusters;
	auto& tb = *batch->taskbuilderControl;
	FindSuperclusterNeighbors(stream, batch->params.cutoff_nm + listBuffer, true);

	int neighborOverflow = 0;
	cudaMemcpyAsync(&neighborOverflow, tb.overflow.Get(), sizeof(int), cudaMemcpyDeviceToHost, stream);
	tb.AllocateEmBuffers();
	tb.Reset(stream); // Synchronizes the stream
	ThrowIfNeighborOverflow(neighborOverflow);
	DistributeEmInteractionsKernel<<<n, 32, 0, stream>>>(tb.neighbors.Get(), tb.nNeighbors.Get(), tb.contents, tb.overflow.Get());
	LIMA_UTILS::genericErrorCheck(stream, "DistributeEmInteractionsKernel");
	SortEmInteractionsKernel<<<n, TaskBuilderControlContents::maxTasksPerSc, 0, stream>>>(tb.contents);
	LIMA_UTILS::genericErrorCheck(stream, "SortEmInteractionsKernel");

	tb.entryStarts.Expand(n + 1, 1.2);
	batch->quarterEntryTasksDevice.Expand(n, 1.2);
	cudaMemsetAsync(tb.contents.nResults + n, 0, sizeof(int), stream);
	cudaMemsetAsync(tb.entryCounts.Get() + n, 0, sizeof(int), stream);
	CubWrappers::ExclusiveScan(tb.contents.nResults, tb.contents.nResults + n + 1, tb.contents.nResultsPrefixsum, stream);
	CubWrappers::ExclusiveScan(tb.entryCounts.Get(), tb.entryCounts.Get() + n + 1, tb.entryStarts.Get(), stream);

	int overflow = 0;
	cudaMemcpyAsync(&overflow, tb.overflow.Get(), sizeof(int), cudaMemcpyDeviceToHost, stream);
	cudaMemcpyAsync(&batch->nResults, tb.contents.nResultsPrefixsum + n, sizeof(int), cudaMemcpyDeviceToHost, stream);
	cudaMemcpyAsync(&batch->nQuarterEntries, tb.entryStarts.Get() + n, sizeof(int), cudaMemcpyDeviceToHost, stream);
	cudaStreamSynchronize(stream);
	if (overflow > 0)
		throw std::runtime_error("A supercluster is the neighbor of " + std::to_string(overflow) + " superclusters, more than the capacity of "
			+ std::to_string(TaskBuilderControlContents::maxTasksPerSc));
	batch->quarterEntriesDevice.Expand(std::max(batch->nQuarterEntries, 1), 1.2);
	batch->quarterEntryResultIndicesDevice.Expand(std::max(batch->nQuarterEntries, 1), 1.2);
	// scResultsDevice is expanded lazily in _deviceMaster

	EmitQuarterEntriesKernel<true><<<n, 32, 0, stream>>>(*batch->superClustersControl, tb.neighbors.Get(), tb.nNeighbors.Get(), tb.scValidMasks.Get(),
		tb.contents.particlesBondedToParticle, tb.entryStarts.Get(), batch->quarterEntryTasksDevice.Get(), batch->quarterEntriesDevice.Get(),
		tb.contents, batch->quarterEntryResultIndicesDevice.Get());
	LIMA_UTILS::genericErrorCheck(stream, "EmitQuarterEntriesKernel");
}

bool Engine::MakeSuperClusterTasksGPU(cudaStream_t stream) {
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
