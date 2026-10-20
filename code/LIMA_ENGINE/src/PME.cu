// PME::Controller and its kernels. A translation unit of its own, compiled in parallel with Engine.cu

#include "PME.cuh"
#include "ChargeBlock.cuh"
#include "DeviceAlgorithmsPrivate.cuh"
#include "EngineBodies.cuh"
#include "PhysicsUtils.cuh"


#include "DeviceAlgorithms.cuh"
#include "BoundaryCondition.cuh"
#include "Filehandling.h"
#include "platform/DeviceIntrinsics.cuh"
#include "Utilities.h"
#include "LimitTesting.cuh"
#include "Virial.cuh"

using namespace ChargeBlock;
using namespace PME;

// --------------------------------------------------------------- Kernel Helpers --------------------------------------------------------------- //	

constexpr int GetGridIndexRealspace(const Int3& gridIndex, const Int3& gridDim) {
	return BoxGrid::Get1dIndex(gridIndex, gridDim);
}
constexpr int GetGridIndexReciprocalspace(const Int3& gridIndex, int gridpointsPerDim, int nGridpointsHalfdim) {
	return gridIndex.x + gridIndex.y * nGridpointsHalfdim + gridIndex.z * nGridpointsHalfdim * gridpointsPerDim;
}
constexpr NodeIndex Get3dIndexReciprocalspace(int index1d, const Int3& gridpointsPerDim, int nGridpointsHalfdim) {
	int z = index1d / (nGridpointsHalfdim * gridpointsPerDim.y);
	index1d -= z * nGridpointsHalfdim * gridpointsPerDim.y;
	int y = index1d / nGridpointsHalfdim;
	index1d -= y * nGridpointsHalfdim;
	int x = index1d;
	return NodeIndex{ x, y, z };
}

constexpr NodeIndex GetParticlesBlockRelativeToCompoundOrigo(const Float3& relpos, bool particleActive) {
	return NodeIndex{
		(relpos.x > 1.f) - (relpos.x < 0.f),
		(relpos.y > 1.f) - (relpos.y < 0.f),
		(relpos.z > 1.f) - (relpos.z < 0.f)
	};
}

struct Direction3 {
	uint8_t data; // store 6 bits: top 2 for x+1, next 2 for y+1, bottom 2 for z+1

	constexpr Direction3() : data(0) {}

	constexpr Direction3(int x, int y, int z)
		: data(static_cast<uint8_t>(((x + 1) << 4) | ((y + 1) << 2) | (z + 1)))
	{}

	constexpr int x() const { return ((data >> 4) & 0x3) - 1; }
	constexpr int y() const { return ((data >> 2) & 0x3) - 1; }
	constexpr int z() const { return (data & 0x3) - 1; }

	constexpr bool operator==(const Direction3& o) const { return data == o.data; }
	constexpr bool operator!=(const Direction3& o) const { return data != o.data; }
	constexpr Direction3 operator-() const { return Direction3{ -x(), -y(), -z() }; }
	constexpr NodeIndex ToNodeIndex() const { return NodeIndex{ x(), y(), z() }; }
	constexpr Float3 ToFloat3() const { return Float3{ static_cast<float>(x()), static_cast<float>(y()), static_cast<float>(z()) }; }
};


namespace device_tables {
	__device__ constexpr Direction3 sIndexToDirection[27] = {
	Direction3(-1,-1,-1), Direction3(0,-1,-1), Direction3(1,-1,-1),
	Direction3(-1, 0,-1), Direction3(0, 0,-1), Direction3(1, 0,-1),
	Direction3(-1, 1,-1), Direction3(0, 1,-1), Direction3(1, 1,-1),

	Direction3(-1,-1, 0), Direction3(0,-1, 0), Direction3(1,-1, 0),
	Direction3(-1, 0, 0), /* (0,0,0) excl */   Direction3(1, 0, 0),
	Direction3(-1, 1, 0), Direction3(0, 1, 0), Direction3(1, 1, 0),

	Direction3(-1,-1, 1), Direction3(0,-1, 1), Direction3(1,-1, 1),
	Direction3(-1, 0, 1), Direction3(0, 0, 1), Direction3(1, 0, 1),
	Direction3(-1, 1, 1), Direction3(0, 1, 1), Direction3(1, 1, 1),
	Direction3(0,0,0)												// 0-Dir is put here, so we can skip it in kernels that doesnt need it
	};

	__device__ constexpr uint8_t sDirectionToIndex[64] = {
		  0,   9,  17, 255,   3,  12,  20, 255,
		  6,  14,  23, 255, 255, 255, 255, 255,
		  1,  10,  18, 255,   4,   26,  21, 255,
		  7,  15,  24, 255, 255, 255, 255, 255,
		  2,  11,  19, 255,   5,  13,  22, 255,
		  8,  16,  25, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255
	};
}

// Bit d+1 is set if a particle at this grid floor index (relative to its chargeblock) spreads to the neighbor block in direction
// d along the axis. The spline reaches 2 gridpoints below and 1 above the floor index
constexpr uint32_t DirectionsAlongAxis(int floorIndex) {
	return uint32_t(floorIndex <= 0)
		| uint32_t(floorIndex >= -2 && floorIndex <= gridpointsPerNm) << 1
		| uint32_t(floorIndex >= gridpointsPerNm - 2) << 2;
}

constexpr Int3 FloorIndex3d(const Float3& relpos) {
	return Int3{
		static_cast<int>(floorf(relpos.x * gridpointsPerNm_f)),
		static_cast<int>(floorf(relpos.y * gridpointsPerNm_f)),
		static_cast<int>(floorf(relpos.z * gridpointsPerNm_f))
	};
}

// --------------------------------------------------------------- PME Kernels --------------------------------------------------------------- //	

// Copies each charged particle to the chargeblocks whose grid it spreads to: its own, and neighbors when it is within the
// spline's reach of the border. One warp per supercluster, a lane per particle. blockDim = (32, distributeScsPerBlock, 1)
constexpr int distributeScsPerBlock = 4;
// The grid works in coordinates scaled by gridScale per axis, to a box of blocksPerDim nm, so it tiles a box of any length with
// 1 nm chargeblocks. These are the fractional coordinates of smooth PME, the reciprocal space uses the real lengths
__global__ void DistributeCompoundchargesToBlocksKernel(const SuperCluster* const superclusters, const ChargeblockBuffers chargeblockBuffers, Int3 blocksPerDim,
	const SuperClusterMeta* metadata, const int* simulationSlots, int nSuperclusters, Float3 gridScale)
{
	const Float3 boxSize = NodeIndex(blocksPerDim).toFloat3();	// Scaled
	const Float3 boxSizeInv = boxSize.Inv();
	const int scId = blockIdx.x * distributeScsPerBlock + threadIdx.y;
	if (scId >= nSuperclusters)
		return;
	const int lane = threadIdx.x;
	const int blockOffset = simulationSlots[metadata[scId].simulationId] * blocksPerDim.InnerProduct();

	bool charged = false;
	Float3 pos{};
	float charge = 0.f;
	if (lane < SuperCluster::maxParticles) {
		charge = superclusters[scId].Charge(lane);
		charged = superclusters[scId].EpsilonSqrt(lane) != -1 && charge != 0.f;
		pos = superclusters[scId].Position(lane) * gridScale;
	}
	// Particles spread relative to a chargeblock, to it and the neighbor blocks within the spline's reach: 2 gridpoints below and 1
	// above a particle's floor index, so floor indices from -gridpointsPerNm + 2 to 2 * gridpointsPerNm - 3 relative to the block.
	// Usually all of a supercluster's particles are within that of one block, per axis the block of its lowest particle or the one
	// above, so they share it and reserve their space together. A wider one, eg. in the thinner last cell of an axis once scaled,
	// would lose charge, so its particles are grouped by their own blocks instead, within which they always reach
	constexpr int minFloorIndex = -gridpointsPerNm + 2, maxFloorIndex = 2 * gridpointsPerNm - 3;
	const Int3 floorIndexAbs = FloorIndex3d(pos);
	Int3 lo = charged ? floorIndexAbs : Int3{ INT_MAX, INT_MAX, INT_MAX };
	Int3 hi = charged ? floorIndexAbs : Int3{ INT_MIN, INT_MIN, INT_MIN };
	lo = Int3{ __reduce_min_sync(LIMA_FULL_WARP_MASK, lo.x), __reduce_min_sync(LIMA_FULL_WARP_MASK, lo.y), __reduce_min_sync(LIMA_FULL_WARP_MASK, lo.z) };
	hi = Int3{ __reduce_max_sync(LIMA_FULL_WARP_MASK, hi.x), __reduce_max_sync(LIMA_FULL_WARP_MASK, hi.y), __reduce_max_sync(LIMA_FULL_WARP_MASK, hi.z) };
	if (lo.x == INT_MAX)	// No charges
		return;
	// The block along an axis whose reach covers the floor indices lo to hi, or INT_MIN if none does
	const auto SharedBlock = [](int lo, int hi) {
		const int below = lo >= 0 ? lo / gridpointsPerNm : -((-lo + gridpointsPerNm - 1) / gridpointsPerNm);	// floor(lo / gridpointsPerNm)
		for (const int block : { below, below + 1 })
			if (lo - block * gridpointsPerNm >= minFloorIndex && hi - block * gridpointsPerNm <= maxFloorIndex)
				return block;
		return INT_MIN;
	};
	const NodeIndex sharedBlock{ SharedBlock(lo.x, hi.x), SharedBlock(lo.y, hi.y), SharedBlock(lo.z, hi.z) };
	const bool shared = sharedBlock.x != INT_MIN && sharedBlock.y != INT_MIN && sharedBlock.z != INT_MIN;
	const NodeIndex ownBlock = shared ? sharedBlock : NodeIndex(pos.Floor().ToInt3());
	const Float3 relPos = pos - ownBlock.toFloat3();

	// List the particle once in the block owning its grid floor index, which interpolates its force
	PeriodicBoundaryCondition::ApplyBC(pos, boxSize, boxSizeInv);
	// Wrapped first, as a position just below 0 can come out of ApplyBC, and integer division would round its floor index -1 up to block 0
	const NodeIndex ownerBlock = PeriodicBoundaryCondition::applyBC(FloorIndex3d(pos), blocksPerDim * gridpointsPerNm) / gridpointsPerNm;
	const int ownerBlockIndex = charged ? blockOffset + BoxGrid::Get1dIndex(ownerBlock, blocksPerDim) : -1;
	const uint32_t sameOwner = __match_any_sync(LIMA_FULL_WARP_MASK, ownerBlockIndex);
	if (charged) {
		const int firstLane = __ffs(sameOwner) - 1;
		int offsetInOwner = 0;
		if (lane == firstLane)
			offsetInOwner = atomicAdd(&chargeblockBuffers.nOwnedInBlock[ownerBlockIndex], __popc(sameOwner));
		const int indexInOwner = __shfl_sync(WarpMask(sameOwner), offsetInOwner, firstLane) + __popc(sameOwner & ((1u << lane) - 1));
		if (indexInOwner >= ChargeBlock::maxParticlesInBlock) {
			chargeblockBuffers.overflow.Report(CapacityOverflow::ChargeBlock, ownerBlockIndex, indexInOwner + 1, ChargeBlock::maxParticlesInBlock);
		}
		else {
			const int index = ownerBlockIndex * ChargeBlock::maxParticlesInBlock + indexInOwner;
			chargeblockBuffers.owned.particles[index].Store(pos, charge);
			chargeblockBuffers.owned.slots[index] = scId * SuperCluster::maxParticles + lane;
		}
	}
	const Int3 floorIndex3d = FloorIndex3d(relPos);
	const uint32_t xDirections = DirectionsAlongAxis(floorIndex3d.x);
	const uint32_t yDirections = DirectionsAlongAxis(floorIndex3d.y);
	const uint32_t zDirections = DirectionsAlongAxis(floorIndex3d.z);

	uint32_t pending = __ballot_sync(LIMA_FULL_WARP_MASK, charged);
	while (pending != 0) {
		const int leader = __ffs(pending) - 1;
		const NodeIndex groupBlock{ __shfl_sync(LIMA_FULL_WARP_MASK, ownBlock.x, leader), __shfl_sync(LIMA_FULL_WARP_MASK, ownBlock.y, leader),
			__shfl_sync(LIMA_FULL_WARP_MASK, ownBlock.z, leader) };
		const bool inGroup = charged && ownBlock.x == groupBlock.x && ownBlock.y == groupBlock.y && ownBlock.z == groupBlock.z;
		pending &= ~__ballot_sync(LIMA_FULL_WARP_MASK, inGroup);

		// Lane d collects which particles of the group go in direction d, and reserves space for them in that block. The reservations
		// are independent, so all 27 atomics are in flight at once
		uint32_t goingMask = 0;
		for (int directionIndex = 0; directionIndex < 27; directionIndex++) {
			const Direction3 direction = device_tables::sIndexToDirection[directionIndex];
			const bool goes = inGroup && ((xDirections >> (direction.x() + 1)) & (yDirections >> (direction.y() + 1)) & (zDirections >> (direction.z() + 1)) & 1);
			const uint32_t mask = __ballot_sync(LIMA_FULL_WARP_MASK, goes);
			if (lane == directionIndex)
				goingMask = mask;
		}
		int targetBlockIndex = 0;
		int offsetInTarget = 0;
		if (goingMask != 0) {
			const Direction3 direction = device_tables::sIndexToDirection[lane];
			targetBlockIndex = blockOffset + BoxGrid::Get1dIndex(PeriodicBoundaryCondition::applyBC(groupBlock + direction.ToNodeIndex(), blocksPerDim), blocksPerDim);
			offsetInTarget = atomicAdd(&chargeblockBuffers.nParticlesInBlock[targetBlockIndex], __popc(goingMask));
		}

		for (int directionIndex = 0; directionIndex < 27; directionIndex++) {
			const uint32_t mask = __shfl_sync(LIMA_FULL_WARP_MASK, goingMask, directionIndex);
			const int target = __shfl_sync(LIMA_FULL_WARP_MASK, targetBlockIndex, directionIndex);
			const int offset = __shfl_sync(LIMA_FULL_WARP_MASK, offsetInTarget, directionIndex);
			if (!((mask >> lane) & 1))
				continue;
			const int indexInTarget = offset + __popc(mask & ((1u << lane) - 1));
			if (indexInTarget >= ChargeBlock::maxParticlesInBlock)
				chargeblockBuffers.overflow.Report(CapacityOverflow::ChargeBlock, target, indexInTarget + 1, ChargeBlock::maxParticlesInBlock);
			else
				ChargeBlock::GetParticles(chargeblockBuffers, target)[indexInTarget].Store(relPos - device_tables::sIndexToDirection[directionIndex].ToFloat3(), charge);
		}
	}
}



// This kernel accumulates charges of it's particle in a local shared buffer. Any charge outside any kernels volume is ignored, and assumed handled elsewhere
// To remain deterministic, the intermediate charges are stored as integers, such that we simply can use atomicAdd to accumulate charges
// gridCellVolumeScale: the volume of a scaled grid cell per real one, see DistributeCompoundchargesToBlocksKernel
__global__ void ChargeblockDistributeToGrid(ChargeblockBuffers chargeblockBuffers, float* const realspaceGrid, Int3 blocksPerDim, Int3 gridpointsPerDim,
	float gridCellVolumeScale) {
	__shared__  int nParticles;
	__shared__ int localGrid[gridpointsPerNm * gridpointsPerNm * gridpointsPerNm];

	// Init shared variables
	float* localGridAsFloat = reinterpret_cast<float*>(localGrid);
	for (int i = threadIdx.x; i < gridpointsPerNm * gridpointsPerNm * gridpointsPerNm; i += blockDim.x)
		localGrid[i] = 0;

	if (threadIdx.x == 0) {
		nParticles = min(chargeblockBuffers.nParticlesInBlock[blockIdx.x], static_cast<uint32_t>(ChargeBlock::maxParticlesInBlock));	
		chargeblockBuffers.nParticlesInBlock[blockIdx.x] = 0; // Reset the particle count for the next round of accumulation	
	}
	__syncthreads();


	// By scaling the charge with the scalar, we can store the charge in the grid as an integer, allowing us to use atomicAdd deterministically	
	constexpr double largestPossibleValue = (4. / 6.) * (5. * elementaryChargeToKiloCoulombPerMole) * 8; // Highest bspline coeff * (highestCharge) * maxExpectedParticleNearNode
	constexpr float scalar = static_cast<double>(INT_MAX - 10) / largestPossibleValue / 10.f; // 10 for safety
	constexpr float invScalar = 1.f / scalar;


	// Loop over all particles, distribute it's charges to the __shared__ grid
	for (int i = threadIdx.x; i < nParticles; i += blockDim.x) {
		const ChargePos* const particlePtr = &ChargeBlock::GetParticles(chargeblockBuffers, blockIdx.x)[i];

		const float charge = particlePtr->charge;	// [kC/mol]
		const Float3 posRelativeToBlock = particlePtr->pos;

		// Map position to fractional local grid coordinates
		Float3 gridPos = posRelativeToBlock * gridpointsPerNm_f;
		int ix = static_cast<int>(floorf(gridPos.x));
		int iy = static_cast<int>(floorf(gridPos.y));
		int iz = static_cast<int>(floorf(gridPos.z));

		float fx = gridPos.x - static_cast<float>(ix);
		float fy = gridPos.y - static_cast<float>(iy);
		float fz = gridPos.z - static_cast<float>(iz);

		float wx[4], wy[4], wz[4];
		LAL::CalcBspline(fx, wx);
		LAL::CalcBspline(fy, wy);
		LAL::CalcBspline(fz, wz);

		// Distribute the charge to the surrounding 4x4x4 cells
		for (int dx = 0; dx < 4; dx++) {
			int X = -1 + ix + dx;
			if (X < 0 || X >= gridpointsPerNm) // Out-of-local-bounds gridpoints are handled by the neighbor blocks
				continue;

			float wxCur = wx[dx];
			for (int dy = 0; dy < 4; dy++) {
				int Y = -1 + iy + dy;
				if (Y < 0 || Y >= gridpointsPerNm)
					continue;

				float wxyCur = wxCur * wy[dy];
				for (int dz = 0; dz < 4; dz++) {
					int Z = -1 + iz + dz;
					if (Z < 0 || Z >= gridpointsPerNm)
						continue;

					const NodeIndex index3d = NodeIndex{ X,Y,Z };
					const int index1D = GetGridIndexRealspace(index3d, Int3(gridpointsPerNm, gridpointsPerNm, gridpointsPerNm));
					const float chargeScaled = charge * wxyCur * wz[dz] * scalar;
					const int clampedChargeDiscretized = static_cast<int>(fminf(fmaxf(chargeScaled, static_cast<float>(INT_MIN)), static_cast<float>(INT_MAX)));
					atomicAdd(&localGrid[index1D], clampedChargeDiscretized);
					
				}
			}
		}
	}
	__syncthreads();

	// Transform the grid back to float values, and multiply with invCellVolume
	for (int i = threadIdx.x; i < gridpointsPerNm * gridpointsPerNm * gridpointsPerNm; i += blockDim.x) {
		localGridAsFloat[i] = static_cast<float>(localGrid[i]) * invScalar * invCellVolume * gridCellVolumeScale;
	}
	__syncthreads();

	// Transform local to global grid coordinates, and push to global memory
	const int localCellCount = gridpointsPerNm * gridpointsPerNm * gridpointsPerNm;
	const NodeIndex blocksFirstIndex3dInRealspacegrid =
		BoxGrid::Get3dIndex(blockIdx.x % blocksPerDim.InnerProduct(), blocksPerDim) * gridpointsPerNm;

	for (int localIndex = threadIdx.x; localIndex < localCellCount; localIndex += blockDim.x) {
		const int x = localIndex % gridpointsPerNm;
		const int y = (localIndex / gridpointsPerNm) % gridpointsPerNm;
		const int z = localIndex / (gridpointsPerNm * gridpointsPerNm);

		const NodeIndex globalIndex3d = blocksFirstIndex3dInRealspacegrid + NodeIndex{ x, y, z };
		const int globalIndex = BoxGrid::Get1dIndex(globalIndex3d, gridpointsPerDim);

		const size_t gridOffset = size_t(blockIdx.x / blocksPerDim.InnerProduct()) * gridpointsPerDim.InnerProduct();
		realspaceGrid[gridOffset + globalIndex] = localGridAsFloat[localIndex];
	}
}

// B-spline interpolation of the potential at gridPos, and of the field from the B-splines' derivatives (standard SPME).
// phi(X, Y, Z) returns the potential at a gridpoint, which may be outside the grid
template <typename Phi>
__device__ ForceEnergy InterpolateForceEnergy(Float3 gridPos, Phi phi) {
	const int ix = static_cast<int>(floorf(gridPos.x));
	const int iy = static_cast<int>(floorf(gridPos.y));
	const int iz = static_cast<int>(floorf(gridPos.z));

	const float fx = gridPos.x - static_cast<float>(ix);
	const float fy = gridPos.y - static_cast<float>(iy);
	const float fz = gridPos.z - static_cast<float>(iz);

	float wx[4], wy[4], wz[4];
	LAL::CalcBspline(fx, wx);
	LAL::CalcBspline(fy, wy);
	LAL::CalcBspline(fz, wz);
	float dwx[4], dwy[4], dwz[4];
	LAL::CalcBsplineDerivative(fx, dwx);
	LAL::CalcBsplineDerivative(fy, dwy);
	LAL::CalcBsplineDerivative(fz, dwz);

	// The gradient is with respect to grid units, gridpointsPerNm converts it to per nm
	Float3 gradient{};
	float potential = 0.f;
	for (int dz = 0; dz < 4; dz++) {
		for (int dy = 0; dy < 4; dy++) {
#pragma unroll
			for (int dx = 0; dx < 4; dx++) {
				const float p = phi(ix - 1 + dx, iy - 1 + dy, iz - 1 + dz);
				gradient.x += dwx[dx] * wy[dy] * wz[dz] * p;
				gradient.y += wx[dx] * dwy[dy] * wz[dz] * p;
				gradient.z += wx[dx] * wy[dy] * dwz[dz] * p;
				potential += wx[dx] * wy[dy] * wz[dz] * p;
			}
		}
	}
	return ForceEnergy{ gradient * -gridpointsPerNm_f, potential };
}

// The grid around a chargeblock that its particles' interpolation reads: the block's gridpoints, plus the spline's reach of
// 1 below and 2 above
constexpr int interpolationTileLen = gridpointsPerNm + 3;

// Interpolates the forces on the particles a chargeblock owns, from a tile of the grid around the block in shared memory.
// One chargeblock per block, blockDim = interpolateThreads
constexpr int interpolateThreads = 128;
__global__ void __launch_bounds__(interpolateThreads) InterpolateForcesKernel(const ChargeblockBuffers chargeblockBuffers, const float* const realspaceGrid,
	Int3 blocksPerDim, Int3 gridDim, const SuperClusterMeta* const scMeta, float ewaldKappa /*[nm^-1]*/,
	ForceEnergy* const forceEnergies /*EM only*/, const ForceAccumulator forceAcc /*MD only, forceAcc.fx == nullptr in EM*/,
	Float3 gridScale /*see DistributeCompoundchargesToBlocksKernel*/)
{
	constexpr int tileLen = interpolationTileLen;
	__shared__ float tile[tileLen][tileLen][tileLen];
	__shared__ int nOwned;

	const int nBlocksPerGrid = blocksPerDim.InnerProduct();
	const float* const grid = realspaceGrid + size_t(blockIdx.x / nBlocksPerGrid) * gridDim.InnerProduct();
	const NodeIndex tileOrigin = BoxGrid::Get3dIndex(blockIdx.x % nBlocksPerGrid, blocksPerDim) * gridpointsPerNm - NodeIndex{ 1, 1, 1 };

	if (threadIdx.x == 0) {
		nOwned = min(chargeblockBuffers.nOwnedInBlock[blockIdx.x], static_cast<uint32_t>(ChargeBlock::maxParticlesInBlock));
		chargeblockBuffers.nOwnedInBlock[blockIdx.x] = 0; // Reset for the next step
	}
	// All of a thread's loads are issued before storing any, so they are in flight together
	constexpr int tileSize = tileLen * tileLen * tileLen;
	constexpr int loadsPerThread = (tileSize + interpolateThreads - 1) / interpolateThreads;
	float loaded[loadsPerThread];
#pragma unroll
	for (int k = 0; k < loadsPerThread; k++) {
		const int i = threadIdx.x + k * interpolateThreads;
		const NodeIndex local{ i % tileLen, (i / tileLen) % tileLen, i / (tileLen * tileLen) };
		if (i < tileSize)
			loaded[k] = grid[GetGridIndexRealspace(PeriodicBoundaryCondition::applyBC(tileOrigin + local, gridDim), gridDim)];
	}
#pragma unroll
	for (int k = 0; k < loadsPerThread; k++) {
		const int i = threadIdx.x + k * interpolateThreads;
		if (i < tileSize)
			(&tile[0][0][0])[i] = loaded[k];
	}
	__syncthreads();

	for (int i = threadIdx.x; i < nOwned; i += interpolateThreads) {
		const int index = blockIdx.x * ChargeBlock::maxParticlesInBlock + i;
		const ChargePos particle = chargeblockBuffers.owned.particles[index];
		const int slot = chargeblockBuffers.owned.slots[index];
		const Float3 gridPos = particle.pos * gridpointsPerNm_f;

		// A floor index on the grid's upper edge is a grid length off its block's tile
		const NodeIndex fromOrigin = FloorIndex3d(particle.pos) - tileOrigin;
		const NodeIndex wrap = PeriodicBoundaryCondition::applyBC(fromOrigin, gridDim) - fromOrigin;
		ForceEnergy fe = InterpolateForceEnergy(gridPos, [&](int X, int Y, int Z) {
			return tile[Z - tileOrigin.z + wrap.z][Y - tileOrigin.y + wrap.y][X - tileOrigin.x + wrap.x];
			});

		// The gradient is per scaled nm
		fe.force = fe.force * gridScale;

		// Now add self charge to calculations
		fe.force *= particle.charge;
		fe.potE *= particle.charge;

		fe.potE *= 0.5f; // Potential is halved because we computing for both this and the other particle's

		// The reciprocal sum includes each charge's interaction with its own Gaussian, which the Ewald self-energy removes
		fe.potE -= ewaldKappa / sqrtf(PI) * PhysicsUtilsDevice::modifiedCoulombConstant * particle.charge * particle.charge;

		const int scId = slot / SuperCluster::maxParticles;

		if (forceAcc.fx) {
			forceAcc.Add(slot, fe);
		}
		else {
			const int indexInSc = slot % SuperCluster::maxParticles;
			forceEnergies[scMeta[scId]._pclusterIds[indexInSc] * PersistentCluster::maxParticles + scMeta[scId].indexInPcluster[indexInSc]] = fe;
		}
	}
}

/// ------------------------------
/// Tiled version
/// ------------------------------
//
//constexpr int tilePadding = 2;
//constexpr int tileSize = 11 + tilePadding * 2;
//__device__ int GetIndexInTile(
//	const NodeIndex& tileStartGlobal,
//	NodeIndex queryIndex)
//{
//	PeriodicBoundaryCondition::applyHyperpos(tileStartGlobal, queryIndex);
//	NodeIndex relativeIndex = queryIndex - tileStartGlobal;
//	
//	// clamp index between 0 and tileSize-1
//	relativeIndex.x = std::clamp(relativeIndex.x, 0, tileSize - 1);
//	relativeIndex.y = std::clamp(relativeIndex.y, 0, tileSize - 1);
//	relativeIndex.z = std::clamp(relativeIndex.z, 0, tileSize - 1);
//	return (relativeIndex.z * tileSize + relativeIndex.y) * tileSize + relativeIndex.x;
//}
//
//__device__ ForceEnergy InterpolateForceEnergyFromGrid1(const float* _, const float* const tile, Float3 gridPos, Int3 gridDim, NodeIndex tileStart) {
//	int ix = static_cast<int>(floorf(gridPos.x));
//	int iy = static_cast<int>(floorf(gridPos.y));
//	int iz = static_cast<int>(floorf(gridPos.z));
//
//	float fx = gridPos.x - static_cast<float>(ix);
//	float fy = gridPos.y - static_cast<float>(iy);
//	float fz = gridPos.z - static_cast<float>(iz);
//
//	float wx[4], wy[4], wz[4];
//	LAL::CalcBspline(fx, wx);
//	LAL::CalcBspline(fy, wy);
//	LAL::CalcBspline(fz, wz);
//
//	ForceEnergy fe{};
//
//	for (int dz = 0; dz < 4; dz++) {
//		int Z = iz - 1 + dz;
//		float wzCur = wz[dz];
//		for (int dy = 0; dy < 4; dy++) {
//			int Y = iy - 1 + dy;
//			float wyzCur = wzCur * wy[dy];
//			for (int dx = 0; dx < 4; dx++) {
//				int X = ix - 1 + dx;
//				float wxyzCur = wyzCur * wx[dx];
//
//				NodeIndex node = PeriodicBoundaryCondition::applyBC(NodeIndex{ X, Y, Z }, gridDim);
//				int plusX = GetIndexInTile(tileStart, PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x + 1, node.y,     node.z }, gridDim));
//				int minusX = GetIndexInTile(tileStart, PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x - 1, node.y,     node.z }, gridDim));
//				int plusY = GetIndexInTile(tileStart, PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y + 1, node.z }, gridDim));
//				int minusY = GetIndexInTile(tileStart, PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y - 1, node.z }, gridDim));
//				int plusZ = GetIndexInTile(tileStart, PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y,     node.z + 1 }, gridDim));
//				int minusZ = GetIndexInTile(tileStart, PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y,     node.z - 1 }, gridDim));
//				int center = GetIndexInTile(tileStart, node);
//
//				float phi_plusX = tile[plusX];
//				float phi_minusX = tile[minusX];
//				float phi_plusY = tile[plusY];
//				float phi_minusY = tile[minusY];
//				float phi_plusZ = tile[plusZ];
//				float phi_minusZ = tile[minusZ];
//				float phi = tile[center];
//
//				float E_x = -(phi_plusX - phi_minusX) * (gridpointsPerNm / 2.0f);
//				float E_y = -(phi_plusY - phi_minusY) * (gridpointsPerNm / 2.0f);
//				float E_z = -(phi_plusZ - phi_minusZ) * (gridpointsPerNm / 2.0f);
//
//				fe.force += Float3{ E_x, E_y, E_z } *wxyzCur;
//				fe.potE += phi * wxyzCur;				
//			}
//		}
//	}
//
//	return fe;
//}
//
//
//// ================================================================
//// Kernel with shared realspace tile per solvent block
//// ================================================================
//__global__ void InterpolateForcesAndPotentialSolventsTiledVersion( // TODO: OPTIM: We could take a dynamic approach where we launch both this and the simple kernel, and this SM heavy one one runs for dense solvetnblocks?
//	const BoxConfig  config,
//	const BoxState   state,
//	const float* __restrict__ realspaceGrid,
//	Int3             gridDim,               // charge grid dimensions
//	ForceEnergy* __restrict__ forceEnergies,
//	float            selfenergyCorrection,  // [J/mol]
//	Int3             blocksPerDim           // solvent blocks grid (1 nm^3 per block)
//)
//{
//	__shared__ float tile[tileSize * tileSize * tileSize];
//	__shared__ NodeIndex tileStart;
//	__shared__ int nParticles;
//	
//	if (threadIdx.x == 0) {
//		nParticles = state.nParticlesInSolventblock[blockIdx.x];
//
//		//const Float3 gridPos = absPos * gridpointsPerNm_f;
//		NodeIndex tileStart = BoxGrid::Get3dIndex(blockIdx.x % blocksPerDim.InnerProduct(), blocksPerDim) * gridpointsPerNm - NodeIndex{tilePadding, tilePadding , tilePadding };
//		tileStart = PeriodicBoundaryCondition::applyBC(tileStart, gridDim);
//	}
//	__syncthreads();
//
//
//	if (nParticles < 96)
//		return;
//
//	// Cooperative load of realspace tile into shared memory
//	//for (int relIndex = threadIdx.x; relIndex < tileSize * tileSize * tileSize; relIndex += blockDim.x) {
//	if (threadIdx.x < tileSize) {
//		for (int z = 0; z < tileSize; z++) {
//			for (int y = 0; y < tileSize; y++) {
//				NodeIndex globalIndex3d = tileStart + NodeIndex{ threadIdx.x, y, z };
//				globalIndex3d = PeriodicBoundaryCondition::applyBC(globalIndex3d, gridDim);
//
//				int globalIndex = BoxGrid::Get1dIndex(globalIndex3d, gridDim);
//				int indexInTile = BoxGrid::Get1dIndex(NodeIndex{ threadIdx.x, y, z }, Int3{ tileSize, tileSize, tileSize });
//
//				tile[indexInTile] = realspaceGrid[globalIndex];
//			}
//		}
//	}
//	__syncthreads();
//
//
//
//	if (threadIdx.x >= nParticles) {
//		return;
//	}
//
//	const ParticleQuickData pqd = state.solventsParticleQuickData[blockIdx.x * SolventBlock::maxParticles + threadIdx.x];
//	const float charge = pqd.params.charge;
//	if (charge == 0.f)
//		return;
//
//	const NodeIndex origo = BoxGrid::Get3dIndex(blockIdx.x, blocksPerDim);
//	const Float3 relpos = pqd.relPos;
//	Float3 absPos = relpos + origo.toFloat3();
//	PeriodicBoundaryCondition::applyBCNM(absPos);
//
//	const Float3 gridPos = absPos * gridpointsPerNm_f;
//	ForceEnergy fe = InterpolateForceEnergyFromGrid1(realspaceGrid + gridOffset, tile, gridPos, gridDim, tileStart);
//
//	// Now add self charge to calculations
//	fe.force *= charge;
//	fe.potE *= charge;
//
//	// Ewald self-energy correction
//	fe.potE += selfenergyCorrections[simulationId];
//
//	fe.potE *= 0.5f; // Potential is halved because we computing for both this and the other particle's
//
//#ifdef FORCE_NAN_CHECK
//	if (force.isNan()) {
//		printf("PME computed NaN force\n");
//		asm("trap;");
//	}
//#endif
//
//	forceEnergies[blockIdx.x * SolventBlock::maxParticles + threadIdx.x] = fe;
//
//}








__global__ void PrecomputeGreensFunctionKernel(float* d_greensFunction, float* reciprocalVirialScalars, Int3 gridpointsPerDim,
	double boxLenX, double boxLenY, double boxLenZ,		// [nm]
	double ewaldKappa, float gridCellVolumeScale		// ewaldKappa [nm^-1]
) {
	const Int3 halfNodes = gridpointsPerDim / 2;
	int nGridpointsHalfdim = gridpointsPerDim.x / 2 + 1;
	const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
	if (index1D >= nGridpointsHalfdim * gridpointsPerDim.y * gridpointsPerDim.z)
		return;

	NodeIndex freqIndex = Get3dIndexReciprocalspace(index1D, gridpointsPerDim, nGridpointsHalfdim);


	// Remap frequencies to negative for indices > N/2
	int kxIndex = freqIndex.x;
	int kyShiftedIndex = (freqIndex.y <= halfNodes.y) ? freqIndex.y : freqIndex.y - gridpointsPerDim.y;
	int kzShiftedIndex = (freqIndex.z <= halfNodes.z) ? freqIndex.z : freqIndex.z - gridpointsPerDim.z;

	// Physical wavevectors
	double kx = (2.0 * PI * (double)kxIndex) / boxLenX;
	double ky = (2.0 * PI * (double)kyShiftedIndex) / boxLenY;
	double kz = (2.0 * PI * (double)kzShiftedIndex) / boxLenZ;

	double kSquared = kx * kx + ky * ky + kz * kz;

	// Spreading the charges and interpolating the forces each smooth by the cubic B-spline, whose discrete Fourier transform
	// is (2 + cos(k h)) / 3 per dimension. Dividing it out twice makes the reciprocal sum exact up to aliasing (Essmann 1995)
	auto splineModulusSquared = [](double kh) {
		const double modulus = (2.0 + cos(kh)) / 3.0;
		return modulus * modulus;
		};
	const double splineCorrection = 1.0 / (splineModulusSquared(kx * boxLenX / gridpointsPerDim.x)
		* splineModulusSquared(ky * boxLenY / gridpointsPerDim.y) * splineModulusSquared(kz * boxLenZ / gridpointsPerDim.z));

	double currentGreensValue = 0.0;
	if (kSquared > 0.0) {
		currentGreensValue = (4.0 * PI / (kSquared))
			* exp(-kSquared / (4.0 * ewaldKappa * ewaldKappa))
			* splineCorrection
			* PhysicsUtilsDevice::modifiedCoulombConstant
			;
	}

	// The inverse FFT is unnormalized, so the 1/N normalization of the realspace grid is folded in here instead of a separate pass
	const double nGridpointsRealspace = double(gridpointsPerDim.x) * gridpointsPerDim.y * gridpointsPerDim.z;
	d_greensFunction[index1D] = static_cast<float>(currentGreensValue / nGridpointsRealspace);
	reciprocalVirialScalars[index1D] = static_cast<float>(currentGreensValue / nGridpointsRealspace
		* (1.0 - kSquared / (2.0 * ewaldKappa * ewaldKappa)) / (invCellVolume * gridCellVolumeScale));
}


template<bool computeVirial>
__global__ void ApplyGreensFunctionKernel(
	cufftComplex* const d_reciprocalFreqData,
	const float* const d_greensFunctionArray,
	const float* const reciprocalVirialScalars, Int3 gridpointsPerDim, int batchCount,
	long long* const virialPartials
)
{
	int nGridpointsHalfdim = gridpointsPerDim.x / 2 + 1;// TODO: add comment here
	const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
	const int gridSize = nGridpointsHalfdim * gridpointsPerDim.y * gridpointsPerDim.z;
	const bool active = index1D < size_t(gridSize) * batchCount;
	if constexpr (!computeVirial)
		if (!active) return;

	float virial = 0.f;
	if constexpr (computeVirial) {
		if (active) {
			const int mode = index1D % gridSize;
			const NodeIndex frequency = Get3dIndexReciprocalspace(mode, gridpointsPerDim, nGridpointsHalfdim);
			const cufftComplex rho = d_reciprocalFreqData[index1D];
			const float weight = frequency.x == 0 || (gridpointsPerDim.x % 2 == 0 && frequency.x == gridpointsPerDim.x / 2) ? 1.f : 2.f;
			virial = 0.5f * weight * (rho.x * rho.x + rho.y * rho.y) * reciprocalVirialScalars[mode];
		}
		for (int offset = 16; offset > 0; offset >>= 1)
			virial += __shfl_down_sync(__activemask(), virial, offset);
		__shared__ float warpVirials[2];
		if ((threadIdx.x & 31) == 0)
			warpVirials[threadIdx.x >> 5] = virial;
		__syncthreads();
		if (threadIdx.x == 0)
			virialPartials[blockIdx.x] = Virial::ToFixed(warpVirials[0] + warpVirials[1]);
	}

	if (active)
		d_reciprocalFreqData[index1D] = cufftComplex{
			d_reciprocalFreqData[index1D].x * d_greensFunctionArray[index1D % gridSize],
			d_reciprocalFreqData[index1D].y * d_greensFunctionArray[index1D % gridSize]
		};

}


// --------------------------------------------------------------- Controller --------------------------------------------------------------- //	


namespace PME {
	inline void CheckFft(cufftResult result) {
		if (result != CUFFT_SUCCESS) throw std::runtime_error("cuFFT failed: " + std::to_string(static_cast<int>(result)));
	}
}

// The chargeblocks per axis: the box length rounded up, so the grid spacing is at most the nominal 1 / gridpointsPerNm
static Int3 ChargeblocksPerDim(Float3 boxlenNm) {
	const auto Blocks = [](float l) { return (std::max)(1, static_cast<int>(std::ceil(l - boxLengthRoundingTolerance))); };
	return Int3{ Blocks(boxlenNm.x), Blocks(boxlenNm.y), Blocks(boxlenNm.z) };
}

PME::Controller::Controller(const std::vector<EngineSimulationData>& simulations, float cutoffNM, cudaStream_t& stream)
	: boxlenNm(simulations.front().simulation->box->boxparams.BoxSizeFloat()),
	  blocksPerDim(ChargeblocksPerDim(boxlenNm)),
	  nChargeblocks(blocksPerDim.InnerProduct()),
	  ewaldKappa(PhysicsUtils::CalcEwaldkappa(cutoffNM)), stream(stream)
{
	gridScale = NodeIndex(blocksPerDim).toFloat3() / boxlenNm;
	gridpointsPerDim = blocksPerDim * gridpointsPerNm;
	nGridpointsRealspace = size_t(gridpointsPerDim.x) * gridpointsPerDim.y * gridpointsPerDim.z;
	nGridpointsReciprocalspace = gridpointsPerDim.z * gridpointsPerDim.y * (gridpointsPerDim.x / 2 + 1);
	if (nGridpointsRealspace > INT_MAX || nGridpointsRealspace * simulations.size() > INT_MAX)
		throw std::runtime_error("PME batch exceeds grid index range");
	cudaMalloc(&greensFunctionScalars, nGridpointsReciprocalspace * sizeof(float));
	cudaMalloc(&reciprocalVirialScalars, nGridpointsReciprocalspace * sizeof(float));
	PrecomputeGreensFunctionKernel<<<(nGridpointsReciprocalspace + 63) / 64, 64, 0, stream>>>(
		greensFunctionScalars, reciprocalVirialScalars, gridpointsPerDim, boxlenNm.x, boxlenNm.y, boxlenNm.z, ewaldKappa,
		gridScale.x * gridScale.y * gridScale.z);
	LIMA_UTILS::genericErrorCheck(stream, "PrecomputeGreensFunctionKernel");
	SetActiveSimulations(simulations);
}

void PME::Controller::SetActiveSimulations(const std::vector<EngineSimulationData>& simulations) {
	std::vector<int> activeIds, slots(simulations.size(), -1);
	for (int id = 0; id < simulations.size(); ++id) {
		if (!simulations[id].device.active) continue;
		slots[id] = static_cast<int>(activeIds.size());
		activeIds.push_back(id);
	}
	if (activeIds == activeSimulationIds) return;
	cudaStreamSynchronize(stream);
	if (planForward) CheckFft(cufftDestroy(planForward));
	if (planInverse) CheckFft(cufftDestroy(planInverse));
	planForward = planInverse = 0;
	cudaFree(realspaceGrid);
	cudaFree(fourierspaceGrid);
	realspaceGrid = nullptr;
	fourierspaceGrid = nullptr;
	if (chargeblockBuffers) chargeblockBuffers->Free();
	chargeblockBuffers.reset();
	activeSimulationIds = std::move(activeIds);
	batchCount = static_cast<int>(activeSimulationIds.size());
	simulationSlots.SetData(slots);
	if (batchCount == 0) return;
	int dimensions[3]{gridpointsPerDim.z, gridpointsPerDim.y, gridpointsPerDim.x};
	int complexDimensions[3]{gridpointsPerDim.z, gridpointsPerDim.y, gridpointsPerDim.x / 2 + 1};
	CheckFft(cufftPlanMany(&planForward, 3, dimensions, dimensions, 1, static_cast<int>(nGridpointsRealspace),
		complexDimensions, 1, nGridpointsReciprocalspace, CUFFT_R2C, batchCount));
	CheckFft(cufftPlanMany(&planInverse, 3, dimensions, complexDimensions, 1, nGridpointsReciprocalspace,
		dimensions, 1, static_cast<int>(nGridpointsRealspace), CUFFT_C2R, batchCount));
	CheckFft(cufftSetStream(planForward, stream));
	CheckFft(cufftSetStream(planInverse, stream));
	cudaMalloc(&realspaceGrid, nGridpointsRealspace * batchCount * sizeof(float));
	cudaMalloc(&fourierspaceGrid, size_t(nGridpointsReciprocalspace) * batchCount * sizeof(cufftComplex));
	chargeblockBuffers = std::make_unique<ChargeBlock::ChargeblockBuffers>(nChargeblocks * batchCount);
}

PME::Controller::~Controller() {
	cudaStreamSynchronize(stream);
	if (planForward) cufftDestroy(planForward);
	if (planInverse) cufftDestroy(planInverse);
	cudaFree(realspaceGrid);
	cudaFree(fourierspaceGrid);
	cudaFree(greensFunctionScalars);
	cudaFree(reciprocalVirialScalars);
	if (chargeblockBuffers) chargeblockBuffers->Free();
}

const CapacityOverflow* PME::Controller::Overflow() const {
	return chargeblockBuffers ? &chargeblockBuffers->overflow : nullptr;
}

bool PME::Controller::RequiresRebuild(Float3 boxLength) const {
	return GridSizeRequiresRebuild(blocksPerDim, ChargeblocksPerDim(boxLength));
}

void PME::Controller::UpdateBoxLength(Float3 boxLength) {
	boxlenNm = boxLength;
	gridScale = NodeIndex(blocksPerDim).toFloat3() / boxlenNm;
	PrecomputeGreensFunctionKernel<<<(nGridpointsReciprocalspace + 63) / 64, 64, 0, stream>>>(
		greensFunctionScalars, reciprocalVirialScalars, gridpointsPerDim, boxlenNm.x, boxlenNm.y, boxlenNm.z, ewaldKappa,
		gridScale.x * gridScale.y * gridScale.z);
}

void PME::Controller::CalcCharges(SuperCluster* scData, SuperClusterMeta* scMeta, int nSuperclusters, ForceEnergy* forceEnergy,
	ForceAccumulator forceAcc, long long* virialPartials) {
	if (nSuperclusters == 0 || batchCount == 0) return;
	DistributeCompoundchargesToBlocksKernel<<<(nSuperclusters + distributeScsPerBlock - 1) / distributeScsPerBlock, dim3(32, distributeScsPerBlock, 1), 0, stream>>>(
		scData, *chargeblockBuffers, blocksPerDim, scMeta, simulationSlots.Get(), nSuperclusters, gridScale);
	ChargeblockDistributeToGrid<<<nChargeblocks * batchCount, 128, 0, stream>>>(
		*chargeblockBuffers, realspaceGrid, blocksPerDim, gridpointsPerDim, gridScale.x * gridScale.y * gridScale.z);
	CheckFft(cufftExecR2C(planForward, realspaceGrid, fourierspaceGrid));
	const int applyBlocks = VirialPartialCount();
	if (virialPartials)
		ApplyGreensFunctionKernel<true><<<applyBlocks, 64, 0, stream>>>(fourierspaceGrid, greensFunctionScalars,
			reciprocalVirialScalars, gridpointsPerDim, batchCount, virialPartials);
	else
		ApplyGreensFunctionKernel<false><<<applyBlocks, 64, 0, stream>>>(fourierspaceGrid, greensFunctionScalars,
			reciprocalVirialScalars, gridpointsPerDim, batchCount, nullptr);
	CheckFft(cufftExecC2R(planInverse, fourierspaceGrid, realspaceGrid)); // Normalization is folded into the greens function
	InterpolateForcesKernel<<<nChargeblocks * batchCount, interpolateThreads, 0, stream>>>(*chargeblockBuffers, realspaceGrid, blocksPerDim, gridpointsPerDim,
		scMeta, ewaldKappa, forceEnergy, forceAcc, gridScale);
	LIMA_UTILS::genericErrorCheckNoSync("Batched PME");
}

void PME::Controller::PlotPotentialSlices() {
	std::vector<float> gridHost;
	GenericCopyToHost(realspaceGrid, gridHost, nGridpointsRealspace);

	int centerSlice = 20;
	int numSlices = 1;
	int spacing = 5;

	std::vector<float> combinedData;
	std::vector<int> sliceIndices;

	for (int i = -numSlices; i <= numSlices; ++i) {
		int sliceIndex = centerSlice + i * spacing;
		if (sliceIndex < 0 || sliceIndex >= gridpointsPerDim.z) {
			throw std::runtime_error("error");
		}

		int firstIndex = gridpointsPerDim.x * gridpointsPerDim.y * sliceIndex;
		int lastIndex = firstIndex + gridpointsPerDim.x * gridpointsPerDim.y;
		combinedData.insert(combinedData.end(), gridHost.begin() + firstIndex, gridHost.begin() + lastIndex);
		sliceIndices.push_back(sliceIndex);
	}

	// Save all slices to one file
	FileUtils::WriteVectorToBinaryFile("C:/Users/Daniel/git_repo/LIMA_data/Pool/PmePot_AllSlices.bin", combinedData);

	// Call the Python script to plot the slices
	std::string pyscriptPath = (FileUtils::GetLimaDir() / "dev" / "PyTools" / "Plot2dVec.py").string();
	std::string command = "python " + pyscriptPath + " " + std::to_string(numSlices * 2 + 1) + " " + gridpointsPerDim.toString() + " " + std::to_string(centerSlice) + " " + std::to_string(spacing);
	std::cout << "Executing command:\n" << command << std::endl;
	std::system(command.c_str());
}

namespace EngineLimitTesting {
	void ChargeBlock(int count) {
		// All positions are inside the central part of charge block zero, away from stencil transfers.
		// Block one is a sentinel: overflowing block zero can corrupt it without leaving the allocation.
		const int nClusters = (count + SuperCluster::maxParticles - 1) / SuperCluster::maxParticles;
		std::vector<SuperCluster> clusters(nClusters);
		std::vector<SuperClusterMeta> metadata(nClusters);
		for (int sc = 0; sc < nClusters; ++sc) {
			for (int lane = 0; lane < SuperCluster::maxParticles; ++lane) {
				const int id = sc * SuperCluster::maxParticles + lane;
				clusters[sc].SetPdata(PData{ Float3{ 0.5f }, NBParams{ 0.f, id < count ? 0.f : -1.f, id < count ? static_cast<float>(id + 1) : 0.f } }, lane);
			}
		}
		CudaBuffer<SuperCluster> clustersDevice;
		CudaBuffer<SuperClusterMeta> metadataDevice;
		CudaBuffer<int> slots;
		clustersDevice.SetData(clusters);
		metadataDevice.SetData(metadata);
		slots.SetData({0});
		std::unique_ptr<ChargeBlock::ChargeblockBuffers, FreeDeviceMembers<ChargeBlock::ChargeblockBuffers>> buffers(new ChargeBlock::ChargeblockBuffers(8));
		CheckCuda();
		DistributeCompoundchargesToBlocksKernel<<<(nClusters + distributeScsPerBlock - 1) / distributeScsPerBlock, dim3(32, distributeScsPerBlock, 1)>>>(
			clustersDevice.Get(), *buffers, Int3{2, 2, 2}, metadataDevice.Get(), slots.Get(), nClusters, Float3{ 1.f });
		CheckCuda();
		const auto actual = GenericCopyToHost(buffers->chargeposBuffer, 8 * ChargeBlock::maxParticlesInBlock);
		for (size_t i = ChargeBlock::maxParticlesInBlock; i < actual.size(); ++i)
			Require(actual[i].charge == 0.f, "PME charge distribution overwrote a neighboring charge block");
		buffers->overflow.Check();
		Require(count <= ChargeBlock::maxParticlesInBlock, "PME accepted more than 384 charge entries without rejecting the overflow");
		const auto counts = GenericCopyToHost(buffers->nParticlesInBlock, 8);
		Require(counts[0] == count, "PME charge distribution lost entries");
		std::vector<float> charges;
		for (int i = 0; i < count; ++i) charges.push_back(actual[i].charge);
		std::ranges::sort(charges);
		for (int i = 0; i < count; ++i) Require(charges[i] == static_cast<float>(i + 1), "PME charge distribution changed charges");
	}

}
