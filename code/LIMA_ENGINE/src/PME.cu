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

constexpr Int3 FloorIndex3d(const Float3& relpos) {
	return Int3{
		static_cast<int>(floorf(relpos.x * gridpointsPerNm_f)),
		static_cast<int>(floorf(relpos.y * gridpointsPerNm_f)),
		static_cast<int>(floorf(relpos.z * gridpointsPerNm_f))
	};
}

// The grid around a chargeblock that its particles' splines reach: the block's gridpoints, plus 1 below and 2 above
constexpr int tileLen = gridpointsPerNm + 3;

// The charges are accumulated as fixed point integers, so the atomics are deterministic. The scale fits the largest charge a
// gridpoint can expect: the highest B-spline coefficient * the highest charge * the particles expected near a gridpoint
constexpr double largestGridCharge = (4. / 6.) * (5. * elementaryChargeToKiloCoulombPerMole) * 8;
constexpr float chargeToFixed = static_cast<double>(INT_MAX - 10) / largestGridCharge / 10.f; // 10 for safety

// --------------------------------------------------------------- PME Kernels --------------------------------------------------------------- //

// Lists each charged particle in the chargeblock owning its grid floor index, which spreads its charge and interpolates its force.
// One warp per supercluster, a lane per particle. blockDim = (32, distributeScsPerBlock, 1)
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

	PeriodicBoundaryCondition::ApplyBC(pos, boxSize, boxSizeInv);
	// Wrapped first, as a position just below 0 can come out of ApplyBC, and integer division would round its floor index -1 up to block 0
	const NodeIndex ownerBlock = PeriodicBoundaryCondition::applyBC(FloorIndex3d(pos), blocksPerDim * gridpointsPerNm) / gridpointsPerNm;
	const int ownerBlockIndex = charged ? blockOffset + BoxGrid::Get1dIndex(ownerBlock, blocksPerDim) : -1;
	const uint32_t sameOwner = __match_any_sync(LIMA_FULL_WARP_MASK, ownerBlockIndex);
	if (!charged)
		return;
	const int firstLane = __ffs(sameOwner) - 1;
	int offsetInOwner = 0;
	if (lane == firstLane)
		offsetInOwner = atomicAdd(&chargeblockBuffers.nOwnedInBlock[ownerBlockIndex], __popc(sameOwner));
	const int indexInOwner = __shfl_sync(WarpMask(sameOwner), offsetInOwner, firstLane) + __popc(sameOwner & ((1u << lane) - 1));
	if (indexInOwner >= ChargeBlock::maxParticlesInBlock) {
		chargeblockBuffers.overflow.Report(CapacityOverflow::ChargeBlock, ownerBlockIndex, indexInOwner + 1, ChargeBlock::maxParticlesInBlock);
		return;
	}
	const int index = ownerBlockIndex * ChargeBlock::maxParticlesInBlock + indexInOwner;
	chargeblockBuffers.owned.particles[index].Store(pos, charge);
	chargeblockBuffers.owned.slots[index] = scId * SuperCluster::maxParticles + lane;
}

// The tile of a chargeblock's particles' splines, and the offset of a particle's floor index from it. A floor index on the
// grid's upper edge is a grid length off its block's tile
struct ChargeblockTile {
	NodeIndex origin;	// Gridpoint of the tile's first point, 1 below the block's
	__device__ ChargeblockTile(int blockIndex, Int3 blocksPerDim)
		: origin(BoxGrid::Get3dIndex(blockIndex % blocksPerDim.InnerProduct(), blocksPerDim) * gridpointsPerNm - NodeIndex{ 1, 1, 1 }) {}
	__device__ NodeIndex Wrap(const Float3& pos, Int3 gridDim) const {
		const NodeIndex fromOrigin = FloorIndex3d(pos) - origin;
		return PeriodicBoundaryCondition::applyBC(fromOrigin, gridDim) - fromOrigin;
	}
};

// Spreads the charges of the particles a chargeblock owns onto a tile of the grid around it in shared memory, then adds the tile
// to the fixed point grid, which the neighbor blocks' tiles overlap. One chargeblock per block, blockDim = spreadThreads
constexpr int spreadThreads = 128;
__global__ void __launch_bounds__(spreadThreads) SpreadChargesKernel(const ChargeblockBuffers chargeblockBuffers, int* const fixedPointGrid, Int3 blocksPerDim, Int3 gridDim)
{
	__shared__ int tile[tileLen][tileLen][tileLen];
	constexpr int tileSize = tileLen * tileLen * tileLen;
	for (int i = threadIdx.x; i < tileSize; i += spreadThreads)
		(&tile[0][0][0])[i] = 0;
	const int nOwned = min(chargeblockBuffers.nOwnedInBlock[blockIdx.x], static_cast<uint32_t>(ChargeBlock::maxParticlesInBlock));
	const ChargeblockTile tileInfo(blockIdx.x, blocksPerDim);
	__syncthreads();

	for (int i = threadIdx.x; i < nOwned; i += spreadThreads) {
		const ChargePos particle = chargeblockBuffers.owned.particles[blockIdx.x * ChargeBlock::maxParticlesInBlock + i];
		const Float3 gridPos = particle.pos * gridpointsPerNm_f;
		const int ix = static_cast<int>(floorf(gridPos.x));
		const int iy = static_cast<int>(floorf(gridPos.y));
		const int iz = static_cast<int>(floorf(gridPos.z));
		float wx[4], wy[4], wz[4];
		LAL::CalcBspline(gridPos.x - static_cast<float>(ix), wx);
		LAL::CalcBspline(gridPos.y - static_cast<float>(iy), wy);
		LAL::CalcBspline(gridPos.z - static_cast<float>(iz), wz);
		// The spline's first gridpoint, 1 below the floor index, in the tile
		const NodeIndex wrap = tileInfo.Wrap(particle.pos, gridDim);
		const int x0 = ix - 1 - tileInfo.origin.x + wrap.x, y0 = iy - 1 - tileInfo.origin.y + wrap.y, z0 = iz - 1 - tileInfo.origin.z + wrap.z;
		const float charge = particle.charge * chargeToFixed;
#pragma unroll
		for (int dz = 0; dz < 4; dz++)
#pragma unroll
			for (int dy = 0; dy < 4; dy++) {
				const float chargeYZ = charge * wy[dy] * wz[dz];
#pragma unroll
				for (int dx = 0; dx < 4; dx++) {
					const float chargeScaled = chargeYZ * wx[dx];
					atomicAdd(&tile[z0 + dz][y0 + dy][x0 + dx], static_cast<int>(fminf(fmaxf(chargeScaled, static_cast<float>(INT_MIN)), static_cast<float>(INT_MAX))));
				}
			}
	}
	__syncthreads();

	int* const grid = fixedPointGrid + size_t(blockIdx.x / blocksPerDim.InnerProduct()) * gridDim.InnerProduct();
	for (int i = threadIdx.x; i < tileSize; i += spreadThreads) {
		const int value = (&tile[0][0][0])[i];
		if (value == 0)
			continue;
		const NodeIndex local{ i % tileLen, (i / tileLen) % tileLen, i / (tileLen * tileLen) };
		atomicAdd(&grid[BoxGrid::Get1dIndex(PeriodicBoundaryCondition::applyBC(tileInfo.origin + local, gridDim), gridDim)], value);
	}
}

// The fixed point charges to the grid of charge densities the FFT takes, and zeroes the fixed point grid for the next step
// gridCellVolumeScale: the volume of a scaled grid cell per real one, see DistributeCompoundchargesToBlocksKernel
__global__ void FixedPointGridToDensityKernel(int* const fixedPointGrid, float* const realspaceGrid, int nGridpoints, float gridCellVolumeScale) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nGridpoints)
		return;
	realspaceGrid[i] = static_cast<float>(fixedPointGrid[i]) * (1.f / chargeToFixed) * invCellVolume * gridCellVolumeScale;
	fixedPointGrid[i] = 0;
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

// Interpolates the forces on the particles a chargeblock owns, from a tile of the grid around the block in shared memory.
// One chargeblock per block, blockDim = interpolateThreads
constexpr int interpolateThreads = 128;
__global__ void __launch_bounds__(interpolateThreads) InterpolateForcesKernel(const ChargeblockBuffers chargeblockBuffers, const float* const realspaceGrid,
	Int3 blocksPerDim, Int3 gridDim, const SuperClusterMeta* const scMeta, float ewaldKappa /*[nm^-1]*/,
	ForceEnergy* const forceEnergies /*EM only*/, const ForceAccumulator forceAcc /*MD only, forceAcc.fx == nullptr in EM*/,
	Float3 gridScale /*see DistributeCompoundchargesToBlocksKernel*/)
{
	__shared__ float tile[tileLen][tileLen][tileLen];
	__shared__ int nOwned;

	const int nBlocksPerGrid = blocksPerDim.InnerProduct();
	const float* const grid = realspaceGrid + size_t(blockIdx.x / nBlocksPerGrid) * gridDim.InnerProduct();
	const ChargeblockTile tileInfo(blockIdx.x, blocksPerDim);
	const NodeIndex tileOrigin = tileInfo.origin;

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

		const NodeIndex wrap = tileInfo.Wrap(particle.pos, gridDim);
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

// Revisit: tiling by density. InterpolateForcesKernel always loads its block's grid tile into shared memory. A removed
// variant (in git history, "Tiled version" in PME.cu) used the tile only for dense blocks (>= 96 particles) and left sparse
// blocks to a kernel reading the grid directly, so sparse blocks skip the tile load. Measure on a system with mostly
// sparse chargeblocks (eg. a membrane or vacuum-padded box) before trying it.

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
	int nGridpointsHalfdim = gridpointsPerDim.x / 2 + 1;	// R2C keeps only x/2+1 frequencies along x, the rest follow by Hermitian symmetry
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
	cudaFree(fixedPointGrid);
	cudaFree(fourierspaceGrid);
	realspaceGrid = nullptr;
	fixedPointGrid = nullptr;
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
	cudaMalloc(&fixedPointGrid, nGridpointsRealspace * batchCount * sizeof(int));
	cudaMemsetAsync(fixedPointGrid, 0, nGridpointsRealspace * batchCount * sizeof(int), stream);
	cudaMalloc(&fourierspaceGrid, size_t(nGridpointsReciprocalspace) * batchCount * sizeof(cufftComplex));
	chargeblockBuffers = std::make_unique<ChargeBlock::ChargeblockBuffers>(nChargeblocks * batchCount);
}

PME::Controller::~Controller() {
	cudaStreamSynchronize(stream);
	if (planForward) cufftDestroy(planForward);
	if (planInverse) cufftDestroy(planInverse);
	cudaFree(realspaceGrid);
	cudaFree(fixedPointGrid);
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
	SpreadChargesKernel<<<nChargeblocks * batchCount, spreadThreads, 0, stream>>>(*chargeblockBuffers, fixedPointGrid, blocksPerDim, gridpointsPerDim);
	const int nGridpoints = static_cast<int>(nGridpointsRealspace) * batchCount;
	FixedPointGridToDensityKernel<<<(nGridpoints + 255) / 256, 256, 0, stream>>>(fixedPointGrid, realspaceGrid, nGridpoints,
		gridScale.x * gridScale.y * gridScale.z);
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
		// All positions are inside charge block zero.
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
		const auto actual = GenericCopyToHost(buffers->owned.particles, 8 * ChargeBlock::maxParticlesInBlock);
		for (size_t i = ChargeBlock::maxParticlesInBlock; i < actual.size(); ++i)
			Require(actual[i].charge == 0.f, "PME charge distribution overwrote a neighboring charge block");
		buffers->overflow.Check();
		Require(count <= ChargeBlock::maxParticlesInBlock, "PME accepted more than 384 charge entries without rejecting the overflow");
		const auto counts = GenericCopyToHost(buffers->nOwnedInBlock, 8);
		Require(counts[0] == count, "PME charge distribution lost entries");
		std::vector<float> charges;
		for (int i = 0; i < count; ++i) charges.push_back(actual[i].charge);
		std::ranges::sort(charges);
		for (int i = 0; i < count; ++i) Require(charges[i] == static_cast<float>(i + 1), "PME charge distribution changed charges");
	}

}
