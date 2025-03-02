//#pragma once // Only allowed to be included by engine.cu

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "SimulationDevice.cuh"
#include "ChargeBlock.cuh"
#include "DeviceAlgorithmsPrivate.cuh"

#include <cufft.h>
#include <memory>

// TODO: Do i need to account for e0, vacuum/spaceial permitivity here? Probably....


namespace PME {
	const int gridpointsPerNm = 10;
	constexpr float gridpointsPerNm_f = static_cast<float>(gridpointsPerNm);
	constexpr float invCellVolume = static_cast<float>(gridpointsPerNm * gridpointsPerNm * gridpointsPerNm);


	class Controller {
		int gridpointsPerDim = -1;
		size_t nGridpointsRealspace = 0;
		int nGridpointsReciprocalspace = -1;
		const float ewaldKappa;
		float boxlenNm{};
		const int nChargeblocks;

		// Always applied constant per particle
		float selfenergyCorrection;

		// FFT
		float* realspaceGrid;
		cufftComplex* fourierspaceGrid;
		float* greensFunctionScalars;

		// Chargeblocks 
		std::unique_ptr<ChargeBlock::ChargeblockBuffers> chargeblockBuffers;

		cufftHandle planForward;
		cufftHandle planInverse;

		// For system with a net charge, we apply to correction to each realspaceGridnode
		//LAL::optional<float> backgroundchargeCorrection;

		void CalcEnergyCorrection(const Box& box);

	public:

		Controller(int boxLenNm, const Box& box, float cutoffNM, cudaStream_t& stream);
		~Controller();

		void CalcCharges(const BoxConfig& config, const BoxState& state, int nCompounds, ForceEnergy* const forceEnergy, cudaStream_t& stream);

	private:
		//Just for debugging
		void PlotPotentialSlices();

	};
}


// ----------------------------------------- IMPLEMENTATION ----------------------------------------- //


#include "DeviceAlgorithms.cuh"
#include "BoundaryCondition.cuh"
#include "Filehandling.h"
#include "Utilities.h"

using namespace ChargeBlock;
using namespace PME;

// --------------------------------------------------------------- Kernel Helpers --------------------------------------------------------------- //	

constexpr int GetGridIndexRealspace(const Int3& gridIndex, int gridpointsPerDim) {
	return BoxGrid::Get1dIndex(gridIndex, gridpointsPerDim);
}
constexpr int GetGridIndexReciprocalspace(const Int3& gridIndex, int gridpointsPerDim, int nGridpointsHalfdim) {
	return gridIndex.x + gridIndex.y * nGridpointsHalfdim + gridIndex.z * nGridpointsHalfdim * gridpointsPerDim;
}
constexpr NodeIndex Get3dIndexReciprocalspace(int index1d, int gridpointsPerDim, int nGridpointsHalfdim) {
	int z = index1d / (nGridpointsHalfdim * gridpointsPerDim);
	index1d -= z * nGridpointsHalfdim * gridpointsPerDim;
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

constexpr bool Floorindex3dShouldBeTransferredThisDirection(const Int3& floorindex3d, const Direction3& queryDirection) {
	return !(
		(queryDirection.x() == -1 && floorindex3d.x > 0) ||
		(queryDirection.x() == 1 && floorindex3d.x < gridpointsPerNm - 2) ||
		(queryDirection.x() == 0 && (floorindex3d.x < -2 || floorindex3d.x > gridpointsPerNm)) ||
		(queryDirection.y() == -1 && floorindex3d.y > 0) ||
		(queryDirection.y() == 1 && floorindex3d.y < gridpointsPerNm - 2) ||
		(queryDirection.y() == 0 && (floorindex3d.y < -2 || floorindex3d.y > gridpointsPerNm)) ||		
		(queryDirection.z() == -1 && floorindex3d.z > 0) ||
		(queryDirection.z() == 1 && floorindex3d.z < gridpointsPerNm - 2) ||
		(queryDirection.z() == 0 && (floorindex3d.z < -2 || floorindex3d.z > gridpointsPerNm))
		);
}

constexpr Int3 FloorIndex3d(const Float3& relpos) {
	return Int3{
		static_cast<int>(floorf(relpos.x * gridpointsPerNm_f)),// TODO should this mul not be in int for precision?
		static_cast<int>(floorf(relpos.y * gridpointsPerNm_f)),
		static_cast<int>(floorf(relpos.z * gridpointsPerNm_f))
	};
}

// --------------------------------------------------------------- PME Kernels --------------------------------------------------------------- //	

// Call with MAX_COMPOUND_PARTICLES threads (atleast 27)
__global__ void DistributeCompoundchargesToBlocksKernel(const BoxConfig boxConfig, const BoxState boxState, const ChargeblockBuffers chargeblockBuffers, int blocksPerDim)
{
	NodeIndex compoundOrigo = boxState.compoundOrigos[blockIdx.x]; // This coincides as the center block, or 0,0,0 direction block we target
	const int nParticles = boxConfig.compounds[blockIdx.x].n_particles;

	__shared__ Float3 relPositions[MAX_COMPOUND_PARTICLES];
	__shared__ float charges[MAX_COMPOUND_PARTICLES];

	__shared__ int outgoingParticlesId[27 * MAX_COMPOUND_PARTICLES];
	__shared__ int offsetsInTarget[27];
	__shared__ int nOutgoingParticles[27];
	relPositions[threadIdx.x] = boxState.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
	charges[threadIdx.x] = boxConfig.compoundsAtomCharges[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];

	if (threadIdx.x < 27) {
		nOutgoingParticles[threadIdx.x] = 0;
	}
	__syncthreads();

	// The first 27 threads are assigned a direction. They then count which particles are in their node, and store the id's
	if (threadIdx.x < 27) {
		const Direction3 myDirection = device_tables::sIndexToDirection[threadIdx.x];
		const int targetBlockIndex = BoxGrid::Get1dIndex(PeriodicBoundaryCondition::applyBC(compoundOrigo + myDirection.ToNodeIndex(), blocksPerDim), blocksPerDim);
		int myCount = 0;

		for (int i = 0; i < nParticles; i++) {
			if (charges[i] == 0.f) // Skip particles with no charge
				continue;

			const Int3 floorIndex3d = FloorIndex3d(relPositions[i]);
			if (Floorindex3dShouldBeTransferredThisDirection(floorIndex3d, myDirection)) {
				outgoingParticlesId[threadIdx.x * MAX_COMPOUND_PARTICLES + myCount] = i; 
				myCount++;
			}
		}

		// Now reserve space for these particles
		offsetsInTarget[threadIdx.x] = ChargeBlock::MakeReservation(myCount, chargeblockBuffers, targetBlockIndex);
		nOutgoingParticles[threadIdx.x] = myCount;
	}

	// Now all threads collaborate in pushing the outbound particles
	for (int directionIndex = 0; directionIndex < 27; directionIndex++) {
		if (threadIdx.x < nOutgoingParticles[directionIndex]) {
			const Direction3 direction = device_tables::sIndexToDirection[directionIndex];
			
			const int designatedParticleId = outgoingParticlesId[directionIndex * MAX_COMPOUND_PARTICLES + threadIdx.x];
			const Float3 relposRelativeToTargetBlock = relPositions[designatedParticleId] - direction.ToFloat3();

			const int targetBlockIndex = BoxGrid::Get1dIndex(PeriodicBoundaryCondition::applyBC(compoundOrigo + direction.ToNodeIndex(), blocksPerDim), blocksPerDim);
			const int indexInTarget = offsetsInTarget[directionIndex] + threadIdx.x;

			ChargeBlock::GetParticles(chargeblockBuffers, targetBlockIndex)[indexInTarget] = ChargePos{ relposRelativeToTargetBlock, charges[designatedParticleId] };
		}
	}
}


// This kernel accumulates charges of it's particle in a local shared buffer. Any charge outside any kernels volume is ignored, and assumed handled elsewhere
// To remain deterministic, the intermediate charges are stored as integers, such that we simply can use atomicAdd to accumulate charges
__global__ void ChargeblockDistributeToGrid(ChargeblockBuffers chargeblockBuffers, float* const realspaceGrid, int blocksPerDim, int gridpointsPerDim) {
	__shared__  int nParticles;
	__shared__ int localGrid[gridpointsPerNm * gridpointsPerNm * gridpointsPerNm];

	// Init shared variables
	float* localGridAsFloat = reinterpret_cast<float*>(localGrid);
	for (int i = threadIdx.x; i < gridpointsPerNm * gridpointsPerNm * gridpointsPerNm; i += blockDim.x)
		localGrid[i] = 0;

	if (threadIdx.x == 0) {
		nParticles = ChargeBlock::nParticlesReserved(chargeblockBuffers, blockIdx.x);
		chargeblockBuffers.reservationKeyBuffer[blockIdx.x] = 0; // Reset the key buffer		
	}
	__syncthreads();


	// By scaling the charge with the scalar, we can store the charge in the grid as an integer, allowing us to use atomicAdd deterministically	
	constexpr double largestPossibleValue = (4. / 6.) * (5. * elementaryChargeToKiloCoulombPerMole) * 8; // Highest bspline coeff * (highestCharge) * maxExpectedParticleNearNode
	constexpr float scalar = static_cast<double>(INT_MAX - 10) / largestPossibleValue; // 10 for safety
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
					const int index1D = GetGridIndexRealspace(index3d, gridpointsPerNm);
					atomicAdd(&localGrid[index1D], static_cast<int>(charge * wxyCur * wz[dz] * scalar));
				}
			}
		}
	}
	__syncthreads();

	// Transform the grid back to float values, and multiply with invCellVolume
	for (int i = threadIdx.x; i < gridpointsPerNm * gridpointsPerNm * gridpointsPerNm; i += blockDim.x) {
		localGridAsFloat[i] = static_cast<float>(localGrid[i]) * invScalar * invCellVolume;
	}
	__syncthreads();

	// Transform local to global grid coordinates, and push to global memory
	const int nRowsAtATime = blockDim.x / gridpointsPerNm;
	const int myRow = threadIdx.x / gridpointsPerNm;
	if (myRow > nRowsAtATime)
		return;
	const int myIndexInRow = threadIdx.x % gridpointsPerNm;

	const NodeIndex blocksFirstIndex3dInRealspacegrid = BoxGrid::Get3dIndex(blockIdx.x, blocksPerDim) * gridpointsPerNm;
	for (int row = myRow; row < gridpointsPerNm * gridpointsPerNm; row += nRowsAtATime) {
		const int zOffset = row / gridpointsPerNm;
		const int yOffset = row % gridpointsPerNm;
		const int xOffset = myIndexInRow;

		const NodeIndex globalIndex3d = blocksFirstIndex3dInRealspacegrid + NodeIndex{ xOffset, yOffset, zOffset };
        //if (globalIndex3d.x < 0 || globalIndex3d.x >= gridpointsPerDim || globalIndex3d.y < 0 || globalIndex3d.y >= gridpointsPerDim || globalIndex3d.z < 0 || globalIndex3d.z >= gridpointsPerDim) {
        //	globalIndex3d.print('g');
        //	printf("index %d %d %d\n", xOffset, yOffset, zOffset);
        //}

		const int globalIndex = BoxGrid::Get1dIndex(globalIndex3d, gridpointsPerDim);
		const int localIndex = row * gridpointsPerNm + myIndexInRow;
		realspaceGrid[globalIndex] = localGridAsFloat[localIndex];
	}
}

__global__ void InterpolateForcesAndPotentialKernel(
	const BoxConfig config,
	const BoxState state,
	const float* realspaceGrid,
	int gridpointsPerDim,
	ForceEnergy* const forceEnergies,
	float selfenergyCorrection			// [J/mol]
)
{
	if (threadIdx.x >= config.compounds[blockIdx.x].n_particles) {
		return;
	}

	const float charge = config.compoundsAtomCharges[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];	// [kC/mol]
	if (charge == 0.f)
		return;

	const NodeIndex origo = state.compoundOrigos[blockIdx.x];
	const Float3 relpos = state.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
	Float3 absPos = relpos + origo.toFloat3();
	PeriodicBoundaryCondition::applyBCNM(absPos);

	const Float3 gridPos = absPos * gridpointsPerNm_f;
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

	Float3 force{};			// [J/mol/nm]
	float potential{};		// [J/mol]

	for (int dx = 0; dx < 4; dx++) {
		int X = ix - 1 + dx;
		float wxCur = wx[dx];
		for (int dy = 0; dy < 4; dy++) {
			int Y = iy - 1 + dy;
			float wxyCur = wxCur * wy[dy];
			for (int dz = 0; dz < 4; dz++) {
				int Z = iz - 1 + dz;
				float wxyzCur = wxyCur * wz[dz];

				const NodeIndex node = PeriodicBoundaryCondition::applyBC(NodeIndex{ X, Y, Z }, gridpointsPerDim);
				const int gridIndex = GetGridIndexRealspace(node, gridpointsPerDim);

				float phi = realspaceGrid[gridIndex];

				NodeIndex plusX = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x + 1, node.y,     node.z }, gridpointsPerDim);
				NodeIndex minusX = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x - 1, node.y,     node.z }, gridpointsPerDim);
				NodeIndex plusY = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y + 1, node.z }, gridpointsPerDim);
				NodeIndex minusY = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y - 1, node.z }, gridpointsPerDim);
				NodeIndex plusZ = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y,     node.z + 1 }, gridpointsPerDim);
				NodeIndex minusZ = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y,     node.z - 1 }, gridpointsPerDim);

				float phi_plusX = realspaceGrid[GetGridIndexRealspace(plusX, gridpointsPerDim)];
				float phi_minusX = realspaceGrid[GetGridIndexRealspace(minusX, gridpointsPerDim)];
				float phi_plusY = realspaceGrid[GetGridIndexRealspace(plusY, gridpointsPerDim)];
				float phi_minusY = realspaceGrid[GetGridIndexRealspace(minusY, gridpointsPerDim)];
				float phi_plusZ = realspaceGrid[GetGridIndexRealspace(plusZ, gridpointsPerDim)];
				float phi_minusZ = realspaceGrid[GetGridIndexRealspace(minusZ, gridpointsPerDim)];

				float E_x = -(phi_plusX - phi_minusX) * (gridpointsPerNm / 2.0f);
				float E_y = -(phi_plusY - phi_minusY) * (gridpointsPerNm / 2.0f);
				float E_z = -(phi_plusZ - phi_minusZ) * (gridpointsPerNm / 2.0f);

				force += Float3{ E_x, E_y, E_z } *wxyzCur;
				potential += phi * wxyzCur;
			}
		}
	}

	// Now add self charge to calculations
	force *= charge;
	potential *= charge;

	// Ewald self-energy correction
	potential += selfenergyCorrection;

	potential *= 0.5f; // Potential is halved because we computing for both this and the other particle's

#ifdef FORCE_NAN_CHECK
	if (force.isNan()) {
		printf("PME computed NaN force\n");
		asm("trap;");
	}
#endif

	forceEnergies[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = ForceEnergy{ force, potential };
}



__global__ void PrecomputeGreensFunctionKernel(float* d_greensFunction, int gridpointsPerDim,
	double boxLen,		// [nm]
	double ewaldKappa	// [nm^-1]
) {
	const int halfNodes = gridpointsPerDim / 2;
	int nGridpointsHalfdim = gridpointsPerDim / 2 + 1;
	const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
	if (index1D >= nGridpointsHalfdim * gridpointsPerDim * gridpointsPerDim)
		return;

	NodeIndex freqIndex = Get3dIndexReciprocalspace(index1D, gridpointsPerDim, nGridpointsHalfdim);


	// Remap frequencies to negative for indices > N/2
	int kxIndex = freqIndex.x;
	int kyShiftedIndex = (freqIndex.y <= halfNodes) ? freqIndex.y : freqIndex.y - gridpointsPerDim;
	int kzShiftedIndex = (freqIndex.z <= halfNodes) ? freqIndex.z : freqIndex.z - gridpointsPerDim;

	// Ewald kappa fixed
	double delta = boxLen / (double)gridpointsPerDim;		// [nm]

	// Physical wavevectors
	double kx = (2.0 * PI * (double)kxIndex) / boxLen;
	double ky = (2.0 * PI * (double)kyShiftedIndex) / boxLen;
	double kz = (2.0 * PI * (double)kzShiftedIndex) / boxLen;

	double kSquared = kx * kx + ky * ky + kz * kz;

	double currentGreensValue = 0.0f;

	// Compute B-spline structure factor (4th order)
	double kHalfX = kx * (delta * 0.5);
	double kHalfY = ky * (delta * 0.5);
	double kHalfZ = kz * (delta * 0.5);

	const double epsilon = 1e-14; // TODO try to change this

	auto splineFactor = [epsilon](double kh) {
		if (fabs(kh) < epsilon) return 1.0;
		double ratio = sin(kh) / kh;
		return pow(ratio, 4);
		};

	double Sx = splineFactor(kHalfX);
	double Sy = splineFactor(kHalfY);
	double Sz = splineFactor(kHalfZ);

	double splineCorrection = (Sx * Sy * Sz);
	splineCorrection = splineCorrection * splineCorrection; // squared for forward+back interpolation

	//splineCorrection = 1;

	if (kSquared > epsilon) {
		currentGreensValue = (4.0 * PI / (kSquared))
			* exp(-kSquared / (4.0 * ewaldKappa * ewaldKappa))
			* splineCorrection
			* PhysicsUtilsDevice::modifiedCoulombConstant
			;
	}

	d_greensFunction[index1D] = static_cast<float>(currentGreensValue);
}


__global__ void ApplyGreensFunctionKernel(
	cufftComplex* const d_reciprocalFreqData,
	const float* const d_greensFunctionArray,
	int gridpointsPerDim,
	float boxLength							// [nm]
)
{
	int nGridpointsHalfdim = gridpointsPerDim / 2 + 1;// TODO: add comment here
	const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
	if (index1D >= nGridpointsHalfdim * gridpointsPerDim * gridpointsPerDim)
		return;

	d_reciprocalFreqData[index1D] = cufftComplex{
		d_reciprocalFreqData[index1D].x * d_greensFunctionArray[index1D],
		d_reciprocalFreqData[index1D].y * d_greensFunctionArray[index1D]
	};

}

__global__ void Normalize(float* realspaceGrid, int nGridpointsRealspace, float normalizationFactor) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= nGridpointsRealspace)
		return;

	realspaceGrid[index] *= normalizationFactor;
}


// --------------------------------------------------------------- Controller --------------------------------------------------------------- //	


PME::Controller::Controller(int boxLenNm, const Box& box, float cutoffNM, cudaStream_t& stream)
	: boxlenNm(boxLenNm), nChargeblocks(boxLenNm* boxLenNm* boxLenNm), ewaldKappa(PhysicsUtils::CalcEwaldkappa(cutoffNM))
{
	gridpointsPerDim = boxLenNm * gridpointsPerNm;
	nGridpointsRealspace = gridpointsPerDim * gridpointsPerDim * gridpointsPerDim;
	nGridpointsReciprocalspace = gridpointsPerDim * gridpointsPerDim * (gridpointsPerDim / 2 + 1);

	if (nGridpointsRealspace > INT32_MAX)
		throw std::runtime_error("Ewald grid too large to index with integers");

	const size_t byteSize = nGridpointsRealspace * sizeof(float) + nGridpointsReciprocalspace * sizeof(float) * 3; // 3= 2 from cufftComplex, 1 from greensfunction
	if (byteSize > 10'000'000'000)
		throw std::runtime_error("Ewald grid too large");

	CalcEnergyCorrection(box);

	cufftPlan3d(&planForward, gridpointsPerDim, gridpointsPerDim, gridpointsPerDim, CUFFT_R2C);
	cufftPlan3d(&planInverse, gridpointsPerDim, gridpointsPerDim, gridpointsPerDim, CUFFT_C2R);
	cufftSetStream(planForward, stream);
	cufftSetStream(planInverse, stream);

	cudaMalloc(&realspaceGrid, nGridpointsRealspace * sizeof(float));
	cudaMalloc(&fourierspaceGrid, nGridpointsReciprocalspace * sizeof(cufftComplex));
	cudaMalloc(&greensFunctionScalars, nGridpointsReciprocalspace * sizeof(float));
	chargeblockBuffers = std::make_unique<ChargeBlock::ChargeblockBuffers>(nChargeblocks);


	const int nBlocks = (nGridpointsReciprocalspace + 63) / 64;
	PrecomputeGreensFunctionKernel << <nBlocks, 64 >> > (greensFunctionScalars, gridpointsPerDim, boxLenNm, ewaldKappa);
	LIMA_UTILS::genericErrorCheck("PrecomputeGreensFunctionKernel failed!");
}

PME::Controller::~Controller() {
	cudaDeviceSynchronize();

	cufftDestroy(planForward);
	cufftDestroy(planInverse);

	cudaFree(realspaceGrid);
	cudaFree(fourierspaceGrid);
	cudaFree(greensFunctionScalars);

	chargeblockBuffers->Free();
}

void PME::Controller::CalcCharges(const BoxConfig& config, const BoxState& state, int nCompounds, ForceEnergy* const forceEnergy, cudaStream_t& stream) {
	if (nCompounds == 0)
		return;

	DistributeCompoundchargesToBlocksKernel << <nCompounds, MAX_COMPOUND_PARTICLES, 0, stream >> > (config, state, *chargeblockBuffers, boxlenNm);
	LIMA_UTILS::genericErrorCheckNoSync("DistributeCompoundchargesToBlocksKernel failed!");
	ChargeblockDistributeToGrid << <boxlenNm * boxlenNm * boxlenNm, 32, 0, stream >> > (*chargeblockBuffers, realspaceGrid, boxlenNm, gridpointsPerDim);
	LIMA_UTILS::genericErrorCheckNoSync("ChargeblockDistributeToGrid failed!");

	// ForwardFFT
	{
		cufftResult result = cufftExecR2C(planForward, realspaceGrid, fourierspaceGrid);
		if (result != CUFFT_SUCCESS) {
			fprintf(stderr, "cufftExecR2C failed with error code %d\n", result);
		}
	}

	ApplyGreensFunctionKernel << <(nGridpointsReciprocalspace + 63) / 64, 64, 0, stream >> > (fourierspaceGrid, greensFunctionScalars, gridpointsPerDim, boxlenNm);
	LIMA_UTILS::genericErrorCheckNoSync("ApplyGreensFunctionKernel failed!");

	// InverseFFT
	{
		cufftResult result = cufftExecC2R(planInverse, fourierspaceGrid, realspaceGrid);
		if (result != CUFFT_SUCCESS) {
			fprintf(stderr, "cufftExecC2R failed with error code %d\n", result);
		}
	}

	Normalize << <(nGridpointsRealspace + 63) / 64, 64, 0, stream >> > (realspaceGrid, nGridpointsRealspace, 1.0 / static_cast<double>(nGridpointsRealspace));

	//PlotPotentialSlices();

	InterpolateForcesAndPotentialKernel << <nCompounds, MAX_COMPOUND_PARTICLES, 0, stream >> > (config, state, realspaceGrid, gridpointsPerDim, forceEnergy, selfenergyCorrection);
	LIMA_UTILS::genericErrorCheckNoSync("InterpolateForcesAndPotentialKernel failed!");
}

void PME::Controller::CalcEnergyCorrection(const Box& box) {
	double chargeSquaredSum = 0;
	double chargeSum = 0;
	for (const Compound& compound : box.compounds) {
		for (int pid = 0; pid < compound.n_particles; ++pid) {
			const double charge = compound.atom_charges[pid];		// [kC/mol]
			chargeSquaredSum += charge * charge;
			chargeSum += charge;
		}
	}

	selfenergyCorrection = static_cast<float>(-ewaldKappa / sqrt(PI) * chargeSquaredSum * PhysicsUtils::modifiedCoulombConstant);

	/*if (std::abs(chargeSum) > 1e-6) {
		backgroundchargeCorrection = chargeSum / static_cast<double>
	}*/
	// Only relevant for systems with a net charge
	//const float volume = static_cast<float>(box.boxparams.boxSize) * static_cast<float>(box.boxparams.boxSize) * static_cast<float>(box.boxparams.boxSize);	
	//const float backgroundEnergyCorrection = static_cast<float>(-PI * chargeSum * chargeSum / (kappa * kappa * volume) * PhysicsUtils::modifiedCoulombConstant);
	//return selfEnergyCorrection + backgroundEnergyCorrection;		// [J/mol]
}

void PME::Controller::PlotPotentialSlices() {
	std::vector<float> gridHost;
	GenericCopyToHost(realspaceGrid, gridHost, nGridpointsRealspace);

	int centerSlice = 40;
	int numSlices = 1;
	int spacing = 10;

	std::vector<float> combinedData;
	std::vector<int> sliceIndices;

	for (int i = -numSlices; i <= numSlices; ++i) {
		int sliceIndex = centerSlice + i * spacing;
		if (sliceIndex < 0 || sliceIndex >= gridpointsPerDim) {
			throw std::runtime_error("error");
		}

		int firstIndex = gridpointsPerDim * gridpointsPerDim * sliceIndex;
		int lastIndex = firstIndex + gridpointsPerDim * gridpointsPerDim;
		combinedData.insert(combinedData.end(), gridHost.begin() + firstIndex, gridHost.begin() + lastIndex);
		sliceIndices.push_back(sliceIndex);
	}

	// Save all slices to one file
	FileUtils::WriteVectorToBinaryFile("C:/Users/Daniel/git_repo/LIMA_data/Pool/PmePot_AllSlices.bin", combinedData);

	// Call the Python script to plot the slices
	std::string pyscriptPath = (FileUtils::GetLimaDir() / "dev" / "PyTools" / "Plot2dVec.py").string();
	std::string command = "python " + pyscriptPath + " " + std::to_string(numSlices * 2 + 1) + " " + std::to_string(gridpointsPerDim) + " " + std::to_string(centerSlice) + " " + std::to_string(spacing);
	std::system(command.c_str());
}
