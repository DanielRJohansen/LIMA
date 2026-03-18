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
		Int3 gridpointsPerDim{};
		size_t nGridpointsRealspace = 0;
		int nGridpointsReciprocalspace = -1;	// TODO: Does this also need to be size_t?
		const float ewaldKappa;
		Float3 boxlenNm{};
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

		Controller(const Box& box, float cutoffNM, cudaStream_t& stream);
		~Controller();

		void CalcCharges(SuperCluster* const scData, SuperClusterMeta* const scMeta, int nSuperclusters, ForceEnergy* const forceEnergy, cudaStream_t& stream);

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
		static_cast<int>(floorf(relpos.x * gridpointsPerNm_f)),
		static_cast<int>(floorf(relpos.y * gridpointsPerNm_f)),
		static_cast<int>(floorf(relpos.z * gridpointsPerNm_f))
	};
}

// --------------------------------------------------------------- PME Kernels --------------------------------------------------------------- //	

// blockDim = (32, 1, 1)
__global__ void DistributeCompoundchargesToBlocksKernel(const SuperCluster* const superclusters, const ChargeblockBuffers chargeblockBuffers, Int3 blocksPerDim)
{
	__shared__ Float3 relPositions[SuperCluster::nParticles ];
	__shared__ float charges[SuperCluster::nParticles];

	__shared__ int outgoingParticlesId[27 * SuperCluster::nParticles];
	__shared__ int offsetsInTarget[27];
	__shared__ int nOutgoingParticles[27];

	NodeIndex nearestGridnode = superclusters[blockIdx.x].pData[0].position.Floor().ToInt3();

	if (threadIdx.x < SuperCluster::nParticles) {
		PData pqd = superclusters[blockIdx.x].pData[threadIdx.x];

		if (pqd.Valid()) {
			Float3 scNodeOrigoPos = nearestGridnode.toFloat3();// superclusters[blockIdx.x].pData[0].position.Floor();
			relPositions[threadIdx.x] = pqd.position - scNodeOrigoPos;// +Float3{ 0.5, 0.5, 0.5 };
			charges[threadIdx.x] = pqd.params.charge;
		}
		else {
			relPositions[threadIdx.x] = Float3{ NAN, NAN, NAN };
			charges[threadIdx.x] = 0.f;
		}
	}
	for (int i = threadIdx.x; i < 27; i+=blockDim.x) {
		nOutgoingParticles[threadIdx.x] = 0;
	}
	__syncthreads();

	

	// The first 27 threads are assigned a direction. They then count which particles are in their node, and store the id's
	if (threadIdx.x < 27) {
		const Direction3 myDirection = device_tables::sIndexToDirection[threadIdx.x];
		const int targetBlockIndex = BoxGrid::Get1dIndex(PeriodicBoundaryCondition::applyBC(nearestGridnode + myDirection.ToNodeIndex(), blocksPerDim), blocksPerDim);
		int myCount = 0;

		for (int i = 0; i < SuperCluster::nParticles; i++) {
			if (charges[i] == 0.f || isnan(charges[i])) // Skip particles with no charge
				continue;

			const Int3 floorIndex3d = FloorIndex3d(relPositions[i]);
			if (Floorindex3dShouldBeTransferredThisDirection(floorIndex3d, myDirection)) {
				outgoingParticlesId[threadIdx.x * SuperCluster::nParticles + myCount] = i; 
				myCount++;
			}
		}

		// Now reserve space for these particles
		offsetsInTarget[threadIdx.x] = ChargeBlock::MakeReservation(myCount, chargeblockBuffers, targetBlockIndex);
		nOutgoingParticles[threadIdx.x] = myCount;
	}
	__syncthreads();

	// Now all threads collaborate in pushing the outbound particles
	for (int directionIndex = 0; directionIndex < 27; directionIndex++) {
		if (threadIdx.x < nOutgoingParticles[directionIndex]) {
			const Direction3 direction = device_tables::sIndexToDirection[directionIndex];
			
			const int designatedParticleId = outgoingParticlesId[directionIndex * SuperCluster::nParticles + threadIdx.x];
			const Float3 relposRelativeToTargetBlock = relPositions[designatedParticleId] - direction.ToFloat3();

			const int targetBlockIndex = BoxGrid::Get1dIndex(PeriodicBoundaryCondition::applyBC(nearestGridnode + direction.ToNodeIndex(), blocksPerDim), blocksPerDim);
			const int indexInTarget = offsetsInTarget[directionIndex] + threadIdx.x;

			ChargeBlock::GetParticles(chargeblockBuffers, targetBlockIndex)[indexInTarget] = ChargePos{ relposRelativeToTargetBlock, charges[designatedParticleId] };
		}
	}
}



// This kernel accumulates charges of it's particle in a local shared buffer. Any charge outside any kernels volume is ignored, and assumed handled elsewhere
// To remain deterministic, the intermediate charges are stored as integers, such that we simply can use atomicAdd to accumulate charges
__global__ void ChargeblockDistributeToGrid(ChargeblockBuffers chargeblockBuffers, float* const realspaceGrid, Int3 blocksPerDim, Int3 gridpointsPerDim) {
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
					const int index1D = GetGridIndexRealspace(index3d, Int3(gridpointsPerNm, gridpointsPerNm, gridpointsPerNm));
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

__device__ ForceEnergy InterpolateForceEnergyFromGrid(const float* realspaceGrid, Float3 gridPos, Int3 gridDim) {
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

	ForceEnergy fe{};

	for (int dz = 0; dz < 4; dz++) {
		int Z = iz - 1 + dz;
		float wzCur = wz[dz];
		for (int dy = 0; dy < 4; dy++) {
			int Y = iy - 1 + dy;
			float wyzCur = wzCur * wy[dy];
			for (int dx = 0; dx < 4; dx++) {
				int X = ix - 1 + dx;
				float wxyzCur = wyzCur * wx[dx];

				const NodeIndex node = PeriodicBoundaryCondition::applyBC(NodeIndex{ X, Y, Z }, gridDim);
				const int gridIndex = GetGridIndexRealspace(node, gridDim);

				float phi = realspaceGrid[gridIndex];

				NodeIndex plusX = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x + 1, node.y,     node.z }, gridDim);
				NodeIndex minusX = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x - 1, node.y,     node.z }, gridDim);
				NodeIndex plusY = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y + 1, node.z }, gridDim);
				NodeIndex minusY = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y - 1, node.z }, gridDim);
				NodeIndex plusZ = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y,     node.z + 1 }, gridDim);
				NodeIndex minusZ = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y,     node.z - 1 }, gridDim);

				float phi_plusX = realspaceGrid[GetGridIndexRealspace(plusX, gridDim)];
				float phi_minusX = realspaceGrid[GetGridIndexRealspace(minusX, gridDim)];
				float phi_plusY = realspaceGrid[GetGridIndexRealspace(plusY, gridDim)];
				float phi_minusY = realspaceGrid[GetGridIndexRealspace(minusY, gridDim)];
				float phi_plusZ = realspaceGrid[GetGridIndexRealspace(plusZ, gridDim)];
				float phi_minusZ = realspaceGrid[GetGridIndexRealspace(minusZ, gridDim)];

				float E_x = -(phi_plusX - phi_minusX) * (gridpointsPerNm / 2.0f);
				float E_y = -(phi_plusY - phi_minusY) * (gridpointsPerNm / 2.0f);
				float E_z = -(phi_plusZ - phi_minusZ) * (gridpointsPerNm / 2.0f);

				fe.force += Float3{ E_x, E_y, E_z } *wxyzCur;
				fe.potE += phi * wxyzCur;
			}
		}
	}

	return fe;
}


__device__ ForceEnergy InterpolateForceEnergyFromGrid1(const float* realspaceGrid, Float3 gridPos, Int3 gridDim) {
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

	ForceEnergy fe{};

	for (int dz = 0; dz < 4; dz++) {
		int Z = iz - 1 + dz;
		float wzCur = wz[dz];
		for (int dy = 0; dy < 4; dy++) {
			int Y = iy - 1 + dy;
			float wyzCur = wzCur * wy[dy];


			// Load all values used in this YZ plane
			float phisPlusY[4];
			float phisMinusY[4];
			float phisPlusZ[4];
			float phisMinusZ[4];
			float phisCenter[6];

#pragma unroll
			for (int dx = 0; dx < 4; dx++) {
				const NodeIndex node = PeriodicBoundaryCondition::applyBC(NodeIndex{ ix - 1 + dx, Y + 1, Z }, gridDim);
				phisPlusY[dx] = realspaceGrid[GetGridIndexRealspace(node, gridDim)];
			}

#pragma unroll
			for (int dx = 0; dx < 4; dx++) {
				const NodeIndex node = PeriodicBoundaryCondition::applyBC(NodeIndex{ ix - 1 + dx, Y - 1, Z }, gridDim);
				phisMinusY[dx] = realspaceGrid[GetGridIndexRealspace(node, gridDim)];
			}
#pragma unroll
			for (int dx = 0; dx < 4; dx++) {
				const NodeIndex node = PeriodicBoundaryCondition::applyBC(NodeIndex{ ix - 1 + dx, Y, Z + 1 }, gridDim);
				phisPlusZ[dx] = realspaceGrid[GetGridIndexRealspace(node, gridDim)];
			}
#pragma unroll
			for (int dx = 0; dx < 4; dx++) {
				const NodeIndex node = PeriodicBoundaryCondition::applyBC(NodeIndex{ ix - 1 + dx, Y, Z - 1 }, gridDim);
				phisMinusZ[dx] = realspaceGrid[GetGridIndexRealspace(node, gridDim)];
			}





//#pragma unroll
//			for (int dx = 0; dx < 4; dx++) {
//				int X = ix - 1 + dx;
//				const NodeIndex node = PeriodicBoundaryCondition::applyBC(NodeIndex{ X, Y, Z }, gridDim);
//				NodeIndex plusY = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y + 1, node.z }, gridDim);
//				NodeIndex minusY = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y - 1, node.z }, gridDim);
//				NodeIndex plusZ = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y,     node.z + 1 }, gridDim);
//				NodeIndex minusZ = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y,     node.z - 1 }, gridDim);
//				phisPlusY[dx] = realspaceGrid[GetGridIndexRealspace(plusY, gridDim)];
//				phisMinusY[dx] = realspaceGrid[GetGridIndexRealspace(minusY, gridDim)];
//				phisPlusZ[dx] = realspaceGrid[GetGridIndexRealspace(plusZ, gridDim)];
//				phisMinusZ[dx] = realspaceGrid[GetGridIndexRealspace(minusZ, gridDim)];
//			}


#pragma unroll
			for (int dx = 0; dx < 6; dx++) {
				int X = ix - 2 + dx;
				const NodeIndex node = PeriodicBoundaryCondition::applyBC(NodeIndex{ X, Y, Z }, gridDim);
				phisCenter[dx] = realspaceGrid[GetGridIndexRealspace(node, gridDim)];
			}

#pragma unroll
			for (int dx = 0; dx < 4; dx++) {
				int X = ix - 1 + dx;
				float wxyzCur = wyzCur * wx[dx];

				float& phi = phisCenter[dx + 1];
				float E_x = -(phisCenter[dx + 2] - phisCenter[dx]) * (gridpointsPerNm / 2.0f);
				float E_y = -(phisPlusY[dx] - phisMinusY[dx]) * (gridpointsPerNm / 2.0f);
				float E_z = -(phisPlusZ[dx] - phisMinusZ[dx]) * (gridpointsPerNm / 2.0f);

				fe.force += Float3{ E_x, E_y, E_z } *wxyzCur;
				fe.potE += phi * wxyzCur;
			}
		}
	}

	return fe;
}


// blockDim = (SuperCluster::nParticles, 1, 1)
__global__ void InterpolateForcesAndPotentialCompounds(
	SuperCluster* const scData,
	SuperClusterMeta* const scMeta,
	const float* realspaceGrid,
	Int3 gridDim,
	ForceEnergy* const forceEnergies,
	float selfenergyCorrection			// [J/mol]
)
{
	PData pqd = scData[blockIdx.x].pData[threadIdx.x];
	if (!pqd.Valid())
		return;

	const float charge = pqd.params.charge;	// [kC/mol]
	if (charge == 0.f)
		return;

	Float3 absPos = pqd.position;
	PeriodicBoundaryCondition::applyBCNM(absPos);

	const Float3 gridPos = absPos * gridpointsPerNm_f;
	ForceEnergy fe = InterpolateForceEnergyFromGrid1(realspaceGrid, gridPos, gridDim);

	// Now add self charge to calculations
	fe.force *= charge;
	fe.potE *= charge;

	// Ewald self-energy correction
	fe.potE += selfenergyCorrection;

	fe.potE *= 0.5f; // Potential is halved because we computing for both this and the other particle's

#ifdef FORCE_NAN_CHECK
	if (force.isNan()) {
		printf("PME computed NaN force\n");
		asm("trap;");
	}
#endif

	int pcId = scMeta[blockIdx.x].pclusterIds[threadIdx.x/4];
	int pid = threadIdx.x % 4;	
	forceEnergies[pcId * PersistentCluster::nParticles + pid] = fe;
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
//		NodeIndex tileStart = BoxGrid::Get3dIndex(blockIdx.x, blocksPerDim) * gridpointsPerNm - NodeIndex{tilePadding, tilePadding , tilePadding };
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
//	const float charge = DeviceConstants::tinymolForcefield.types[pqd.atomType].charge;
//	if (charge == 0.f)
//		return;
//
//	const NodeIndex origo = BoxGrid::Get3dIndex(blockIdx.x, blocksPerDim);
//	const Float3 relpos = pqd.relPos;
//	Float3 absPos = relpos + origo.toFloat3();
//	PeriodicBoundaryCondition::applyBCNM(absPos);
//
//	const Float3 gridPos = absPos * gridpointsPerNm_f;
//	ForceEnergy fe = InterpolateForceEnergyFromGrid1(realspaceGrid, tile, gridPos, gridDim, tileStart);
//
//	// Now add self charge to calculations
//	fe.force *= charge;
//	fe.potE *= charge;
//
//	// Ewald self-energy correction
//	fe.potE += selfenergyCorrection;
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








__global__ void PrecomputeGreensFunctionKernel(float* d_greensFunction, Int3 gridpointsPerDim,
	Double3 boxLen,		// [nm]
	double ewaldKappa	// [nm^-1]
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

	// Ewald kappa fixed
	//double delta = boxLen / (double)gridpointsPerDim;		// [nm]
	double delta = std::min(std::min(	// TODO: Compute this on host instead
		boxLen.x / (double)gridpointsPerDim.x,
		boxLen.y / (double)gridpointsPerDim.y),
		boxLen.z / (double)gridpointsPerDim.z);	// [nm]


	// Physical wavevectors
	double kx = (2.0 * PI * (double)kxIndex) / boxLen.x;
	double ky = (2.0 * PI * (double)kyShiftedIndex) / boxLen.y;
	double kz = (2.0 * PI * (double)kzShiftedIndex) / boxLen.z;

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
	Int3 gridpointsPerDim
)
{
	int nGridpointsHalfdim = gridpointsPerDim.x / 2 + 1;// TODO: add comment here
	const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
	if (index1D >= nGridpointsHalfdim * gridpointsPerDim.y * gridpointsPerDim.z)
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


PME::Controller::Controller(const Box& box, float cutoffNM, cudaStream_t& stream)
	: boxlenNm(box.boxparams.BoxSizeFloat()), nChargeblocks(box.boxparams.boxSize.InnerProduct()), ewaldKappa(PhysicsUtils::CalcEwaldkappa(cutoffNM))
{
	gridpointsPerDim = box.boxparams.boxSize * gridpointsPerNm;
	nGridpointsRealspace = gridpointsPerDim.x * gridpointsPerDim.y * gridpointsPerDim.z;
	nGridpointsReciprocalspace = gridpointsPerDim.z * gridpointsPerDim.y * (gridpointsPerDim.x / 2 + 1);

	if (nGridpointsRealspace > INT32_MAX)
		throw std::runtime_error("Ewald grid too large to index with integers");

	const size_t byteSize = nGridpointsRealspace * sizeof(float) + nGridpointsReciprocalspace * sizeof(float) * 3; // 3= 2 from cufftComplex, 1 from greensfunction
	if (byteSize > 10'000'000'000)
		throw std::runtime_error("Ewald grid too large");

	CalcEnergyCorrection(box);

	cufftPlan3d(&planForward, gridpointsPerDim.z, gridpointsPerDim.y, gridpointsPerDim.x, CUFFT_R2C);	// TODO I think this should be the opposite way around
	cufftPlan3d(&planInverse, gridpointsPerDim.z, gridpointsPerDim.y, gridpointsPerDim.x, CUFFT_C2R);
	cufftSetStream(planForward, stream);
	cufftSetStream(planInverse, stream);

	cudaMalloc(&realspaceGrid, nGridpointsRealspace * sizeof(float));
	cudaMalloc(&fourierspaceGrid, nGridpointsReciprocalspace * sizeof(cufftComplex));
	cudaMalloc(&greensFunctionScalars, nGridpointsReciprocalspace * sizeof(float));
	chargeblockBuffers = std::make_unique<ChargeBlock::ChargeblockBuffers>(nChargeblocks);


	const int nBlocks = (nGridpointsReciprocalspace + 63) / 64;
	PrecomputeGreensFunctionKernel << <nBlocks, 64 >> > (greensFunctionScalars, gridpointsPerDim, Double3{ boxlenNm }, ewaldKappa);
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

void PME::Controller::CalcCharges(SuperCluster* const scData, SuperClusterMeta* const scMeta, int nSuperclusters, ForceEnergy* const forceEnergy, cudaStream_t& stream) {
	if (nSuperclusters == 0)
		return;

	Int3 bpd = boxlenNm.ToInt3();
	DistributeCompoundchargesToBlocksKernel << <nSuperclusters, 32, 0, stream >> > (scData, *chargeblockBuffers, bpd);
	LIMA_UTILS::genericErrorCheckNoSync("DistributeCompoundchargesToBlocksKernel failed!");

	ChargeblockDistributeToGrid<<<bpd.InnerProduct(), 32, 0, stream >> > (*chargeblockBuffers, realspaceGrid, bpd, gridpointsPerDim);
	LIMA_UTILS::genericErrorCheckNoSync("ChargeblockDistributeToGrid failed!");

	// ForwardFFT
	{
		cufftResult result = cufftExecR2C(planForward, realspaceGrid, fourierspaceGrid);
		if (result != CUFFT_SUCCESS) {
			fprintf(stderr, "cufftExecR2C failed with error code %d\n", result);
		}
	}

	ApplyGreensFunctionKernel << <(nGridpointsReciprocalspace + 63) / 64, 64, 0, stream >> > (fourierspaceGrid, greensFunctionScalars, gridpointsPerDim);
	LIMA_UTILS::genericErrorCheckNoSync("ApplyGreensFunctionKernel failed!");

	// InverseFFT
	{
		cufftResult result = cufftExecC2R(planInverse, fourierspaceGrid, realspaceGrid);
		if (result != CUFFT_SUCCESS) {
			fprintf(stderr, "cufftExecC2R failed with error code %d\n", result);
		}
	}

	Normalize << <(nGridpointsRealspace + 63) / 64, 64, 0, stream >> > (realspaceGrid, nGridpointsRealspace, 1.0 / static_cast<double>(nGridpointsRealspace));	

	InterpolateForcesAndPotentialCompounds << <nSuperclusters, SuperCluster::nParticles, 0, stream >> > (scData, scMeta, realspaceGrid, gridpointsPerDim, forceEnergy, selfenergyCorrection);
	LIMA_UTILS::genericErrorCheckNoSync("InterpolateForcesAndPotentialCompounds failed!");

	//PlotPotentialSlices();
}

void PME::Controller::CalcEnergyCorrection(const Box& box) {
	double chargeSquaredSum = 0;
	double chargeSum = 0;


	for (const auto& pc : box.persistentClusters) {
		for (const auto& pqd : pc.pqd) {
			if (!pqd.Valid())
				continue;
			chargeSquaredSum += pqd.params.charge * pqd.params.charge;
			chargeSum += pqd.params.charge;
		}
	}

    selfenergyCorrection = static_cast<float>(-ewaldKappa / std::sqrt(PI) * chargeSquaredSum * PhysicsUtils::modifiedCoulombConstant);

	 //Only relevant for systems with a net charge
	/*const float volume = static_cast<float>(box.boxparams.boxSize.InnerProduct());	
	const float backgroundEnergyCorrection = static_cast<float>(-PI * chargeSum * chargeSum / (ewaldKappa * ewaldKappa * volume) * PhysicsUtils::modifiedCoulombConstant);*/
	///return selfEnergyCorrection + backgroundEnergyCorrection;		// [J/mol]
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
