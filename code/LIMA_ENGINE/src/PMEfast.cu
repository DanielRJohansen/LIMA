#include "PMEfast.cuh"
#include "PhysicsUtilsDevice.cuh"
#include "LimaPositionSystem.cuh"
#include "DeviceAlgorithms.cuh"

using namespace PME;

__device__ __host__ int GetGridIndexRealspace(const Int3& gridIndex, int gridpointsPerDim) {
	return BoxGrid::Get1dIndex(gridIndex, gridpointsPerDim);
	//return gridIndex.x + gridIndex.y * gridpointsPerDim + gridIndex.z * gridpointsPerDim * gridpointsPerDim;
}
__device__ __host__ int GetGridIndexReciprocalspace(const Int3& gridIndex, int gridpointsPerDim, int nGridpointsHalfdim) {
	return gridIndex.x + gridIndex.y * nGridpointsHalfdim + gridIndex.z * nGridpointsHalfdim * gridpointsPerDim;
}
//__device__ NodeIndex Get3dIndexRealspace(int index1d, int gridpointsPerDim) {
//	int z = index1d / (gridpointsPerDim * gridpointsPerDim);
//	index1d -= z * gridpointsPerDim * gridpointsPerDim;
//	int y = index1d / gridpointsPerDim;
//	index1d -= y * gridpointsPerDim;
//	int x = index1d;
//	return NodeIndex{ x, y, z };
//}
__device__ NodeIndex Get3dIndexReciprocalspace(int index1d, int gridpointsPerDim, int nGridpointsHalfdim) {
	int z = index1d / (nGridpointsHalfdim * gridpointsPerDim);
	index1d -= z * nGridpointsHalfdim * gridpointsPerDim;
	int y = index1d / nGridpointsHalfdim;
	index1d -= y * nGridpointsHalfdim;
	int x = index1d;
	return NodeIndex{ x, y, z };
}

__device__ inline int convertRelative3dIndexToAbsolute1d(const NodeIndex& nodeIndex) {
	return (nodeIndex.x + 1) + (nodeIndex.y + 1) * 3 + (nodeIndex.z + 1) * 3 * 3;
}

__device__ inline NodeIndex ConvertAbsolute1dToRelative3d(int index) {
	NodeIndex nodeIndex;

	nodeIndex.z = index / 9 - 1;
	index = index % 9;
	nodeIndex.y = index / 3 - 1;
	nodeIndex.x = index % 3 - 1;

	return nodeIndex;
}

__device__ inline NodeIndex GetParticlesBlockRelativeToCompoundOrigo(const Float3& relpos, bool particleActive) {
	if (particleActive) {
		NodeIndex relativeLocalIndex;
		relativeLocalIndex = LIMAPOSITIONSYSTEM::PositionToNodeIndex(relpos);
		relativeLocalIndex.x = std::max(relativeLocalIndex.x, -1);
		relativeLocalIndex.x = std::min(relativeLocalIndex.x, 1);
		relativeLocalIndex.y = std::max(relativeLocalIndex.y, -1);
		relativeLocalIndex.y = std::min(relativeLocalIndex.y, 1);
		relativeLocalIndex.z = std::max(relativeLocalIndex.z, -1);
		relativeLocalIndex.z = std::min(relativeLocalIndex.z, 1);
		return relativeLocalIndex;
	}
	return NodeIndex{};
}



struct ChargePos {
	Float3 pos;		// [nm] Relative to a chargeBlock
	float charge;	// [kC/mol]
};

namespace ChargeBlock {
	const int maxParticlesInBlock = 64; // TODO boxGrid should be able to compute this, based on a size
	const int maxParticlesFromNeighborBlock = 16;
	const int nNeighborBlocks = 3 * 3 * 3 - 1;

	const int maxReservations = 32; // This is probably too low.. especially for small compounds
	//CompoundReservation compoundReservations[maxReservations];

	const int maxParticlesInNode = 256 + 128;

	// To remain deterministic AND fast, each compound will push their particles data to a block, 
// in the order in which atomicAdd decides.
// This means, when computing the chargeSum in a block, we need to sort wrt compoundIds
// To do this, each compound needs to leave a reservation
	struct CompoundReservation {
		uint32_t compoundId = UINT32_MAX;
		uint32_t firstParticlesOffsetInThisNode = UINT32_MAX;
		uint32_t nParticlesInThisNode = UINT32_MAX;

		__device__ static uint32_t MakeKey(uint32_t nParticles) {
			return uint32_t{ nParticles | 0x0001'0000 };
		}
		__device__ static uint32_t GetReservationIndex(uint32_t key) {
			return key >> 16;
		}
		__device__ static uint32_t GetOffset(uint32_t key) {
			return key & 0xFFFF;
		}
	};


	struct ChargeblockBuffers {
		ChargeblockBuffers(int nChargeblocks) {
			//ChargeBlock::MallocParticles(&chargeposBuffer, &chargeposFromNearbyBlockBuffer, &compoundReservationsBuffer, nChargeblocks);
			cudaMalloc(&reservationKeyBuffer, nChargeblocks * sizeof(uint32_t));
			cudaMalloc(&compoundReservationsBuffer, ChargeBlock::maxReservations * sizeof(ChargeBlock::CompoundReservation) * nChargeblocks);
			cudaMalloc(&chargeposBuffer, ChargeBlock::maxParticlesInNode * sizeof(ChargePos) * nChargeblocks);
			cudaMalloc(&chargeposFromNearbyBlockBuffer, ChargeBlock::maxParticlesFromNeighborBlock * ChargeBlock::nNeighborBlocks * sizeof(ChargePos) * nChargeblocks);
		}
		/*~ChargeblockBuffers() {
			cudaFree(reservationKeyBuffer);
			cudaFree(compoundReservationsBuffer);
			cudaFree(chargeposBuffer);
			cudaFree(chargeposFromNearbyBlockBuffer);
		}*/
		// Overwrite the reservation buffer with 0's, such that prev reservations and chargepos'es are ignored
		void Refresh(int nChargeblocks) {
			cudaMemset(reservationKeyBuffer, 0, nChargeblocks * sizeof(uint32_t));
		}
		void Free() {
			cudaFree(reservationKeyBuffer);
			cudaFree(compoundReservationsBuffer);
			cudaFree(chargeposBuffer);
			cudaFree(chargeposFromNearbyBlockBuffer);
		}
		uint32_t* reservationKeyBuffer;
		ChargeBlock::CompoundReservation* compoundReservationsBuffer;
		ChargePos* chargeposBuffer;
		ChargePos* chargeposFromNearbyBlockBuffer;
	};


	//float charges[maxParticlesInNode];	// This could just be a set of uint8 vals referencing a charge index, as we wont have many different charges


	/*ChargePos particles[maxParticlesInBlock];
	ChargePos particlesFromNearbyBlocks[maxParticlesFromNeighborBlock * nNeighborBlocks];*/

	__device__ ChargePos* const GetParticles(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.chargeposBuffer[blockIndex * maxParticlesInNode];
	}
	/*__device__ const ChargePos* const GetParticles(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.chargeposBuffer[blockIndex * maxParticlesInNode];
	}*/
	__device__ ChargePos* const GetParticlesFromNearbyBlocks(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.chargeposFromNearbyBlockBuffer[blockIndex * maxParticlesFromNeighborBlock * nNeighborBlocks];
	}
	/*__device__ const ChargePos* const GetParticlesFromNearbyBlocks(const ChargeblockBuffers chargeposFromNearbyBlockBuffer, int blockIndex) {
		return &buffers.chargeposFromNearbyBlockBuffer[blockIndex * maxParticlesFromNeighborBlock * nNeighborBlocks];
	}*/
	__device__ CompoundReservation* const GetCompoundReservations(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.compoundReservationsBuffer[blockIndex * maxReservations];
	}
	//const CompoundReservation* const GetCompoundReservations(const CompoundReservation* compoundReservationsBuffer, int blockIndex) {
	//	return &compoundReservationsBuffer[blockIndex * maxReservations];
	//}


	__device__ int nParticlesReserved(const ChargeblockBuffers buffers, int blockIndex) {
		return CompoundReservation::GetOffset(buffers.reservationKeyBuffer[blockIndex]);
	}
	__device__ int nReservations(const ChargeblockBuffers buffers, int blockIndex) {
		return CompoundReservation::GetReservationIndex(buffers.reservationKeyBuffer[blockIndex]);
	}

	// Returns the index where the compounds first particle should be inserted
	__device__ uint32_t MakeReservation(uint32_t compoundId, uint32_t nParticles, ChargeblockBuffers chargeblockBuffers, int blockIndex) {
		const uint32_t key = ChargeBlock::CompoundReservation::MakeKey(nParticles);
		const uint32_t prevKey = atomicAdd(&chargeblockBuffers.reservationKeyBuffer[blockIndex], key);


		const uint32_t reservationIndex = ChargeBlock::CompoundReservation::GetReservationIndex(prevKey);
		const uint32_t offset = ChargeBlock::CompoundReservation::GetOffset(prevKey);

		if constexpr (!LIMA_PUSH) {
			if (reservationIndex >= maxReservations || offset + nParticles >= maxParticlesInNode) {
				printf("Illegal reservation ri: %d offset: %d nP: %d \n", reservationIndex, offset, nParticles);
			}
		}
		//compoundReservations[reservationIndex] = CompoundReservation{ compoundId, offset, nParticles };
		chargeblockBuffers.compoundReservationsBuffer[blockIndex * maxReservations + reservationIndex] = CompoundReservation{ compoundId, offset, nParticles };
		return offset;
	}
	//__device__ void Reset() {
	//	reservationKey = 0;
	//}

};
using namespace ChargeBlock;



// The number of threads per block (blockDim.x) must match arraySize.
// ArraySize must be a power of 2.
__device__ void BitonicSort(ChargeBlock::CompoundReservation* reservations) { // TODO: move to LAL
	static const int arraySize = ChargeBlock::maxReservations;
	for (int k = 2; k <= arraySize; k *= 2) {
		for (int j = k / 2; j > 0; j /= 2) {
			int ixj = threadIdx.x ^ j;
			if (ixj > threadIdx.x) {
				if ((threadIdx.x & k) == 0) {
					if (reservations[threadIdx.x].compoundId > reservations[ixj].compoundId) {
						// Swap
						ChargeBlock::CompoundReservation temp = reservations[threadIdx.x];
						reservations[threadIdx.x] = reservations[ixj];
						reservations[ixj] = temp;
					}
				}
				else {
					if (reservations[threadIdx.x].compoundId < reservations[ixj].compoundId) {
						// Swap
						ChargeBlock::CompoundReservation temp = reservations[threadIdx.x];
						reservations[threadIdx.x] = reservations[ixj];
						reservations[ixj] = temp;
					}
				}
			}
			__syncthreads();
		}
	}
}

__global__ void DistributeCompoundchargesToBlocksKernel(const BoxConfig boxConfig, const BoxState boxState, const ChargeblockBuffers chargeblockBuffers, int blocksPerDim)
{
	__shared__ int numParticlesInNodeLocal[27];
	__shared__ int numParticlesInNodeGlobal[27];

	NodeIndex compoundOrigo = boxState.compoundOrigos[blockIdx.x];
	const Float3 relPos = boxState.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES]; // I want this to be relpos to the block, no?
	const int nParticles = boxConfig.compounds[blockIdx.x].n_particles;
	const float myCharge = boxConfig.compoundsAtomCharges[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];

	// todo: threads with a charge of 0 should NOT push that charge (duh)
	if (threadIdx.x < 27) {
		numParticlesInNodeLocal[threadIdx.x] = 0;
		numParticlesInNodeGlobal[threadIdx.x] = 0;
	}
	__syncthreads();

	// First each particle figure out which node it belongs to relative to it's compounds origo. Then it decide the offset to put values, wrt other particles in this compound
	int myOffsetInLocalNode = 0;																				// Between 0 and maxParticlesInNode
	const NodeIndex relativeNode = GetParticlesBlockRelativeToCompoundOrigo(relPos, threadIdx.x < nParticles);	// Between {-1,-1,-1} and {1,1,1}
	const int localNodeId = convertRelative3dIndexToAbsolute1d(relativeNode);									// Between 0 and 26	
	for (int i = 0; i < nParticles; i++) {
		if (threadIdx.x == i)
			myOffsetInLocalNode = numParticlesInNodeLocal[localNodeId]++;
		__syncthreads();
	}

	// Fetch the current particles in the global node counts, and figure out where to push our counts
	if (threadIdx.x < 27) {
		const NodeIndex absIndex3d = PeriodicBoundaryCondition::applyBC(compoundOrigo + ConvertAbsolute1dToRelative3d(threadIdx.x), blocksPerDim);
		
		// No need to fetch this data, if no particles goes there...
		if (numParticlesInNodeLocal[threadIdx.x] > 0) {
			//numParticlesInNodeGlobal[threadIdx.x] = chargeBlocks[BoxGrid::Get1dIndex(absIndex3d, blocksPerDim)].MakeReservation(blockIdx.x, numParticlesInNodeLocal[threadIdx.x]);
			const int targetBlockIndex = BoxGrid::Get1dIndex(absIndex3d, blocksPerDim);
			//numParticlesInNodeGlobal[threadIdx.x] = chargeblockBuffers.  chargeBlocks[BoxGrid::Get1dIndex(absIndex3d, blocksPerDim)].MakeReservation(blockIdx.x, numParticlesInNodeLocal[threadIdx.x]);
			numParticlesInNodeGlobal[threadIdx.x] = ChargeBlock::MakeReservation(blockIdx.x, numParticlesInNodeLocal[threadIdx.x], chargeblockBuffers, targetBlockIndex);
		}
	}
	__syncthreads();

	// Finally compute the correct index to insert our data, and push that data
	if (threadIdx.x < nParticles) {
		const NodeIndex absoluteTargetblockIndex3D = PeriodicBoundaryCondition::applyBC(compoundOrigo + relativeNode, blocksPerDim);
		//PeriodicBoundaryCondition::applyBC(absoluteTargetblockIndex3D);
		const int absoluteTargetblockIndex = BoxGrid::Get1dIndex(absoluteTargetblockIndex3D, blocksPerDim);

		const int offsetInTargetsBuffer = numParticlesInNodeGlobal[convertRelative3dIndexToAbsolute1d(relativeNode)] + myOffsetInLocalNode;

		const Float3 relposRelativeToTargetBlock = relPos - relativeNode.toFloat3();

		if (absoluteTargetblockIndex < 0 || absoluteTargetblockIndex >= blocksPerDim * blocksPerDim * blocksPerDim)
			printf("Error %d\n", absoluteTargetblockIndex);

		ChargeBlock::GetParticles(chargeblockBuffers, absoluteTargetblockIndex)[offsetInTargetsBuffer] = ChargePos{ relposRelativeToTargetBlock, myCharge };
		//chargeBlocks[absoluteTargetblockIndex].particles[offsetInTargetsBuffer] = ChargePos{ relposRelativeToTargetBlock, myCharge };
	}
}


struct Direction3 {
	uint8_t data; // store 6 bits: top 2 for x+1, next 2 for y+1, bottom 2 for z+1

	__device__ constexpr Direction3() : data(0) {}

	__device__ constexpr Direction3(int x, int y, int z)
		: data(static_cast<uint8_t>(((x + 1) << 4) | ((y + 1) << 2) | (z + 1)))
	{}

	__device__ inline int x() const { return ((data >> 4) & 0x3) - 1; }
	__device__ inline int y() const { return ((data >> 2) & 0x3) - 1; }
	__device__ inline int z() const { return (data & 0x3) - 1; }

	__device__ bool operator==(const Direction3& o) const { return data == o.data; }
	__device__ bool operator!=(const Direction3& o) const { return data != o.data; }
	__device__ Direction3 operator-() const { return Direction3{ -x(), -y(), -z() }; }
};



namespace device_tables {
	__device__ constexpr Direction3 sIndexToDirection[26] = {
	Direction3(-1,-1,-1), Direction3(0,-1,-1), Direction3(1,-1,-1),
	Direction3(-1, 0,-1), Direction3(0, 0,-1), Direction3(1, 0,-1),
	Direction3(-1, 1,-1), Direction3(0, 1,-1), Direction3(1, 1,-1),

	Direction3(-1,-1, 0), Direction3(0,-1, 0), Direction3(1,-1, 0),
	Direction3(-1, 0, 0), /* (0,0,0) excl */   Direction3(1, 0, 0),
	Direction3(-1, 1, 0), Direction3(0, 1, 0), Direction3(1, 1, 0),

	Direction3(-1,-1, 1), Direction3(0,-1, 1), Direction3(1,-1, 1),
	Direction3(-1, 0, 1), Direction3(0, 0, 1), Direction3(1, 0, 1),
	Direction3(-1, 1, 1), Direction3(0, 1, 1), Direction3(1, 1, 1)
	};

	//__device__ constexpr uint8_t sDirectionToIndex[64] = {
	//	  0,   9,  17, 255,   3,  12,  20, 255,
	//	  6,  14,  23, 255, 255, 255, 255, 255,
	//	  1,  10,  18, 255,   4, 255,  21, 255,
	//	  7,  15,  24, 255, 255, 255, 255, 255,
	//	  2,  11,  19, 255,   5,  13,  22, 255,
	//	  8,  16,  25, 255, 255, 255, 255, 255,
	//	255, 255, 255, 255, 255, 255, 255, 255,
	//	255, 255, 255, 255, 255, 255, 255, 255
	//};
}


__device__ bool Floorindex3dShouldBeTransferredThisDirection(const Int3& floorindex3d, const Direction3& queryDirection) {
	if (queryDirection.x() == -1 && floorindex3d.x > 0)
		return false;
	if (queryDirection.x() == 1 && floorindex3d.x < gridpointsPerNm - 2)
		return false;
	if (queryDirection.y() == -1 && floorindex3d.y > 0)
		return false;
	if (queryDirection.y() == 1 && floorindex3d.y < gridpointsPerNm - 2)
		return false;
	if (queryDirection.z() == -1 && floorindex3d.z > 0)
		return false;
	if (queryDirection.z() == 1 && floorindex3d.z < gridpointsPerNm - 2)
		return false;

	return true;
}

__device__ Int3 FloorIndex3d(const Float3& relpos) {
	return Int3{
		static_cast<int>(floorf(relpos.x * gridpointsPerNm_f)),// TODO should this mul not be in int for precision?
		static_cast<int>(floorf(relpos.y * gridpointsPerNm_f)),
		static_cast<int>(floorf(relpos.z * gridpointsPerNm_f))
	};
}


static_assert(MAX_COMPOUND_PARTICLES == ChargeBlock::maxReservations, "Expecting just 1 particle/thread in a reservation");
// Spawn with nthreads = ChargeBlock::maxReservations
__global__ void DistributeOverlappingParticlesToNeighborBlocks(ChargeblockBuffers chargeblockBuffers, int blocksPerDim)
{
	// First we must sort the particles we own, wrt to the reservations to make this deterministic
	// We apply the sort both on local and global mem

	__shared__ ChargeBlock::CompoundReservation reservations[ChargeBlock::maxReservations];
	reservations[threadIdx.x] = ChargeBlock::CompoundReservation{};
	NodeIndex blockIndex3D = BoxGrid::Get3dIndex(blockIdx.x, blocksPerDim);

	const uint32_t nReservations = ChargeBlock::nReservations(chargeblockBuffers, blockIdx.x);
	if (threadIdx.x < nReservations)
		reservations[threadIdx.x] = ChargeBlock::GetCompoundReservations(chargeblockBuffers, blockIdx.x)[threadIdx.x];
	__syncthreads();

	BitonicSort(reservations);
	__syncthreads();

	__shared__ ChargePos myParticles[ChargeBlock::maxParticlesInBlock];
	int nParticlesLoaded = 0;
	for (int reservationIndex = 0; reservationIndex < nReservations; reservationIndex++) {
		if (threadIdx.x < reservations[reservationIndex].nParticlesInThisNode) {
			const int particleIndexUnsorted = reservations[reservationIndex].firstParticlesOffsetInThisNode + threadIdx.x;
			const int particleIndexSorted = nParticlesLoaded + threadIdx.x;
			//myParticles[particleIndexSorted] = myBlock_GlobalMem->particles[particleIndexUnsorted];
			myParticles[particleIndexSorted] = ChargeBlock::GetParticles(chargeblockBuffers, blockIdx.x)[particleIndexUnsorted];
		}
		nParticlesLoaded += reservations[reservationIndex].nParticlesInThisNode;
		__syncthreads(); // Dont think this is needed
	}

	// Now the particles are sorted in local mem, but we also want a them sorted in global mem
	for (int i = threadIdx.x; i < nParticlesLoaded; i+=blockDim.x)
		ChargeBlock::GetParticles(chargeblockBuffers, blockIdx.x)[i] = myParticles[i];

	// Now we have sorted particles in myParticles, and we can start pushing them to the neighbor blocks locally
	// 
	// Temp storage buffer for outgoing particles
	const int maxOutgoingParticles = ChargeBlock::maxParticlesFromNeighborBlock * ChargeBlock::nNeighborBlocks;
	__shared__ ChargePos particlesToNearbyBlocks[maxOutgoingParticles];
	for (int i = threadIdx.x; i < maxOutgoingParticles; i += blockDim.x)
		particlesToNearbyBlocks[i] = ChargePos{ Float3{}, 0.f };
	__syncthreads();

	if (threadIdx.x < 26) {// only first 26 threads have work
		const Float3 blockAbsPos = BoxGrid::Get3dIndex(blockIdx.x, blocksPerDim).toFloat3();

		const int myOffset = threadIdx.x * ChargeBlock::maxParticlesFromNeighborBlock;
		int myOutgoingParticlesCount = 0;
		const Direction3 myDirection = device_tables::sIndexToDirection[threadIdx.x];

		for (int i = 0; i < nParticlesLoaded; i++) {
			const Int3 floorIndex3d = FloorIndex3d(myParticles[i].pos);

			if (Floorindex3dShouldBeTransferredThisDirection(floorIndex3d, myDirection)) 
			{
				if (myOffset + myOutgoingParticlesCount >= maxOutgoingParticles || i >= maxOutgoingParticles) {
					printf("Illegal index %d %d %d mydir %d %d %d\n Floorindex %d %d %d\n", myOffset, myOutgoingParticlesCount, i, myDirection.x(), myDirection.y(), myDirection.z(), 
						floorIndex3d.x, floorIndex3d.y, floorIndex3d.z);		
				}

				particlesToNearbyBlocks[myOffset + myOutgoingParticlesCount] = myParticles[i];
				myOutgoingParticlesCount++;
				//printf("Sending a particle along direction %d %d %d fi %d %d %d\n", myDirection.x(), myDirection.y(), myDirection.z(), floorIndex3d.x, floorIndex3d.y,floorIndex3d.z );
			}
		}
	}
	__syncthreads();

	// Finally push the particles to the neighbor blocks
	for (int i = threadIdx.x; i < maxOutgoingParticles; i += blockDim.x) {
		// Start by figuring out which incomingBuffer to target in the targetBlock
		const Direction3 direction = device_tables::sIndexToDirection[i / ChargeBlock::maxParticlesFromNeighborBlock];

		// Transform current particle, so it's position is relative to the target block
		particlesToNearbyBlocks[i].pos -= Float3{ static_cast<float>(direction.x()), static_cast<float>(direction.y()), static_cast<float>(direction.z()) };

		const NodeIndex targetBlockIndex3D = PeriodicBoundaryCondition::applyBC(blockIndex3D + NodeIndex{ direction.x(), direction.y(), direction.z() }, blocksPerDim);		
		const int targetBlockIndex = BoxGrid::Get1dIndex(targetBlockIndex3D, blocksPerDim);

		//chargeBlocks[targetBlockIndex].particlesFromNearbyBlocks[i] = particlesToNearbyBlocks[i]; // Being a little cheeky with the indices here, as the target block will have the inverse layout of the source block.. But it works, right?
		ChargeBlock::GetParticlesFromNearbyBlocks(chargeblockBuffers, targetBlockIndex)[i] = particlesToNearbyBlocks[i];
	}
}

__global__ void ChargeblockDistributeToGrid(ChargeblockBuffers chargeblockBuffers, float* const realspaceGrid, int blocksPerDim, int gridpointsPerDim) {
	__shared__ ChargePos particles[ChargeBlock::maxParticlesInBlock + ChargeBlock::nNeighborBlocks * ChargeBlock::maxParticlesFromNeighborBlock];

	for (int i = threadIdx.x; i < ChargeBlock::maxParticlesInBlock + ChargeBlock::nNeighborBlocks * ChargeBlock::maxParticlesFromNeighborBlock; i += blockDim.x)
		particles[i].charge = 0.f;


	for (int i = threadIdx.x; i < ChargeBlock::nParticlesReserved(chargeblockBuffers, blockIdx.x); i += blockDim.x) {
		particles[i] = ChargeBlock::GetParticles(chargeblockBuffers, blockIdx.x)[i];
	}
	for (int i = threadIdx.x; i < ChargeBlock::nNeighborBlocks * ChargeBlock::maxParticlesFromNeighborBlock; i += blockDim.x) {
		particles[ChargeBlock::maxParticlesInBlock + i] = ChargeBlock::GetParticlesFromNearbyBlocks(chargeblockBuffers, blockIdx.x)[i];
	}
	__syncthreads();

	// By scaling the charge with the scalar, we can store the charge in the grid as an integer, allowing us to use atomicAdd atomically	
	constexpr double largestPossibleValue = (4. / 6.) * (5. * elementaryChargeToKiloCoulombPerMole) * 8; // Highest bspline coeff * (highestCharge) * maxExpectedParticleNearNode
	constexpr float scalar = static_cast<double>(INT_MAX - 10) / largestPossibleValue; // 10 for safety
	constexpr float invScalar = 1.f / scalar;

	__shared__ int localGrid[gridpointsPerNm * gridpointsPerNm * gridpointsPerNm];
	float* localGridAsFloat = reinterpret_cast<float*>(localGrid);
	for (int i = threadIdx.x; i < gridpointsPerNm * gridpointsPerNm * gridpointsPerNm; i += blockDim.x)
		localGrid[i] = 0;


	for (int i = threadIdx.x; i < ChargeBlock::maxParticlesInBlock + ChargeBlock::nNeighborBlocks * ChargeBlock::maxParticlesFromNeighborBlock; i+=blockDim.x) {
		const float charge = particles[i].charge;	// [kC/mol]
		if (charge == 0.f)
			continue;;

		const Float3 posRelativeToBlock = particles[i].pos;

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
			if (X < 0 || X >= gridpointsPerNm) // Out-of-local-bounds gridpoints are assumed to be handled by the neighbor blocks
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
		if (globalIndex3d.x < 0 || globalIndex3d.x >= gridpointsPerDim || globalIndex3d.y < 0 || globalIndex3d.y >= gridpointsPerDim || globalIndex3d.z < 0 || globalIndex3d.z >= gridpointsPerDim) {
			globalIndex3d.print('g');
			printf("index %d %d %d\n", xOffset, yOffset, zOffset);
		}

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

	float wx[4];
	wx[0] = (1.f - fx) * (1.f - fx) * (1.f - fx) / 6.f;
	wx[1] = (4.f - 6.f * fx * fx + 3.f * fx * fx * fx) / 6.f;
	wx[2] = (1.f + 3.f * fx + 3.f * fx * fx - 3.f * fx * fx * fx) / 6.f;
	wx[3] = (fx * fx * fx) / 6.f;

	float wy[4];
	wy[0] = (1.f - fy) * (1.f - fy) * (1.f - fy) / 6.f;
	wy[1] = (4.f - 6.f * fy * fy + 3.f * fy * fy * fy) / 6.f;
	wy[2] = (1.f + 3.f * fy + 3.f * fy * fy - 3.f * fy * fy * fy) / 6.f;
	wy[3] = (fy * fy * fy) / 6.f;

	float wz[4];
	wz[0] = (1.f - fz) * (1.f - fz) * (1.f - fz) / 6.f;
	wz[1] = (4.f - 6.f * fz * fz + 3.f * fz * fz * fz) / 6.f;
	wz[2] = (1.f + 3.f * fz + 3.f * fz * fz - 3.f * fz * fz * fz) / 6.f;
	wz[3] = (fz * fz * fz) / 6.f;

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

	forceEnergies[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = ForceEnergy{ force, potential };
}



__global__ void PrecomputeGreensFunctionKernel(float* d_greensFunction, int gridpointsPerDim,
	double boxLen,		// [nm]
	double ewaldKappa	// [nm^-1]
) {
	const int halfNodes = gridpointsPerDim / 2;
	int nGridpointsHalfdim = gridpointsPerDim / 2 + 1;
	const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
	/*if (index1D > nGridpointsHalfdim * gridpointsPerDim * gridpointsPerDim)
		return;*/
	NodeIndex freqIndex = Get3dIndexReciprocalspace(index1D, gridpointsPerDim, nGridpointsHalfdim);

	if (freqIndex.x >= nGridpointsHalfdim || freqIndex.y >= gridpointsPerDim || freqIndex.z >= gridpointsPerDim)
		return;

	// Remap frequencies to negative for indices > N/2
	int kxIndex = freqIndex.x;
	int kyShiftedIndex = (freqIndex.y <= halfNodes) ? freqIndex.y : freqIndex.y - gridpointsPerDim;
	int kzShiftedIndex = (freqIndex.z <= halfNodes) ? freqIndex.z : freqIndex.z - gridpointsPerDim;

	// Ewald kappa fixed
	double volume = boxLen * boxLen * boxLen;				// [nm^3]
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
			* PhysicsUtilsDevice::modifiedCoulombConstant_Force
			;
	}

	//if (freqIndex.y == halfNodes || freqIndex.z == halfNodes) {
	//	currentGreensValue = 0.0;
	//}

	d_greensFunction[index1D] = static_cast<float>(currentGreensValue);
}


__global__ void ApplyGreensFunctionKernel(
	cufftComplex* d_reciprocalFreqData,
	const float* d_greensFunctionArray,
	int gridpointsPerDim,
	float boxLength							// [nm]
)
{
	int nGridpointsHalfdim = gridpointsPerDim / 2 + 1;// TODO: add comment here
	const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
	if (index1D >= nGridpointsHalfdim * gridpointsPerDim * gridpointsPerDim)
		return;


	// Retrieve the precomputed Green’s function value at this frequency point.
	const float currentGreenValue = d_greensFunctionArray[index1D];

	// The frequency-domain data is complex, and we multiply by the real Green’s function to apply the long-range correction.
	cufftComplex currentFreqValue = d_reciprocalFreqData[index1D];

	// Scaling the complex frequency data by the Green’s function imposes the correct electrostatic filtering.
	currentFreqValue.x *= currentGreenValue;
	currentFreqValue.y *= currentGreenValue;

	d_reciprocalFreqData[index1D] = currentFreqValue;
}

__global__ void Normalize(float* realspaceGrid, int nGridpointsRealspace, float normalizationFactor) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= nGridpointsRealspace)
		return;

	realspaceGrid[index] *= normalizationFactor;
}




PME::Controller::Controller(int boxLenNm, const Box& box, float cutoffNM) 
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

void PME::Controller::CalcCharges(const BoxConfig& config, const BoxState& state, int nCompounds, ForceEnergy* const forceEnergy) {
	if (nCompounds == 0)
		return;

	cudaMemset(realspaceGrid, 0, nGridpointsRealspace * sizeof(float));// maybe do this async, after the last kernel?
	cudaMemset(fourierspaceGrid, 0, nGridpointsReciprocalspace * sizeof(cufftComplex));
	chargeblockBuffers->Refresh(nChargeblocks);
	cudaDeviceSynchronize();

	DistributeCompoundchargesToBlocksKernel << <nCompounds, MAX_COMPOUND_PARTICLES >> > (config, state, *chargeblockBuffers, boxlenNm);
	LIMA_UTILS::genericErrorCheck("DistributeCompoundchargesToBlocksKernel failed!");
	DistributeOverlappingParticlesToNeighborBlocks << <boxlenNm * boxlenNm * boxlenNm, ChargeBlock::maxReservations >> > (*chargeblockBuffers, boxlenNm);
	LIMA_UTILS::genericErrorCheck("DistributeOverlappingParticlesToNeighborBlocks failed!");
	ChargeblockDistributeToGrid << <boxlenNm * boxlenNm * boxlenNm, 32 >> > (*chargeblockBuffers, realspaceGrid, boxlenNm, gridpointsPerDim);
	LIMA_UTILS::genericErrorCheck("ChargeblockDistributeToGrid failed!");
	cudaDeviceSynchronize();

	ForwardFFT();
	cudaDeviceSynchronize();

	ApplyGreensFunctionKernel << <(nGridpointsReciprocalspace + 63) / 64, 64 >> > (fourierspaceGrid, greensFunctionScalars, gridpointsPerDim, boxlenNm);
	LIMA_UTILS::genericErrorCheckNoSync("ApplyGreensFunctionKernel failed!");

	InverseFFT();
	cudaDeviceSynchronize();

	Normalize << <(nGridpointsRealspace + 63) / 64, 64 >> > (realspaceGrid, nGridpointsRealspace, 1.0 / static_cast<double>(nGridpointsRealspace));

	//PlotPotentialSlices();

	InterpolateForcesAndPotentialKernel << <nCompounds, MAX_COMPOUND_PARTICLES >> > (config, state, realspaceGrid, gridpointsPerDim, forceEnergy, selfenergyCorrection);
	LIMA_UTILS::genericErrorCheckNoSync("InterpolateForcesAndPotentialKernel failed!");

	cudaDeviceSynchronize();
}

