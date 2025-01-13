#pragma once

#include "BoxGrid.cuh"

#include <cuda_runtime.h>


namespace ChargeBlock {
	const int maxParticlesInBlock = 256; // TODO boxGrid should be able to compute this, based on a size
	const int maxParticlesFromNeighborBlock = 64;
	const int nNeighborBlocks = 3 * 3 * 3 - 1;

	const int maxReservations = 32; // This is probably too low.. especially for small compounds
	//CompoundReservation compoundReservations[maxReservations];

	//const int maxParticlesInNode = 256 + 128;

	struct ChargePos {
		Float3 pos;		// [nm] Relative to a chargeBlock
		float charge;	// [kC/mol]
	};

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
			cudaMalloc(&reservationKeyBuffer, nChargeblocks * sizeof(uint32_t));
			cudaMalloc(&compoundReservationsBuffer, ChargeBlock::maxReservations * sizeof(ChargeBlock::CompoundReservation) * nChargeblocks);
			cudaMalloc(&chargeposBuffer, ChargeBlock::maxParticlesInBlock * sizeof(ChargePos) * nChargeblocks);
			cudaMalloc(&chargeposFromNearbyBlockBuffer, ChargeBlock::maxParticlesFromNeighborBlock * ChargeBlock::nNeighborBlocks * sizeof(ChargePos) * nChargeblocks);

			cudaMemset(reservationKeyBuffer, 0, nChargeblocks * sizeof(uint32_t));
		}

		void Free() const {
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

	__device__ ChargePos* GetParticles(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.chargeposBuffer[blockIndex * maxParticlesInBlock];
	}
	/*__device__ const ChargePos* const GetParticles(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.chargeposBuffer[blockIndex * maxParticlesInNode];
	}*/
	__device__ ChargePos* GetParticlesFromNearbyBlocks(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.chargeposFromNearbyBlockBuffer[blockIndex * maxParticlesFromNeighborBlock * nNeighborBlocks];
	}
	/*__device__ const ChargePos* const GetParticlesFromNearbyBlocks(const ChargeblockBuffers chargeposFromNearbyBlockBuffer, int blockIndex) {
		return &buffers.chargeposFromNearbyBlockBuffer[blockIndex * maxParticlesFromNeighborBlock * nNeighborBlocks];
	}*/
	__device__ CompoundReservation* GetCompoundReservations(const ChargeblockBuffers buffers, int blockIndex) {
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
			if (reservationIndex >= maxReservations || offset + nParticles >= maxParticlesInBlock) {
				printf("Illegal reservation ri: %d offset: %d nP: %d \n", reservationIndex, offset, nParticles);
			}
		}
		//compoundReservations[reservationIndex] = CompoundReservation{ compoundId, offset, nParticles };
		chargeblockBuffers.compoundReservationsBuffer[blockIndex * maxReservations + reservationIndex] = CompoundReservation{ compoundId, offset, nParticles };
		return offset;
	}

	// The number of threads per block (blockDim.x) must match arraySize.
	// ArraySize must be a power of 2.
	//__device__ void BitonicSort(CompoundReservation* reservations) { // TODO: move to LAL
	//	static const int arraySize = maxReservations;
	//	for (int k = 2; k <= arraySize; k *= 2) {
	//		for (int j = k / 2; j > 0; j /= 2) {
	//			int ixj = threadIdx.x ^ j;
	//			if (ixj > threadIdx.x) {
	//				if ((threadIdx.x & k) == 0) {
	//					if (reservations[threadIdx.x].compoundId > reservations[ixj].compoundId) {
	//						// Swap
	//						CompoundReservation temp = reservations[threadIdx.x];
	//						reservations[threadIdx.x] = reservations[ixj];
	//						reservations[ixj] = temp;
	//					}
	//				}
	//				else {
	//					if (reservations[threadIdx.x].compoundId < reservations[ixj].compoundId) {
	//						// Swap
	//						CompoundReservation temp = reservations[threadIdx.x];
	//						reservations[threadIdx.x] = reservations[ixj];
	//						reservations[ixj] = temp;
	//					}
	//				}
	//			}
	//			__syncthreads();
	//		}
	//	}
	//}
};