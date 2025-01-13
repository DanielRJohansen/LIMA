#pragma once

#include <cuda_runtime.h>


namespace ChargeBlock {
	const int maxParticlesInBlock = 256; // TODO boxGrid should be able to compute this, based on a size
	const int maxParticlesFromNeighborBlock = 64;
	const int nNeighborBlocks = 3 * 3 * 3 - 1;

	const int maxReservations = 32; // This is probably too low.. especially for small compounds

	struct ChargePos {
		Float3 pos;		// [nm] Relative to a chargeBlock
		float charge;	// [kC/mol]
	};

	struct ChargeblockBuffers {
		ChargeblockBuffers(int nChargeblocks) {
			cudaMalloc(&reservationKeyBuffer, nChargeblocks * sizeof(uint32_t));
			cudaMalloc(&chargeposBuffer, ChargeBlock::maxParticlesInBlock * sizeof(ChargePos) * nChargeblocks);
			cudaMalloc(&chargeposFromNearbyBlockBuffer, ChargeBlock::maxParticlesFromNeighborBlock * ChargeBlock::nNeighborBlocks * sizeof(ChargePos) * nChargeblocks);

			cudaMemset(reservationKeyBuffer, 0, nChargeblocks * sizeof(uint32_t));
		}

		void Free() const {
			cudaFree(reservationKeyBuffer);
			cudaFree(chargeposBuffer);
			cudaFree(chargeposFromNearbyBlockBuffer);
		}
		uint32_t* reservationKeyBuffer;
		ChargePos* chargeposBuffer;
		ChargePos* chargeposFromNearbyBlockBuffer;
	};

	namespace { // private functions
		__device__ static uint32_t MakeKey(uint32_t nParticles) {
			return uint32_t{ nParticles | 0x0001'0000 };
		}
		__device__ static uint32_t GetOffset(uint32_t key) {
			return key & 0xFFFF;
		}
	}

	__device__ ChargePos* GetParticles(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.chargeposBuffer[blockIndex * maxParticlesInBlock];
	}
	__device__ ChargePos* GetParticlesFromNearbyBlocks(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.chargeposFromNearbyBlockBuffer[blockIndex * maxParticlesFromNeighborBlock * nNeighborBlocks];
	}
	__device__ int nParticlesReserved(const ChargeblockBuffers buffers, int blockIndex) {
		return GetOffset(buffers.reservationKeyBuffer[blockIndex]);
	}
	// Returns the index where the compounds first particle should be inserted
	__device__ uint32_t MakeReservation(uint32_t compoundId, uint32_t nParticles, ChargeblockBuffers chargeblockBuffers, int blockIndex) {
		const uint32_t key = MakeKey(nParticles);
		const uint32_t prevKey = atomicAdd(&chargeblockBuffers.reservationKeyBuffer[blockIndex], key);
		const uint32_t offset = GetOffset(prevKey);
		return offset;
	}
};