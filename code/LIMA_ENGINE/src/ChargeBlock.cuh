#pragma once

#include <cuda_runtime.h>


namespace ChargeBlock {
	const int maxParticlesInBlock = 256; // TODO boxGrid should be able to compute this, based on a size

	struct ChargePos {
		Float3 pos;		// [nm] Relative to a chargeBlock
		float charge;	// [kC/mol]
	};

	struct ChargeblockBuffers {
		ChargeblockBuffers(int nChargeblocks) {
			cudaMalloc(&reservationKeyBuffer, nChargeblocks * sizeof(uint32_t));
			cudaMalloc(&chargeposBuffer, ChargeBlock::maxParticlesInBlock * sizeof(ChargePos) * nChargeblocks);

			cudaMemset(reservationKeyBuffer, 0, nChargeblocks * sizeof(uint32_t));
		}

		void Free() const {
			cudaFree(reservationKeyBuffer);
			cudaFree(chargeposBuffer);
		}
		uint32_t* reservationKeyBuffer;
		ChargePos* chargeposBuffer;
	};

	namespace { // private functions
		constexpr uint32_t MakeKey(uint32_t nParticles) {
			return uint32_t{ nParticles | 0x0001'0000 };
		}
		constexpr uint32_t GetOffset(uint32_t key) {
			return key & 0xFFFF;
		}
	}

	__device__ ChargePos* GetParticles(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.chargeposBuffer[blockIndex * maxParticlesInBlock];
	}
	__device__ int nParticlesReserved(const ChargeblockBuffers buffers, int blockIndex) {
		return GetOffset(buffers.reservationKeyBuffer[blockIndex]);
	}
	// Returns the index where the compounds first particle should be inserted
	__device__ uint32_t MakeReservation(uint32_t nParticles, ChargeblockBuffers chargeblockBuffers, int blockIndex) {
		const uint32_t key = MakeKey(nParticles);
		const uint32_t prevKey = atomicAdd(&chargeblockBuffers.reservationKeyBuffer[blockIndex], key);
		const uint32_t offset = GetOffset(prevKey);
		return offset;
	}
};