#pragma once

#include <cuda_runtime.h>


namespace ChargeBlock {
	const int maxParticlesInBlock = 256 + 128; // TODO boxGrid should be able to compute this, based on a size

	struct ChargePos {
		Float3 pos;		// [nm] Relative to a chargeBlock
		float charge;	// [kC/mol]

		constexpr bool operator!= (const ChargePos& a) const { return (a.pos != pos || a.charge != charge); }
	};

	struct ChargeblockBuffers {
		ChargeblockBuffers(int nChargeblocks) {
			cudaMalloc(&nParticlesInBlock, nChargeblocks * sizeof(uint32_t));
			cudaMalloc(&chargeposBuffer, ChargeBlock::maxParticlesInBlock * sizeof(ChargePos) * nChargeblocks);

			cudaMemset(nParticlesInBlock, 0, nChargeblocks * sizeof(uint32_t));
			cudaMemset(chargeposBuffer, 0, nChargeblocks * maxParticlesInBlock * sizeof(ChargePos));
		}

		void Free() const {
			cudaFree(nParticlesInBlock);
			cudaFree(chargeposBuffer);
		}
		uint32_t* nParticlesInBlock;
		ChargePos* chargeposBuffer;
	};


	__device__ ChargePos* GetParticles(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.chargeposBuffer[blockIndex * maxParticlesInBlock];
	}
};