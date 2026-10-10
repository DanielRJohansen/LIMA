#pragma once

#include <cuda_runtime.h>
#include "CapacityOverflow.cuh"


namespace ChargeBlock {
	const int maxParticlesInBlock = 256 + 128; // TODO boxGrid should be able to compute this, based on a size

	struct ChargePos {
		Float3 pos;		// [nm] Relative to a chargeBlock
		float charge;	// [kC/mol]

		// One 16 byte store. The buffer is 16 byte aligned
		__device__ void Store(Float3 position, float charge) { *reinterpret_cast<float4*>(this) = make_float4(position.x, position.y, position.z, charge); }

		constexpr bool operator!= (const ChargePos& a) const { return (a.pos != pos || a.charge != charge); }
	};

	// The particles a chargeblock owns, whose grid floor index is in the block. Their forces are interpolated by the block
	struct OwnedParticles {
		ChargePos* particles;	// pos is absolute, inside the box [nm]
		int* slots;				// Supercluster slot (scId * 16 + index)
	};

	struct ChargeblockBuffers {
		ChargeblockBuffers(int nChargeblocks) {
			cudaMalloc(&nParticlesInBlock, nChargeblocks * sizeof(uint32_t));
			cudaMalloc(&chargeposBuffer, ChargeBlock::maxParticlesInBlock * sizeof(ChargePos) * nChargeblocks);
			cudaMalloc(&nOwnedInBlock, nChargeblocks * sizeof(uint32_t));
			cudaMalloc(&owned.particles, ChargeBlock::maxParticlesInBlock * sizeof(ChargePos) * nChargeblocks);
			cudaMalloc(&owned.slots, ChargeBlock::maxParticlesInBlock * sizeof(int) * nChargeblocks);

			cudaMemset(nParticlesInBlock, 0, nChargeblocks * sizeof(uint32_t));
			cudaMemset(chargeposBuffer, 0, nChargeblocks * maxParticlesInBlock * sizeof(ChargePos));
			cudaMemset(nOwnedInBlock, 0, nChargeblocks * sizeof(uint32_t));
			overflow = CapacityOverflow::Create();
		}

		void Free() const {
			cudaFree(nParticlesInBlock);
			cudaFree(chargeposBuffer);
			cudaFree(nOwnedInBlock);
			cudaFree(owned.particles);
			cudaFree(owned.slots);
			CapacityOverflow{ overflow }.Free();
		}
		uint32_t* nParticlesInBlock;
		ChargePos* chargeposBuffer;		// Particles spreading charge to the block's grid, including from neighbor blocks
		uint32_t* nOwnedInBlock;
		OwnedParticles owned;
		CapacityOverflow overflow;
	};


	__device__ ChargePos* GetParticles(const ChargeblockBuffers buffers, int blockIndex) {
		return &buffers.chargeposBuffer[blockIndex * maxParticlesInBlock];
	}
};
