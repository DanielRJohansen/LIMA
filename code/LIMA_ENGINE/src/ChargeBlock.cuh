#pragma once

#include "GpuRuntime.h"
#include "CapacityOverflow.cuh"


namespace ChargeBlock {
	const int maxParticlesInBlock = 256 + 128;

	struct ChargePos {
		Float3 pos;		// [nm] Relative to a chargeBlock
		float charge;	// [kC/mol]

		// One 16 byte store. The buffer is 16 byte aligned
		__device__ void Store(Float3 position, float charge) { *reinterpret_cast<float4*>(this) = make_float4(position.x, position.y, position.z, charge); }

		constexpr bool operator!= (const ChargePos& a) const { return (a.pos != pos || a.charge != charge); }
	};

	// The particles a chargeblock owns, whose grid floor index is in the block. The block spreads their charges and interpolates their forces
	struct OwnedParticles {
		ChargePos* particles;	// pos is absolute, inside the box [nm]
		int* slots;				// Supercluster slot (scId * 16 + index)
	};

	struct ChargeblockBuffers {
		ChargeblockBuffers(int nChargeblocks) {
			cudaMalloc(&nOwnedInBlock, nChargeblocks * sizeof(uint32_t));
			cudaMalloc(&owned.particles, ChargeBlock::maxParticlesInBlock * sizeof(ChargePos) * nChargeblocks);
			cudaMalloc(&owned.slots, ChargeBlock::maxParticlesInBlock * sizeof(int) * nChargeblocks);

			cudaMemset(nOwnedInBlock, 0, nChargeblocks * sizeof(uint32_t));
			cudaMemset(owned.particles, 0, nChargeblocks * maxParticlesInBlock * sizeof(ChargePos));
			overflow = CapacityOverflow::Create();
		}

		void Free() const {
			cudaFree(nOwnedInBlock);
			cudaFree(owned.particles);
			cudaFree(owned.slots);
			CapacityOverflow{ overflow }.Free();
		}
		uint32_t* nOwnedInBlock;
		OwnedParticles owned;
		CapacityOverflow overflow;
	};

};
