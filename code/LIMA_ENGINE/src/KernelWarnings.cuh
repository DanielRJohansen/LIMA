// This file is for print the various warnings the kernels will display
// Each of these functions will be omitted based on the preprocessor flags.
// There will be a separate namespace for each kernel!

// Each function determines if its priority is under LIMASAFEMODE or LIMAPUSH

#include <iostream>

#include "Constants.h"
#include "Bodies.cuh"
#include "EngineBodies.cuh"

//
//namespace EngineUtilsWarnings {
//
//	__device__  inline void verifyNodeIndexShiftIsSafe(const NodeIndex& nodeshift_right_to_left) {
//#ifdef LIMASAFEMODE
//		if (nodeshift_right_to_left.manhattanLen() > MAX_SAFE_SHIFT) {
//			printf("Shifting compound further than what is safe! Block %d Thread %d Shift %d\n", blockIdx.x, threadIdx.x, nodeshift_right_to_left.manhattanLen());
//		}
//#endif
//	}
//}








namespace KernelHelpersWarnings {
	__device__ void assertHyperorigoIsValid(const NodeIndex& querycompound_hyperorigo, const NodeIndex& origo_self) {
#if defined LIMASAFEMODE
		if ((querycompound_hyperorigo - origo_self).maxElement() > 10) {
			printf("Here:: %d\n", blockIdx.x);
			origo_self.print('s');
			querycompound_hyperorigo.print('q');
		}
#endif
	}

	__device__ void verifyOnehotRemaindersIsValid(uint8_t* onehot_remainers, int i) {
#if defined LIMASAFEMODE
		if (onehot_remainers[i] > 230) { printf("Sequential-Prefix-Sum algo is about to crash!"); }
#endif
	}

	__device__ void transferoutVerifyQueueIndex(int queue_index, const NodeIndex& transfer_dir) {
#if defined LIMASAFEMODE
		if (queue_index < 0 || queue_index > 5) { printf("\nGot unexpected queue index %d\n", queue_index); transfer_dir.print(); }
#endif
	}

	__device__ void transferoutVerifyInsertion(bool success) {
#if defined LIMASAFEMODE
		if (!success)
			printf("\nTried to add too many solvents in outgoing transferqueue\n");
#endif
	}

	__device__ void assertValidBlockId(const int blockid) {
#if defined LIMASAFEMODE
		if (blockid < 0 || blockid >= SolventBlocksCircularQueue::blocks_per_grid) {
			printf("\nGot unexpected Block id index %d\n", blockid);
		}
#endif		
	}

	__device__ static void transferOutDebug(STransferQueue* queue_global, 
		const STransferQueue& queue_local, const NodeIndex& transferdir_queue, const int queue_index) {
#if defined LIMASAFEMODE
		if (queue_global->rel_positions[threadIdx.x].x < -2 * static_cast<int32_t>(NANO_TO_LIMA) || queue_global->rel_positions[threadIdx.x].x > 2 * static_cast<int32_t>(NANO_TO_LIMA)
			|| queue_global->rel_positions[threadIdx.x].y < -2 * static_cast<int32_t>(NANO_TO_LIMA) || queue_global->rel_positions[threadIdx.x].y > 2 * static_cast<int32_t>(NANO_TO_LIMA)
			|| queue_global->rel_positions[threadIdx.x].z < -2 * static_cast<int32_t>(NANO_TO_LIMA) || queue_global->rel_positions[threadIdx.x].z > 2 * static_cast<int32_t>(NANO_TO_LIMA)
			) {
			printf("\n");
			transferdir_queue.print('t');
			queue_local.rel_positions[threadIdx.x].print('q');
			queue_global->rel_positions[threadIdx.x].print('Q');
		}

		if (threadIdx.x == 0 && queue_global->n_elements != 0) {
			printf("\nN elements was: %d in queue %d\n", queue_global->n_elements, queue_index);
			transferdir_queue.print('d');
		}
#endif
#ifndef LIMAPUSH
		if (threadIdx.x == 0 && queue_local.n_elements > 15) {
			printf("\nTransferring %d elements\n", queue_local.n_elements);
		}
#endif
	}
}


namespace SolventWarnings {


}






namespace SolventTransferWarnings {
	__device__ void assertSolventsEqualNRemain(const SolventBlock& solventblock_next, const SolventBlockTransfermodule& transfermodule) {
#if defined LIMASAFEMODE
		if (solventblock_next.n_solvents != transfermodule.n_remain) {
			printf("Solventblock_next size doesn't match remain-size %d %d\n", solventblock_next.n_solvents, transfermodule.n_remain);
		}
#endif
	};
	
	__device__ void assertMaxPlacedSolventsIsWithinLimits(int n_solvents_next, bool& critical_error_encountered) {
#if defined LIMASAFEMODE
		if (threadIdx.x == 0 && n_solvents_next >= SolventBlock::MAX_SOLVENTS_IN_BLOCK) {
			printf("Tried to put %d solvents in a single block\n", n_solvents_next);
			critical_error_encountered = true;
		}
#endif
	}
}

namespace BridgeWarnings {
	__device__ void verifyPRefValid(const ParticleReference& p_ref, const CompoundBridge& bridge) {
#if defined LIMASAFEMODE
		if (p_ref.compoundid_local_to_bridge >= bridge.n_compounds) {
			printf("What the fuck! %d %d\n", p_ref.compoundid_local_to_bridge, bridge.n_compounds);
			printf("CIDs %d %d %d %d\n", bridge.compound_ids[0], bridge.compound_ids[1], bridge.compound_ids[2], bridge.compound_ids[3]);
		}
#endif
	}


}
