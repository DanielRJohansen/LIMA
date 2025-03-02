// This file is for print the various warnings the kernels will display
// Each of these functions will be omitted based on the preprocessor flags.
// There will be a separate namespace for each kernel!

// Each function determines if its priority is under LIMASAFEMODE or LIMAPUSH

#pragma once


#include "Bodies.cuh"
#include "EngineBodies.cuh"

#include "Constants.h"

#include <iostream>





namespace KernelHelpersWarnings {
	__device__ inline void assertHyperorigoIsValid(const NodeIndex& querycompound_hyperorigo, const NodeIndex& origo_self) {
#if defined LIMASAFEMODE
		if ((querycompound_hyperorigo - origo_self).maxElement() > 10) {
			printf("Here:: %d\n", blockIdx.x);
			origo_self.print('s');
			querycompound_hyperorigo.print('q');
		}
#endif
	}

	__device__ inline void verifyOnehotRemaindersIsValid(uint8_t* onehot_remainers, int i) {
#if defined LIMASAFEMODE
		if (onehot_remainers[i] > 230) { printf("Sequential-Prefix-Sum algo is about to crash!"); }
#endif
	}

	__device__ inline void transferoutVerifyQueueIndex(int queue_index, const NodeIndex& transfer_dir) {
#if defined LIMASAFEMODE
		if (queue_index < 0 || queue_index > 5) { printf("\nGot unexpected queue index %d\n", queue_index); transfer_dir.print(); }
#endif
	}

	__device__ inline void transferoutVerifyInsertion(bool success) {
#if defined LIMASAFEMODE
		if (!success)
			printf("\nTried to add too many solvents in outgoing transferqueue\n");
#endif
	}

	__device__ inline void assertValidBlockId(const int blockid) {
#if defined LIMASAFEMODE
		if (blockid < 0 || blockid >= SolventBlocksCircularQueue::blocks_per_grid) {
			printf("\nGot unexpected Block id index %d\n", blockid);
		}
#endif		
	}
}


namespace SolventWarnings {


}






namespace SolventTransferWarnings {
//	__device__ static void assertSolventsEqualNRemain(const SolventBlock& solventblock_next, const SolventBlockTransfermodule& transfermodule) {
//#if defined LIMASAFEMODE
//		if (solventblock_next.n_solvents != transfermodule.n_remain) {
//			printf("Solventblock_next size doesn't match remain-size %d %d\n", solventblock_next.n_solvents, transfermodule.n_remain);
//		}
//#endif
//	};
	
	__device__ inline void assertMaxPlacedSolventsIsWithinLimits(int n_solvents_next, bool& critical_error_encountered) {
#if defined LIMASAFEMODE
		if (threadIdx.x == 0 && n_solvents_next >= SolventBlock::MAX_SOLVENTS_IN_BLOCK) {
			printf("Tried to put %d solvents in a single block\n", n_solvents_next);
			critical_error_encountered = true;
		}
#endif
	}
}