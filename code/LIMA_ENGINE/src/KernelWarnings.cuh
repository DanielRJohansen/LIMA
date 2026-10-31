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

	__device__ inline void ForceCheck(Float3 force/*, const std::string message*/) {
		if constexpr (FORCE_CHECKS) {		
			if (isnan(force.x) || isnan(force.y) || isnan(force.z) || isinf(force.x) || isinf(force.y) || isinf(force.z)) {
				//printf("Force check failed: %s\n", message.c_str());
				printf("Force check failed:\n");
			}
		}
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
	
}