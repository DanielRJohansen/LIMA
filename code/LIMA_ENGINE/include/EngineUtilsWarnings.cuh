// This file is for print the various warnings the kernels will display
// Each of these functions will be omitted based on the preprocessor flags.
// There will be a separate namespace for each kernel!

// Each function determines if its priority is under LIMASAFEMODE or LIMAPUSH

#pragma once

#include <iostream>

#include "Constants.cuh"
#include "Bodies.cuh"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// TODO HARD: make this a namespace
class EngineUtilsWarnings {
public:
	__device__ static void verifyNodeIndexShiftIsSafe(const NodeIndex& nodeshift_right_to_left) {
#ifdef LIMASAFEMODE
		if (nodeshift_right_to_left.manhattanLen() > MAX_SAFE_SHIFT) {
			printf("Shifting compound further than what is safe! Block %d Thread %d Shift %d\n", blockIdx.x, threadIdx.x, nodeshift_right_to_left.manhattanLen());
		}
#endif	
	}


	__device__ static void logcompoundVerifyVelocity(const Compound& compound, 
		SimParams& simparams, const CompoundCoords& compound_coords, const Float3& force) {
#ifdef LIMASAFEMODE
		if (compound.vels_prev[threadIdx.x] * simparams.constparams.dt > BOXGRID_NODE_LEN_i / 20) {	// Do we move more than 1/20 of a box per step?
			printf("\nParticle %d in compound %d is moving too fast\n", threadIdx.x, blockIdx.x);
			(compound.vels_prev[threadIdx.x] * simparams.constparams.dt).print('V');
			force.print('F');
			//LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(compound_coords.origo).print('O');
			//LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compound_coords.origo, compound_coords.rel_positions[threadIdx.x]).print('P');

			simparams.critical_error_encountered = true;
		}
#endif
	}

	__device__ static void verifyValidRelpos(const Coord& relpos) {
#ifdef LIMASAFEMODE
		const int32_t blocklen_half = BOXGRID_NODE_LEN_i / 2;
		const Coord rel_blockcenter{ blocklen_half };
		if (relpos.x < INT32_MIN + blocklen_half || relpos.y < INT32_MIN + blocklen_half || relpos.z < INT32_MIN + blocklen_half) {
			printf("\nWe have underflow!\n");
			relpos.print('R');
		}
		if (relpos.x > INT32_MAX - blocklen_half || relpos.y > INT32_MAX - blocklen_half || relpos.z > INT32_MAX - blocklen_half) {
			printf("\nWe have overflow!\n");
			relpos.print('R');
		}
#endif
	}

	__device__ static void verifyOrigoShiftIsValid(const NodeIndex& from, const NodeIndex& to) {
#ifdef LIMASAFEMODE
		const NodeIndex origo_shift = from - to;
		if (abs(origo_shift.x) > 10 || abs(origo_shift.y) > 10 || abs(origo_shift.z) > 10) {
			printf("block %d thread %d\n", blockIdx.x, threadIdx.x);
			from.print('f');
			to.print('t');
		}
#endif
	}
};
