#pragma once

/*
// Precomputed LUTs are intentionally disabled, but retained for possible reuse.
namespace DeviceConstants {
	static constexpr int BSPLINE_LUT_SIZE = 128;
	// We store only w0, w1 => total size = 2*N.
	__constant__ float bsplineTable[2 * BSPLINE_LUT_SIZE]; // Precomputed fourth-order B-splines [0, 1].

	static constexpr int ERFC_LUT_SIZE = 32;
	__constant__ float erfcForcescalarTable[ERFC_LUT_SIZE]; // Precomputed scalars [0, 1], where 1 = cutoffNM.
	__constant__ float erfcPotentialscalarTable[ERFC_LUT_SIZE];
}
*/
