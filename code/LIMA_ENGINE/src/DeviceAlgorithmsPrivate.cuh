#pragma once

#include "KernelConstants.cuh"

#include <cuda_runtime.h>

namespace LAL {
	__device__ void CalcBspline(float f, float* w) {
		if constexpr (!USE_PRECOMPUTED_BSPLINES) {
			w[0] = (1.f - f) * (1.f - f) * (1.f - f) / 6.f;
			w[1] = (4.f - 6.f * f * f + 3.f * f * f * f) / 6.f;
			w[2] = (1.f + 3.f * f + 3.f * f * f - 3.f * f * f * f) / 6.f;
			w[3] = (f * f * f) / 6.f;
		}
		else {
			const float alpha = f * (BSPLINE_LUT_SIZE - 1);
			const int index = static_cast<int>(floor(alpha));
			const int idxUp = min(index + 1, BSPLINE_LUT_SIZE - 1);
			const float frac = alpha - index;   // in [0,1)
			w[0] = bsplineTable_device[index] * (1 - frac)
				+ bsplineTable_device[idxUp] * frac;
			w[1] = bsplineTable_device[BSPLINE_LUT_SIZE + index] * (1 - frac)
				+ bsplineTable_device[BSPLINE_LUT_SIZE + idxUp] * frac;

			float alphaInv = (1.f - f) * (BSPLINE_LUT_SIZE - 1);
			int idxInv = (int)floor(alphaInv);
			int idxInvUp = min(idxInv + 1, BSPLINE_LUT_SIZE - 1);
			float fracInv = alphaInv - idxInv;

			w[2] = bsplineTable_device[BSPLINE_LUT_SIZE + idxInv] * (1 - fracInv)
				+ bsplineTable_device[BSPLINE_LUT_SIZE + idxInvUp] * fracInv;
			w[3] = bsplineTable_device[idxInv] * (1 - fracInv)
				+ bsplineTable_device[idxInvUp] * fracInv;
		}		
	}
}