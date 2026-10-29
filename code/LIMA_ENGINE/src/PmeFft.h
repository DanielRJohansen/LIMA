#pragma once

// The host side of PME's z pass, see PmeFft.cuh

#include <cmath>
#include <vector>
#include "GpuRuntime.h"

namespace PmeFft {
	constexpr int maxLength = 256;
	constexpr int maxStages = 8;
	constexpr int threads = 256;

	// An axis length's radices, 4s first
	struct Axis {
		int length = 0;
		int nStages = 0;
		int radices[maxStages]{};
	};

	// Whether the length factors into radices 2, 3, 4 and 5, and is at most maxLength
	inline bool Factor(int length, Axis& axis) {
		axis = Axis{ length };
		for (const int radix : { 4, 2, 3, 5 })
			while (length % radix == 0 && axis.nStages < maxStages) {
				axis.radices[axis.nStages++] = radix;
				length /= radix;
			}
		return length == 1 && axis.length <= maxLength;
	}

	// exp(-2 pi i m / length) for m < length
	inline std::vector<float2> Twiddles(int length) {
		std::vector<float2> twiddles(length);
		for (int m = 0; m < length; m++) {
			const double angle = -2.0 * 3.14159265358979323846 * m / length;
			twiddles[m] = float2{ static_cast<float>(std::cos(angle)), static_cast<float>(std::sin(angle)) };
		}
		return twiddles;
	}

	// x frequencies per block: 16, or 8 for axes longer than 96, so blocks fit 4 per SM
	__host__ __device__ inline int PencilWidth(int length) { return length > 96 ? 8 : 16; }
	inline size_t SharedBytes(int length) { return sizeof(float2) * length * (1 + 2 * PencilWidth(length)); }
}
