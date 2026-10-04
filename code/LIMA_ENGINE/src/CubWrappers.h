#pragma once

// The device-wide algorithms the engine uses, in their own translation unit to keep CUB's headers out of Engine.cu.
// Calls CUB directly rather than through thrust, which dispatches to the same CUB algorithms but takes twice as
// long to compile.
// Like the thrust calls they replace, each call synchronizes its stream and throws on CUDA errors.

#include <cuda_runtime.h>
#include <cstddef>

namespace CubWrappers {
	void ExclusiveScan(const int* first, const int* last, int* result, cudaStream_t stream);
	void FillN(int* first, size_t count, int value, cudaStream_t stream);
	// Sums in double precision
	double Sum(const float* first, const float* last, cudaStream_t stream);
	// On the default stream
	void SqrtInPlace(float* data, size_t count);
}
