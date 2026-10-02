#pragma once

// Thin wrappers around the thrust algorithms the engine uses. Thrust's headers take ~16 s to compile, so they are
// kept out of Engine.cu (the slowest file in the build) and only ThrustWrappers.cu includes them.
// The calls are the same as before, so the generated code is unchanged.

#include <cuda_runtime.h>
#include <cstddef>

namespace ThrustWrappers {
	void ExclusiveScan(const int* first, const int* last, int* result, cudaStream_t stream);
	void FillN(int* first, size_t count, int value, cudaStream_t stream);
	// Sums in double precision
	double Sum(const float* first, const float* last, cudaStream_t stream);
	// On the default stream
	void SqrtInPlace(float* data, size_t count);
}
