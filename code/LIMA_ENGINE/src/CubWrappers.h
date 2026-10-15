#pragma once

// The device-wide algorithms the engine uses, in their own translation unit to keep CUB's headers out of Engine.cu.
// Calls CUB directly rather than through thrust, which dispatches to the same CUB algorithms but takes twice as
// long to compile.
// Like the thrust calls they replace, each call synchronizes its stream and throws on CUDA errors, except ExclusiveScan and GatherInts.

#include "GpuRuntime.h"
#include <cstddef>

#include "CudaBuffer.h"

namespace CubWrappers {
	// Stream ordered, does not synchronize. tempStorage is only used by work queued on the stream, so it may be reused by
	// later calls on the same stream; growing it frees the old buffer, which synchronizes the device.
	void ExclusiveScan(const int* first, const int* last, int* result, CudaBuffer<char>& tempStorage, cudaStream_t stream);
	// dst[i] = src[indices[i]], or src[i] without indices. Stream ordered, does not synchronize. Pointers may be mapped host memory
	void GatherInts(const int* src, const int* indices, int count, int* dst, cudaStream_t stream);
	void FillN(int* first, size_t count, int value, cudaStream_t stream);
	// Sums in double precision
	double Sum(const float* first, const float* last, cudaStream_t stream);
	// On the default stream
	void SqrtInPlace(float* data, size_t count);
}
