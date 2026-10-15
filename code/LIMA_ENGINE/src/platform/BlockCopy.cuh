#pragma once

// Block-wide async copy (CUDA: cooperative_groups::memcpy_async). Separate from DeviceIntrinsics.cuh, as the
// cooperative_groups headers add device globals to every file including them, which would shift their SASS
#include "GpuRuntime.h"

#if defined(LIMA_GPU_COMPILER)
#if !defined(LIMA_HIP)
#include <cooperative_groups.h>
#include <cooperative_groups/memcpy_async.h>
#endif

namespace DeviceIntrinsics {
	// Starts copying bytes from src to dst with the whole block. Call BlockCopyWait before reading dst.
	// Typed, as memcpy_async picks its vector width from the element alignment.
	// HIP has no memcpy_async: there the block copies directly, 4 bytes per thread when aligned
	template <typename T>
	__device__ __forceinline__ void BlockCopyAsync(T* dst, const T* src, size_t bytes) {
#if defined(LIMA_HIP)
		const size_t nThreads = size_t(blockDim.x) * blockDim.y * blockDim.z;
		const size_t tid = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
		if constexpr (alignof(T) >= 4) {
			if (bytes % 4 == 0) {
				for (size_t i = tid; i < bytes / 4; i += nThreads)
					reinterpret_cast<unsigned int*>(dst)[i] = reinterpret_cast<const unsigned int*>(src)[i];
				return;
			}
		}
		for (size_t i = tid; i < bytes; i += nThreads)
			reinterpret_cast<char*>(dst)[i] = reinterpret_cast<const char*>(src)[i];
#else
		cooperative_groups::memcpy_async(cooperative_groups::this_thread_block(), dst, src, bytes);
#endif
	}
	// Waits for the block's BlockCopyAsync copies, and makes them visible to the whole block
	__device__ __forceinline__ void BlockCopyWait() {
#if defined(LIMA_HIP)
		__syncthreads();
#else
		cooperative_groups::wait(cooperative_groups::this_thread_block());
#endif
	}
}
#endif
