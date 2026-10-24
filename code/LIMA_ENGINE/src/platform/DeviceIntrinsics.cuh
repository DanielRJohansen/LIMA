#pragma once

// The few device-side constructs that differ between CUDA and HIP. Under CUDA each expands to exactly the code it
// replaced, so the CUDA build's SASS is unchanged. Kernels use these instead of #if on the vendor.
// Block-wide async copies are in BlockCopy.cuh, as cooperative_groups adds device globals to every file including it.
// Safe to include from host code: the device helpers only exist when compiling for the GPU (LIMA_GPU_COMPILER)

#include "GpuRuntime.h"

// Lane masks for the warp intrinsics (__shfl_sync, __ballot_sync, ...). HIP requires 64-bit masks; LIMA only targets
// 32-wide warps, so the upper half is always 0 there
#if defined(LIMA_HIP)
using WarpMask = unsigned long long;
#else
using WarpMask = unsigned int;
#endif
#define LIMA_FULL_WARP_MASK WarpMask(0xFFFFFFFFu)

// A __shared__ array of a type with default member initializers (Float3), left uninitialized as all shared memory is.
// nvcc ignores the initializers (warning 20054); clang rejects them, so under HIP it is raw storage viewed as the type.
// Usage: LIMA_SHARED_UNINIT(Float3, positions, [4][n]);
#if defined(LIMA_HIP)
#define LIMA_SHARED_UNINIT(Type, name, dims) 	alignas(Type) __shared__ unsigned char name##Storage[sizeof(Type dims)]; 	auto& name = *reinterpret_cast<Type(*) dims>(name##Storage)
#else
#define LIMA_SHARED_UNINIT(Type, name, dims) __shared__ Type name dims
#endif

#if defined(LIMA_GPU_COMPILER)
namespace DeviceIntrinsics {
	// Crashes the kernel, for unrecoverable device errors
	__device__ __forceinline__ void Trap() {
#if defined(LIMA_HIP)
		__builtin_trap();
#else
		asm("trap;");
#endif
	}

	// 2^a, approximate and flushing denormals to zero
	__device__ __forceinline__ float Exp2Approx(float a) {
#if defined(LIMA_HIP)
		return __builtin_amdgcn_exp2f(a);
#else
		float r;
		asm("ex2.approx.ftz.f32 %0,%1;" : "=f"(r) : "f"(a));
		return r;
#endif
	}

	// A load that bypasses the SM's/CU's own cache, for data another block wrote and fenced during this kernel.
	// On AMD a plain load may hit a stale per-CU L1 line, so it is an agent-scope atomic load there
	template <typename T>
	__device__ __forceinline__ T LoadCoherent(const T* ptr) {
#if defined(LIMA_HIP)
		return __hip_atomic_load(ptr, __ATOMIC_RELAXED, __HIP_MEMORY_SCOPE_AGENT);
#else
		return __ldcg(ptr);
#endif
	}
	#if defined(LIMA_HIP)
	// HIP's atomic load takes scalars only, so an int2 is loaded as one 64-bit word
	__device__ __forceinline__ int2 LoadCoherent(const int2* ptr) {
		const unsigned long long word = __hip_atomic_load(reinterpret_cast<const unsigned long long*>(ptr), __ATOMIC_RELAXED, __HIP_MEMORY_SCOPE_AGENT);
		return int2{ static_cast<int>(static_cast<unsigned int>(word)), static_cast<int>(static_cast<unsigned int>(word >> 32)) };
	}
#endif
}
#endif
