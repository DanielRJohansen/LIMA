#pragma once

// The GPU runtime, included everywhere instead of cuda_runtime.h. All projects use the CUDA runtime names (cudaMalloc,
// cudaStream_t, ...); a HIP build (LIMA_GPU_BACKEND=HIP, which defines LIMA_HIP) maps them to HIP here, so the code
// using them stays the same. Only the names LIMA uses are mapped: add one here when the HIP build reports it missing
// Compiling a .cu file (nvcc or hipcc), so __device__ code and device intrinsics are available
#if defined(__CUDACC__) || defined(__HIPCC__)
#define LIMA_GPU_COMPILER
#endif

#if defined(LIMA_HIP)
#include <hip/hip_runtime.h>

// Types and constants
#define cudaError_t hipError_t
#define cudaSuccess hipSuccess
#define cudaErrorNotReady hipErrorNotReady
#define cudaStream_t hipStream_t
#define cudaEvent_t hipEvent_t
#define cudaDeviceProp hipDeviceProp_t
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaEventDisableTiming hipEventDisableTiming
#define cudaStreamNonBlocking hipStreamNonBlocking
#define cudaHostAllocMapped hipHostAllocMapped

// Errors
#define cudaGetLastError hipGetLastError
#define cudaPeekAtLastError hipPeekAtLastError
#define cudaGetErrorString hipGetErrorString

// Device
#define cudaSetDevice hipSetDevice
#define cudaGetDevice hipGetDevice
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaDriverGetVersion hipDriverGetVersion
#define cudaRuntimeGetVersion hipRuntimeGetVersion
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaDeviceSetCacheConfig hipDeviceSetCacheConfig
#define cudaFuncCachePreferShared hipFuncCachePreferShared
// hipProfilerStart/Stop return hipErrorNotSupported on Windows, which then sticks as the last error and fails the next
// error check. Profilers attach without them, so they do nothing under HIP
inline hipError_t cudaProfilerStart() { return hipSuccess; }
inline hipError_t cudaProfilerStop() { return hipSuccess; }

// Memory
#define cudaMalloc hipMalloc
#define cudaMallocManaged hipMallocManaged
#define cudaFree hipFree
#define cudaMallocHost(ptr, size) hipHostMalloc(ptr, size, 0)
#define cudaHostAlloc hipHostAlloc
#define cudaFreeHost hipHostFree
#define cudaMemcpy hipMemcpy
#define cudaMemcpyAsync hipMemcpyAsync
#define cudaMemset hipMemset
#define cudaMemsetAsync hipMemsetAsync

// CUDA's has a typed overload, HIP's takes void** only
template <typename T>
inline hipError_t cudaHostGetDevicePointer(T** devicePtr, void* hostPtr, unsigned int flags) {
	return hipHostGetDevicePointer(reinterpret_cast<void**>(devicePtr), hostPtr, flags);
}

// Streams and events
#define cudaStreamCreate hipStreamCreate
#define cudaStreamCreateWithFlags hipStreamCreateWithFlags
#define cudaStreamDestroy hipStreamDestroy
#define cudaStreamSynchronize hipStreamSynchronize
#define cudaStreamQuery hipStreamQuery
#define cudaStreamWaitEvent hipStreamWaitEvent
#define cudaEventCreateWithFlags hipEventCreateWithFlags
#define cudaEventDestroy hipEventDestroy
#define cudaEventRecord hipEventRecord
#define cudaEventSynchronize hipEventSynchronize
#define cudaEventQuery hipEventQuery

// Graphs: HIP builds launch steps directly (Engine.cu), but Engine keeps the members
#define cudaGraphExec_t hipGraphExec_t
#define cudaGraphExecDestroy hipGraphExecDestroy

// Graphs. Engine captures steps into graphs on CUDA only, HIP just needs the type
#define cudaGraphExec_t hipGraphExec_t
#define cudaGraphExecDestroy hipGraphExecDestroy

#else
#include <cuda_runtime.h>
#include <device_launch_parameters.h>	// threadIdx etc. for host compilers, which also parse __device__ functions
#include <cuda_profiler_api.h>
#endif
