#pragma once

// The GPU runtime, included everywhere instead of cuda_runtime.h. All projects use the CUDA runtime names (cudaMalloc,
// cudaStream_t, ...); a HIP build will map them to HIP here, so the code using them stays the same
#if defined(LIMA_HIP)
#error "The HIP backend is not implemented yet"
#else
#include <cuda_runtime.h>
#include <device_launch_parameters.h>	// threadIdx etc. for host compilers, which also parse __device__ functions
#include <cuda_profiler_api.h>
#endif
