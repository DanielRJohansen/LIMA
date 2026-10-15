#pragma once

// The spike's code uses the CUDA names. Under hipcc they are mapped to HIP here, the way GpuRuntime.h's HIP branch will.
// Only the names the spike (and LIMA's PME) use
#if defined(__HIP_PLATFORM_AMD__)
#include <hip/hip_runtime.h>
#include <hipfft/hipfft.h>

#define cudaError_t hipError_t
#define cudaSuccess hipSuccess
#define cudaGetErrorString hipGetErrorString
#define cudaGetLastError hipGetLastError
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaDeviceProp hipDeviceProp_t
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaDriverGetVersion hipDriverGetVersion
#define cudaRuntimeGetVersion hipRuntimeGetVersion
#define cudaMalloc hipMalloc
#define cudaFree hipFree
#define cudaMemcpy hipMemcpy
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaStream_t hipStream_t
#define cudaStreamCreate hipStreamCreate
#define cudaStreamDestroy hipStreamDestroy
#define cudaStreamSynchronize hipStreamSynchronize

#define cufftHandle hipfftHandle
#define cufftResult hipfftResult
#define cufftComplex hipfftComplex
#define cufftReal hipfftReal
#define CUFFT_SUCCESS HIPFFT_SUCCESS
#define CUFFT_R2C HIPFFT_R2C
#define CUFFT_C2R HIPFFT_C2R
#define cufftPlanMany hipfftPlanMany
#define cufftSetStream hipfftSetStream
#define cufftExecR2C hipfftExecR2C
#define cufftExecC2R hipfftExecC2R
#define cufftDestroy hipfftDestroy

#define LIMA_GPU_VENDOR "HIP"
#else
#include <cuda_runtime.h>
#include <cufft.h>

#define LIMA_GPU_VENDOR "CUDA"
#endif
