#pragma once

// The FFT library: cuFFT, or hipFFT under HIP with the cuFFT names PME uses mapped to it
#if defined(LIMA_HIP)
#include <hipfft/hipfft.h>

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
#else
#include <cufft.h>
#endif
