#include "ThrustWrappers.h"

#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/transform.h>

namespace {
	struct SqrtFloat {
		__device__ float operator()(float x) const { return sqrtf(x); }
	};
}

void ThrustWrappers::ExclusiveScan(const int* first, const int* last, int* result, cudaStream_t stream) {
	thrust::exclusive_scan(thrust::cuda::par.on(stream), first, last, result);
}

void ThrustWrappers::FillN(int* first, size_t count, int value, cudaStream_t stream) {
	thrust::fill_n(thrust::cuda::par.on(stream), first, count, value);
}

double ThrustWrappers::Sum(const float* first, const float* last, cudaStream_t stream) {
	return thrust::reduce(thrust::cuda::par.on(stream), first, last, 0.0);
}

void ThrustWrappers::SqrtInPlace(float* data, size_t count) {
	thrust::device_ptr<float> begin(data);
	thrust::transform(thrust::device, begin, begin + count, begin, SqrtFloat{});
}
