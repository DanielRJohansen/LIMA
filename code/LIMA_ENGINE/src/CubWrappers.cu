#include "CubWrappers.h"

#include <cub/device/device_reduce.cuh>
#include <cub/device/device_scan.cuh>
#include <cuda/std/functional>

#include <stdexcept>
#include <string>

namespace {
	void Check(cudaError_t error, const char* what) {
		if (error != cudaSuccess)
			throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(error));
	}
	void SyncAndCheck(cudaStream_t stream, const char* what) {
		Check(cudaGetLastError(), what);
		Check(cudaStreamSynchronize(stream), what);
	}

	__global__ void FillKernel(int* data, size_t count, int value) {
		for (size_t i = blockIdx.x * size_t(blockDim.x) + threadIdx.x; i < count; i += size_t(gridDim.x) * blockDim.x)
			data[i] = value;
	}
	__global__ void SqrtKernel(float* data, size_t count) {
		for (size_t i = blockIdx.x * size_t(blockDim.x) + threadIdx.x; i < count; i += size_t(gridDim.x) * blockDim.x)
			data[i] = sqrtf(data[i]);
	}
	unsigned int GridSize(size_t count, int blockSize) {
		const size_t nBlocks = (count + blockSize - 1) / blockSize;
		return static_cast<unsigned int>(nBlocks < 65535 ? nBlocks : 65535);
	}
}

void CubWrappers::ExclusiveScan(const int* first, const int* last, int* result, cudaStream_t stream) {
	const auto count = static_cast<int>(last - first);
	if (count <= 0) return;
	size_t tempBytes = 0;
	Check(cub::DeviceScan::ExclusiveSum(nullptr, tempBytes, first, result, count, stream), "ExclusiveScan");
	void* temp = nullptr;
	Check(cudaMallocAsync(&temp, tempBytes, stream), "ExclusiveScan");
	Check(cub::DeviceScan::ExclusiveSum(temp, tempBytes, first, result, count, stream), "ExclusiveScan");
	Check(cudaFreeAsync(temp, stream), "ExclusiveScan");
	SyncAndCheck(stream, "ExclusiveScan");
}

void CubWrappers::FillN(int* first, size_t count, int value, cudaStream_t stream) {
	if (count == 0) return;
	FillKernel<<<GridSize(count, 256), 256, 0, stream>>>(first, count, value);
	SyncAndCheck(stream, "FillN");
}

double CubWrappers::Sum(const float* first, const float* last, cudaStream_t stream) {
	const auto count = static_cast<int>(last - first);
	if (count <= 0) return 0.0;
	size_t tempBytes = 0;
	Check(cub::DeviceReduce::Reduce(nullptr, tempBytes, first, static_cast<double*>(nullptr), count,
		cuda::std::plus<double>{}, 0.0, stream), "Sum");
	// The result is stored after the temp storage, so there is one allocation
	const size_t resultOffset = (tempBytes + alignof(double) - 1) / alignof(double) * alignof(double);
	void* temp = nullptr;
	Check(cudaMallocAsync(&temp, resultOffset + sizeof(double), stream), "Sum");
	double* resultDevice = reinterpret_cast<double*>(static_cast<char*>(temp) + resultOffset);
	Check(cub::DeviceReduce::Reduce(temp, tempBytes, first, resultDevice, count,
		cuda::std::plus<double>{}, 0.0, stream), "Sum");
	double result = 0.0;
	Check(cudaMemcpyAsync(&result, resultDevice, sizeof(double), cudaMemcpyDeviceToHost, stream), "Sum");
	Check(cudaFreeAsync(temp, stream), "Sum");
	SyncAndCheck(stream, "Sum");
	return result;
}

void CubWrappers::SqrtInPlace(float* data, size_t count) {
	if (count == 0) return;
	SqrtKernel<<<GridSize(count, 256), 256>>>(data, count);
	SyncAndCheck(0, "SqrtInPlace");
}
