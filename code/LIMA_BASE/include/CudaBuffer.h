#pragma once
#include "GpuRuntime.h"
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

// Throws on a failed allocation, and clears the error so it isn't reported again by a later error check
inline void ThrowIfAllocationFailed(cudaError_t status, const char* what) {
	if (status == cudaSuccess)
		return;
	cudaGetLastError();
	throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(status));
}

template <typename T>
class CudaBuffer {
	T* devicePtr = nullptr;
	size_t size = 0;

public:
	CudaBuffer() {}
	CudaBuffer(const CudaBuffer& other) = delete;
	~CudaBuffer() {
		if (devicePtr)
			cudaFree(devicePtr);
	}
	T* Get() const {
		return devicePtr;
	}
	void Expand(size_t requiredSize, std::optional<double> margin=std::nullopt) {
		if (requiredSize <= size)
			return;

		size_t newSize = static_cast<size_t>((double)requiredSize * margin.value_or(1.));		
		if (devicePtr)
			cudaFree(devicePtr);
		devicePtr = nullptr;
		size = 0;
		ThrowIfAllocationFailed(cudaMalloc(&devicePtr, newSize * sizeof(T)), "CudaBuffer allocation failed");
		size = newSize;
	}
	void SetData(const std::vector<T>& v){
		Expand(v.size());
		cudaMemcpy(devicePtr, v.data(), v.size() * sizeof(T), cudaMemcpyHostToDevice);
	}
	std::vector<T> GetData() const { 
		std::vector<T> hostData(size);
		cudaMemcpy(hostData.data(), devicePtr, sizeof(T) * size, cudaMemcpyDeviceToHost);
		return hostData;
	}
	void GetData(std::vector<T>& dst) const {
		dst.resize(size);
		cudaMemcpy(dst.data(), devicePtr, sizeof(T) * size, cudaMemcpyDeviceToHost);
	}
	size_t Size() const {
		return size;
	}
};

// Pinned host ints that kernels read and write directly (zero-copy), for the few values the host needs mid-step. Unlike
// a cudaMemcpy to or from the host, this needs no copy engine, whose round trip costs ~60 us under WDDM.
// The host reads the values after Wait, which waits for the kernels queued before Written. Work queued after
// Written keeps the GPU busy while the host wakes up and queues the work that depends on the values
class MappedHostBuffer {
	int* host = nullptr;
	int* device = nullptr;
	size_t size = 0;
	cudaEvent_t written = nullptr;

public:
	MappedHostBuffer() {}
	MappedHostBuffer(const MappedHostBuffer&) = delete;
	~MappedHostBuffer() {
		cudaFreeHost(host);
		if (written) cudaEventDestroy(written);
	}
	// Discards the contents when growing
	void Reserve(size_t requiredSize) {
		if (requiredSize <= size)
			return;
		cudaFreeHost(host);
		host = nullptr;
		device = nullptr;
		size = 0;
		ThrowIfAllocationFailed(cudaHostAlloc(&host, requiredSize * sizeof(int), cudaHostAllocMapped), "MappedHostBuffer allocation failed");
		cudaHostGetDevicePointer(&device, host, 0);
		size = requiredSize;
		if (!written) cudaEventCreateWithFlags(&written, cudaEventDisableTiming);
	}
	int* Host() const { return host; }
	int* Device() const { return device; }
	// Marks the kernels writing the values as queued
	void Written(cudaStream_t stream) const { cudaEventRecord(written, stream); }
	// Submits the work queued since Written, which WDDM may otherwise hold back, then waits for the values
	void Wait(cudaStream_t stream) const {
		cudaStreamQuery(stream);
		cudaEventSynchronize(written);
	}
};
