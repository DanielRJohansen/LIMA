#include <cuda_runtime.h>
#include <optional>

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
		cudaMalloc(&devicePtr, newSize * sizeof(T));
		size = newSize;
	}
	void SetData(const std::vector<T>& v){
		Expand(v.size());
		cudaMemcpy(devicePtr, v.data(), v.size() * sizeof(T), cudaMemcpyHostToDevice);
	}
	std::vector<T> GetData() {
		std::vector<T> hostData(size);
		cudaMemcpy(hostData.data(), devicePtr, sizeof(T) * size, cudaMemcpyDeviceToHost);
		return hostData;
	}
	size_t Size() const {
		return size;
	}
};
