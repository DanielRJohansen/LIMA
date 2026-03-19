#include <cuda_runtime.h>
#include <optional>

template <typename T>
class CudaBuffer {
	T* devicePtr = nullptr;
	size_t size = 0;

public:
	CudaBuffer() {}
	~CudaBuffer() {
		if (devicePtr)
			cudaFree(devicePtr);
	}
	T* Get() const {
		return devicePtr;
	}
	void Expand(size_t requiredSize, std::optional<double> margin) {
		if (requiredSize <= size)
			return;

		size_t newSize = static_cast<size_t>((double)requiredSize * margin.value_or(1.));		
		if (devicePtr)
			cudaFree(devicePtr);		
		cudaMalloc(&devicePtr, newSize * sizeof(T));
		size = newSize;
	}
};
