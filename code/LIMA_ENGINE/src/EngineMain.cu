#include <iostream>

#include <cuda_runtime.h>

__global__ void addArrays(int* a, int* b, int* c, int size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        c[idx] = a[idx] + b[idx];
    }
}


template<typename T>
T* genericMoveToDevice(T* data_ptr, int n_elements) {	// Currently uses MallocManaged, switch to unmanaged for safer operation
    if (n_elements == 0) { return nullptr; }

    T* gpu_ptr = nullptr;
    size_t bytesize = n_elements * sizeof(T);

    cudaMallocManaged(&gpu_ptr, bytesize);

    auto cuda_status = cudaMemcpy(gpu_ptr, data_ptr, bytesize, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess) {
        std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
        throw std::runtime_error("Move to device failed");
    }

    if (n_elements == 1)
        delete data_ptr;
    else
        delete[] data_ptr;

    return gpu_ptr;
}

int main() {
    const int size = 256;
    int* a = new int[size];
    int* b = new int[size];
    int* c = new int[size];

    for (int i = 0; i < size; i++) {
        a[i] = i;
        b[i] = i;
        c[i] = 0;
    }

    genericMoveToDevice(a, size);
    genericMoveToDevice(b, size);
    genericMoveToDevice(c, size);

    addArrays << <256 / 32, 32 >> > (a, b, c, size);

    for (int i = 0; i < size; i++) {
        if (c[i] != i * 2) {
            std::printf("Failed self-test\n");
            return 1;
        }
    }

    cudaFree(a);
    cudaFree(b);
    cudaFree(c);
    std::printf("\tLIMA_ENGINE self-test success\n");

    return 0;
}
