#include <iostream>

#include "EngineUtils.cuh"

#include <cuda_runtime.h>

__global__ void addArrays(int* a, int* b, int* c, int size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        c[idx] = a[idx] + b[idx];
    }
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
