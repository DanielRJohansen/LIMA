#pragma once

#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

// Kernels with fixed-capacity buffers skip writes past capacity and report the first overflow here,
// instead of corrupting neighboring data. The status is sticky: the host reads it at an existing sync
// point and aborts the run, since dropped entries mean the step's physics is already wrong.
struct CapacityOverflow {
	enum Code : int { None = 0, ClusterTransfer, ClusterOccupancy, ChargeBlock };
	static constexpr int nValues = 4; // code, block, count, capacity

	int* status = nullptr;

	__host__ static CapacityOverflow Create() {
		CapacityOverflow overflow;
		cudaMalloc(&overflow.status, sizeof(int) * nValues);
		cudaMemset(overflow.status, 0, sizeof(int) * nValues);
		return overflow;
	}
	__host__ void Free() { cudaFree(status); status = nullptr; }

#ifdef __CUDACC__
	__device__ void Report(Code code, int block, int count, int capacity) const {
		if (atomicCAS(status, static_cast<int>(None), static_cast<int>(code)) == None) {
			status[1] = block;
			status[2] = count;
			status[3] = capacity;
		}
	}
#endif

	// Blocking read, for use after the stream has been synchronized anyway
	__host__ void Check() const {
		int values[nValues]{};
		cudaMemcpy(values, status, sizeof(values), cudaMemcpyDeviceToHost);
		ThrowIfSet(values);
	}

	__host__ static void ThrowIfSet(const int values[nValues]) {
		if (values[0] == None) return;
		const char* what = values[0] == ClusterTransfer ? "outgoing cluster transfers per direction"
			: values[0] == ClusterOccupancy ? "persistent clusters per grid block"
			: "PME charge entries per charge block";
		throw std::length_error("Engine capacity exceeded: " + std::to_string(values[2]) + " " + what + " in block "
			+ std::to_string(values[1]) + " (max " + std::to_string(values[3]) + "). The system is locally too dense for the engine's fixed buffers");
	}
};
