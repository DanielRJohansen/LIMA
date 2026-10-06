#pragma once

#include "EngineBodies.cuh"
#include "Utilities.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>

// Shared helpers for fixtures defined alongside the production kernels they exercise.
// Each fixture runs in a disposable process, including the within-capacity controls.
namespace EngineLimitTesting {

	inline void Require(bool condition, const char* message) {
		if (!condition) throw std::runtime_error(message);
	}

	inline void CheckCuda() {
		LIMA_UTILS::genericErrorCheck(cudaDeviceSynchronize());
		LIMA_UTILS::genericErrorCheck(cudaGetLastError());
	}

	template<typename T>
	struct FreeDeviceMembers {
		void operator()(T* value) const { value->Free(); delete value; }
	};

	inline auto MakeTransferModule(int blocks) {
		return std::unique_ptr<PClusterTransfermodule, FreeDeviceMembers<PClusterTransfermodule>>(
			new PClusterTransfermodule(PClusterTransfermodule::Create(blocks)));
	}

	void ClusterTransfer(int count);
	void ClusterOccupancy(int count);
	void ChargeBlock(int count);
}
