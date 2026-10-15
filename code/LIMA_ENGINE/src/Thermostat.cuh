#pragma once

#include "LimaTypes.cuh"
#include "PhysicsUtilsDevice.cuh"
#include "BatchLayout.cuh"

#include "CubWrappers.h"

#include "GpuRuntime.h"

namespace _Thermostat {

	__global__ void ComputeKineticEnergyKernel(const PersistentclusterInterimState* states, const PersistentClusterMeta* metadata,
		int nPclusters, float* intermediate) {
		const int index = blockIdx.x * blockDim.x + threadIdx.x;
		if (index >= nPclusters * PersistentCluster::maxParticles) return;
		const int pc = index / PersistentCluster::maxParticles;
		const int lane = index % PersistentCluster::maxParticles;
		intermediate[pc * PersistentCluster::maxParticles + lane] = PhysicsUtils::calcKineticEnergy(
			states[pc].vels_prev[lane].len(), metadata[pc].mass[lane]);
	}


} // namespace Thermostat

class Thermostat {
	float* intermediate = nullptr;
	int nPclusters = 0;

public:
	Thermostat(int nPclusters)
		: nPclusters(nPclusters)
	{
		cudaMalloc(&intermediate, sizeof(float) * nPclusters * PersistentCluster::maxParticles);
		cudaMemset(intermediate, 0, sizeof(float) * nPclusters * PersistentCluster::maxParticles);
	}

	void ComputeKineticEnergy(const PersistentclusterInterimState* states, const PersistentClusterMeta* metadata,
		int count, cudaStream_t stream) {
		if (count == 0) return;
		_Thermostat::ComputeKineticEnergyKernel<<<(count * PersistentCluster::maxParticles + 127) / 128, 128, 0, stream>>>(
			states, metadata, count, intermediate);
	}

	// Keep each simulation's reduction order and degrees of freedom independent.
	// [K]. The thermostat's factor is computed from it by VRescale
	float Temperature(const BoxParams& boxparams, BatchRange pclusters, cudaStream_t stream) {
		const float* begin = intermediate + pclusters.offset * PersistentCluster::maxParticles;
		const double totalKineticEnergy = CubWrappers::Sum(begin, begin + pclusters.count * PersistentCluster::maxParticles, stream);
		return PhysicsUtils::kineticEnergyToTemperature(totalKineticEnergy, boxparams.degreesOfFreedom);
	}

	~Thermostat() {
		cudaFree(intermediate);
	}
};

