#pragma once

#include "LimaTypes.cuh"
#include "PhysicsUtilsDevice.cuh"
#include "BatchLayout.cuh"

#include "CubWrappers.h"

#include <cuda_runtime.h>

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
	// [J/mol]. The same reduction feeds temperature and pressure when their coupling steps coincide.
	double KineticEnergy(BatchRange pclusters, cudaStream_t stream) {
		const float* begin = intermediate + pclusters.offset * PersistentCluster::maxParticles;
		return CubWrappers::Sum(begin, begin + pclusters.count * PersistentCluster::maxParticles, stream);
	}

	static float Temperature(const BoxParams& boxparams, double kineticEnergy) {
		return PhysicsUtils::kineticEnergyToTemperature(kineticEnergy, boxparams.degreesOfFreedom);
	}

	const float* KineticEnergyValues(BatchRange pclusters) const {
		return intermediate + pclusters.offset * PersistentCluster::maxParticles;
	}

	~Thermostat() {
		cudaFree(intermediate);
	}
};

