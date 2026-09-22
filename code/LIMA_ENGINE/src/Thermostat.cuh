#pragma once

#include "LimaTypes.cuh"
#include "PhysicsUtilsDevice.cuh"
#include "BatchLayout.cuh"

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

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


	float ComputeThermostatScalar(float temperature, const SimParams& simparams) {
		const float target_temp = 310.f;  // Target temperature in [K]

		// Avoid division by zero
		const float temp_safe = (temperature == 0.f) ? 1.0f : temperature;

		// Compute the temperature scalar
		float temp_scalar = target_temp / temp_safe;

		// Clamp the temperature scalar to avoid rapid temperature changes
		const float max_scalar = 0.001f / static_cast<float>(simparams.steps_per_temperature_measurement);  // Change velocity by 0.1% over NSTEPS
		return std::clamp(temp_scalar, 1.f - max_scalar, 1.f + max_scalar);
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
	std::pair<float, float> Temperature(const BoxParams& boxparams, const SimParams& simparams, BatchRange pclusters, cudaStream_t stream) {
		const float* begin = intermediate + pclusters.offset * PersistentCluster::maxParticles;
		const double totalKineticEnergy = thrust::reduce(thrust::cuda::par.on(stream),
			begin, begin + pclusters.count * PersistentCluster::maxParticles, 0.0);
		const float temperature = PhysicsUtils::kineticEnergyToTemperature(totalKineticEnergy, boxparams.degreesOfFreedom);
		return {temperature, _Thermostat::ComputeThermostatScalar(temperature, simparams)};
	}

	~Thermostat() {
		cudaFree(intermediate);
	}
};

