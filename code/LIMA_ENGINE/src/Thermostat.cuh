#pragma once

#include "LimaTypes.cuh"
#include "PhysicsUtilsDevice.cuh"
#include "KernelConstants.cuh"
#include "SimulationDevice.cuh"

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

#include <cuda_runtime.h>

namespace _Thermostat {

	struct TotalKineticEnergyCompounds {
		const PersistentclusterInterimState* const states;
		const PersistentClusterMeta* const pcMeta;

		__host__ __device__
			TotalKineticEnergyCompounds(const PersistentclusterInterimState* const _states, const PersistentClusterMeta* const _pcMeta)
			: states(_states), pcMeta(_pcMeta){}
		__host__ __device__
			float operator()(int idx) const {
			int pcId = idx / PersistentCluster::maxParticles;
			int pId = idx % PersistentCluster::maxParticles;
			const float mass = pcMeta[pcId].mass[pId];

			const Float3& velocity = states[pcId].vels_prev[pId];
			return PhysicsUtils::calcKineticEnergy(velocity.len(), mass); // TODO OPTIM: calcKineticEnergy can use lenSquared instead, save a sqrtf!!		
		}
	};


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

	// {temp,thermostatScalar}
	std::pair<float, float> Temperature(SimulationDevice* simDev, const BoxParams& boxparams, const SimParams& simparams, int step, const PersistentClusterMeta* const pcMetaDevice) {
		// Step 1: Calculate kinetic energy for each Pcluster and store in the intermediate buffer
		thrust::transform(thrust::device, thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(nPclusters * PersistentCluster::maxParticles),
			intermediate, _Thermostat::TotalKineticEnergyCompounds(simDev->boxState.pclusterInterimStates, pcMetaDevice));
		LIMA_UTILS::genericErrorCheckNoSync("TotalKineticEnergyCompounds");
		cudaDeviceSynchronize();

		// Step 3: Sum up all kinetic energy values (compounds + solvents)
		double totalKineticEnergy = thrust::reduce(thrust::device, intermediate, intermediate + nPclusters * PersistentCluster::maxParticles, 0.0);

		//printf("Total kinetic energy: %f\n", totalKineticEnergy); 
		const float temperature = PhysicsUtils::kineticEnergyToTemperature(totalKineticEnergy, boxparams.degreesOfFreedom);
		const float scalar = _Thermostat::ComputeThermostatScalar(temperature, simparams);
		return { temperature, scalar };		
	}

	~Thermostat() {
		//cudaFree(intermediate);
	}
};

