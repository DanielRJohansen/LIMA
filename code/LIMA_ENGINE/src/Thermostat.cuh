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
		const CompoundInterimState* const states;
		const Compound* const compounds;

		__host__ __device__
			TotalKineticEnergyCompounds(const CompoundInterimState* const _states, const Compound* const compounds)
			: states(_states), compounds(compounds){}

		__host__ __device__
			float operator()(int idx) const {
			int compoundIdx = idx / MAX_COMPOUND_PARTICLES;
			int particleIdx = idx % MAX_COMPOUND_PARTICLES;
			const float mass = compounds[compoundIdx].atomMasses[particleIdx];

			const Float3& velocity = states[compoundIdx].vels_prev[particleIdx];
			return PhysicsUtils::calcKineticEnergy(velocity.len(), mass); // TODO: calcKineticEnergy can use lenSquared instead, save a sqrtf!!		
		}
	};

	struct TotalKineticEnergySolvents {
		const TinyMol* const states;

		__host__ __device__
			TotalKineticEnergySolvents(const TinyMol* const _states)
			: states(_states){}

		__host__ __device__
			float operator()(int idx) const {
			const float mass = SOLVENT_MASS;

			const Float3& velocity = states[idx].vel_prev;
			return PhysicsUtils::calcKineticEnergy(velocity.len(), mass); // TODO: calcKineticEnergy can use lenSquared instead, save a sqrtf!!
		}
	};

	__global__ void ComputeTemperatureScalarKernel(float temperature, SimulationDevice* simDev) {

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
	float* intermediate;
	int totalParticlesUpperbound;
	int nCompounds;
	int nSolvents;

public:
	Thermostat(int nCompounds, int nSolvents, int totalParticlesUpperbound) : 
		totalParticlesUpperbound(totalParticlesUpperbound),
		nCompounds(nCompounds),
		nSolvents(nSolvents)
	{
		cudaMalloc(&intermediate, sizeof(float) * totalParticlesUpperbound); // totalParticlesUpperbound = MAX_COMPOUND_PARTICLES * nCompounds + nSolvents
	}

	// {temp,thermostatScalar}
	std::pair<float, float> Temperature(SimulationDevice* simDev, const BoxParams& boxparams, const SimParams& simparams) {
		// Step 1: Calculate kinetic energy for each compound particle and store in the intermediate buffer
		thrust::transform(thrust::device, thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(nCompounds * MAX_COMPOUND_PARTICLES),
			intermediate, _Thermostat::TotalKineticEnergyCompounds(simDev->boxState->compoundsInterimState, simDev->boxConfig.compounds));

		// Step 2: Calculate kinetic energy for each solvent particle and store in the next segment of the intermediate buffer
		thrust::transform(thrust::device, thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(nSolvents),
			intermediate + (nCompounds * MAX_COMPOUND_PARTICLES), _Thermostat::TotalKineticEnergySolvents(simDev->boxState->tinyMols));

		// Step 3: Sum up all kinetic energy values (compounds + solvents)
		double totalKineticEnergy = thrust::reduce(thrust::device, intermediate, intermediate + totalParticlesUpperbound, 0.0);

		//printf("Total kinetic energy: %f\n", totalKineticEnergy); 
		const float temperature = PhysicsUtils::kineticEnergyToTemperature(totalKineticEnergy, boxparams.degreesOfFreedom);
		const float scalar = _Thermostat::ComputeThermostatScalar(temperature, simparams);
		return { temperature, scalar };
	}

	~Thermostat() {
		cudaFree(intermediate);
	}
};

