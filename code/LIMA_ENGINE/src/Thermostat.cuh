#pragma once

#include "LimaTypes.cuh"
#include "PhysicsUtilsDevice.cuh"
#include "BatchLayout.cuh"
#include "EngineBodies.cuh"
#include "VRescale.h"
#include "Utilities.h"

#include "CubWrappers.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

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

	float* KineticEnergyValues() { return intermediate; }

	// Set once an MD integration kernel is queued that writes each particle's kinetic energy after its step to KineticEnergyValues,
	// so the device thermostat's measurements need not compute them from the slots. A repacked batch gets a new Thermostat
	bool integrationWritesValues = false;

	~Thermostat() {
		cudaFree(intermediate);
	}
};

namespace _Thermostat {
	// As ComputeKineticEnergyKernel, from the integration states in the supercluster slots, which StoreIntegrationStates copies to the
	// pclusters' vels_prev. Each particle's value goes to its pcluster lane, so the sums are the same. Unused lanes stay 0
	__global__ void ComputeKineticEnergyFromSlotsKernel(const ParticleIntegrationState* states, int nSlots, const PersistentClusterMeta* metadata,
		float* intermediate) {
		const int slot = blockIdx.x * blockDim.x + threadIdx.x;
		if (slot >= nSlots) return;
		const ParticleIntegrationState state = states[slot];
		if (state.mass == 0.f) return;
		const int pc = state.pclusterParticle / PersistentCluster::maxParticles;
		const int lane = state.pclusterParticle % PersistentCluster::maxParticles;
		intermediate[state.pclusterParticle] = PhysicsUtils::calcKineticEnergy(state.velocity.len(), metadata[pc].mass[lane]);
	}

	struct DeviceThermostatSimulation {
		int64_t degreesOfFreedom = 0;
		double referenceTemperature = 0.0;	// [K]
		double tauInCouplings = 0.0;
		bool applyThermostat = false;
		bool measure = false;	// Active and with particles
	};

	// One thread per simulation: the temperature after the previous step and, with the thermostat, the v-rescale factor for this step, as
	// Engine::hostMaster did between them. step is this step, or -1 for the step counter's. onlySimulation, if not -1, measures only that
	// simulation, and without coupling: its final temperature, after its last step
	__global__ void VRescaleKernel(const DeviceThermostatSimulation* simulations, int nSimulations, const double* kineticEnergies,
		VRescale* rngs, ThermostatImpulse* impulses, float* temperatureRing, int ringSize, const int64_t* stepCounter, int64_t step,
		int interval, int onlySimulation) {
		const int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= nSimulations || !simulations[id].measure || (onlySimulation != -1 && id != onlySimulation)) return;
		const DeviceThermostatSimulation simulation = simulations[id];
		if (step == -1)
			step = *stepCounter;
		const float temperature = PhysicsUtils::kineticEnergyToTemperature(kineticEnergies[id], simulation.degreesOfFreedom);
		temperatureRing[size_t(id) * ringSize + (step / interval) % ringSize] = temperature;
		if (simulation.applyThermostat && onlySimulation == -1) {
			const double factor = rngs[id].Factor(temperature, simulation.referenceTemperature, simulation.degreesOfFreedom, simulation.tauInCouplings);
			impulses[id] = ThermostatImpulse{ static_cast<float>(factor), step };
		}
	}
} // namespace _Thermostat

// The temperature measurement and v-rescale thermostat on the GPU, inside the captured step, so a coupling step needs no host roundtrip
// and does not invalidate the step graph. A coupling step, the one after a measurement interval, starts with the kinetic energies after
// the previous step (written by its integration, or from the slots after a batch is (re)packed), their sum per simulation and
// VRescaleKernel, on a side stream overlapping the force kernels, which
// do not touch the velocities. The factor goes to this step's integration as a ThermostatImpulse, and the temperature to a ring in
// mapped host memory, which the host reads without waiting (Harvest). Each simulation's random stream is the one the host
// thermostat would draw. Not with the barostat, whose pressure measurement reads the kinetic energy on the host anyway
class DeviceThermostat {
	CudaBuffer<_Thermostat::DeviceThermostatSimulation> simulationsDevice;
	CudaBuffer<VRescale> rngs;
	CudaBuffer<double> kineticEnergies;
	CudaBuffer<char> sumTempStorage;
	std::vector<BatchRange> pclusters;
	int interval = 1;
	int64_t firstStep = 0;
	float* ringHost = nullptr;
	float* ringDevice = nullptr;
	std::vector<int64_t> nextMeasurement;	// Per simulation, the index (step / interval) of the next temperature to harvest

public:
	// Holds the temperatures not yet harvested. The host harvests at least every Engine::maxStepsPerCall steps
	static constexpr int ringSize = 4096;
	CudaBuffer<ThermostatImpulse> impulses;	// Always allocated, as the integration kernels read it

	// step is the engine's first: the later steps that are a multiple of interval measure
	DeviceThermostat(const std::vector<VRescale>& initialRngs, int interval, int64_t step) : interval(interval), firstStep(step) {
		const size_t n = initialRngs.size();
		rngs.SetData(initialRngs);
		impulses.SetData(std::vector<ThermostatImpulse>(n));
		kineticEnergies.Expand(n);
		LIMA_UTILS::genericErrorCheck(cudaHostAlloc(&ringHost, sizeof(float) * n * ringSize, cudaHostAllocMapped));
		LIMA_UTILS::genericErrorCheck(cudaHostGetDevicePointer(&ringDevice, ringHost, 0));
		std::fill(ringHost, ringHost + n * ringSize, NAN);
		nextMeasurement.assign(n, step / interval + 1);
	}
	DeviceThermostat(const DeviceThermostat&) = delete;
	~DeviceThermostat() { cudaFreeHost(ringHost); }

	// Must be called before the first coupling step and whenever the batch is repacked. Synchronizes
	void SetSimulations(const std::vector<_Thermostat::DeviceThermostatSimulation>& simulations, const std::vector<BatchRange>& simulationPclusters) {
		simulationsDevice.SetData(simulations);
		pclusters = simulationPclusters;
		size_t tempBytes = 0;
		for (const BatchRange range : pclusters)
			if (range.count > 0)
				tempBytes = std::max(tempBytes, CubWrappers::SumTempBytes(range.count * PersistentCluster::maxParticles));
		sumTempStorage.Expand(std::max<size_t>(tempBytes, 1));
	}

	bool MeasuresAt(int64_t step) const { return step > firstStep && step % interval == 0; }

	// Capturable. kineticEnergyValues is the batch's Thermostat buffer, one value per pcluster lane. If valuesCurrent, the previous
	// step's integration wrote them (Thermostat::integrationWritesValues), else they are computed from the slots. See VRescaleKernel
	// for step and onlySimulation
	void Queue(const ParticleIntegrationState* states, int nSlots, const PersistentClusterMeta* metadata, float* kineticEnergyValues,
		bool valuesCurrent, const int64_t* stepCounter, cudaStream_t stream, int64_t step = -1, int onlySimulation = -1) {
		if (nSlots > 0 && !valuesCurrent)
			_Thermostat::ComputeKineticEnergyFromSlotsKernel<<<(nSlots + 127) / 128, 128, 0, stream>>>(states, nSlots, metadata, kineticEnergyValues);
		for (size_t id = 0; id < pclusters.size(); id++)
			if (pclusters[id].count > 0 && (onlySimulation == -1 || id == size_t(onlySimulation)))
				CubWrappers::SumInto(kineticEnergyValues + pclusters[id].offset * PersistentCluster::maxParticles,
					pclusters[id].count * PersistentCluster::maxParticles, kineticEnergies.Get() + id, sumTempStorage, stream);
		const int n = static_cast<int>(pclusters.size());
		_Thermostat::VRescaleKernel<<<(n + 31) / 32, 32, 0, stream>>>(simulationsDevice.Get(), n, kineticEnergies.Get(), rngs.Get(),
			impulses.Get(), ringDevice, ringSize, stepCounter, step, interval, onlySimulation);
	}

	// Appends the temperatures the GPU has written since the last harvest, in order, and returns the latest. Does not wait: a temperature
	// is either fully written or still NaN. After a synchronization it gets all of them
	std::optional<float> Harvest(size_t simulationId, std::vector<float>& temperatures) {
		volatile float* ring = ringHost + simulationId * ringSize;
		std::optional<float> latest;
		while (true) {
			const size_t slot = nextMeasurement[simulationId] % ringSize;
			const float temperature = ring[slot];
			if (std::isnan(temperature))
				return latest;
			ring[slot] = NAN;
			temperatures.push_back(temperature);
			latest = temperature;
			nextMeasurement[simulationId]++;
		}
	}
};

