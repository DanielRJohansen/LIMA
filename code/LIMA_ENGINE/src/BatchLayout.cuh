#pragma once

#include "Bodies.cuh"

// Offsets refer to the batch allocations. Particle slots (pc * 4 + lane) and
// dense particle IDs deliberately use different ranges.
struct BatchRange {
	int offset = 0;
	int count = 0;
};

struct SimulationDeviceData {
	BatchRange particles;
	BatchRange pclusters;
	BatchRange bondgroups;
	BatchRange gridnodes;
	size_t logOffset = 0;
	float dt = 0.f;
	float thermostatScalar = 1.f;
	bool active = true;
	UniformElectricField uniformElectricField;
};

// Compact per-simulation data consumed by the batch-wide integration kernel.
// Device copy of per-simulation data needed by the integration kernel. Built from EngineSimulationData::device,
// so after changing any of these fields there, call Engine::UploadIntegrationSimulationData()
struct IntegrationSimulationData {
	BatchRange particles;
	BatchRange pclusters;
	size_t logOffset = 0;
	float dt = 0.f;
	float thermostatScalar = 1.f;
};

// A v-rescale factor the device thermostat computed after a coupling step, see DeviceThermostat. The integration kernels apply it to
// the velocities of step, in place of IntegrationSimulationData::thermostatScalar, so no step resets it
struct ThermostatImpulse {
	float scalar = 1.f;
	int64_t step = -1;
};
