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
	bool hasFixedMovement = false;
	bool hasFixedRotation = false;
	bool hasForceMask = false;
	bool hasElasticPositions = false;
	UniformElectricField uniformElectricField;
};
