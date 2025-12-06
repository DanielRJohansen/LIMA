#pragma once
#include "Bodies.cuh"
#include "Simulation.cuh"
#include "CompoundBuilder.h"

#include <vector>

namespace BoxBuilder
{
	std::unique_ptr<Box> BuildBox(const SimParams& params, BoxImage& boxImage);

	// This function expects all ptr's of simulation->box to be pre-allocated on host
	void copyBoxState(Simulation& simulation, std::unique_ptr<Box> boxsrc, uint32_t boxsrc_current_step);

	bool verifyAllParticlesIsInsideBox(Simulation& sim, float padding = 0.f, bool verbose = true);
};