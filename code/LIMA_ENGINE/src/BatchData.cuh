#pragma once

#include "SimulationData.h"

#include <vector>

// Host-only, defined in EngineBatch.cpp
namespace EngineBatch {
	void Validate(const std::vector<Simulation*>& simulations);
	// Packs the simulations into one batch. With active, the inactive simulations are kept as retired members
	void Pack(EngineBatchData& batch, const std::vector<Simulation*>& simulations, const std::vector<bool>* active = nullptr);
}
