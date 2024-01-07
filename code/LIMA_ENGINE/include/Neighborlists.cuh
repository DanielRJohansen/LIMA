#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "SimulationDevice.cuh"
#include "EngineUtils.cuh"

#include <chrono>
#include <thread>
#include <mutex>
#include <functional>

namespace NListUtils {
	__device__  __host__ static bool neighborWithinCutoff(const Float3* pos_a, const Float3* pos_b, const float cutoff_nm) {		// This is used for compounds with a confining_particle_sphere from key_particle BEFORE CUTOFF begins
		const float dist = EngineUtils::calcHyperDistNM(pos_a, pos_b);
		return dist < cutoff_nm;
	}
};



namespace NeighborLists {

	void updateNlists(SimulationDevice<PeriodicBoundaryCondition>* sim_dev, int n_compounds, int& timing);


};
