#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
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


struct NListDataCollection {
	NListDataCollection(Simulation* simulation);

	void preparePositionData(const Simulation& simulation, uint32_t step_at_update);

	Float3 compound_key_positions[MAX_COMPOUNDS];	// [nm] absolute position
	//NodeIndex compound_origos[MAX_COMPOUNDS];		// compound's corresponding gridnode

	// These are loaded before simulaiton start. Kept on host, and copied to device each update.
	std::vector<NeighborList> compound_neighborlists;
	std::unique_ptr<CompoundGrid> compoundgrid;
};

class NListManager {
public:
	NListManager(Simulation* simulation);


	void handleNLISTS(Simulation* simulation, bool async, bool force_update, int* timing);
	void handleNlistGPU(Simulation* simulation, int& timing);


	void pushNlistsToDevice(Simulation* simulation);


	int stepsSinceUpdate(uint64_t currentSimStep) const { 
		return static_cast<int>(currentSimStep - prev_update_step); 
	}

private:

	std::mutex m_mutex{};

	volatile bool updated_neighborlists_ready = 0;
	std::unique_ptr<NListDataCollection> nlist_data_collection;

	uint64_t prev_update_step = 0;
};
