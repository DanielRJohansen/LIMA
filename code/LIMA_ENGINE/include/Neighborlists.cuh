#pragma once

#include "LIMA_BASE/include/LimaTypes.cuh"
#include "LIMA_BASE/include/Simulation.cuh"
#include "EngineUtils.cuh"
#include <chrono>
#include <thread>
#include <mutex>
#include <functional>


struct NListDataCollection {
	NListDataCollection(Simulation* simulation);

	void preparePositionData(const Simulation& simulation, uint32_t step_at_update);


	int n_compounds = -1;

	// I guess this is not critical but temp, needed to load pos device->host
	CompoundState* compoundstates = nullptr;
	//Solvent* solvents;

	Float3 compound_key_positions[MAX_COMPOUNDS];	// [nm] absolute position
	//NodeIndex compound_origos[MAX_COMPOUNDS];		// compound's corresponding gridnode

	// These are loaded before simulaiton start. Kept on host, and copied to device each update.
	NeighborList* compound_neighborlists = nullptr;
	CompoundGrid* compoundgrid = nullptr;
};

class NListManager {
public:
	NListManager(){}	// This is messy, necessary to 
	NListManager(Simulation* simulation);


	void handleNLISTS(Simulation* simulation, bool async, bool force_update, int* timing);

	void pushNlistsToDevice(Simulation* simulation);


	int stepsSinceUpdate(uint64_t currentSimStep) const { 
		return static_cast<int>(currentSimStep - prev_update_step); 
	}

private:

	std::mutex m_mutex{};

	volatile bool updated_neighborlists_ready = 0;
	NListDataCollection* nlist_data_collection = nullptr;

	uint64_t prev_update_step = 0;
};
