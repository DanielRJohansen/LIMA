#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "EngineUtils.cuh"
#include <chrono>
#include <thread>

struct NListDataCollection {
	NListDataCollection(Simulation* simulation);

	void preparePositionData(const Simulation& simulation, uint32_t step_at_update);


	int n_compounds;
	int n_solvents;

	// I guess this is not critical but temp, needed to load pos device->host
	CompoundState* compoundstates;
	//Solvent* solvents;

	Float3 compound_key_positions[MAX_COMPOUNDS];
	Float3 solvent_positions[MAX_SOLVENTS];

	// These are loaded before simulaiton start. Kept on host, and copied to device each update.
	NeighborList* compound_neighborlists;
	NeighborList* solvent_neighborlists;
};


namespace NListUtils {
	//extern bool neighborWithinCutoff(const Float3* pos_a, const Float3* pos_b, float cutoff_offset);
	extern void cullDistantNeighbors(Simulation* simulation, NListDataCollection* nlist_data_collection);

	void static updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection,
		volatile bool* finished, int* timing, bool* mutex_lock, uint32_t step_at_update);
}

class NListManager {
public:
	NListManager(){}	// This is messy, necessary to 
	NListManager(Simulation* simulation);

	void updateNeighborLists(Simulation* simulation, bool* unused_lock, bool force_update, 
		bool async, int* timings, bool* critical_error);
	void pushNlistsToDevice(Simulation* simulation);


	int stepsSinceUpdate(uint64_t currentSimStep) const { 
		return static_cast<int>(currentSimStep - prev_update_step); 
	}
	uint64_t getPrevUpdateStep() { return prev_update_step; }


	// The new stuff
	void updateCompoundGrid(Simulation* simulation);
	void bootstrapCompoundgrid(Simulation* simulation);
	













	volatile bool updated_neighborlists_ready = 0;
	NListDataCollection* nlist_data_collection = nullptr;

private:
	// This is used for compounds with a confining_particle_sphere from key_particle BEFORE CUTOFF begins
	uint64_t prev_update_step = 0;


	void distributeCompoundsInGrid(Simulation* simulation);
	void transferCompoundgridToDevice(Simulation* simulation);
	CompoundGrid* compoundgrid_host = nullptr;
};
