#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "EngineUtils.cuh"
#include <chrono>
#include <thread>

struct NListDataCollection {
	//NListDataCollection() {}
	NListDataCollection(Simulation* simulation) {
		n_compounds = simulation->n_compounds;
		n_solvents = simulation->n_solvents;
		compoundstates = new CompoundState[n_compounds];
		solvents = new Solvent[simulation->n_solvents];
		compound_neighborlists = new NeighborList[MAX_COMPOUNDS];
		solvent_neighborlists = new NeighborList[MAX_SOLVENTS];
		cudaMemcpy(compound_neighborlists, simulation->box->compound_neighborlists, sizeof(NeighborList) * n_compounds, cudaMemcpyDeviceToHost);
		cudaMemcpy(solvent_neighborlists, simulation->box->solvent_neighborlists, sizeof(NeighborList) * n_solvents, cudaMemcpyDeviceToHost);
	}
	/*
	void preparePositionData() {
		for (int i = 0; i < n_compounds; i++) {
			compound_key_positions[i] = compoundstates[i].positions[0];
		}
		for (int i = 0; i < n_solvents; i++) {
			solvent_positions[i] = solvents[i].pos;
		}
	}*/
	void preparePositionData(Compound* compounds) {
		for (int i = 0; i < n_compounds; i++) {
			compound_key_positions[i] = compoundstates[i].positions[compounds[i].key_particle_index];
		}
		for (int i = 0; i < n_solvents; i++) {
			solvent_positions[i] = solvents[i].pos;
		}
	}
	int n_compounds;
	int n_solvents;

	// I guess this is not critical but temp, needed to load pos device->host
	CompoundState* compoundstates;
	Solvent* solvents;

	Float3 compound_key_positions[MAX_COMPOUNDS];
	Float3 solvent_positions[MAX_SOLVENTS];

	// These are loaded before simulaiton start. Kept on host, and copied to device each update.
	NeighborList* compound_neighborlists;
	NeighborList* solvent_neighborlists;


};


namespace NListUtils {
	 extern void updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection,
		volatile bool* finished, int* timing, bool* mutex_lock);

	extern bool neighborWithinCutoff(Float3* pos_a, Float3* pos_b, float cutoff_offset);
	extern void cullDistantNeighbors(Simulation* simulation, NListDataCollection* nlist_data_collection);
}

class NListManager {
public:
	NListManager(){}	// This is messy, necessary to 
	NListManager(Simulation* simulation);

	void updateNeighborLists(Simulation* simulation, bool* unused_lock, bool force_update, bool async, int* timings);
	void offloadPositionDataNLIST(Simulation* simulation);
	void pushNlistsToDevice(Simulation* simulation);


	int stepsSinceUpdate(uint64_t currentSimStep) { return static_cast<int>(currentSimStep - prev_update_step); }
	uint64_t getPrevUpdateStep() { return prev_update_step; }


	volatile bool updated_neighborlists_ready = 0;
	NListDataCollection* nlist_data_collection = nullptr;

private:
	// This is used for compounds with a confining_particle_sphere from key_particle BEFORE CUTOFF begins
	uint64_t prev_update_step = 0;
};