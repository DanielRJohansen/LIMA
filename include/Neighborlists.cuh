#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "EngineUtils.cuh"
#include <chrono>

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

	//Static as it must be able to run in a separate thread
	//static void updataNeighborLists(bool async);
	void NListManager::offloadPositionDataNLIST(Simulation* simulation);
	void NListManager::pushNlistsToDevice(Simulation* simulation);


	bool updated_neighborlists_ready = 0;
	NListDataCollection* nlist_data_collection = nullptr;

private:
	// This is used for compounds with a confining_particle_sphere from key_particle BEFORE CUTOFF begins
};