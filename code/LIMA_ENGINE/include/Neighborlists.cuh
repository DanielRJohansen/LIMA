#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "EngineUtils.cuh"
#include <chrono>
#include <thread>

struct CandidateList {
	int n_candidates = 0;
	std::array<uint32_t, NEIGHBORLIST_MAX_SOLVENTS> candidates{};
};

struct NListDataCollection {
	//NListDataCollection() {}
	NListDataCollection(Simulation* simulation) {
		n_compounds = simulation->n_compounds;
		n_solvents = simulation->n_solvents;
		compoundstates = new CompoundState[n_compounds];
		//solvents = new Solvent[simulation->n_solvents];
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
			//solvent_positions[i] = solvents[i].pos;
		}
	}

	void preparePositionData(const Simulation& simulation) {
		// Data for the current step has not yet been generated so we need to use the previous step.
		// For the very first step, engine has cheated and already written the traj from the initial setup.
		auto step = simulation.getStep();
		if (step != 0) { step--; }
		
		for (int compound_id = 0; compound_id < n_compounds; compound_id++) {
			const size_t index = EngineUtils::getAlltimeIndexOfParticle(step, simulation.total_particles_upperbound, compound_id, 0);
			compound_key_positions[compound_id] = simulation.traj_buffer[index];
		}
		for (int solvent_id = 0; solvent_id < n_solvents; solvent_id++) {
			const size_t index = EngineUtils::getAlltimeIndexOfParticle(step, simulation.total_particles_upperbound, simulation.n_compounds, solvent_id);
			solvent_positions[solvent_id] = simulation.traj_buffer[index];
		}
	}


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
		volatile bool* finished, int* timing, bool* mutex_lock);
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

#include <array>

class SolventBlockCollection {
public:
	SolventBlockCollection(const Float3* positions, int n_solvents);
	void addSolventId(uint32_t id, const Float3& pos);

	// Returns a vector with 1 element per solvent in sim. Each element consists of a vector
	// of the Ids of nearby solvents of this solvent
	std::vector<CandidateList> getNeighborSolventForAllSolvents(uint32_t n_solvents);




private:
	static constexpr int blocks_per_dim = BOX_LEN_NM / CUTOFF_NM;
	static constexpr int blocks_total = blocks_per_dim * blocks_per_dim * blocks_per_dim;
	static constexpr float blocks_per_dim_float = static_cast<float>(blocks_per_dim);
	static constexpr float uncovered_len = (BOX_LEN - CUTOFF_LM * blocks_per_dim_float);
	static constexpr float block_len = CUTOFF_LM + uncovered_len / blocks_per_dim_float;


	struct SolventBlock {
		static const int max_solvents_in_block = 512;	// i dunno
		std::array<uint32_t, max_solvents_in_block> solvent_ids;
		int n_elements{ 0 };
		void insert(uint32_t solvent_id) {
			solvent_ids[n_elements++] = solvent_id;
		}
	};



	std::array<std::array< std::array<SolventBlock, blocks_per_dim>, blocks_per_dim>, blocks_per_dim> m_blocks;

	static Int3 getSolventblockIndex(const Float3& pos);
	static const int a = blocks_per_dim;
	static constexpr std::array<Int3, blocks_per_dim> getAllIndices();
	static constexpr std::array<Int3, 2*2*2> getAdjacentIndicesThatAreGreater(Int3 index);

	SolventBlock& getBlock(Int3 index);

	static void addAllInsideBlock(std::vector<CandidateList>& neighborCandidates,
		const SolventBlock& block);
	static void addAllBetweenBlocks(std::vector<CandidateList>& neighborCandidates,
		const SolventBlock& blocka, const SolventBlock& blockb);
};

