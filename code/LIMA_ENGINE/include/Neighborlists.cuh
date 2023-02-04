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
	static_assert(blocks_per_dim > 1, "Too few solvents blocks in box");
	static constexpr int blocks_total = blocks_per_dim * blocks_per_dim * blocks_per_dim;
	static constexpr float blocks_per_dim_float = static_cast<float>(blocks_per_dim);
	static constexpr float uncovered_len = (BOX_LEN - CUTOFF_LM * blocks_per_dim_float);
	// Make the block len a little longer than it needs to be, so we include particles on the box-boundary
	static constexpr float block_len = CUTOFF_LM + uncovered_len / (blocks_per_dim_float - 1.f);


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
	
	static constexpr std::array<Int3, blocks_total> getAllIndices();
	static constexpr std::array<Int3, 2*2*2 - 1> getAdjacentIndicesThatAreGreater(Int3 index);
	static constexpr std::array<std::array<std::array<std::array<Int3, 2 * 2 * 2 - 1>,
		blocks_per_dim>, blocks_per_dim>, blocks_per_dim> precalcGreaterIndices();


	const std::array<std::array<std::array<std::array<Int3, 2 * 2 * 2 - 1>,
		blocks_per_dim>, blocks_per_dim>, blocks_per_dim> precalcedGreaterIndices = precalcGreaterIndices();

	SolventBlock& getBlock(const Int3& index);

	static void addAllInsideBlock(std::vector<CandidateList>& neighborCandidates,
		const SolventBlock& block);
	static void addAllBetweenBlocks(std::vector<CandidateList>& neighborCandidates,
		const SolventBlock& blocka, const SolventBlock& blockb);
};