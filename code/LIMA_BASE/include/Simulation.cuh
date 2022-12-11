#pragma once

//#include "QuantomTypes.cuh"
#include "Constants.cuh"
#include "Bodies.cuh"
#include <map>




#ifndef __linux__
//const string MOL_FOLDER = "C:\\PROJECTS\\Quantom\\molecules\\t4lys_full\\";
//const std::string MOL_FOLDER = "C:\\PROJECTS\\Quantom\\Simulation\\Molecule\\";
//const std::string OUT_DIR = "C:\\PROJECTS\\Quantom\\Simulation\\";
#else
//const string MOL_FOLDER = "../Compounds/t4lys/";
const string MOL_FOLDER = "../../Simulation/Molecule/";
//const string OUT_DIR = "/home/lima/Desktop/LIMA";
const string OUT_DIR = "../../Simulation/";
#endif




//const int BLOCKS_PER_SOLVENTKERNEL = ceil((float)N_SOLVATE_MOLECULES / (float)THREADS_PER_SOLVENTBLOCK);




constexpr float SOLVENT_MASS = 18.01528f * 1e-3f;	// kg/mol		// TODO: Remove this constant from the program!!
//constexpr double SOLVENT_MASS = 12.0107 * 1e-3;	// kg/mol
//constexpr double COMPOUNDPARTICLE_MASS = 12.0107 * 1e-3;

// This goes on Device
class Box {
public:
	Box() {}
	//~Box();
	void moveToDevice();				// Loses pointer to RAM location!


	Compound* compounds = nullptr;

	uint32_t n_compounds = 0;
	uint16_t n_solvents = 0;
	uint32_t total_particles_upperbound = 0;


	// These are shared for all compounds, MUST be allocated before adding any compounds to box, so not in moveToDevice //
	CompoundState* compound_state_array = nullptr;
	CompoundState* compound_state_array_prev = nullptr;
	CompoundState* compound_state_array_next = nullptr;

	CompoundCoords* compound_coord_array = nullptr;
	CompoundCoords* compound_coord_array_prev = nullptr;
	CompoundCoords* compound_coord_array_next = nullptr;

	NeighborList* compound_neighborlists = nullptr;
	NeighborList* solvent_neighborlists = nullptr;
	//------------------------------------//

	Solvent* solvents = nullptr;
	Solvent* solvents_next = nullptr;

	CompoundBridgeBundleCompact* bridge_bundle = nullptr;

	BondedParticlesLUTManager* bonded_particles_lut_manager = nullptr;

	uint32_t step = 0;
	float dt = 0;
	bool critical_error_encountered = 0;

	float* potE_buffer = nullptr;		// For total energy summation
	Float3* traj_buffer = nullptr;

	float* outdata = nullptr;			// Temp, for longging values to whatever
	Float3* data_GAN = nullptr;			// Only works if theres 1 compounds right now.


	float thermostat_scalar = 1.f;
};

// All members of this struct are double's so they can be parsed easily by std::map, without using variant
struct SimulationParams {
	void overloadParams(std::map<std::string, double>& dict);

	float dt = 1.f * 1e-6f;	// [ns]
	int n_steps = 1000;
private:
	template <typename T> void overloadParam(std::map <std::string, double>& dict, 
		T* param, std::string key, float scalar = 1.f) {
		if (dict.count(key)) { *param = static_cast<T>(dict[key]) * scalar; }
	}
};

// This stays on host
class Simulation {
public:
	Simulation(SimulationParams& sim_params);
	~Simulation();
	void deleteBoxMembers();
	void moveToDevice();
	void copyBoxVariables();
	void incStep() {
		step++;
		box->step++;
	}
	
	inline uint64_t getStep() { return step; }
	

	bool finished = false;



	float* potE_buffer = nullptr;	// Not really a buffer yet, just one large array that holds full simulation data
	Float3* traj_buffer = nullptr;
	float* temperature_buffer = nullptr;
	int n_temp_values = 0;
	Float3* traindata_buffer = nullptr;		// Position and force data for all particles, for NN training
	float* logging_data = nullptr;				// Used for debugging/logging any values. 10 floats per step!

	uint32_t total_particles_upperbound = 0;
	uint32_t total_compound_particles = 0;			// Precise number, but DO NOT EVER USE IN INDEXING!!
	uint32_t total_particles = 0;				// Precise number, but DO NOT EVER USE IN INDEXING!!
	
	uint64_t n_steps = 0;

	const float dt = 1.f;					// [ns] first val corresponds to fs
	const float dt_pico = dt * 1000.f;
	const int steps_per_render = STEPS_PER_RENDER;
	//int n_bodies = N_BODIES_START;

	Box* box;
	bool box_is_on_device = false;

	Compound* compounds_host = nullptr;				// For reading static data, for example during nlist-search

	// Box variable copies, here for ease of access.
	int n_compounds = 0;
	int n_bridges = 0; 
	int n_solvents = 0;



	//std::string out_dir = OUT_DIR;
	


	//int blocks_per_solventkernel = ceil((float)n_solvents / (float)THREADS_PER_SOLVENTBLOCK);
	int blocks_per_solventkernel = 0;
private:
	uint64_t step = 0;

	
	
};

