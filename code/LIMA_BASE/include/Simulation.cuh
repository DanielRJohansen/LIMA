#pragma once

//#include "QuantomTypes.cuh"
#include "Constants.cuh"
#include "Bodies.cuh"
#include <map>
#include <assert.h>


struct ForceField_NB;




constexpr float SOLVENT_MASS = 18.01528f * 1e-3f;	// kg/mol		// TODO: Remove this constant from the program!!

const int DEBUGDATAF3_NVARS = 4;


// All members of this struct are double's so they can be parsed easily by std::map, without using variant
// TODO: Change these to optionals, so we can easily overload only those which values exist
struct InputSimParams {
	void overloadParams(std::map<std::string, double>& dict);

	float dt = 100.f;			// [ls]
	uint32_t n_steps = 1000;
private:
	template <typename T> void overloadParam(std::map <std::string, double>& dict, T* param, std::string key, float scalar = 1.f) {
		if (dict.count(key)) { *param = static_cast<T>(dict[key] * scalar); }
	}
};

struct SimParamsConst {
	uint64_t n_steps;
	float dt;									// [ls]
};

struct SimParams {
	SimParams(const InputSimParams& ip);
	SimParams(const SimParamsConst& spc) : constparams(spc) {}


	uint64_t step = 0;
	bool critical_error_encountered = false;	// Move into struct SimFlags, so SimParams can be const inside kernels
	float thermostat_scalar = 1.f;

	const SimParamsConst constparams;


};

struct DatabuffersDevice {
	DatabuffersDevice(const DatabuffersDevice&) = delete;
	DatabuffersDevice(size_t total_particles_upperbound);
	//~DatabuffersDevice();


	float* potE_buffer = nullptr;				// For total energy summation
	Float3* traj_buffer = nullptr;				// Absolute positions [nm]
};


// This goes on Device
struct Box {
	Box() {}
	//~Box();
	void moveToDevice();				// Loses pointer to RAM location!


	Compound* compounds = nullptr;

	uint32_t n_compounds = 0;
	uint16_t n_solvents = 0;
	uint32_t total_particles_upperbound = 0;


	CompoundCoords* coordarray_circular_queue = nullptr;
	static constexpr size_t coordarray_circular_queue_n_elements = MAX_COMPOUNDS * STEPS_PER_LOGTRANSFER;
	SolventBlockGrid* solventblockgrid_circular_queue = nullptr;
	SolventBlockTransfermodule* transfermodule_array = nullptr;

	CompoundGrid* compound_grid = nullptr;

	NeighborList* compound_neighborlists = nullptr;
	//------------------------------------//




	ForceField_NB* forcefield = nullptr;	// a replika is made available as __constant__ memory to the simulation kernels only

	CompoundBridgeBundleCompact* bridge_bundle = nullptr;
	BondedParticlesLUTManager* bonded_particles_lut_manager = nullptr;
	

	//float* potE_buffer = nullptr;				// For total energy summation
	//Float3* traj_buffer = nullptr;				// Absolute positions [nm]

	float* outdata = nullptr;					// Temp, for longging values to whatever
	Float3* data_GAN = nullptr;					// Only works if theres 1 compounds right now.
	Float3* debugdataf3 = nullptr;


	//float thermostat_scalar = 1.f;
};

struct SimulationDevice {
	SimulationDevice(const SimulationDevice&) = delete;
	SimulationDevice(const SimParams& params_host, Box* box);
	void deleteMembers();

	SimParams* params;
	Box* box;
	DatabuffersDevice* databuffers;
};

// This stays on host
class Simulation {
public:
	Simulation(const SimParams& sim_params);

	~Simulation();
	void deleteBoxMembers();
	void moveToDevice();
	void copyBoxVariables();
	void incStep() {
		assert(sim_dev);
		simparams_host.step++;
		//simparams_device->step++;
		sim_dev->params->step++;
	}
	
	inline uint64_t getStep() const { return simparams_host.step; }
	
	bool ready_to_run = false;
	bool finished = false;


	std::vector<Float3> traj_buffer;
	std::vector<float> potE_buffer;
	std::vector<float> temperature_buffer;

	Float3* traindata_buffer = nullptr;		// LimaPosition and force data for all particles, for NN training
	float* logging_data = nullptr;				// Used for debugging/logging any values. 10 floats per step!

	uint32_t total_particles_upperbound = 0;
	uint32_t total_compound_particles = 0;			// Precise number, but DO NOT EVER USE IN INDEXING!!
	uint32_t total_particles = 0;				// Precise number, but DO NOT EVER USE IN INDEXING!!
	

	//const float dt = 100.f;					// [ls]
	const int steps_per_render = STEPS_PER_RENDER;

	float temperature = -1.f;			// Current temperature [k]

	Box* box = nullptr;

	SimParams simparams_host;
	//SimParams* simparams_device = nullptr;


	//Box* box;
	bool box_is_on_device = false;

	Compound* compounds_host = nullptr;				// For reading static data, for example during nlist-search

	// Box variable copies, here for ease of access.
	int n_compounds = 0;
	int n_bridges = 0; 
	int n_solvents = 0;

	SimulationDevice* sim_dev = nullptr;

	//std::string out_dir = OUT_DIR;
	


	//int blocks_per_solventkernel = ceil((float)n_solvents / (float)THREADS_PER_SOLVENTBLOCK);
	int blocks_per_solventkernel = 0;
};

namespace SimUtils {
	Box copyToHost(const Box* box_dev);
};