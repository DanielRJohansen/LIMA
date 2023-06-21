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

struct BoxParams {
	uint32_t n_compounds = 0;
	uint16_t n_solvents = 0;
	uint32_t total_particles_upperbound = 0;
};

// Params in simulation host-side only
struct SimparamsExtra {
	uint32_t total_compound_particles = 0;			// Precise number. DO NOT USE IN INDEXING!!
	uint32_t total_particles = 0;					// Precise number. DO NOT USE IN INDEXING!!
};

struct DatabuffersDevice {
	DatabuffersDevice(const DatabuffersDevice&) = delete;
	DatabuffersDevice(size_t total_particles_upperbound, int n_compounds);
	void freeMembers();


	float* potE_buffer = nullptr;				// For total energy summation
	Float3* traj_buffer = nullptr;				// Absolute positions [nm]

	float* outdata = nullptr;					// Temp, for longging values to whatever
	Float3* data_GAN = nullptr;					// Only works if theres 1 compounds right now.
	//Float3* debugdataf3 = nullptr;
};

template <typename T>
class ParticleDataBuffer {
public:
	ParticleDataBuffer(size_t n_particles_upperbound, size_t n_compounds, size_t n_steps) 
		: n_particles_upperbound(n_particles_upperbound), n_compounds(n_compounds)
	{
		buffer.resize(n_particles_upperbound * n_steps);
	}

	T* data() { return buffer.data(); }	// temporary: DO NOT USE IN NEW CODE

	T* getBufferAtStep(size_t step) {
		return &buffer[n_particles_upperbound * step];
	}

	T& getCompoundparticleDatapoint(int compound_id, int particle_id_compound, size_t step) {
		const uint32_t step_offset = static_cast<uint32_t>(step) * n_particles_upperbound;
		const uint32_t compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
		return buffer[step_offset + compound_offset + particle_id_compound];
	}

	T& getSolventparticleDatapoint(int solvent_id, size_t step) {
		const uint32_t step_offset = static_cast<uint32_t>(step) * n_particles_upperbound;
		const uint32_t firstsolvent_offset = n_compounds * MAX_COMPOUND_PARTICLES;
		return buffer[step_offset + firstsolvent_offset + solvent_id];
	}

private:
	const size_t n_particles_upperbound = 0;
	const size_t n_compounds = 0;
	std::vector<T> buffer;
};

// This goes on Device
struct Box {
	Box() {}
	~Box();
	void moveToDevice();				// Loses pointer to RAM location!
	void deleteMembers();

	BoxParams boxparams;

	static constexpr size_t coordarray_circular_queue_n_elements = MAX_COMPOUNDS * STEPS_PER_LOGTRANSFER;

	// flags used for destructing only!
	bool is_on_device = false;
	bool owns_members = false;
	//bool marked_for_delete = true;

	Compound* compounds = nullptr;
	CompoundCoords* coordarray_circular_queue = nullptr;
	SolventBlockGrid* solventblockgrid_circular_queue = nullptr;

	SolventBlockTransfermodule* transfermodule_array = nullptr;
	CompoundGrid* compound_grid = nullptr;

	NeighborList* compound_neighborlists = nullptr;


	// TODO: this should be removed from box, i dont think it is used in the engine kernels
	ForceField_NB* forcefield = nullptr;	// a replika is made available as __constant__ memory to the simulation kernels only

	CompoundBridgeBundleCompact* bridge_bundle = nullptr;
	BondedParticlesLUTManager* bonded_particles_lut_manager = nullptr;

};

struct SimulationDevice {
	SimulationDevice(const SimulationDevice&) = delete;
	SimulationDevice(const SimParams& params_host, std::unique_ptr<Box> box);
	
	// Recursively free members
	void deleteMembers();  // Use cudaFree on *this immediately after


	SimParams* params;
	Box* box;
	DatabuffersDevice* databuffers;
};

// This stays on host
class Simulation {
public:
	Simulation(const SimParams& sim_params);

	~Simulation();
	void moveToDevice();
	void copyBoxVariables();
	void incStep() {
		assert(sim_dev);
		simparams_host.step++;
		sim_dev->params->step++;
	}
	
	inline uint64_t getStep() const { return simparams_host.step; }
	
	bool ready_to_run = false;
	bool finished = false;


	//std::vector<Float3> traj_buffer;
	std::unique_ptr<ParticleDataBuffer<Float3>> traj_buffer;
	std::vector<float> potE_buffer;
	std::vector<float> temperature_buffer;

	// TODO: Make these vectors instead
	//Float3* traindata_buffer = nullptr;		// LimaPosition and force data for all particles, for NN training
	//float* logging_data = nullptr;				// Used for debugging/logging any values. 10 floats per step!
	std::vector<Float3> trainingdata;
	std::vector<float> loggingdata;




	//const float dt = 100.f;					// [ls]
	const int steps_per_render = STEPS_PER_RENDER;

	float temperature = -1.f;			// Current temperature [k]

	std::unique_ptr<Box> box_host;

	SimParams simparams_host;
	BoxParams boxparams_host;	// only available after box_device has been created
	SimparamsExtra extraparams;	// only available after box_device has been created

	std::vector<Compound> compounds_host;

	// Box variable copies, here for ease of access.
	//int n_compounds = 0;
	int n_bridges = 0; 
	//int n_solvents = 0;

	SimulationDevice* sim_dev = nullptr;
	//int blocks_per_solventkernel = 0;
};

namespace SimUtils {
	std::unique_ptr<Box> copyToHost(Box* box_dev);
};