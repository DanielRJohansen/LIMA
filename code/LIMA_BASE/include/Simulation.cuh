#pragma once

#include "Constants.h"
#include "Bodies.cuh"
#include "Forcefield.cuh"
#include <memory>


//struct ForceField_NB;
//class Forcefield;

enum ColoringMethod { Atomname, Charge };


constexpr float SOLVENT_MASS = 18.01528f * 1e-3f;	// kg/mol		// TODO: Remove this constant from the program!!

const int DEBUGDATAF3_NVARS = 4;

enum BoundaryConditionSelect{NoBC, PBC};

enum SupernaturalForcesSelect{None, HorizontalSqueeze, HorizontalChargeField};

struct SimParams {
	SimParams() {}
	SimParams(const std::string& path);
	SimParams(uint64_t ns, float dt, bool ev, BoundaryConditionSelect bc) 
		: n_steps(ns), dt(dt), em_variant(ev), bc_select(bc) {}

	void dumpToFile(const std::string& filename = "sim_params.txt");

	uint64_t n_steps = 1000;
	float dt = 100.f;									// [ls]
	bool em_variant = false;
	BoundaryConditionSelect bc_select{ PBC };
	SupernaturalForcesSelect snf_select{ None };	// This should probably be a bitmask instead
	float box_size = 7.f;								// [nm]
	bool enable_electrostatics = false;


	ColoringMethod coloring_method = ColoringMethod::Atomname;

	int data_logging_interval = 5;

	static std::string defaultPath() { return "sim_params.txt"; };
};

struct SimSignals {
	int64_t step = 0;
	bool critical_error_encountered = false;	// Move into struct SimFlags, so SimParams can be const inside kernels
	float thermostat_scalar = 1.f;

	//const SimParams constparams;
};

struct BoxParams {
	Float3 dims{};	 // [nm]

	int n_compounds = 0;
	int n_bridges = 0;
	int64_t n_solvents = 0;
	int64_t total_particles_upperbound = 0;
	uint32_t total_particles = 0;					// Precise number. DO NOT USE IN INDEXING!!
};

// Params in simulation host-side only
struct SimparamsExtra {
	uint32_t total_compound_particles = 0;			// Precise number. DO NOT USE IN INDEXING!!
};

struct DatabuffersDevice {
	DatabuffersDevice(const DatabuffersDevice&) = delete;
	DatabuffersDevice(int total_particles_upperbound, int n_compounds, int loggingInterval);
	void freeMembers();

	static const int nStepsInBuffer = 10;

	static bool IsBufferFull(size_t step, int loggingInterval) {
		return step % (nStepsInBuffer * loggingInterval) == 0;
	}
	static int StepsReadyToTransfer(size_t step, int loggingInterval) {
		const int stepsSinceTransfer = step % (nStepsInBuffer * loggingInterval);
		return stepsSinceTransfer / loggingInterval;
	}

	__device__ int GetLogIndexOfParticle(int particleIdLocal, int compound_id, int step, int loggingInterval) const {
		const int steps_since_transfer = step % (nStepsInBuffer * loggingInterval);
		
		//const int64_t step_offset = LIMALOGSYSTEM::getDataentryIndex(steps_since_transfer, loggingInterval) * total_particles_upperbound;
		const int stepOffset = steps_since_transfer / loggingInterval * total_particles_upperbound;
		const int compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
		return stepOffset + compound_offset + particleIdLocal;
	}

	float* potE_buffer = nullptr;				// For total energy summation
	Float3* traj_buffer = nullptr;				// Absolute positions [nm]
	float* vel_buffer = nullptr;				// Dont need direciton here, so could be a float

	const int total_particles_upperbound;
#ifdef GENERATETRAINDATA
	float* outdata = nullptr;					// Temp, for longging values to whatever
	Float3* data_GAN = nullptr;					// Only works if theres 1 compounds right now.
#endif
	//Float3* debugdataf3 = nullptr;
};

template <typename T>
class ParticleDataBuffer {
public:
	ParticleDataBuffer(size_t n_particles_upperbound, size_t n_compounds, size_t n_steps, int loggingInterval) 
		: n_particles_upperbound(n_particles_upperbound), n_compounds(n_compounds), 
		n_indices(std::max(n_steps/ loggingInterval,1ull)), buffer(n_particles_upperbound* n_indices, T{})
	{
		//buffer.resize(n_particles_upperbound * n_indices);
		//for (size_t i = 0; i < n_particles_upperbound * n_indices; i++) {
		//	buffer[i] = T{};
		//}
	}

	T* data() { return buffer.data(); }	// temporary: DO NOT USE IN NEW CODE

	// Get entryindex from LIMALOGSYSTEM
	T* getBufferAtIndex(size_t entryindex) {
		return &buffer[n_particles_upperbound * entryindex];
	}

	const T* getBufferAtIndexConst(size_t entryindex) const {
		return &buffer[n_particles_upperbound * entryindex];
	}

	T& getCompoundparticleDatapointAtIndex(int compound_id, int particle_id_compound, size_t entryindex) {
		const size_t index_offset = entryindex * n_particles_upperbound;
		const size_t compound_offset = static_cast<size_t>(compound_id) * MAX_COMPOUND_PARTICLES;
		return buffer[index_offset + compound_offset + particle_id_compound];
	}

	T& getSolventparticleDatapointAtIndex(int solvent_id, size_t entryindex) {
		const size_t index_offset = entryindex * n_particles_upperbound;
		const size_t firstsolvent_offset = n_compounds * MAX_COMPOUND_PARTICLES;
		return buffer[index_offset + firstsolvent_offset + solvent_id];
	}

private:
	const size_t n_particles_upperbound;
	const size_t n_compounds;
	const size_t n_indices;
	std::vector<T> buffer;
};

// TODO: Temp, i dont like this being here
namespace LIMALOGSYSTEM {
	// Same as below, but we dont expect to be an even interval
	static constexpr int64_t getMostRecentDataentryIndex(int64_t step, int loggingInterval) {
		return step / loggingInterval;
	}

	static constexpr int64_t getNIndicesBetweenSteps(int64_t from, int64_t to, int loggingInterval) {
		return getMostRecentDataentryIndex(to, loggingInterval) - getMostRecentDataentryIndex(from, loggingInterval);
	}

	__device__ __host__ static constexpr int64_t getDataentryIndex(int64_t step, int loggingInterval) {
		if (step % loggingInterval != 0) {
			//throw std::runtime_error("This step was not expected, there is no equivalent entryindex for it.");	// TODO Maybe then return previous valid entryindex? FIXFIX DANGER
			return 0;
		}
		assert(step % loggingInterval == 0);
		return step / loggingInterval;
	}
}

// This goes on Device
struct Box {
	Box() {}
	~Box();
	void moveToDevice();				// Loses pointer to RAM location!
	void deleteMembers();

	BoxParams boxparams;


	static constexpr size_t coordarray_circular_queue_n_elements = MAX_COMPOUNDS * 3;

	// flags used for destructing only!
	bool is_on_device = false;
	bool owns_members = false;
	//bool marked_for_delete = true;

	Compound* compounds = nullptr;
	//CompoundCoords* coordarray_circular_queue = nullptr;
	CompoundcoordsCircularQueue* compoundcoordsCircularQueue = nullptr;


	Solvent* solvents = nullptr;

	//SolventCoord

	// Positions and solvent_ids
	SolventBlocksCircularQueue* solventblockgrid_circularqueue = nullptr;

	// TODO: this should be removed from box, i dont think it is used in the engine kernels
	ForceField_NB* forcefield = nullptr;	// a replika is made available as __constant__ memory to the simulation kernels only

	CompoundBridgeBundleCompact* bridge_bundle = nullptr;
	BondedParticlesLUTManager* bonded_particles_lut_manager = nullptr;


	UniformElectricField uniformElectricField;
};



// This stays on host
class Simulation {
public:
	Simulation(const SimParams& sim_params, const std::string& molecule_path, EnvMode envmode);

	void copyBoxVariables();
	
	inline int64_t getStep() const { return simsignals_host.step; }
	
	bool ready_to_run = false;
	bool finished = false;


	std::unique_ptr<ParticleDataBuffer<Float3>> traj_buffer;
	std::unique_ptr<ParticleDataBuffer<float>> potE_buffer;
	std::unique_ptr<ParticleDataBuffer<float>> vel_buffer;	// We dont need direction here, so could simply be a float

	std::vector<float> temperature_buffer;

#ifdef GENERATETRAINDATA
	std::vector<Float3> trainingdata;
	std::vector<float> loggingdata;
#endif



	float temperature = -1.f;			// Current temperature [k]

	std::unique_ptr<Box> box_host;


	SimSignals simsignals_host;	// I think this is a mistake, there should be no copy, only a pipeline to access
	SimParams simparams_host;
	BoxParams boxparams_host;	// only available after box_device has been created
	SimparamsExtra extraparams;	// only available after box_device has been created

	std::vector<Compound> compounds_host;
	std::unique_ptr<Forcefield> forcefield;
};

namespace SimUtils {
	std::unique_ptr<Box> copyToHost(Box* box_dev);
};




