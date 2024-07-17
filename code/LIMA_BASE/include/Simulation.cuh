#pragma once

#include "Constants.h"
#include "Bodies.cuh"
#include <memory>
#include <filesystem>
#include "BoxGrid.cuh"


enum ColoringMethod { Atomname, Charge };



enum BoundaryConditionSelect{NoBC, PBC};

enum SupernaturalForcesSelect{None, HorizontalSqueeze, HorizontalChargeField};

struct SimParams {
	SimParams() {}
	SimParams(const std::filesystem::path& path);
	SimParams(uint64_t ns, float dt, bool ev, BoundaryConditionSelect bc) 
		: n_steps(ns), dt(dt), em_variant(ev), bc_select(bc) {}

	void dumpToFile(const std::filesystem::path& filename = "sim_params.txt");

	uint64_t n_steps = 1000;
	float dt = 100.f;									// [ls]
	bool em_variant = false;
	BoundaryConditionSelect bc_select{ PBC };
	SupernaturalForcesSelect snf_select{ None };	// This should probably be a bitmask instead
	float box_size = 7.f;								// [nm]
	bool enable_electrostatics = false;


	ColoringMethod coloring_method = ColoringMethod::Atomname;

	int data_logging_interval = 5;

	float cutoff_nm = 1.2f;

	static std::string defaultPath() { return "sim_params.txt"; };
};

struct SimSignals {
	int64_t step = 0;
	bool critical_error_encountered = false;	// Move into struct SimFlags, so SimParams can be const inside kernels
	float thermostat_scalar = 1.f;
};



struct BoxParams {
	//Float3 dims{};	 // [nm]
	//BoxSize boxSize{ 0 };
	int boxSize;	// [nm]

	int n_compounds = 0;
	int n_bridges = 0;
	int64_t n_solvents = 0;
	int64_t total_particles_upperbound = 0;
	uint32_t total_particles = 0;					// Precise number. DO NOT USE IN INDEXING!!
	uint32_t total_compound_particles = 0;			// Precise number. DO NOT USE IN INDEXING!!
};

// Params in simulation host-side only
struct SimparamsExtra {
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
		n_indices(std::max(n_steps/ loggingInterval,static_cast<size_t>(1))), buffer(n_particles_upperbound* n_indices, T{}),
		loggingInterval(loggingInterval)
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

	const T* GetBufferAtStep(size_t step) const {
		const size_t entryIndex = step / loggingInterval;
		return &buffer[n_particles_upperbound * entryIndex];
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

	T& GetMostRecentCompoundparticleDatapoint(int compound_id, int particle_id_compound, size_t step) {
		const size_t entryIndex = step / loggingInterval;		
		return getCompoundparticleDatapointAtIndex(compound_id, particle_id_compound, entryIndex);
	}

	T& GetMostRecentSolventparticleDatapointAtIndex(int solvent_id, size_t step) {
		const size_t entryIndex = step / loggingInterval;
		return getSolventparticleDatapointAtIndex(solvent_id, entryIndex);
	}

private:
	const size_t loggingInterval;
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
	//ForceField_NB* forcefield = nullptr;	// a replika is made available as __constant__ memory to the simulation kernels only

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


	std::unique_ptr<ParticleDataBuffer<Float3>> traj_buffer;// [nm]
	std::unique_ptr<ParticleDataBuffer<float>> potE_buffer;	// [J/mol]
	std::unique_ptr<ParticleDataBuffer<float>> vel_buffer;	// [m/s]

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
	//std::unique_ptr<Forcefield> forcefield;
	ForceField_NB forcefield;
};

namespace SimUtils {
	std::unique_ptr<Box> copyToHost(Box* box_dev);
};




