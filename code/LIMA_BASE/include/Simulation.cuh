#pragma once

#include "Constants.h"
#include "Bodies.cuh"
#include <memory>
#include <filesystem>
#include "BoxGrid.cuh"

namespace MDFiles { struct TrrFile; }

enum ColoringMethod { Atomname, Charge, GradientFromAtomid, GradientFromCompoundId };

enum BoundaryConditionSelect{NoBC, PBC};

enum SupernaturalForcesSelect{None, HorizontalSqueeze, HorizontalChargeField, BoxEdgePotential};

struct SimParams {
	SimParams() {}
	SimParams(const std::filesystem::path& path);
	SimParams(std::initializer_list<int>) = delete;

	void dumpToFile(const std::filesystem::path& filename = "sim_params.txt");

	// Main params
	uint64_t n_steps = 1000;
	float dt = 100.f;									// [ls]
	bool em_variant = false;
	float em_force_tolerance = 1000; // [kJ/mol/nm]

	// Physics params
	BoundaryConditionSelect bc_select{ PBC };
	bool enable_electrostatics = false;
	float cutoff_nm = 1.2f;
	SupernaturalForcesSelect snf_select{ None };	// This should probably be a bitmask instead

	// Output params
	int data_logging_interval = 5;
	bool save_trajectory = false;
	bool save_energy = false;
	ColoringMethod coloring_method = ColoringMethod::Atomname;

	// Thermostat
	int steps_per_temperature_measurement = 200;
	bool apply_thermostat = false;

};

struct SimSignals {
	int64_t step = 0;
	bool critical_error_encountered = false;	// Move into struct SimFlags, so SimParams can be const inside kernels
	float thermostat_scalar = 1.f;
};



struct BoxParams {
	int boxSize = 0;	// [nm]

	int n_compounds = 0;
	int n_bridges = 0;
	int64_t n_solvents = 0;
	int64_t total_particles_upperbound = 0;
	uint32_t total_particles = 0;					// Precise number. DO NOT USE IN INDEXING!!
	uint32_t total_compound_particles = 0;			// Precise number. DO NOT USE IN INDEXING!!
};

// TODO: Move to Simulation device?
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
	Float3* forceBuffer = nullptr;				// [1/l N/mol] // For debug only
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
	ParticleDataBuffer(size_t n_particles_upperbound, size_t n_compounds, size_t n_steps, int loggingInterval) : 
		n_particles_upperbound(n_particles_upperbound), n_compounds(n_compounds), 
		n_indices(std::max(n_steps/ loggingInterval,static_cast<size_t>(1))), 
		buffer(n_particles_upperbound* n_indices, T{}),
		loggingInterval(loggingInterval)
	{}

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

	size_t GetLoggingInterval() const { return loggingInterval; }
	size_t EntriesPerStep() const { return n_particles_upperbound; }
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
	Box(int boxSize);

	BoxParams boxparams;


	std::vector<Compound> compounds;
	std::unique_ptr<CompoundcoordsCircularQueue> compoundcoordsCircularQueue = nullptr;

	std::vector<Solvent> solvents;
	std::unique_ptr<SolventBlocksCircularQueue> solventblockgrid_circularqueue = nullptr;

	std::unique_ptr<CompoundBridgeBundleCompact> bridge_bundle = nullptr;
	std::vector<BondedParticlesLUT> bpLutCollection;

	UniformElectricField uniformElectricField;
};



// This stays on host
class Simulation {
public:
	// Empty simulation, i dont like this very much
	Simulation(const SimParams& simparams);
	Simulation(const SimParams& simparams, std::unique_ptr<Box> box);
	
	void PrepareDataBuffers();

	inline int64_t getStep() const { return simsignals_host.step; }
	
	std::unique_ptr<MDFiles::TrrFile> ToTracjectoryFile() const;

	bool ready_to_run = false;
	bool finished = false;


	std::unique_ptr<ParticleDataBuffer<Float3>> traj_buffer;// [nm]
	std::unique_ptr<ParticleDataBuffer<float>> potE_buffer;	// [J/mol]
	std::unique_ptr<ParticleDataBuffer<float>> vel_buffer;	// [m/s]
	std::unique_ptr<ParticleDataBuffer<Float3>> forceBuffer;	// [1/l N/mol] // For debug only

	std::vector<float> temperature_buffer;
	std::vector<float> maxForceBuffer; // The maximum force experienced by any particle in the system

#ifdef GENERATETRAINDATA
	std::vector<Float3> trainingdata;
	std::vector<float> loggingdata;
#endif

	std::unique_ptr<Box> box_host = nullptr;


	SimSignals simsignals_host;	// I think this is a mistake, there should be no copy, only a pipeline to access
	SimParams simparams_host;

	ForceField_NB forcefield;
};

namespace SimUtils {
	std::unique_ptr<Box> copyToHost(Box* box_dev);
};




