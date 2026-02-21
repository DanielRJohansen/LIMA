#pragma once

#include "Constants.h"
#include "Bodies.cuh"
#include <memory>
#include <filesystem>
#include "BoxGrid.cuh"
#include "SimParams.h"
#include <set>

namespace MDFiles { struct TrrFile; }



struct SimSignals {
	bool critical_error_encountered = false;	// Move into struct SimFlags, so SimParams can be const inside kernels
};



struct BoxParams {
	Int3 boxSize{};	// [nm]
	int n_compounds = 0;
	int n_bridges = 0;
	int nTinymols = 0;
	int nTinymolParticles = 0;
	int total_particles_upperbound = 0;
	int total_particles = 0;					// Precise number. DO NOT USE IN INDEXING!!
	int total_compound_particles = 0;			// Precise number. DO NOT USE IN INDEXING!!
	int totalParticles = 0;
	int64_t degreesOfFreedom=0;

	__host__ Float3 BoxSizeFloat() const {
		return Float3{ static_cast<float>(boxSize.x), static_cast<float>(boxSize.y), static_cast<float>(boxSize.z) };
	}
};


template <typename T>
class ParticleDataBuffer {
public:
	ParticleDataBuffer(size_t n_particles_upperbound, size_t n_compounds, size_t n_steps, 
		int loggingInterval, int nPclusters	) :
		n_particles_upperbound(nPclusters * PersistentCluster::nParticles), n_compounds(n_compounds), 
		n_indices(std::max(n_steps/ loggingInterval,static_cast<size_t>(1))), 
		buffer(nPclusters* PersistentCluster::nParticles* n_indices, T{}),
		loggingInterval(loggingInterval)
		,nPclusters(nPclusters)
	{}

	T* data() { return buffer.data(); }	// temporary: DO NOT USE IN NEW CODE

	const std::vector<T>& GetBuffer() {
		return buffer;
	}

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

	T& GetDatapoint(int pcid, int pid, size_t entryindex) {
		const size_t indexOffset = entryindex * nPclusters * PersistentCluster::nParticles;
		const size_t pcOffset = static_cast<size_t>(pcid) * PersistentCluster::nParticles;
		return buffer[indexOffset + pcOffset + pid];
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
	const size_t n_particles_upperbound;

private:
	const size_t loggingInterval;
	const size_t n_compounds;
	const size_t n_indices;
	const size_t nPclusters;
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

struct Box {
	Box() {}
	Box(Float3 boxSize);

	BoxParams boxparams;


	std::vector<Compound> compounds;
	//std::vector<CompoundInterimState> compoundInterimStates;
	std::vector<PersistentclusterInterimState> pclusterInterimStates;
	std::vector<CompoundCoords> compoundCoordsBuffer;
	std::vector<TinyMolParticleState> tinyMolParticlesState;
	std::vector<SolventBlock> solventblockgrid_circularqueue;
	std::vector<BondedParticlesLUT> bpLutCollection;

	std::vector<BondGroup> bondgroups;

	UniformElectricField uniformElectricField;

	// Clusters
	std::vector<PersistentCluster> persistentClusters;
	std::vector<PersistentClusterMeta> persistentClustersMetadata;
	std::vector<std::set<int>> particleBondedToParticle;
	std::vector<std::set<int>> pclusterBondedToPcluster;
	std::vector<ParticleToCompoundOrSolventMapping> particleToCompoundOrSolventMapping;

};



// This stays on host
class Simulation {
	int64_t step=0;
public:
	// Empty simulation, i dont like this very much
	Simulation(const SimParams& simparams);
	Simulation(const SimParams& simparams, std::unique_ptr<Box> box);
	
	void PrepareDataBuffers();

	inline int64_t getStep() const { return step; }
	
	std::unique_ptr<MDFiles::TrrFile> ToTracjectoryFile() const;

	
	bool ready_to_run = false;
	bool finished = false;


	std::unique_ptr<ParticleDataBuffer<Float3>> traj_buffer;	// [nm]
	std::unique_ptr<ParticleDataBuffer<float>> potE_buffer;		// [J/mol]
	std::unique_ptr<ParticleDataBuffer<float>> vel_buffer;		// [m/s]
	std::unique_ptr<ParticleDataBuffer<Float3>> forceBuffer;	// [J/mol/nm] // For debug only

	std::vector<float> temperature_buffer;	
	std::vector<std::pair<int64_t,float>> maxForceBuffer; // {step,force} The maximum force experienced by any particle in the system

#ifdef GENERATETRAINDATA
	std::vector<Float3> trainingdata;
	std::vector<float> loggingdata;
#endif

	std::unique_ptr<Box> box_host = nullptr;


	SimSignals simsignals_host;	// I think this is a mistake, there should be no copy, only a pipeline to access
	SimParams simparams_host;

	ForceField_NB forcefield;
	ForcefieldTinymol forcefieldTinymol;
	std::vector<NonbondedInteractionParams> forcefieldTest;


	friend class Engine;
};


