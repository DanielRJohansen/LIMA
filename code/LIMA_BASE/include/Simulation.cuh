#pragma once

#include "Constants.h"
#include "Bodies.cuh"
#include <memory>
#include <filesystem>
#include "BoxGrid.cuh"
#include "SimParams.h"
#include "Backbone.h"
#include "LipidModel.h"
#include "ParticleCategory.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>

namespace MDFiles { struct TrrFile; }
struct BoxImage;



// A box length within this of an integer counts as that integer when rounding it to whole cells or chargeblocks. Above the
// precision of a .gro box (1e-5) and of a float length of up to ~100 nm, far below anything physical [nm]
constexpr float boxLengthRoundingTolerance = 1e-4f;

// A rectangular box of any size. The engine bins particles in a grid of equal cells, boxSize of them per axis. A newly built
// grid rounds the box length down; pressure coupling may retain that count across small size changes to avoid reallocating it.
struct BoxParams {
	Int3 boxSize{};		// Grid cells per axis
	Float3 boxLength{};	// [nm]
	int n_bridges = 0;
	int totalParticles = 0;
	int64_t degreesOfFreedom=0;

	__host__ static Int3 GridSizeForLength(Float3 length) {
		const auto Cells = [](float l) { return (std::max)(1, static_cast<int>(std::floor(l + boxLengthRoundingTolerance))); };
		return Int3{ Cells(length.x), Cells(length.y), Cells(length.z) };
	}
	__host__ void SetBoxLength(Float3 length) {
		boxLength = length;
		boxSize = GridSizeForLength(length);
	}
	__host__ Float3 BoxSizeFloat() const { return boxLength; }
};

inline bool GridSizeRequiresRebuild(Int3 current, Int3 preferred) {
	const auto Outside = [](int existing, int desired) {
		const float ratio = static_cast<float>(existing) / static_cast<float>(desired);
		return ratio < 0.8f || ratio > 1.2f;
	};
	return Outside(current.x, preferred.x) || Outside(current.y, preferred.y) || Outside(current.z, preferred.z);
}


template <typename T>
class ParticleDataBuffer {
public:
	ParticleDataBuffer(size_t n_particles_upperbound, size_t n_steps, 
		int loggingInterval, int nPclusters	) :
		n_particles_upperbound(nPclusters * PersistentCluster::maxParticles),
		n_indices(loggingInterval > 0 ? std::max(n_steps / loggingInterval + (n_steps % loggingInterval != 0), size_t{1}) : size_t{1}),
		buffer(nPclusters* PersistentCluster::maxParticles* n_indices, T{}),
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

	T& GetDatapoint(int pcid, int pid, size_t entryindex) {
		const size_t indexOffset = entryindex * nPclusters * PersistentCluster::maxParticles;
		const size_t pcOffset = static_cast<size_t>(pcid) * PersistentCluster::maxParticles;
		return buffer[indexOffset + pcOffset + pid];
	}

	T GetDatapoint(int pcid, int pid, size_t entryindex) const {
		const size_t indexOffset = entryindex * nPclusters * PersistentCluster::maxParticles;
		const size_t pcOffset = static_cast<size_t>(pcid) * PersistentCluster::maxParticles;
		return buffer[indexOffset + pcOffset + pid];
	}

	T& GetDatapointAtStep(int pcid, int pid, size_t step) {
		const size_t entryIndex = step / loggingInterval;
		return GetDatapoint(pcid, pid, entryIndex);
	}
	T GetDatapointAtStep(int pcid, int pid, size_t step) const {
		const size_t entryIndex = step / loggingInterval;
		return GetDatapoint(pcid, pid, entryIndex);
	}
	size_t GetLoggingInterval() const { return loggingInterval; }
	size_t EntriesPerStep() const { return n_particles_upperbound; }
	const size_t n_particles_upperbound;

private:
	const size_t loggingInterval;
	const size_t n_indices;
	const size_t nPclusters;
	std::vector<T> buffer;
};

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

	std::vector<PersistentclusterInterimState> pclusterInterimStates;

	BondGroups bondgroups;
	std::vector<PclusterConstraints> pclusterConstraints;	// Indexed by pcluster, empty without constraints
	ConstraintsSelect builtConstraints = ConstraintsSelect::None;	// The constraints setting it was built with. Settled water is constrained also with None

	UniformElectricField uniformElectricField;

	// Clusters
	std::vector<PersistentCluster> persistentClusters;
	std::vector<PersistentClusterMeta> persistentClustersMetadata;

	std::vector<ParticlesBondedToParticle> particlesBondedToParticle;
	std::vector<PclustersBondedToPcluster> pclustersBondedToPcluster;

	// Backbone chains
	BackboneChains backboneChains; // Only for rendering
	LipidModels lipidModels; // Only for rendering
	ParticleCategories particleCategories; // Only for rendering
};



// The per-particle data logged every nstenergy steps, kept on host for the whole run. Only what a run's
// outputs read is allocated, as it is O(particles * nsteps / interval)
struct LoggedData {
	bool trajectory = true;	// traj_buffer: trajectory files, self_diffusion, density_profile
	bool energy = true;		// potE_buffer, vel_buffer: energy.xvg, SimAnalysis::analyzeEnergy
	bool forces = true;		// forceBuffer: SimAnalysis::GetForces
	static constexpr LoggedData None() { return { false, false, false }; }
};

// This stays on host
struct TrajectoryStream;	// Simulation.cpp

class Simulation {
	int64_t step=0;
	std::unique_ptr<TrajectoryStream> trajectoryStream;
public:
	// Empty simulation, i dont like this very much
	Simulation(const SimParams& simparams);
	Simulation(const SimParams& simparams, std::unique_ptr<Box> box);
	
	~Simulation();

	// Allocates the buffers selected by 'logged'. Throws if they would not fit in the available RAM
	void PrepareDataBuffers();
	LoggedData logged;

	// A trajectory file written during the run, a frame every nstxout_compressed steps and the final one, instead of
	// keeping the frames in traj_buffer. The engine's log drain calls WriteTrajectoryFrame from its own thread, and
	// FinalizeSimulation the final frame once the drains are joined
	void StartTrajectoryStream(const std::filesystem::path& path);
	bool StreamsTrajectory() const { return trajectoryStream != nullptr; }
	bool TrajectoryFrameAt(int64_t step) const;
	// Positions by logging slot (pcluster * maxParticles + lane), as the engine logs them
	void WriteTrajectoryFrame(int64_t step, Float3 box, const Float3* slots);
	void FinishTrajectoryStream();

	inline int64_t getStep() const { return step; }
	
	std::unique_ptr<MDFiles::TrrFile> ToTracjectoryFile() const;

	
	bool finished = false;
	std::string name;


	std::unique_ptr<ParticleDataBuffer<Float3>> traj_buffer;	// [nm]
	std::unique_ptr<ParticleDataBuffer<float>> potE_buffer;		// [J/mol]
	std::unique_ptr<ParticleDataBuffer<float>> vel_buffer;		// [m/s]
	std::unique_ptr<ParticleDataBuffer<Float3>> forceBuffer;	// [J/mol/nm] // For debug only

	std::vector<float> temperature_buffer;	
	std::vector<float> pressure_buffer;		// [bar], one entry per pressure-coupling step
	std::vector<std::pair<int64_t,float>> maxForceBuffer; // {step,force} The maximum force experienced by any particle in the system
	struct EmLogEntry {
		int64_t step;
		float maxForce;		// [kJ/mol/nm]
		float dt;			// Unitless, see EnergyMinimizationTypes.h
	};
	std::vector<EmLogEntry> emLog; // One entry per step of energy minimization

#ifdef GENERATETRAINDATA
	std::vector<Float3> trainingdata;
	std::vector<float> loggingdata;
#endif

	std::unique_ptr<Box> box = nullptr;
	// BoxImage is defined in LIMA_TOOLS. shared_ptr keeps LIMA_BASE independent of LIMA_TOOLS.
	std::shared_ptr<BoxImage> boxImage = nullptr;
	SimParams simParams;

private:
	// Energy-minimization preconditioner made ahead of time by Engine::PrepareEnergyMinimization, so the host work
	// can overlap another running engine. Moved into the first Engine that runs this simulation
	std::vector<float> emInverseStiffness;
	std::vector<uint8_t> emWholeMolecule;

	friend class Engine;
};


struct SimStatus {
	// SimulationStatus
	std::optional<size_t> step = 0;	
	std::optional<float> progress = std::nullopt;				// [0,1]. For bounded simulations, step/nsteps. For EM the log reduction in maximum force from its initial value to emtol.
	std::optional<float> temperature = std::nullopt;			// [K]
	std::optional<float> maxForce = std::nullopt;				// [kJ/mol/nm]
	std::optional<double> expectedTimeToFinish = std::nullopt;	// [s]

	// Engine Performance
	std::optional<float> avgStepTime = std::nullopt;			// [ms]
	std::optional<float> simulationPerformance = std::nullopt;  // [ns/day]
};
