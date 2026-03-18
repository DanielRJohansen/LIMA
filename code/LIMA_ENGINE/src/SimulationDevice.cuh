#pragma once

#include "EngineBodies.cuh"
#include "Bodies.cuh"
#include "Simulation.cuh"


struct BoxConfig {
	BoxConfig() {};
	static BoxConfig Create(const Box& boxHost);
	void FreeMembers() const;
};

struct BoxState {
	BoxState() {};
	BoxState(PersistentclusterInterimState*);
	static BoxState Create(const Box& boxHost);
	void CopyDataToHost(Box& boxDev) const;
	void FreeMembers() const;

	PersistentclusterInterimState* const pclusterInterimStates = nullptr;	
};

struct AdamState {
	Float3 firstMoment;
	Float3 secondMoment;
};


//struct alignas(128) CompoundQuickData {
//	Float3 relPos[MAX_COMPOUND_PARTICLES];
//	ForceField_NB::ParticleParameters ljParams[MAX_COMPOUND_PARTICLES];
//	float charges[MAX_COMPOUND_PARTICLES];
//
//	// Returns ptr to device buffer
//	__host__ static CompoundQuickData* CreateBuffer(const Simulation& sim);
//};


struct DatabuffersDeviceController {
	DatabuffersDeviceController(const DatabuffersDeviceController&) = delete;
	DatabuffersDeviceController(int nPclusters, int loggingInterval);
	~DatabuffersDeviceController();

	static const int nStepsInBuffer = 5; // TODO: I want this to be dynamic.

	static bool IsBufferFull(size_t step, int loggingInterval) {
		return step % (nStepsInBuffer * loggingInterval) == 0;
	}
	static int StepsReadyToTransfer(size_t step, int loggingInterval) {
		const int64_t stepsSinceTransfer = step % (nStepsInBuffer * loggingInterval);
		return stepsSinceTransfer / loggingInterval;
	}

	__device__ static int GetLogIndexOfParticle(int pidInPclusters, int pcId, int step,
		int loggingInterval, const int totalParticleUpperbound) {
		const int steps_since_transfer = step % (nStepsInBuffer * loggingInterval);

		const int stepOffset = steps_since_transfer / loggingInterval * totalParticleUpperbound;
		const int pclusterOffset = pcId * PersistentCluster::nParticles;
		return stepOffset + pclusterOffset + pidInPclusters;
	}

	float* potE_buffer = nullptr;				// For total energy summation
	Float3* traj_buffer = nullptr;				// Absolute positions [nm]
	float* vel_buffer = nullptr;				// Dont need direciton here, so could be a float
	Float3* forceBuffer = nullptr;				// [J/mol/nm] // For debug only

	const int nParticlesUpperbound;
};


/// <summary>
/// All members of this will only ever exist on device. Immediately after creating this class
/// it must also be moved to device.
/// </summary>
struct SimulationDevice {
	SimulationDevice(const SimulationDevice&) = delete;

	SimulationDevice(const SimParams& params_host, Box* box_host, const BoxConfig& boxConfig,
	const BoxState& boxState, const DatabuffersDeviceController&);

	// Recursively free members. Use cudaFree on *this immediately after
	void FreeMembers();

	
	

	// Module used to move solvents to a new block, in parallel
	//SolventBlockTransfermodule* transfermodule_array = nullptr;

	const SimParams params;
	SimSignals* signals = nullptr;

	const BoxConfig boxConfig;
	const BoxState boxState;
	const BoxParams boxparams;

	uint8_t* nParticlesInCompoundsBuffer = nullptr;

	// Databuffers, NOT owned by this class, so dont free them
	float* potE_buffer = nullptr;
	Float3 * traj_buffer = nullptr;
	float* vel_buffer = nullptr;
	Float3 * forceBuffer = nullptr;

	// Only used in EM
	AdamState* adamState = nullptr;
};

struct ForceEnergyInterims {
	ForceEnergyInterims(int nBondgroups, int nParticles, int nPclusters);
	void Free() const;

	// These are temp, pushed into fromSuperclusters*
	ForceEnergy* bonded = nullptr; // TODO: Also temp, not sure how i wanna proceed here..
	ForceEnergy* pme = nullptr;

	// Pushed to from pclusters
	ForceEnergy* snf = nullptr;

	// Bondgroups
	ForceEnergy* forceEnergiesBondgroups = nullptr;
};
