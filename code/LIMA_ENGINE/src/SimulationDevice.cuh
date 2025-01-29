#pragma once

#include "EngineBodies.cuh"
#include "Bodies.cuh"
#include "Simulation.cuh"


struct BoxConfig {	
	BoxConfig(Compound* compounds, uint8_t* compoundsAtomTypes, float* compoundsAtomCharges, BondedParticlesLUT* bpLUTs,
	const BoxGrid::TinymolBlockAdjacency::BlockRef* tinymolNearbyBlockIds);
	static BoxConfig Create(const Box& boxHost); // Returns a ptr to device
	void FreeMembers() const;

	// CompoundData used ALOT, kept here for memory locality
	const uint8_t* const compoundsAtomtypes;
	const float* const compoundsAtomCharges;	// [kC/mol]
	const Compound* const compounds;

	// BondedParticlesLUT data - NEVER access directly, use the bpLUTHelpers namespace
	const BondedParticlesLUT* const bpLUTs;

	const BoxGrid::TinymolBlockAdjacency::BlockRef* tinymolNearbyBlockIds;
};

struct BoxState {
	BoxState(NodeIndex* compoundsOrigos, Float3* compoundsRelpos, TinyMolState* tinyMols,
		SolventBlock* solventblockgrid_circularqueue, CompoundInterimState* compoundInterimState);
	static BoxState Create(const Box& boxHost);
	void CopyDataToHost(Box& boxDev) const;
	void FreeMembers() const;

	CompoundInterimState* const compoundsInterimState;
	NodeIndex* const compoundOrigos;
	Float3* const compoundsRelposNm;

	TinyMolState* const tinyMols;
	SolventBlock* const solventblockgrid_circularqueue;
};

struct AdamState {
	Float3 firstMoment;
	Float3 secondMoment;
};


struct alignas(128) CompoundQuickData {
	Float3 relPos[MAX_COMPOUND_PARTICLES];
	ForceField_NB::ParticleParameters ljParams[MAX_COMPOUND_PARTICLES];
	float charges[MAX_COMPOUND_PARTICLES];

	// Returns ptr to device buffer
	__host__ static CompoundQuickData* CreateBuffer(const Simulation& sim);
};


struct DatabuffersDeviceController {
	DatabuffersDeviceController(const DatabuffersDeviceController&) = delete;
	DatabuffersDeviceController(int total_particles_upperbound, int n_compounds, int loggingInterval);
	~DatabuffersDeviceController();

	static const int nStepsInBuffer = 5;

	static bool IsBufferFull(size_t step, int loggingInterval) {
		return step % (nStepsInBuffer * loggingInterval) == 0;
	}
	static int StepsReadyToTransfer(size_t step, int loggingInterval) {
		const int64_t stepsSinceTransfer = step % (nStepsInBuffer * loggingInterval);
		return stepsSinceTransfer / loggingInterval;
	}

	__device__ static int GetLogIndexOfParticle(int particleIdLocal, int compound_id, int64_t step,
		int loggingInterval, int totalParticleUpperbound) {
		const int64_t steps_since_transfer = step % (nStepsInBuffer * loggingInterval);

		const int64_t stepOffset = steps_since_transfer / loggingInterval * totalParticleUpperbound;
		const int compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
		return stepOffset + compound_offset + particleIdLocal;
	}

	float* potE_buffer = nullptr;				// For total energy summation
	Float3* traj_buffer = nullptr;				// Absolute positions [nm]
	float* vel_buffer = nullptr;				// Dont need direciton here, so could be a float
	Float3* forceBuffer = nullptr;				// [J/mol/nm] // For debug only

	const int total_particles_upperbound;
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
	SolventBlockTransfermodule* transfermodule_array = nullptr;

	const SimParams params;
	SimSignals* signals = nullptr;

	const BoxConfig boxConfig;
	const BoxState boxState;
	const BoxParams boxparams;

	uint8_t* nParticlesInCompoundsBuffer = nullptr;
	CompoundInteractionBoundary* compoundsInteractionBoundaryBuffer = nullptr;

	// Databuffers, NOT owned by this class, so dont free them
	float* potE_buffer = nullptr;
	Float3 * traj_buffer = nullptr;
	float* vel_buffer = nullptr;
	Float3 * forceBuffer = nullptr;

	// Only used in EM
	AdamState* adamState = nullptr;
};

struct ForceEnergyInterims {
	ForceEnergyInterims(int nCompounds, int nSolvents, int nSolventblocks, int nBondgroups);
	void Free() const;

	__device__ ForceEnergy SumCompound(int compoundId, int particleId) const {
		ForceEnergy pmeFE = {};
		if constexpr (ENABLE_ES_LR) {
			pmeFE = forceEnergiesPME[compoundId * MAX_COMPOUND_PARTICLES + particleId];
		}

		return forceEnergyFarneighborShortrange[compoundId * MAX_COMPOUND_PARTICLES + particleId]
			+ forceEnergyImmediateneighborShortrange[compoundId * MAX_COMPOUND_PARTICLES + particleId]
			+ forceEnergyBonds[compoundId * MAX_COMPOUND_PARTICLES + particleId]
			+ pmeFE;
	}

	// Compounds
	ForceEnergy* forceEnergyFarneighborShortrange = nullptr;
	ForceEnergy* forceEnergyImmediateneighborShortrange = nullptr;
	ForceEnergy* forceEnergyBonds = nullptr;
	ForceEnergy* forceEnergiesPME = nullptr;

	// Bondgroups
	ForceEnergy* forceEnergiesBondgroups = nullptr;

	// Tinymol
	ForceEnergy* forceEnergiesCompoundinteractions = nullptr;
	ForceEnergy* forceEnergiesTinymolinteractions = nullptr;
};
