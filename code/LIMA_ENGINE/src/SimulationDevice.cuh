#pragma once

#include "EngineBodies.cuh"
#include "Bodies.cuh"
#include "Simulation.cuh"


struct BoxConfig {	
	BoxConfig(Compound* compounds, uint8_t* compoundsAtomTypes, float* compoundsAtomCharges, BondedParticlesLUT* bpLUTs,
	const BoxGrid::TinymolBlockAdjacency::BlockRef* tinymolNearbyBlockIds);
	static BoxConfig Create(const Box& boxHost); // Returns a ptr to device
	void FreeMembers() const;// Free *this immediately after calling this function

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
	static BoxState* Create(const Box& boxHost); // Returns a ptr to device
	void CopyDataToHost(Box& boxDev) const;
	void FreeMembers();// Free *this immediately after calling this function

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
	BoxState* boxState, const DatabuffersDeviceController&);

	// Recursively free members. Use cudaFree on *this immediately after
	void FreeMembers();

	// Compounds signal where they are on a grid, handled by NLists. Used by solvents to load nearby compounds.
	CompoundGridNode* compound_grid = nullptr;

	// Compounds can see which compounds are near them
	NeighborList* compound_neighborlists = nullptr;

	// Module used to move solvents to a new block, in parallel
	SolventBlockTransfermodule* transfermodule_array = nullptr;

	const SimParams params;
	SimSignals* signals = nullptr;

	const BoxConfig boxConfig;
	BoxState* const boxState;
	const BoxParams boxparams;

	// potE should be divided equally between all the particles in the node
	ForceAndPotential* chargeGridOutputForceAndPot = nullptr; // {Float3 force [J/mol/nm], float potE [J/mol]}

	// Databuffers, NOT owned by this class, so dont free them
	float* potE_buffer = nullptr;
	Float3 * traj_buffer = nullptr;
	float* vel_buffer = nullptr;
	Float3 * forceBuffer = nullptr;

	// Only used in EM
	AdamState* adamState = nullptr;
};