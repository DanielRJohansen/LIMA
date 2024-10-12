#pragma once

#include "Bodies.cuh"
#include "Simulation.cuh"
#include "EngineBodies.cuh"
#include "ChargeOcttree.cuh"


struct BoxDevice {
	BoxDevice() {}
	void CopyDataToHost(Box& boxDev) const;
	void DeleteBox();

	BoxParams boxparams;

	Compound* compounds = nullptr;
	CompoundCoords* compoundcoordsCircularQueue = nullptr;

	Solvent* solvents = nullptr;
	SolventBlocksCircularQueue* solventblockgrid_circularqueue = nullptr;

	CompoundBridgeBundleCompact* bridge_bundle = nullptr;

	// BondedParticlesLUT data - NEVER access directly, use the bpLUTHelpers namespace
	BondedParticlesLUT* bpLUTs = nullptr;

	UniformElectricField uniformElectricField;	
};


struct DatabuffersDeviceController {
	DatabuffersDeviceController(const DatabuffersDeviceController&) = delete;
	DatabuffersDeviceController(int total_particles_upperbound, int n_compounds, int loggingInterval);
	~DatabuffersDeviceController();

	static const int nStepsInBuffer = 10;

	static bool IsBufferFull(size_t step, int loggingInterval) {
		return step % (nStepsInBuffer * loggingInterval) == 0;
	}
	static int StepsReadyToTransfer(size_t step, int loggingInterval) {
		const int stepsSinceTransfer = step % (nStepsInBuffer * loggingInterval);
		return stepsSinceTransfer / loggingInterval;
	}

	__device__ static int GetLogIndexOfParticle(int particleIdLocal, int compound_id, int step,
		int loggingInterval, int totalParticleUpperbound) {
		const int steps_since_transfer = step % (nStepsInBuffer * loggingInterval);

		const int stepOffset = steps_since_transfer / loggingInterval * totalParticleUpperbound;
		const int compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
		return stepOffset + compound_offset + particleIdLocal;
	}

	float* potE_buffer = nullptr;				// For total energy summation
	Float3* traj_buffer = nullptr;				// Absolute positions [nm]
	float* vel_buffer = nullptr;				// Dont need direciton here, so could be a float
	Float3* forceBuffer = nullptr;				// [1/l N/mol] // For debug only
	const int total_particles_upperbound;
};


/// <summary>
/// All members of this will only ever exist on device. Immediately after creating this class
/// it must also be moved to device.
/// </summary>
struct SimulationDevice {
	SimulationDevice(const SimulationDevice&) = delete;

	SimulationDevice(const SimParams& params_host, Box* box_host, const DatabuffersDeviceController&);

	// Recursively free members. Use cudaFree on *this immediately after
	void deleteMembers();

	// Compounds signal where they are on a grid, handled by NLists. Used by solvents to load nearby compounds.
	CompoundGridNode* compound_grid = nullptr;

	// Compounds can see which compounds are near them
	NeighborList* compound_neighborlists = nullptr;

	// Module used to move solvents to a new block, in parallel
	SolventBlockTransfermodule* transfermodule_array = nullptr;

	SimParams* params;
	SimSignals* signals;

	BoxDevice* box;

	//ChargeOctTree* charge_octtree;
	Electrostatics::ChargeNode* chargeGrid = nullptr;
	float* chargeGridChargeSums = nullptr;

	// potE should be divided equally between all the particles in the node
	ForceAndPotential* chargeGridOutputForceAndPot = nullptr; // {Float3 force [1/l N/mol], float potE [J/mol]}

	// Databuffers, NOT owned by this class, so dont free them
	float* potE_buffer = nullptr;
	Float3 * traj_buffer = nullptr;
	float* vel_buffer = nullptr;
	Float3 * forceBuffer = nullptr;
};