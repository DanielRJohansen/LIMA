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
	ForceEnergyInterims(int nBondgroups, int nParticles, int nPclusters);
	void Free() const;

	__device__ ForceEnergy SumCompound(int compoundId, int particleId) const {
		ForceEnergy pmeFE = {};
		if constexpr (ENABLE_ES_LR) {
			pmeFE = forceEnergiesPME[compoundId * MAX_COMPOUND_PARTICLES + particleId];
		}

		Float3 fOld = 
			forceEnergyImmediateneighborShortrange[compoundId * MAX_COMPOUND_PARTICLES + particleId].force +
			forceEnergyFarneighborShortrange[compoundId * MAX_COMPOUND_PARTICLES + particleId].force;
		Float3 nNew = fromSuperclusters[compoundId * MAX_COMPOUND_PARTICLES + particleId].force;
		float vecErr = (fOld - nNew).len() / fOld.len();
		float magDiff = std::abs(fOld.len() - nNew.len());
		float magErr = magDiff / fOld.len();
		float threshold = 4000;

		//if (vecErr > 0.1 && magDiff > threshold) {
		//	printf("\nCompound %5d Particle %2d: relative error %.6f Old force %10.1f %10.1f %10.1f, New force %10.1f %10.1f %10.1f\n",
		//		compoundId, particleId, vecErr, fOld.x, fOld.y, fOld.z, nNew.x, nNew.y, nNew.z);
		//}
		

		return/* forceEnergyFarneighborShortrange[compoundId * MAX_COMPOUND_PARTICLES + particleId]
			+ forceEnergyImmediateneighborShortrange[compoundId * MAX_COMPOUND_PARTICLES + particleId]
			+ */forceEnergySNF[compoundId * MAX_COMPOUND_PARTICLES + particleId]
			+ fromSuperclusters[compoundId * MAX_COMPOUND_PARTICLES + particleId]
			+ pmeFE;
	}

	// These are temp, pushed into fromSuperclusters*
	ForceEnergy* nbNonlocal = nullptr;// Currently 1 per particle, i guess i want them in pclustergroups lateron
	ForceEnergy* bonded = nullptr; // TODO: Also temp, not sure how i wanna proceed here..


	// Compounds
	ForceEnergy* forceEnergyFarneighborShortrange = nullptr;
	ForceEnergy* forceEnergyImmediateneighborShortrange = nullptr;
	ForceEnergy* forceEnergySNF = nullptr;
	ForceEnergy* forceEnergiesPME = nullptr;
	ForceEnergy* fromSuperclusters = nullptr;

	// Bondgroups
	ForceEnergy* forceEnergiesBondgroups = nullptr;

	// Tinymol
	struct {
		ForceEnergy* compoundsInteractions = nullptr;
		ForceEnergy* solventsInteractions = nullptr;
		ForceEnergy* bondgroupsInteractions = nullptr;
		ForceEnergy* pmeInteraction = nullptr;
		ForceEnergy* fromSuperclusters = nullptr;
	} solvents;
};
