#pragma once

#include <optional>
#include <memory>
#include <vector>

#include "EngineBodies.cuh"
#include "Bodies.cuh"
#include "Simulation.cuh"
#include "CudaBuffer.h"

class Thermostat;
struct SuperClustersControl;
struct PClusterTransfermodule;
class SuperclusterStagingControl;
class TaskBuilderControl;

namespace PME { class Controller; }



struct BoxState {
	BoxState() {};
	BoxState(PersistentclusterInterimState*);
	static BoxState Create(const Box& boxHost);
	void CopyDataToHost(Box& boxDev) const;
	void FreeMembers() const;

	PersistentclusterInterimState* pclusterInterimStates = nullptr;
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
		if (loggingInterval == 0)
			return false;
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
		const int pclusterOffset = pcId * PersistentCluster::maxParticles;
		return stepOffset + pclusterOffset + pidInPclusters;
	}

	float* potE_buffer = nullptr;				// For total energy summation
	Float3* traj_buffer = nullptr;				// Absolute positions [nm]
	float* vel_buffer = nullptr;				// Dont need direciton here, so could be a float
	Float3* forceBuffer = nullptr;				// [J/mol/nm] // For debug only

	const int nParticlesUpperbound;
};


struct ForceEnergyInterims {
	ForceEnergyInterims(int nBondgroupParticles, int nParticles, int nPclusters);
	void Free() const;

	// These are temp, pushed into fromSuperclusters*
	ForceEnergy* bonded = nullptr; // TODO: Also temp, not sure how i wanna proceed here..
	ForceEnergy* pme = nullptr;

	// Pushed to from pclusters
	ForceEnergy* snf = nullptr;

	// Bondgroups
	ForceEnergy* forceEnergiesBondgroups = nullptr;
};

// All state owned by one simulation. Engine deliberately owns only execution
// resources (streams) and orchestration; a future batch executor can own many
// of these objects and concatenate the buffer categories.
struct EngineSimulationData {
	EngineSimulationData(Simulation* simulation);
	~EngineSimulationData();

	EngineSimulationData(const EngineSimulationData&) = delete;
	EngineSimulationData& operator=(const EngineSimulationData&) = delete;

	uint64_t step_at_last_traj_transfer = 0;
	int64_t stepAtLastEarlystopCheck = INT_MIN;
	Simulation* simulation = nullptr; // nonowning

	int nSuperclusters = 0;
	float ewaldKappa = 0.f;
	float thermostatScalar = 1.f;
	size_t nResults = 0;

	CudaBuffer<PersistentCluster> pClusterDevice;
	CudaBuffer<PersistentClusterMeta> pClusterMetaDevice;
	CudaBuffer<float> forcesMagnitudeSquareDevice;
	BoxState boxState;
	AdamState* adamState = nullptr;

	CudaBuffer<BondGroup> bondgroupDescriptors;
	CudaBuffer<BondGroup::ParticleRef> bondgroupParticles;
	CudaBuffer<SingleBond> bondgroupSinglebonds;
	CudaBuffer<PairBond> bondgroupPairbonds;
	CudaBuffer<AngleUreyBradleyBond> bondgroupAnglebonds;
	CudaBuffer<DihedralBond> bondgroupDihedralbonds;
	CudaBuffer<ImproperDihedralBond> bondgroupImproperdihedralbonds;

	CudaBuffer<ScScTask> scscTasksDevice;
	CudaBuffer<int> idsOfQuerySuperclustersDevice;
	CudaBuffer<int> resultIndicesDevice;
	CudaBuffer<BoolMatrix16x16> noInteractionMatricesDevice;
	CudaBuffer<SCResult> scResultsDevice;

	std::unique_ptr<SuperClustersControl> superClustersControl;
	std::unique_ptr<PClusterTransfermodule> pclusterTransfermodule;
	std::unique_ptr<SuperclusterStagingControl> superclusterStagingControl;
	std::unique_ptr<TaskBuilderControl> taskbuilderControl;
	std::unique_ptr<PME::Controller> pmeController;
	std::unique_ptr<DatabuffersDeviceController> dataBuffersDevice;
	std::unique_ptr<Thermostat> thermostat;
	std::unique_ptr<ForceEnergyInterims> forceEnergyInterims;

	std::vector<ParticlesBondedToParticle> particlesBondedToParticle;
	std::vector<PclustersBondedToPcluster> pclustersBondedToPcluster;

	CudaBuffer<PersistentCluster> pdataCopyBuffer;
	CudaBuffer<float> forcesMagnitudeCopyBuffer;
	std::optional<CudaBuffer<Float3>> fixedParticleMovementBuffer;
	std::optional<CudaBuffer<Rotation>> fixedParticleRotationBuffer;
	std::optional<CudaBuffer<Float3>> forceMaskBuffer;
	std::optional<CudaBuffer<Float3>> elasticPositionsBuffer;
};
