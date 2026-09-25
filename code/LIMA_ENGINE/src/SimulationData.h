#pragma once

#include <optional>
#include <memory>
#include <vector>

#include "EngineBodies.cuh"
#include "Bodies.cuh"
#include "Simulation.cuh"
#include "CudaBuffer.h"
#include "Engine.cuh"
#include "BatchLayout.cuh"

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
		return step % (int64_t{nStepsInBuffer} * loggingInterval) == 0;
	}
	static int StepsReadyToTransfer(size_t step, int loggingInterval) {
		if (loggingInterval == 0) return 0;
		const int64_t stepsSinceTransfer = step % (int64_t{nStepsInBuffer} * loggingInterval);
		return stepsSinceTransfer / loggingInterval;
	}

	__device__ static size_t GetLogIndexOfParticle(int pidInPclusters, int pcId, int64_t step,
		int loggingInterval, const int totalParticleUpperbound) {
		const int64_t steps_since_transfer = step % (int64_t{nStepsInBuffer} * loggingInterval);

		const size_t stepOffset = size_t(steps_since_transfer / loggingInterval) * totalParticleUpperbound;
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

// Host bookkeeping for one member; GPU allocations belong to EngineBatchData.
struct EngineSimulationData {
	Simulation* simulation = nullptr; // nonowning; caller outlives Engine
	SimulationDeviceData device;
	RunStatus runstatus;
	BatchRange superclusters;
	int64_t step = 0;
	int64_t stepAtLastEarlystopCheck = INT_MIN;
	size_t nLogEntriesTransferred = 0;
	bool finalized = false;
	std::vector<float> finalForcesMagnitudeSquared;
};

struct EngineBatchData {
	EngineBatchData() = default;
	~EngineBatchData();
	EngineBatchData(const EngineBatchData&) = delete;
	EngineBatchData& operator=(const EngineBatchData&) = delete;

	std::vector<EngineSimulationData> simulations;
	SimParams params;
	Int3 boxSize;
	int64_t step = 0;
	int nPclusters = 0;
	int nBondgroups = 0;
	int nParticles = 0;
	int nGridnodes = 0;
	int nSuperclusters = 0;
	int nResults = 0;
	float ewaldKappa = 0.f;

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
	CudaBuffer<IntegrationSimulationData> integrationSimulationDataDevice;

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
