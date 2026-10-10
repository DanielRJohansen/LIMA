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
#include "EnergyMinimizationTypes.h"
#include "CompactBondTable.h"

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
	size_t nLogEntriesTransferred = 0;
	bool finalized = false;
	std::vector<float> finalForcesMagnitudeSquared;
	EM::Preconditioner emPreconditioner; // Empty unless the engine uses energy minimization
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
	int nBondgroupParticles = 0;
	int nParticles = 0;
	int nGridnodes = 0;
	int nSuperclusters = 0;
	int nResults = 0;
	float ewaldKappa = 0.f;

	CudaBuffer<PersistentCluster> pClusterDevice;
	CudaBuffer<PersistentClusterMeta> pClusterMetaDevice;
	CudaBuffer<float> forcesMagnitudeSquareDevice;
	BoxState boxState;

	EM::Config emConfig;
	CudaBuffer<EM::ParticleState> emParticles;		// Indexed by pcluster slot
	CudaBuffer<Float3> emForces;					// Indexed by pcluster slot [J/mol/nm]
	CudaBuffer<Float3> emPreconditionedForce;		// Indexed by pcluster slot
	CudaBuffer<float> emInverseStiffness;			// Indexed by pcluster slot
	CudaBuffer<uint8_t> emWholeMolecule;			// Indexed by pcluster
	CudaBuffer<EM::Sums> emBlockSums;				// Indexed by simulation * nBlocks + block
	CudaBuffer<unsigned int> emBlocksDone;			// Indexed by simulation
	CudaBuffer<EM::SimState> emStates;				// Indexed by simulation
	std::vector<EM::SimState> emStatesHost;

	CudaBuffer<BondGroup> bondgroupDescriptors;
	CudaBuffer<BondGroup::ParticleRef> bondgroupParticles;
	CompactBondTable<SingleBond> bondgroupSinglebonds;
	CompactBondTable<PairBond> bondgroupPairbonds;
	CompactBondTable<AngleUreyBradleyBond> bondgroupAnglebonds;
	CompactBondTable<DihedralBond> bondgroupDihedralbonds;
	CompactBondTable<ImproperDihedralBond> bondgroupImproperdihedralbonds;

	CudaBuffer<QuarterEntryTask> quarterEntryTasksDevice;	// See NbNonlocalKernel
	CudaBuffer<QuarterEntry> quarterEntriesDevice;
	CudaBuffer<int> quarterEntryResultIndicesDevice;		// EM only, the SCResult each entry's j forces are stored in
	int nQuarterEntries = 0;
	CudaBuffer<SCResult> scResultsDevice;					// EM only
	CudaBuffer<unsigned long long> forceAccumulatorDevice;	// MD only: four atomic planes, then four primary bonded planes, each nSuperclusters*16
	CudaBuffer<int> pclusterParticleSlots;					// Slot (scId * 16 + index) of each pcluster particle in the current superclusters, -1 if none
	CudaBuffer<ulonglong4> extraBondForceResults; // Only secondary bondgroup appearances
	CudaBuffer<int> bondgroupExtraResultIndices; // -1 for primary appearances, otherwise an extra result index
	CudaBuffer<int> slotExtraBondReferences; // Three planes of secondary indices, rebuilt with superclusters
	CudaBuffer<int> bondgroupParticleSlots;					// Slot of each of bondgroupParticles, so BondgroupsKernel finds them directly
	CudaBuffer<ParticleIntegrationState> integrationStates;	// Indexed by supercluster slot
	bool integrationStatesLoaded = false;					// integrationStates hold the current state, which pcluster states may lag
	bool forceMagnitudesInStates = false;					// The last step left forcesMagnitudeSquareDevice to StoreIntegrationStates
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
