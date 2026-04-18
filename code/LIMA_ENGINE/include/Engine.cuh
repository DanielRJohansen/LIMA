#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "CudaBuffer.h"

#include "Constants.h"
#include "Utilities.h"

#include <iostream>
#include <memory>
#include <thread>
#include <optional>



class SimulationDevice;
class DatabuffersDeviceController;
class Thermostat;
class BoxState;
class BoxConfig;
class CompoundGridNode;
struct CompoundQuickData;
struct ForceEnergyInterims;
class TinymolTransferModule;
struct SuperClustersControl;
struct PClusterTransfermodule;
struct PersistentCluster;
class SuperclusterStagingControl;
class TaskBuilderControl;

namespace NeighborList { class Controller; }

namespace PME { class Controller; };
namespace NeighborList{struct IdAndRelshift;}


struct RunStatus {
	Float3* most_recent_positions = nullptr; // TODO: Refactor this out
	int64_t stepForMostRecentData = -1;
	int current_step = 0;
	float current_temperature = NAN;
	float greatestForce = NAN; // measured in a single particle

	bool simulation_finished = false;
	bool critical_error_occured = false;
};


class Engine {
public:
	Engine(Simulation*, BoundaryConditionSelect, std::unique_ptr<LimaLogger>);
	~Engine();

	void step();

	/// <summary>
	/// Engine takes ownedship of sim. Noone else is allowed to access
	/// </summary>
	void runAsync(std::unique_ptr<Simulation>, RunStatus& runstatus);


	void CopySimulationToHost();


	volatile RunStatus runstatus;

	void terminateSimulation();

	SimulationDevice* getSimDev() { return sim_dev; }

	static bool TestAlgorithms();

	// Offloads current pcluster state to another (existing) device buffer, and returns a reference to that
	// 1. This ensure that this funciton is rather quick, and the caller can continue sim immediately after this
	// kernel, and do copytohost async afterwards
	// 2. We return a reference to an existing buffer, so we wont have to allocate mem each time!
	CudaBuffer<PersistentCluster>& OffloadPclusterState();
	// TODO: Make another version of the func above, that does the copy-to-host-part async, and can reuse
	// the host memory..

	// Overwrites force in IntegrationKernel during EM if present
	void SetFixedParticleMovementBuffer(const std::vector<Float3>& velocities);
	void SetFixedParticleRotationBuffer(const std::vector<Rotation>& rotations);
	// Is multiplied with forces in integration kernel, for partial fixing of particles
	void SetForceMask(const std::vector<Float3>& mask); 

private:


	void hostMaster();
	void deviceMaster();
	template <typename BoundaryCondition, bool emvariant, bool computePotE>
	void _deviceMaster();

	template <typename BoundaryCondition, bool emvariant>
	void SnfHandler(cudaStream_t& stream);

	// -------------------------------------- CPU LOAD -------------------------------------- //
	void setDeviceConstantMemory();
	void verifyEngine();

	// streams every n steps
	void offloadLoggingData(const int64_t steps_to_transfer);
	void offloadTrainData();

	// Needed to get positions before initial kernel call. Necessary in order to get positions for first NList call
	void bootstrapTrajbufferWithCoords();

	void BootstrapSolventblockDistributeFromDensity();

	void HandleEarlyStoppingInEM();
	int64_t stepAtLastEarlystopCheck = 0;

	std::unique_ptr<LimaLogger> m_logger;

	std::array<cudaStream_t, 5> cudaStreams;
	cudaStream_t pmeStream;
	// ################################# VARIABLES AND ARRAYS ################################# //

	uint64_t step_at_last_traj_transfer = 0;
	Simulation* simulation;

	// Owned
	SimulationDevice* sim_dev = nullptr;
	BondGroup* bondgroups = nullptr;

	//SuperClusterControl// TODO: Handle lifetimes!
	std::unique_ptr<SuperClustersControl> superClustersControl;
	std::unique_ptr<PClusterTransfermodule> pclusterTransfermodule;

	size_t nTasks = 0;
	int nSuperclusters = 0;
	PersistentCluster* pClusterDevice = nullptr; // TODO: Handle lifetime somethwere
	PersistentClusterMeta* pClusterMetaDevice = nullptr;
	size_t nResults = 0;

	CudaBuffer<ScScTask> scscTasksDevice;
	CudaBuffer<BoolMatrix16x16> noInteractionMatricesDevice;
	CudaBuffer<SCResult> scResultsDevice;

	std::unique_ptr<SuperclusterStagingControl> superclusterStagingControl;
	std::unique_ptr<TaskBuilderControl> taskbuilderControl;

	std::vector<ParticlesBondedToParticle> particlesBondedToParticle;
	std::vector<PclustersBondedToPcluster> pclustersBondedToPcluster;

	// Copies of device ptrs kept here for performance. The data array data is NOT owned here, so dont clean that up!
	std::unique_ptr<BoxState> boxStateCopy;
	std::unique_ptr<BoxConfig> boxConfigCopy;

	uint8_t* nParticlesInCompoundsBufferPtr = nullptr;// dont own data!

	std::unique_ptr<PME::Controller> pmeController;
	std::unique_ptr<DatabuffersDeviceController> dataBuffersDevice;
	std::unique_ptr<Thermostat> thermostat;
	std::unique_ptr<ForceEnergyInterims> forceEnergyInterims;
	std::unique_ptr<NeighborList::Controller> nlistController;

	//CudaBuffer<ForceEnergy> nbGatherForceenergy;

	const BoundaryConditionSelect bc_select;

	// Available to be copied to, while sim is running
	CudaBuffer<PersistentCluster> pdataCopyBuffer; 

	// For EM only, overwrites forces in integration kernel. 
	std::optional<CudaBuffer<Float3>> fixedParticleMovementBuffer; 
	std::optional<CudaBuffer<Rotation>> fixedParticleRotationBuffer;
	std::optional<CudaBuffer<Float3>> forceMaskBuffer;	// Multiplied with forces in integration kernel, for partial fixing of particles


	// Temp
	bool MakeSuperClusterTasksCPU();
	bool MakeSuperClusterTasksGPU();
	void RunClustering(bool runPclustering = true);
	void BootstrapClustering();
};

 
