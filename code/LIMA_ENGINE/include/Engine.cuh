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



class Thermostat;
class CompoundGridNode;
struct CompoundQuickData;
struct ForceEnergyInterims;
class TinymolTransferModule;
struct SuperClustersControl;
struct PClusterTransfermodule;
struct PersistentCluster;
class SuperclusterStagingControl;
class TaskBuilderControl;
struct EngineSimulationData;
struct EngineBatchData;

namespace NeighborList { class Controller; }

namespace PME { class Controller; };
namespace NeighborList{struct IdAndRelshift;}




struct RunStatus {
	Float3* most_recent_positions = nullptr; // TODO: Refactor this out
	int64_t stepForMostRecentData = -1;
	int64_t current_step = 0;
	float current_temperature = NAN;
	float greatestForce = NAN; // measured in a single particle

	bool simulation_finished = false;
	bool critical_error_occured = false;
};


enum class EngineRunMode { Simulation, Interactive };

class Engine {
public:
	// Simulations are nonowning and must outlive this engine. Interactive mode
	// supports the existing single-simulation live editor without a step limit.
	explicit Engine(const std::vector<Simulation*>& simulations, EngineRunMode mode = EngineRunMode::Simulation);
	~Engine();

	void step();

	void CopySimulationToHost(size_t simulationId = 0);
	const RunStatus& GetRunStatus(size_t simulationId = 0) const;
	bool IsFinished() const;
	void terminateSimulation();

	static bool TestAlgorithms();

	// Offloads current pcluster state to another (existing) device buffer, and returns a reference to that
	// 1. This ensure that this funciton is rather quick, and the caller can continue sim immediately after this
	// kernel, and do copytohost async afterwards
	// 2. We return a reference to an existing buffer, so we wont have to allocate mem each time!
	CudaBuffer<PersistentCluster>& OffloadPclusterState(size_t simulationId = 0);
	// TODO: Make another version of the func above, that does the copy-to-host-part async, and can reuse
	// the host memory..

	// Similarly to above, returns a ref to a copy buffer
	CudaBuffer<float>& OffloadForcesMagnitudeBuffer(size_t simulationId = 0);

	// Overwrites force in IntegrationKernel during EM if present
	void SetFixedParticleMovementBuffer(const std::vector<Float3>& velocities, size_t simulationId = 0);
	void SetFixedParticleRotationBuffer(const std::vector<Rotation>& rotations, size_t simulationId = 0);
	// Is multiplied with forces in integration kernel, for partial fixing of particles
	void SetForceMask(const std::vector<Float3>& mask, size_t simulationId = 0);
	void SetElasticPositions(const std::vector<Float3>& mask, size_t simulationId = 0);

private:


	bool hostMaster();
	void deviceMaster();
	template <typename BoundaryCondition, bool emvariant, bool computePotE>
	void _deviceMaster();

	template <typename BoundaryCondition, bool emvariant>
	void SnfHandler(cudaStream_t& stream);

	// -------------------------------------- CPU LOAD -------------------------------------- //
	void verifyEngine();

	// streams every n steps
	void OffloadLoggingData(EngineSimulationData& simData);


	// Needed to get positions before initial kernel call. Necessary in order to get positions for first NList call
	void BootstrapTrajbufferWithCoords(EngineSimulationData& simData);
	void Synchronize();


	void HandleEarlyStoppingInEM(EngineSimulationData& simData);
	void RebuildActiveBatch();
	void InitializePME();
	void FinalizeSimulation(EngineSimulationData& simData);


	std::array<cudaStream_t, 5> cudaStreams{};
	cudaStream_t pmeStream = nullptr;
	EngineRunMode mode;
	// ################################# VARIABLES AND ARRAYS ################################# //

	std::unique_ptr<EngineBatchData> batch;


	// Temp
	bool MakeSuperClusterTasksGPU(cudaStream_t stream);
	void RunClustering(cudaStream_t stream, bool runPclustering = true);
	void BootstrapClustering(cudaStream_t stream);
};

 
