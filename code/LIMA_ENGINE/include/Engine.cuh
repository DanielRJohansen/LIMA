#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "CudaBuffer.h"

#include "Constants.h"

#include <memory>
#include <optional>



class Thermostat;
struct ForceEnergyInterims;
struct SuperClustersControl;
struct ForceAccumulator;
struct PClusterTransfermodule;
struct PersistentCluster;
class SuperclusterStagingControl;
class TaskBuilderControl;
struct EngineSimulationData;
struct EngineBatchData;
class RenderDataPipe;

namespace NeighborList { class Controller; }

namespace PME { class Controller; };
namespace NeighborList{struct IdAndRelshift;}




struct RunStatus {
	int64_t current_step = 0;
	float current_temperature = NAN;
	float greatestForce = NAN; // measured in a single particle

	bool simulation_finished = false;
	bool critical_error_occured = false;
};


enum class EngineRunMode { Simulation, Interactive };

// Direct kernel fixtures used by the isolated limit tests. Never call DeviceFailure in a live engine process.
enum class EngineLimitProbe { ClusterTransfer, ClusterOccupancy, ChargeBlock, DeviceFailure };

class Engine {
public:
	// Simulations are nonowning and must outlive this engine. Interactive mode
	// supports the existing single-simulation live editor without a step limit.
	explicit Engine(const std::vector<Simulation*>& simulations, EngineRunMode mode = EngineRunMode::Simulation,
		const std::vector<RenderDataPipe*>& renderDataPipes = {});
	~Engine();

	void step();

	void CopySimulationToHost(size_t simulationId = 0);
	const RunStatus& GetRunStatus(size_t simulationId = 0) const;
	void StopSimulation(size_t simulationId = 0);
	bool IsFinished() const;
	void terminateSimulation();

	static bool TestAlgorithms();
	static void TestLimit(EngineLimitProbe probe, int count = 0);

	// Optional. Makes the energy-minimization preconditioner of an EM simulation before an Engine is constructed for it,
	// so this host work does not delay the GPU. The simulation must not be modified before it is run
	static void PrepareEnergyMinimization(Simulation& simulation);

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
	void SnfHandler(cudaStream_t& stream, const ForceAccumulator& forceAcc);

	// -------------------------------------- CPU LOAD -------------------------------------- //
	void verifyEngine();

	// streams every n steps
	void OffloadLoggingData(EngineSimulationData& simData);
	void PublishRenderData();
	void PublishRenderData(size_t simulationIndex);


	// Needed to get positions before initial kernel call. Necessary in order to get positions for first NList call
	void BootstrapTrajbufferWithCoords(EngineSimulationData& simData);
	void Synchronize();
	void StopRenderDataPipes();


	void HandleEarlyStoppingInEM(EngineSimulationData& simData);
	template <typename BoundaryCondition>
	void UpdateEnergyMinimization(Float3 boxSize);
	void ResetEnergyMinimization();
	bool UsesEnergyMinimization() const;
	void UploadEnergyMinimizationPreconditioner();
	void RebuildActiveBatch();
	void InitializePME();
	void FinalizeSimulation(EngineSimulationData& simData);


	std::array<cudaStream_t, 5> cudaStreams{};
	cudaStream_t pmeStream = nullptr;
	// Used to make cudaStreams[0] wait for the other streams on the GPU, without a host roundtrip. [0] is for pmeStream, [i] for cudaStreams[i]
	std::array<cudaEvent_t, 5> streamJoinEvents{};
	void JoinStreamsIntoMainStream();
	// Recorded on cudaStreams[0] at the start of each step, so the other streams wait for the previous step on the GPU
	cudaEvent_t stepStartEvent = nullptr;
	void ForkStreamsFromMainStream();
	// Must be called after changing any field that IntegrationSimulationData is built from
	void UploadIntegrationSimulationData();
	EngineRunMode mode;
	// ################################# VARIABLES AND ARRAYS ################################# //

	std::unique_ptr<EngineBatchData> batch;
	std::vector<RenderDataPipe*> renderDataPipes;
	static constexpr int StepsPerRender = 20;


	// Temp
	bool MakeSuperClusterTasksGPU(cudaStream_t stream);
	void FindSuperclusterNeighbors(cudaStream_t stream, float listRadius, bool allQuarters);
	void MakeNbTasksMD(cudaStream_t stream);
	void MakeNbTasksEM(cudaStream_t stream);
	bool tasksBuiltForEm = false; // The MD and EM nonbonded kernels use different tasks
	void RunClustering(cudaStream_t stream, bool runPclustering = true);
	void BootstrapClustering(cudaStream_t stream);
	// MD keeps the integration states in supercluster slots, see ParticleIntegrationState
	void LoadIntegrationStates(cudaStream_t stream);
	void StoreIntegrationStates(cudaStream_t stream);
};

 
