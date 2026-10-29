#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "CudaBuffer.h"

#include "Constants.h"

#include <memory>
#include <optional>



class Thermostat;
class DeviceThermostat;
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
	float current_pressure = NAN;	// [bar]
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
	// Runs a step, then the following steps that need no host work (see IsPlainReplayStep) as relaunches of the step graph,
	// queued in a tight loop. Returns the number of steps run, at most maxStepsPerCall
	int StepUntilHostWork();

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
	template <typename BoundaryCondition, bool emvariant, bool computePotE, bool computeVirial>
	void _deviceMaster(bool measureTemperature);

	template <typename BoundaryCondition, bool emvariant>
	void SnfHandler(cudaStream_t& stream, const ForceAccumulator& forceAcc);

	// -------------------------------------- CPU LOAD -------------------------------------- //
	void verifyEngine();

	// streams every n steps
	void OffloadLoggingData(EngineSimulationData& simData);
	// Logging data is drained without stalling the steps: the ring is snapshotted on the GPU, copied to pinned memory on
	// logCopyStream while the simulation continues, and unpacked into the host buffers by a thread. One drain per simulation index
	struct LogDrain;
	std::vector<std::unique_ptr<LogDrain>> logDrains;
	cudaStream_t logCopyStream = nullptr;
	// Waits until every drained log entry is in the host buffers
	void JoinLogDrains();
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
	// Retired members are finalized at once but stay in the batch until the next nlist update step, where the survivors
	// are re-clustered anyway. Rebuilding on any other step changes their neighbor lists and so their trajectories.
	bool rebuildPending = false;
	void InitializePME();
	void UploadNbChunks();
	void FinalizeSimulation(EngineSimulationData& simData);


	std::array<cudaStream_t, 5> cudaStreams{};
	cudaStream_t pmeStream = nullptr;
	// Used to make cudaStreams[0] wait for the other streams on the GPU, without a host roundtrip. [0] is for pmeStream, [i] for cudaStreams[i]
	std::array<cudaEvent_t, 5> streamJoinEvents{};
	void JoinStreamsIntoMainStream();
	// Recorded on cudaStreams[0] at the start of each step, so the other streams wait for the previous step on the GPU
	cudaEvent_t stepStartEvent = nullptr;
	void ForkStreamsFromMainStream();
	// MD steps are captured into a graph and launched as one unit, so the GPU sees the whole step at once rather than as launches
	// the driver submits in batches. Captured every step and applied to the instantiated graph with cudaGraphExecUpdate, so changed
	// arguments and buffers need no tracking. Anything launched in a captured step must be capturable: no syncs, no allocations,
	// and only on the forked streams. Ignored on HIP
	static constexpr bool useStepGraphs = true;
	// [logData][computeVirial][measureTemperature], the last for the DeviceThermostat's coupling steps
	std::array<std::array<std::array<cudaGraphExec_t, 2>, 2>, 2> stepGraphs{};
	void BeginStepCapture();
	void LaunchCapturedStep(bool logData, bool computeVirial, bool measureTemperature);
	// Steps without logging or virial relaunch stepGraphs[0][0][measureTemperature] without capturing while it is valid: the host then queues a step in
	// a few us rather than the ~30 us of capturing and updating it, and stays ahead of the GPU on small systems. Its arguments change
	// only on a clustering, a task build, an upload of IntegrationSimulationData or a changed live-edit buffer, which invalidate it.
	// The step index the kernels read is advanced on the GPU, see stepCounterDevice. Not in Interactive mode, where params change freely
	std::array<bool, 2> replayableStepGraph{};	// [measureTemperature]
	void InvalidateStepGraph() { replayableStepGraph = {}; }
	// The current step's index, advanced by the step's first kernel. Read by the integration kernels for error reports
	int64_t* stepCounterDevice = nullptr;
	void LaunchStepGraph(cudaGraphExec_t exec, bool captured);
	// Bounds the time between the caller's status and display updates
	static constexpr int maxStepsPerCall = 100;
	void FullStep();
	// Whether the next step would only replay the step graph: no logging, virial, measurement, render, log drain, rebuild or
	// thermostat reset, and no simulation finishing. Must match the conditions of FullStep, deviceMaster and hostMaster
	bool IsPlainReplayStep() const;
	void AdvanceStepCounters();
	// Set after the host waited on the GPU, see LaunchStepGraph
	bool flushNextStepLaunch = false;
	// Must be called after changing any field that IntegrationSimulationData is built from
	void UploadIntegrationSimulationData();
	// The temperature measurement and thermostat run on the GPU in MD without the barostat, see DeviceThermostat. Else on the host, in hostMaster
	std::unique_ptr<DeviceThermostat> deviceThermostat;
	bool UsesDeviceThermostat() const;
	// Whether the step batch->step starts with the DeviceThermostat's measurement and coupling
	bool DeviceThermostatStep() const;
	// Must be called whenever the batch is packed
	void UploadDeviceThermostatSimulations();
	void HarvestTemperatures(EngineSimulationData& sim);
	bool ApplyBarostat(float coordinateScale);
	EngineRunMode mode;
	// ################################# VARIABLES AND ARRAYS ################################# //

	std::unique_ptr<EngineBatchData> batch;
	std::vector<RenderDataPipe*> renderDataPipes;
	static constexpr int StepsPerRender = 20;


	// Temp
	bool MakeSuperClusterTasksGPU(cudaStream_t stream);
	void FindSuperclusterNeighbors(cudaStream_t stream, float listRadius, bool allQuarters);
	void MakeNbTasksMD(cudaStream_t stream);
	// Throws if the last MD task build's neighbor search overflowed and it was not checked yet. The GPU must have finished that build
	void CheckDeferredNeighborOverflow();
	void MakeNbTasksEM(cudaStream_t stream);
	bool tasksBuiltForEm = false; // The MD and EM nonbonded kernels use different tasks
	void RunClustering(cudaStream_t stream, bool runPclustering = true);
	void BootstrapClustering(cudaStream_t stream);
	// MD keeps the integration states in supercluster slots, see ParticleIntegrationState
	void LoadIntegrationStates(cudaStream_t stream);
	void StoreIntegrationStates(cudaStream_t stream);
};

 
