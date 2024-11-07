#pragma once

#include <iostream>
#include <chrono>
#include <thread>
#include <memory>

#include "Constants.h"
#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "Utilities.h"






class SimulationDevice;
class DatabuffersDeviceController;
class Thermostat;
class BoxState;
class BoxConfig;
class NeighborList;

const int cbkernel_utilitybuffer_size = sizeof(DihedralBond) * MAX_DIHEDRALBONDS_IN_COMPOUND;
constexpr int clj_utilitybuffer_bytes = sizeof(CompoundCoords); // TODO: Make obsolete and remove
static_assert(sizeof(int) * 3 * 3 * 3 <= cbkernel_utilitybuffer_size,
	"Not enough space for Electrostatics::DistributeChargesToChargegrid local offsets buffer");
static_assert(sizeof(int) * 3 * 3 * 3 * 2 <= cbkernel_utilitybuffer_size,
	"Not enough space for Electrostatics::DistributeChargesToChargegrid global offsets buffer");

struct EngineTimings {
	int compound_kernels{};
	int solvent_kernels{};
	int cpu_master{};
	int nlist{};
	int electrostatics{};

	void reset() {
		compound_kernels = 0;
		solvent_kernels = 0;
		cpu_master = 0;
		nlist = 0;
		electrostatics = 0;
	}
};

struct RunStatus {
	Float3* most_recent_positions = nullptr;
	int64_t stepForMostRecentData = 0;
	int current_step = 0;
	float current_temperature = NAN;

	//int64_t stepsSinceEnergycheck = 0;
	//float highestEnergy = 0.f; // measured in a single particle
	float greatestForce = NAN; // measured in a single particle

	bool simulation_finished = false;
	bool critical_error_occured = false;
};


class Engine {
public:
	Engine(std::unique_ptr<Simulation>, BoundaryConditionSelect, std::unique_ptr<LimaLogger>);
	~Engine();

	// Todo: Make env run in another thread, so engine has it's own thread entirely
	// I'm sure that will help the branch predictor alot! Actually, probably no.
	void step();

	/// <summary>
	/// Engine takes ownedship of sim. Noone else is allowed to access
	/// </summary>
	void runAsync(std::unique_ptr<Simulation>, RunStatus& runstatus);


	std::unique_ptr<Simulation> takeBackSim();


	EngineTimings timings{};
	volatile RunStatus runstatus;

	void terminateSimulation();

	SimulationDevice* getSimDev() { return sim_dev; }

private:


	void hostMaster();
	void deviceMaster();
	template <typename BoundaryCondition, bool emvariant, bool computePotE>
	void _deviceMaster();

	// -------------------------------------- CPU LOAD -------------------------------------- //
	void setDeviceConstantMemory();
	void verifyEngine();

	// streams every n steps
	void offloadLoggingData(const int64_t steps_to_transfer);
	void offloadTrainData();

	// Needed to get positions before initial kernel call. Necessary in order to get positions for first NList call
	void bootstrapTrajbufferWithCoords();

	void HandleEarlyStoppingInEM();
	int64_t stepAtLastEarlystopCheck = 0;

	std::unique_ptr<LimaLogger> m_logger;

	bool updatenlists_mutexlock = 0;


	// ################################# VARIABLES AND ARRAYS ################################# //

	int testval = 0;

	//ForceField_NB forcefield_host;
	uint64_t step_at_last_traj_transfer = 0;
	std::unique_ptr<Simulation> simulation = nullptr;

	SimulationDevice* sim_dev = nullptr;
	// Copies of device ptrs kept here for performance. The data array data is NOT owned here, so dont clean that up!
	std::unique_ptr<BoxState> boxStateCopy;
	std::unique_ptr<BoxConfig> boxConfigCopy;
	NeighborList* neighborlistsPtr = nullptr; // dont own data!
	//std::unique_ptr<NeighborList> neighborlistsCopy = nullptr;

	std::unique_ptr<DatabuffersDeviceController> dataBuffersDevice;

	std::unique_ptr<Thermostat> thermostat;

	const BoundaryConditionSelect bc_select;
};

 
