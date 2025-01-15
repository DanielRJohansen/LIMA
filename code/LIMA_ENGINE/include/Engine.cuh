#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"

#include "Constants.h"
#include "Utilities.h"

#include <iostream>
#include <chrono>
#include <thread>
#include <memory>




class SimulationDevice;
class DatabuffersDeviceController;
class Thermostat;
class BoxState;
class BoxConfig;
class NeighborList;
class CompoundGridNode;

namespace PME { class Controller; };

struct CompoundForceEnergyInterims {
	CompoundForceEnergyInterims(int nCompounds);
	void Free() const;

	__device__ ForceEnergy Sum(int compoundId, int particleId) const {
		return forceEnergyFarneighborShortrange[compoundId * MAX_COMPOUND_PARTICLES + particleId]
			+ forceEnergyImmediateneighborShortrange[compoundId * MAX_COMPOUND_PARTICLES + particleId]
			+ forceEnergyBonds[compoundId * MAX_COMPOUND_PARTICLES + particleId];
	}

	ForceEnergy* forceEnergyFarneighborShortrange;
	ForceEnergy* forceEnergyImmediateneighborShortrange;
	ForceEnergy* forceEnergyBonds;
};

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

	void HandleEarlyStoppingInEM();
	int64_t stepAtLastEarlystopCheck = 0;

	std::unique_ptr<LimaLogger> m_logger;

	bool updatenlists_mutexlock = 0;

	std::array<cudaStream_t, 5> cudaStreams;
	cudaStream_t pmeStream;
	// ################################# VARIABLES AND ARRAYS ################################# //

	int testval = 0;

	//ForceField_NB forcefield_host;
	uint64_t step_at_last_traj_transfer = 0;
	std::unique_ptr<Simulation> simulation;

	// Owned
	SimulationDevice* sim_dev = nullptr;
	// These are owned, but this is temporary place to store them
	ForceEnergy* forceEnergiesBondgroups = nullptr;
	ForceEnergy* forceEnergiesPME = nullptr;
	BondGroup* bondgroups = nullptr;
	ForceField_NB::ParticleParameters* compoundLjParameters = nullptr;

	// Copies of device ptrs kept here for performance. The data array data is NOT owned here, so dont clean that up!
	std::unique_ptr<BoxState> boxStateCopy;
	std::unique_ptr<BoxConfig> boxConfigCopy;
	NeighborList* neighborlistsPtr = nullptr; // dont own data!
	CompoundGridNode* compoundgridPtr = nullptr;// dont own data!




	std::unique_ptr<PME::Controller> pmeController;


	CompoundForceEnergyInterims compoundForceEnergyInterims;


	//std::unique_ptr<NeighborList> neighborlistsCopy = nullptr;

	std::unique_ptr<DatabuffersDeviceController> dataBuffersDevice;

	std::unique_ptr<Thermostat> thermostat;

	const BoundaryConditionSelect bc_select;
};

 
