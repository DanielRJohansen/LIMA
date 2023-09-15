#pragma once

#include <iostream>
#include <chrono>
#include <thread>

#include "Constants.h"
#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "Forcefield.cuh"
#include "Utilities.h"

#include "Neighborlists.cuh"




#include <memory>
#include <vector>
//#include <algorithm>

__global__ void compoundKernel(SimulationDevice* sim);
__global__ void solventForceKernel(SimulationDevice* sim);

__global__ void compoundBridgeKernel(SimulationDevice* sim);
__global__ void solventTransferKernel(SimulationDevice* sim);

struct EngineTimings {
	int compound_kernels{};
	int solvent_kernels{};
	int cpu_master{};
	int nlist{};

	void reset() {
		compound_kernels = 0;
		solvent_kernels = 0;
		cpu_master = 0;
		nlist = 0;
	}
};

struct RunStatus {
	Float3* most_recent_positions = nullptr;
	int current_step = 0;
	float current_temperature = 0.f;



	bool simulation_finished = false;
	bool critical_error_occured = false;
};

class Engine {
public:
	Engine(std::unique_ptr<Simulation>, ForceField_NB, std::unique_ptr<LimaLogger>);
	~Engine();

	// Todo: Make env run in another thread, so engine has it's own thread entirely
	// I'm sure that will help the branch predictor alot!
	void step();

	/// <summary>
	/// Engine takes ownedship of sim. Noone else is allowed to access
	/// </summary>
	void runAsync(std::unique_ptr<Simulation>, RunStatus& runstatus);


	std::unique_ptr<Simulation> takeBackSim() { return std::move(simulation); }


	EngineTimings timings{};
	RunStatus runstatus;

	ForceField_NB getForcefield() { return forcefield_host; }
	void terminateSimulation();


private:


	void hostMaster();
	void deviceMaster();

	// -------------------------------------- CPU LOAD -------------------------------------- //
	std::unique_ptr<NListManager> nlist_manager;
	void setDeviceConstantMemory();


	// streams every n steps
	void offloadLoggingData(const int steps_to_transfer = STEPS_PER_LOGTRANSFER);
	void offloadTrajectory(const int steps_to_transfer = STEPS_PER_LOGTRANSFER);
	void offloadTrainData();

	// Needed to get positions before initial kernel call. Necessary in order to get positions for first NList call
	void bootstrapTrajbufferWithCoords();


	void handleBoxtemp();

	std::unique_ptr<LimaLogger> m_logger;

	bool updatenlists_mutexlock = 0;


	// ################################# VARIABLES AND ARRAYS ################################# //

	int testval = 0;

	ForceField_NB forcefield_host;
	uint64_t step_at_last_traj_transfer = 0;
	std::unique_ptr<Simulation> simulation;
};



struct TemperaturPackage {	// kinE is for a single particle in compound, not sum of particles in said compound. Temp in [k], energy in [J]
	float temperature = 0;			// [k]
	float avg_kinE_compound = 0;	// [J/mol]
	float max_kinE_compound = 0;	// [J/mol]
	float avg_kinE_solvent = 0;		// [J/mol]
	float max_kinE_solvent = 0;		// [J/mol]
};
