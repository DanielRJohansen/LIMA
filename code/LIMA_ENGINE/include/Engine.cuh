#pragma once

#include <iostream>
#include <chrono>
#include <thread>

#include "LIMA_BASE/include/Constants.cuh"
#include "LIMA_BASE/include/LimaTypes.cuh"
#include "LIMA_BASE/include/Simulation.cuh"
#include "LIMA_BASE/include/Forcefield.cuh"


#include "Neighborlists.cuh"
#include "EngineUtils.cuh"


#include "LIMA_BASE/include/utilities.h"

#include <memory>
#include <vector>
//#include <algorithm>

template <bool em_variant>
__global__ void compoundKernel(SimulationDevice* sim);
template __global__ void compoundKernel<true>(SimulationDevice* sim);
template __global__ void compoundKernel<false>(SimulationDevice* sim);

template <bool em_variant>
__global__ void solventForceKernel(SimulationDevice* sim);
template __global__ void solventForceKernel<true>(SimulationDevice* sim);
template __global__ void solventForceKernel<false>(SimulationDevice* sim);

__global__ void compoundBridgeKernel(SimulationDevice* sim);
__global__ void solventTransferKernel(SimulationDevice* sim);


class Engine {
public:
	Engine(Simulation*, ForceField_NB, std::unique_ptr<LimaLogger>);

	// Todo: Make env run in another thread, so engine has it's own thread entirely
	// I'm sure that will help the branch predictor alot!
	template <bool em_variant> void step();



	Int3 timings = Int3(0, 0, 0);


	ForceField_NB getForcefield() { return forcefield_host; }
	void terminateSimulation();


private:


	template <bool em_variant> void hostMaster();
	template <bool em_variant> void deviceMaster();

	// -------------------------------------- CPU LOAD -------------------------------------- //
	std::unique_ptr<NListManager> nlist_manager;
	void setDeviceConstantMemory();


	// streams every n steps
	void offloadLoggingData(const int steps_to_transfer = STEPS_PER_LOGTRANSFER);
	void offloadTrajectory(const int steps_to_transfer = STEPS_PER_LOGTRANSFER);
	void offloadTrainData();

	// Needed to get positions before initial kernel call. Necessary in order to get positions for first NList call
	void bootstrapTrajbufferWithCoords();


	void handleBoxtemp(bool em_variant);

	std::unique_ptr<LimaLogger> m_logger;

	bool updatenlists_mutexlock = 0;


	// ################################# VARIABLES AND ARRAYS ################################# //

	int testval = 0;

	ForceField_NB forcefield_host;
	uint64_t step_at_last_traj_transfer = 0;
	Simulation* simulation;
};

template void Engine::deviceMaster<false>();
template void Engine::deviceMaster<true>();
template void Engine::hostMaster<false>();
template void Engine::hostMaster<true>();
template void Engine::step<false>();
template void Engine::step<true>();



struct TemperaturPackage {	// kinE is for a single particle in compound, not sum of particles in said compound. Temp in [k], energy in [J]
	float temperature = 0;			// [k]
	float avg_kinE_compound = 0;	// [J/mol]
	float max_kinE_compound = 0;	// [J/mol]
	float avg_kinE_solvent = 0;		// [J/mol]
	float max_kinE_solvent = 0;		// [J/mol]
};