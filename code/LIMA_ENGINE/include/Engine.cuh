#pragma once

#include <iostream>
#include <chrono>
#include <thread>

#include "LIMA_BASE/include/Constants.cuh"
#include "LIMA_BASE/include/LimaTypes.cuh"
#include "LIMA_BASE/include/Simulation.cuh"
#include "Forcefield.cuh"


#include "Neighborlists.cuh"
#include "EngineUtils.cuh"

#include "LIMA_BASE/include/Printer.h"

#include <memory>
#include <vector>
//#include <algorithm>

template <bool em_variant>
__global__ void compoundKernel(Box* box, SimulationDevice* sim);
template __global__ void compoundKernel<true>(Box* box, SimulationDevice* sim);
template __global__ void compoundKernel<false>(Box* box, SimulationDevice* sim);

template <bool em_variant>
__global__ void solventForceKernel(Box* box, SimulationDevice* sim);
template __global__ void solventForceKernel<true>(Box* box, SimulationDevice* sim);
template __global__ void solventForceKernel<false>(Box* box, SimulationDevice* sim);

__global__ void compoundBridgeKernel(Box* box, SimulationDevice* sim);
__global__ void solventTransferKernel(Box* box, SimulationDevice* sim);


class Engine {
public:
	Engine();
	Engine(Simulation* simulation, ForceField_NB forcefield);

	// Todo: Make env run in another thread, so engine has it's own thread entirely
	// I'm sure that will help the branch predictor alot!
	template <bool em_variant> void step();



	Int3 timings = Int3(0, 0, 0);


	ForceField_NB getForcefield() { return forcefield_host; }
	void terminateSimulation();


private:


	void hostMaster();
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


	void handleBoxtemp();


	bool updatenlists_mutexlock = 0;


	// ################################# VARIABLES AND ARRAYS ################################# //

	int testval = 0;

	ForceField_NB forcefield_host;
	uint32_t step_at_last_traj_transfer = 0;
	Simulation* simulation;

	// Simulation variables
	//cudaStream_t stream[N_STREAMS];
	//dim3 gridblock_size;
	//int threads_per_gridblock;

	//bool finished = false;
	//int sim_blocks;

	//double block_dist;
	//int bpd;
	//double box_base;				// Of box (left, back, down-most), is negative!
	//double block_center_base;	// Including that edge blocks focus area is halfway outside the box
	//cudaError_t cuda_status;
	//bool critical_error = false;	// Not used yet
};

template void Engine::deviceMaster<false>();
template void Engine::deviceMaster<true>();
template void Engine::step<false>();
template void Engine::step<true>();



struct TemperaturPackage {	// kinE is for a single particle in compound, not sum of particles in said compound. Temp in [k], energy in [J]
	float temperature = 0;			// [k]
	float avg_kinE_compound = 0;	// [J/mol]
	float max_kinE_compound = 0;	// [J/mol]
	float avg_kinE_solvent = 0;		// [J/mol]
	float max_kinE_solvent = 0;		// [J/mol]
};