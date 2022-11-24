#pragma once

#include <iostream>
#include <chrono>
#include <thread>

#include "Constants.cuh"
#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "Forcefield.cuh"


#include "Neighborlists.cuh"
#include "EngineUtils.cuh"

#include "Printer.h"

//#include <algorithm>

__global__ void compoundKernel(Box* box);
__global__ void solventForceKernel(Box* box);
__global__ void compoundBridgeKernel(Box* box);



class Engine {
public:
	Engine();
	Engine(Simulation* simulation, ForceField_NB forcefield);

	void deviceMaster();
	void hostMaster();
	void terminateSimulation();

	Int3 timings = Int3(0, 0, 0);


	ForceField_NB getForcefield() { return forcefield_host; }

private:
	Simulation* simulation;
	//ForceFieldMaker FFM;

	// -------------------------------------- GPU LOAD -------------------------------------- //
	void step();

	// -------------------------------------- CPU LOAD -------------------------------------- //
	NListManager* nlist_manager = nullptr;
	void handleNLISTS(Simulation* simulation, bool async = true, bool force_update=false);
	void setDeviceConstantMemory();


	// streams every n steps
	void offloadLoggingData(const int steps_to_transfer = STEPS_PER_LOGTRANSFER);
	void offloadPositionData();
	void offloadTrainData();

	void handleBoxtemp();


	bool updatenlists_mutexlock = 0;


	// ################################# VARIABLES AND ARRAYS ################################# //

	int testval = 0;

	ForceField_NB forcefield_host;


	// Simulation variables
	//cudaStream_t stream[N_STREAMS];
	dim3 gridblock_size;
	int threads_per_gridblock;

	bool finished = false;
	int sim_blocks;

	double block_dist;
	int bpd;
	double box_base;				// Of box (left, back, down-most), is negative!
	double block_center_base;	// Including that edge blocks focus area is halfway outside the box
	cudaError_t cuda_status;
};

struct TemperaturPackage {	// kinE is for a single particle in compound, not sum of particles in said compound. Temp in [k], energy in [J]
	float temperature{}, avg_kinE_compound{}, max_kinE_compound{}, avg_kinE_solvent{}, max_kinE_solvent{};
};