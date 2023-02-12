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

#include <memory>
#include <vector>
//#include <algorithm>

__global__ void compoundKernel(Box* box);
__global__ void solventForceKernel(Box* box);
__global__ void compoundBridgeKernel(Box* box);
__global__ void solventTransferKernel(Box* box);


class Engine {
public:
	Engine();
	Engine(Simulation* simulation, ForceField_NB forcefield);

	// Todo: Make env run in another thread, so engine has it's own thread entirely
	// I'm sure that will help the branch predictor alot!
	void runOnce();

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
	void offloadTrajectory(const int steps_to_transfer = STEPS_PER_LOGTRANSFER);
	void offloadTrainData();

	bool neighborlistUpdateRequired() const;

	// Needed to get positions before initial kernel call. Necessary in order to get positions for first NList call
	void bootstrapTrajbufferWithCoords();


	void handleBoxtemp();


	bool updatenlists_mutexlock = 0;


	// ################################# VARIABLES AND ARRAYS ################################# //

	int testval = 0;

	ForceField_NB forcefield_host;
	uint32_t step_at_last_traj_transfer = 0;

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
	bool critical_error = false;	// Not used yet
};

struct TemperaturPackage {	// kinE is for a single particle in compound, not sum of particles in said compound. Temp in [k], energy in [J]
	float temperature{}, avg_kinE_compound{}, max_kinE_compound{}, avg_kinE_solvent{}, max_kinE_solvent{};
};