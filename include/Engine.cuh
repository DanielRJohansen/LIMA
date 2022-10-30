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

__global__ void compoundKernel(Box* box);
__global__ void solventForceKernel(Box* box);
__global__ void compoundBridgeKernel(Box* box);



class Engine {
public:
	Engine();
	Engine(Simulation* simulation, ForceField forcefield);

	void deviceMaster();
	void hostMaster();


	Int3 timings = Int3(0, 0, 0);


	ForceField getForcefield() { return forcefield_host; }

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
	void offloadLoggingData();
	void offloadPositionData();
	void offloadTrainData();

	void handleBoxtemp();


	bool updatenlists_mutexlock = 0;


	// ################################# VARIABLES AND ARRAYS ################################# //

	int testval = 0;

	ForceField forcefield_host;


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



// FUNCTIONS IN THIS CLASS CAN BE CALLED FROM OTHER FILES
/*
class DeviceFunctionCollection {
public:
	__device__ __host__ DeviceFunctionCollection(){}

	__device__ static double getAngle(Float3 v1, Float3 v2);
	__device__ static void determineMoleculerHyperposOffset(Float3* utility_buffer, Float3* compund_center, 
		Float3* neighbor_center);	

	__device__ static Float3 calcLJForce(Float3* pos0, Float3* pos1);

	__device__ static Float3 calcPairbondForce(Float3* self_pos, Float3* other_pos, double* reference_dist);
	__device__ static Float3 calcAngleForce(CompoundState* statebuffer, AngleBond* anglebond);

	__device__ static Float3 computeLJForces(Box* box, Compound* compound, CompoundNeighborList* neighborlist,
		CompoundState* self_state, CompoundState* neighborstate_buffer, Float3* utility_buffer);
	__device__ static Float3 computePairbondForces(Compound* compound, CompoundState* self_state);
	__device__ static Float3 computeAnglebondForces(Compound* compound, CompoundState* self_state);

};
*/