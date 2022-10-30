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

	//double* getDatabuffer();
	bool testFunction();
	double* analyzeEnergy();


	Int3 timings = Int3(0, 0, 0);



	// ----- Functions used by analyzer aswell ----- //
	//__device__ Float3 computeLJForces(Box* box, Compound* compound, CompoundNeighborList* neighborlist, 
		//CompoundState* self_state, CompoundState* neighborstate_buffer, Float3* utility_buffer);


	ForceField getForcefield() { return forcefield_host; }


	~Engine() {
		//for (int i = 0; i < N_STREAMS; i++)
			//cudaStreamDestroy(stream[i]);
	}
private:
	Simulation* simulation;
	//ForceFieldMaker FFM;

	// -------------------------------------- GPU LOAD -------------------------------------- //
	void step();

	// -------------------------------------- CPU LOAD -------------------------------------- //
	NListManager* nlist_manager = nullptr;
	void handleNLISTS(Simulation* simulation, bool async = true, bool force_update=false);
//	void offloadPositionDataNLIST(Simulation* simulation);	// all at once
//	void pushNlistsToDevice();
//	static void updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection, 
//		volatile bool* finished, int* timing, bool* mutex_lock);	// thread worker, can't own engine object, thus pass ref
////	static bool neighborWithinCutoff(Float3* pos_a, Float3* pos_b);
//	static bool neighborWithinCutoff(Float3* pos_a, Float3* pos_b, float cutoff_offset);
//	/*static bool removeFromNeighborlists(NeighborList* nlist_self, NeighborList* nlist_neighbor,
//		NeighborList::NEIGHBOR_TYPE type_self, NeighborList::NEIGHBOR_TYPE type_other);*/
//	static void cullDistantNeighbors(Simulation* simulation, NListDataCollection* nlist_data_collection);
//	NListDataCollection* nlist_data_collection;

	void setDeviceConstantMemory();


	// streams every n steps
	void offloadLoggingData();
	void offloadPositionData();
	void offloadTrainData();

	void handleBoxtemp();


	//int prev_nlist_update_step = 0;
	bool updatenlists_mutexlock = 0;
	//volatile bool updated_neighborlists_ready = 0;

	// -------------------------------------- HELPERS -------------------------------------- //






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





class EngineCtrl {
public:
	EngineCtrl(){}

	static int get();
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