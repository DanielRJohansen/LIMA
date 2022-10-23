#pragma once

#include <iostream>
#include <chrono>
#include <thread>

#include "Constants.cuh"
#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "Forcefield.cuh"


	



__global__ void compoundKernel(Box* box);
__global__ void solventForceKernel(Box* box);
__global__ void compoundBridgeKernel(Box* box);



class LIMAENG {
public:
	static void __device__ __host__ applyHyperpos(Float3* static_particle, Float3* movable_particle) {
		//#pragma unroll
		for (int i = 0; i < 3; i++) {
			*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
			*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		}	
	}

	static float __device__ __host__ calcKineticEnergy(Float3* pos1, Float3* pos2, float mass, double dt) {	// pos1/2 MUST be 2 steps apart!!!!
		LIMAENG::applyHyperpos(pos1, pos2);

		if ((*pos1 - *pos2).len() > 1) {
			//printf("KinE Dist over 1 nm!\n");
			//pos1->print('1');
			//pos2->print('2');
		}
			

		float vel = (*pos1 - *pos2).len() * (float) (0.5f / dt);
		float kinE = 0.5f * mass * vel * vel;
		return kinE;
	}

	static void __host__ genericErrorCheck(const char* text) {
		cudaDeviceSynchronize();
		cudaError_t cuda_status = cudaGetLastError();
		if (cuda_status != cudaSuccess) {
			cout << "\nCuda error code: " <<  cuda_status<< endl;
			fprintf(stderr, text);
			exit(1);
		}
	}
};




struct NListDataCollection {
	//NListDataCollection() {}
	NListDataCollection(Simulation* simulation) {
		n_compounds = simulation->n_compounds;
		n_solvents = simulation->n_solvents;
		compoundstates = new CompoundState[n_compounds];
		solvents = new Solvent[simulation->n_solvents];
		compound_neighborlists = new NeighborList[MAX_COMPOUNDS];
		solvent_neighborlists = new NeighborList[MAX_SOLVENTS];
		cudaMemcpy(compound_neighborlists, simulation->box->compound_neighborlists, sizeof(NeighborList) * n_compounds, cudaMemcpyDeviceToHost);
		cudaMemcpy(solvent_neighborlists, simulation->box->solvent_neighborlists, sizeof(NeighborList) * n_solvents, cudaMemcpyDeviceToHost);
	}
	/*
	void preparePositionData() {
		for (int i = 0; i < n_compounds; i++) {
			compound_key_positions[i] = compoundstates[i].positions[0];
		}
		for (int i = 0; i < n_solvents; i++) {
			solvent_positions[i] = solvents[i].pos;
		}
	}*/
	void preparePositionData(Compound* compounds) {
		for (int i = 0; i < n_compounds; i++) {
			compound_key_positions[i] = compoundstates[i].positions[compounds[i].key_particle_index];
		}
		for (int i = 0; i < n_solvents; i++) {
			solvent_positions[i] = solvents[i].pos;
		}
	}
	int n_compounds;
	int n_solvents;

	// I guess this is not critical but temp, needed to load pos device->host
	CompoundState* compoundstates;
	Solvent* solvents;

	Float3 compound_key_positions[MAX_COMPOUNDS];
	Float3 solvent_positions[MAX_SOLVENTS];

	// These are loaded before simulaiton start. Kept on host, and copied to device each update.
	NeighborList* compound_neighborlists;
	NeighborList* solvent_neighborlists;


};



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
	void handleNLISTS(Simulation* simulation, bool async = true, bool force_update=false);
	void offloadPositionDataNLIST(Simulation* simulation);	// all at once
	void pushNlistsToDevice();
	static void updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection, 
		volatile bool* finished, int* timing, bool* mutex_lock);	// thread worker, can't own engine object, thus pass ref
//	static bool neighborWithinCutoff(Float3* pos_a, Float3* pos_b);
	static bool neighborWithinCutoff(Float3* pos_a, Float3* pos_b, float cutoff_offset);
	/*static bool removeFromNeighborlists(NeighborList* nlist_self, NeighborList* nlist_neighbor,
		NeighborList::NEIGHBOR_TYPE type_self, NeighborList::NEIGHBOR_TYPE type_other);*/
	static void cullDistantNeighbors(Simulation* simulation, NListDataCollection* nlist_data_collection);
	NListDataCollection* nlist_data_collection;

	// streams every n steps
	void offloadLoggingData();
	void offloadPositionData();
	void offloadTrainData();

	Float3 getBoxTemperature();
	void handleBoxtemp();


	int prev_nlist_update_step = 0;
	bool updatenlists_mutexlock = 0;
	volatile bool updated_neighborlists_ready = 0;

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