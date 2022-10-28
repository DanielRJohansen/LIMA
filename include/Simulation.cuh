#pragma once

//#include "QuantomTypes.cuh"
#include "Constants.cuh"
#include "Bodies.cuh"





#ifndef __linux__
//const string MOL_FOLDER = "C:\\PROJECTS\\Quantom\\molecules\\t4lys_full\\";
const std::string MOL_FOLDER = "C:\\PROJECTS\\Quantom\\Simulation\\Molecule\\";
const std::string OUT_DIR = "C:\\PROJECTS\\Quantom\\Simulation\\";
#else
//const string MOL_FOLDER = "../Compounds/t4lys/";
const string MOL_FOLDER = "../../Simulation/Molecule/";
//const string OUT_DIR = "/home/lima/Desktop/LIMA";
const string OUT_DIR = "../../Simulation/";
#endif




//const int BLOCKS_PER_SOLVENTKERNEL = ceil((float)N_SOLVATE_MOLECULES / (float)THREADS_PER_SOLVENTBLOCK);




constexpr double SOLVENT_MASS = 18.01528f * 1e-3;	// kg/mol
//constexpr double SOLVENT_MASS = 12.0107 * 1e-3;	// kg/mol
//constexpr double COMPOUNDPARTICLE_MASS = 12.0107 * 1e-3;

// This goes on Device
class Box {
public:
	Box() {}
	Compound* compounds = nullptr;

	uint32_t n_compounds = 0;
	uint16_t n_solvents = 0;
	uint32_t total_particles_upperbound = 0;


	// These are shared for all compounds, MUST be allocated before adding any compounds to box, so not in moveToDevice //
	CompoundState* compound_state_array = nullptr;
	CompoundState* compound_state_array_next = nullptr;

	NeighborList* compound_neighborlists = nullptr;
	NeighborList* solvent_neighborlists = nullptr;
	//------------------------------------//

	Solvent* solvents = nullptr;
	Solvent* solvents_next = nullptr;

	CompoundBridgeBundleCompact* bridge_bundle = nullptr;

	BondedParticlesLUTManager* bonded_particles_lut_manager = nullptr;

	uint32_t step = 0;
	double dt = 0;
	bool critical_error_encountered = 0;

	float* potE_buffer = nullptr;		// For total energy summation
	Float3* traj_buffer = nullptr;

	float* outdata = nullptr;			// Temp, for longging values to whatever
	Float3* data_GAN = nullptr;			// Only works if theres 1 compounds right now.


	float thermostat_scalar = 1.f;


	void moveToDevice() {	// Loses pointer to RAM location!
		int bytes_total = sizeof(Compound) * n_compounds
			+ sizeof(Solvent) * MAX_SOLVENTS * 2
			+ sizeof(CompoundState) * MAX_COMPOUNDS * 2
			+ sizeof(NeighborList) * (MAX_SOLVENTS + MAX_COMPOUNDS);
		printf("BOX: moving %.2f MB to device\n", (float)bytes_total * 1e-6);

		compounds = genericMoveToDevice(compounds, n_compounds);
		solvents = genericMoveToDevice(solvents, MAX_SOLVENTS);
		bridge_bundle = genericMoveToDevice(bridge_bundle, 1);
		compound_state_array = genericMoveToDevice(compound_state_array, MAX_COMPOUNDS);


		compound_neighborlists = genericMoveToDevice(compound_neighborlists, MAX_COMPOUNDS);
		solvent_neighborlists = genericMoveToDevice(solvent_neighborlists, MAX_SOLVENTS);

		bonded_particles_lut_manager = genericMoveToDevice(bonded_particles_lut_manager, 1);

		cudaDeviceSynchronize();
		printf("Box transferred to device\n");
	}
};

	

// This stays on host
class Simulation {
public:
	Simulation() {
		box = new Box();
	}

	__host__ void moveToDevice() {
		box = genericMoveToDevice(box, 1);
		cudaDeviceSynchronize();
		if (cudaGetLastError() != cudaSuccess) {
			fprintf(stderr, "Error during Simulation Host->Device transfer\n");
			exit(1);
		}
		printf("Simulation ready for device\n\n");
	}
	__host__ void copyBoxVariables() {
		n_compounds = box->n_compounds;
		n_bridges = box->bridge_bundle->n_bridges;


		n_solvents = box->n_solvents;
		blocks_per_solventkernel = (int) ceil((float)n_solvents / (float)THREADS_PER_SOLVENTBLOCK);


		compounds_host = new Compound[n_compounds];
		for (int i = 0; i < n_compounds; i++)
			compounds_host[i] = box->compounds[i];

	}
	__host__ void incStep() {
		step++;
		box->step++;
	}
	__host__ inline uint64_t getStep() {
		return step;
	}
	~Simulation() {}

	bool finished = false;
	//int step = 0;


	float* potE_buffer = NULL;	// Not really a buffer yet, just one large array that holds full simulation data
	Float3* traj_buffer = nullptr;
	float* temperature_buffer = nullptr;
	int n_temp_values = 0;
	Float3* traindata_buffer = nullptr;		// Position and force data for all particles, for NN training
	float* logging_data = nullptr;				// Used for debugging/logging any values. 10 floats per step!

	uint32_t total_particles_upperbound = 0;
	uint32_t total_compound_particles = 0;			// Precise number, but DO NOT EVER USE IN INDEXING!!
	uint32_t total_particles = 0;				// Precise number, but DO NOT EVER USE IN INDEXING!!
	
	int n_steps = SIMULATION_STEPS;

	const double dt = 1 * 1e-6;					// [ns] first val corresponds to fs
	const double dt_pico = dt * 1000.f;
	int steps_per_render = STEPS_PER_RENDER;
	//int n_bodies = N_BODIES_START;
	Box* box;


	Compound* compounds_host = nullptr;				// For reading static data, for example during nlist-search

	// Box variable copies, here for ease of access.
	int n_compounds = 0;
	int n_bridges = 0; 
	int n_solvents = 0;



	std::string out_dir = OUT_DIR;
	


	//int blocks_per_solventkernel = ceil((float)n_solvents / (float)THREADS_PER_SOLVENTBLOCK);
	int blocks_per_solventkernel = 0;
private:
	uint64_t step = 0;
};

