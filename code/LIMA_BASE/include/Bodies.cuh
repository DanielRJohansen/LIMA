#pragma once
#include <iostream>

#include "Constants.h"
#include "LimaTypes.cuh"
#include <vector>
#include <array>

#include <cuda_fp16.h>

// Hyper-fast objects for kernel, so we turn all safety off here!
#pragma warning (push)
#pragma warning (disable : 26495)







//--------------------------- THE FOLLOWING IS FOR HANDLING INTRAMOLECULAR FORCES ---------------------------//

struct NBAtomtype {
	NBAtomtype(){}
	NBAtomtype(float m, float s, float e) : mass(m), sigma(s), epsilon(e) {}
	float mass = 0.f;			// kg / mol
	float sigma = 0.f;			// nm
	float epsilon = 0.f;		// J/mol
};


struct SingleBond {	// IDS and indexes are used interchangeably here!
	SingleBond(){}
	SingleBond(std::array<uint8_t, 2> ids, float b0, float kb);

	//SingleBond(uint32_t particleindex_a, uint32_t particleindex_b) {
	//	atom_indexes[0] = particleindex_a;
	//	atom_indexes[1] = particleindex_b;
	//}

	float b0 = 0.f;
	float kb = 0.f;
	uint8_t atom_indexes[2] = {0,0};	// Relative to the compund - NOT ABSOLUTE INDEX. Used in global table with compunds start-index
	const static int n_atoms = 2;
};
struct SingleBondFactory : public SingleBond {
	SingleBondFactory(std::array<uint32_t, 2> ids, float b0, float kb);

	uint32_t global_atom_indexes[2] = { 0,0 };
};

struct AngleBond {
	AngleBond() {}
	AngleBond(std::array<uint8_t, 3> ids, float theta_0, float k_theta);

	float theta_0 = 0.f;
	float k_theta = 0.f;
	uint8_t atom_indexes[3] = {0,0,0}; // i,j,k angle between i and k
	const static int n_atoms = 3;
};
struct AngleBondFactory : public AngleBond {
	AngleBondFactory(std::array<uint32_t, 3> ids, float theta_0, float k_theta);

	uint32_t global_atom_indexes[3] = { 0,0, 0 };
};

struct DihedralBond {
	const static int n_atoms = 4;
	DihedralBond() {}
	DihedralBond(std::array<uint8_t, 4> ids, float phi0, float kphi, float n);
	
	half phi_0 = 0.f;
	half k_phi = 0.f;
	half n = 0.f;		// n parameter, how many energy equilibriums does the dihedral have // OPTIMIZE: maybe float makes more sense, to avoid conversion in kernels?
	uint8_t atom_indexes[4] = {0,0,0,0};
};
struct DihedralBondFactory : public DihedralBond {
	DihedralBondFactory(std::array<uint32_t, 4> ids, float phi0, float kphi, float n);

	uint32_t global_atom_indexes[4] = { 0,0, 0, 0 };
};

struct ImproperDihedralBond {
	ImproperDihedralBond() {}
	ImproperDihedralBond(std::array<uint8_t, 4> ids, float psi_0, float k_psi);

	//TODO: make const?
	float psi_0 = 0.f;
	float k_psi = 0.f;

	uint8_t atom_indexes[4] = { 0,0,0,0 };
	const static int n_atoms = 4;
};
struct ImproperDihedralBondFactory : public ImproperDihedralBond {
	ImproperDihedralBondFactory(std::array<uint32_t, 4> ids, float psi_0, float k_psi);

	uint32_t global_atom_indexes[4] = { 0,0, 0, 0 };
};




// ------------------------------------------------- COMPOUNDS ------------------------------------------------- //







struct CompoundCoords {
	__device__ void loadData(const CompoundCoords& coords) {
		if (threadIdx.x == 0) { origo = coords.origo; };
		rel_positions[threadIdx.x] = coords.rel_positions[threadIdx.x];
	}

	NodeIndex origo{};								// [nm]
	Coord rel_positions[MAX_COMPOUND_PARTICLES];	// [lm]


	__host__ void static copyInitialCoordConfiguration(CompoundCoords* coords,
		CompoundCoords* coords_prev, CompoundCoords* coordarray_circular_queue);
	//__host__ Float3 getAbsolutePositionLM(int particle_id);
};

const int MAX_SINGLEBONDS_IN_COMPOUND = MAX_COMPOUND_PARTICLES+4;	// Due to AA such a TRP, thhere might be more bonds than atoms
const int MAX_ANGLEBONDS_IN_COMPOUND = 128;
const int MAX_DIHEDRALBONDS_IN_COMPOUND = 128+64;
const int MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND = 32;
struct CompoundState {							// Maybe delete this soon?
	__device__ void setMeta(int n_p) {
		n_particles = n_p;
	}
	__device__ void loadData(const CompoundCoords& coords) {
		if (threadIdx.x < n_particles)
			positions[threadIdx.x] = coords.rel_positions[threadIdx.x].toFloat3();
	}



	Float3 positions[MAX_COMPOUND_PARTICLES];
	uint8_t n_particles = 0;
};

struct SolventCoord {
	__device__ __host__ SolventCoord() {}
	__device__ __host__ SolventCoord(NodeIndex ori , Coord rel) : origo(ori ), rel_position(rel) {}
	NodeIndex origo{};									// [nm]
	Coord rel_position{};							// [lm]

	__host__ void static copyInitialCoordConfiguration(SolventCoord* coords,
		SolventCoord* coords_prev, SolventCoord* coordarray_circular_queue);
};


// struct with data that only the solvent itself needs
struct Solvent {
	// Needed for VelocityVS
	Float3 vel_prev{};
	Float3 force_prev{};
};








// blocks are notcentered 
struct SolventBlock {
	static const int MAX_SOLVENTS_IN_BLOCK = 64;

	__device__ __host__ SolventBlock() {};
	__device__ __host__ void loadMeta(const SolventBlock& block) {
		origo = block.origo; // Not necessary, is given by the blockIdx.x
		n_solvents = block.n_solvents;
		if (n_solvents >= MAX_SOLVENTS_IN_BLOCK) {
			printf("Too many solvents in block!\n");
		}
	}
	__device__ __host__ void loadData(const SolventBlock& block) {
		rel_pos[threadIdx.x] = Coord{};	// temp
		if (threadIdx.x < n_solvents) {

			if (block.rel_pos[threadIdx.x] == Coord{ 0 }) {
				printf("Loading zeroes blockid %d nsol %d\n", blockIdx.x, n_solvents);
			}

			rel_pos[threadIdx.x] = block.rel_pos[threadIdx.x];
			ids[threadIdx.x] = block.ids[threadIdx.x];
		}
	}

	__host__ bool addSolvent(const Coord& rel_position, uint32_t id) {
		if (n_solvents == MAX_SOLVENTS_IN_BLOCK) {
			return false;
		}
		ids[n_solvents] = id;
		rel_pos[n_solvents++] = rel_position;
		return true;
	}

	NodeIndex origo{};
	int n_solvents = 0;
	Coord rel_pos[MAX_SOLVENTS_IN_BLOCK];	// Pos rel to lower left forward side of block, or floor() of pos
	uint32_t ids[MAX_SOLVENTS_IN_BLOCK];
};

static constexpr int BOXGRID_NODE_LEN_pico = 1000;
static constexpr float BOXGRID_NODE_LEN = static_cast<float>(BOXGRID_NODE_LEN_pico * PICO_TO_LIMA);	// [lm]
static_assert((static_cast<int64_t>(BOX_LEN) % static_cast<int64_t>(BOXGRID_NODE_LEN)) == 0, "Illegal box dimension");
static constexpr int32_t BOXGRID_NODE_LEN_i = BOXGRID_NODE_LEN_pico * PICO_TO_LIMA;

static const int BOXGRID_N_NODES = static_cast<int>(BOX_LEN / BOXGRID_NODE_LEN);	// Only reference to compile time BOXLEN in this file!

class SolventBlocksCircularQueue {
	static const int queue_len = STEPS_PER_SOLVENTBLOCKTRANSFER;


	// Please dont add other non-static vars to this class without checking movetodevice is not broken
	// {Step,z,y,x}
	SolventBlock* blocks = nullptr;
	bool has_allocated_data = false;
	bool is_on_device = false;


public:
	static const int blocks_per_grid = BOXGRID_N_NODES * BOXGRID_N_NODES * BOXGRID_N_NODES;
	static const int grid_bytesize = sizeof(SolventBlock) * blocks_per_grid;
	static const int blocks_total = blocks_per_grid * queue_len;
	static const int gridqueue_bytesize = sizeof(SolventBlock) * blocks_total;
	static const int first_step_prev = queue_len - 1;

	SolventBlocksCircularQueue() {};	// C

	__host__ static SolventBlocksCircularQueue* createQueue() {
		auto queue = new SolventBlocksCircularQueue();
		queue->allocateData();
		queue->initializeBlocks();
		return queue;
	}

	__host__ bool addSolventToGrid(const SolventCoord& coord, uint32_t solvent_id, int step) {
		// TODO: Implement safety feature checking and failing if PBC is not met!
		return getBlockPtr(coord.origo, step)->addSolvent(coord.rel_position, solvent_id);
	}

	__host__ void allocateData() {
		if (has_allocated_data) {
			throw std::runtime_error("BoxGridCircularQueue may not allocate data multiple times");
		}
		blocks = new SolventBlock[blocks_total]();
		has_allocated_data = true;
	}

	__host__ void initializeBlocks() {
		for (int step = 0; step < queue_len; step++) {
			for (int i = 0; i < blocks_per_grid; i++) {
				getBlockPtr(i, step)->origo = get3dIndex(i);
			}
		}
	}


	__host__ SolventBlocksCircularQueue* moveToDevice() {
		SolventBlock* blocks_dev_ptr = nullptr;
		cudaMallocManaged(&blocks_dev_ptr, gridqueue_bytesize);
		cudaMemcpy(blocks_dev_ptr, blocks, gridqueue_bytesize, cudaMemcpyHostToDevice);
		delete[] blocks;
		blocks = blocks_dev_ptr;
		is_on_device = true;


		SolventBlocksCircularQueue* this_dev_ptr;
		cudaMallocManaged(&this_dev_ptr, sizeof(SolventBlocksCircularQueue));
		cudaMemcpy(this_dev_ptr, this, sizeof(SolventBlocksCircularQueue), cudaMemcpyHostToDevice);	
		// LEAK: I think we leak *this here
		return this_dev_ptr;
	}

	// Fuck me this seems overcomplicated
	__host__ SolventBlocksCircularQueue* copyToHost() {
		SolventBlocksCircularQueue* this_host = new SolventBlocksCircularQueue();
		this_host->allocateData();
		cudaMemcpy(this_host->blocks, blocks, gridqueue_bytesize, cudaMemcpyDeviceToHost);
		return this_host;
	}




	// This function assumes the user has used PBC
	__host__ SolventBlock* getBlockPtr(const NodeIndex& index3d, const int step) {
#if defined LIMASAFEMODE
		if (index3d.x >= BOXGRID_N_NODES || index3d.y >= BOXGRID_N_NODES || index3d.z >= BOXGRID_N_NODES
			|| index3d.x < 0 || index3d.y < 0 || index3d.z < 0) {
			throw std::runtime_error("Bad 3d index for blockptr\n");
		}
#endif

		return getBlockPtr(get1dIndex(index3d), step);
	}

	__device__ __host__ bool static isTransferStep(int step) {
		return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == SOLVENTBLOCK_TRANSFERSTEP;
	}
	__device__ bool static isFirstStepAfterTransfer(int step) {
		return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == 0;
	}

	// This function assumes the user has used PBC
	__device__ __host__ SolventBlock* getBlockPtr(const size_t index1d, const size_t step) {
		const size_t step_offset = (step % queue_len) * blocks_per_grid;
		return &blocks[index1d + step_offset];
	}

	__device__ __host__ static int get1dIndex(const NodeIndex& index3d) {
		static const int bpd = BOXGRID_N_NODES;
		return index3d.x + index3d.y * bpd + index3d.z * bpd * bpd;
	}
	__device__ __host__ static NodeIndex get3dIndex(int index1d) {
		static const int bpd = BOXGRID_N_NODES;
		auto z = index1d / (bpd * bpd);
		index1d -= z * bpd * bpd;
		auto y = index1d / bpd;
		index1d -= y * bpd;
		auto x = index1d;
		return NodeIndex{ x, y, z };
	}
};









namespace CoordArrayQueueHelpers {
	__host__ __device__ static CompoundCoords* getCoordarrayRef(CompoundCoords* coordarray_circular_queue,
		uint32_t step, uint32_t compound_index) {
		const int index0_of_currentstep_coordarray = (step % STEPS_PER_LOGTRANSFER) * MAX_COMPOUNDS;
		return &coordarray_circular_queue[index0_of_currentstep_coordarray + compound_index];
	}
}






// Instead of having a single key_particle and an single radius, we now have multiple
struct CompoundInteractionBoundary {
	static const int k = 2;

	float radii[k];	// [nm]
	int key_particle_indices[k];
};

struct CompoundCompact {
	__host__ __device__ CompoundCompact() {}

	uint8_t n_particles = 0;			
	uint8_t atom_types[MAX_COMPOUND_PARTICLES];

#ifdef LIMAKERNELDEBUGMODE
	uint32_t particle_global_ids[MAX_COMPOUND_PARTICLES];
#endif


	int n_singlebonds = 0;
	int n_anglebonds = 0;
	int n_dihedrals = 0;
	int n_improperdihedrals = 0;

	// Use this to quickly lookup wheter a bondedparticleslut exists with another compound
	static const int max_bonded_compounds = MAX_COMPOUNDS_IN_BRIDGE * 2 - 2;
	int n_bonded_compounds = 0;

	__device__ void loadMeta(CompoundCompact* compound) {
		n_particles = compound->n_particles;
		n_singlebonds = compound->n_singlebonds;
		n_anglebonds = compound->n_anglebonds;
		n_dihedrals = compound->n_dihedrals;
		n_improperdihedrals = compound->n_improperdihedrals;
		n_bonded_compounds = compound->n_bonded_compounds;
	}

	__device__ void loadData(CompoundCompact* compound) {
		if (threadIdx.x < n_particles) {
			atom_types[threadIdx.x] = compound->atom_types[threadIdx.x];

			#ifdef LIMAKERNELDEBUGMODE
			particle_global_ids[threadIdx.x] = compound->particle_global_ids[threadIdx.x];
			#endif
		}
	}
};

// Rather large unique structures in global memory, that can be partly loaded when needed
struct Compound : public CompoundCompact {
	__host__ __device__ Compound() {}

	// For drawing pretty spheres :)	//TODO Move somewhere else
	uint8_t atom_color_types[MAX_COMPOUND_PARTICLES];

	// Interims from the bridgekernel to compoundkernel
	//bool is_in_bridge[MAX_COMPOUND_PARTICLES];	// TODO: implement this?
	float potE_interim[MAX_COMPOUND_PARTICLES];
	Float3 forces_interim[MAX_COMPOUND_PARTICLES];

	// Used specifically for Velocity Verlet stormer, and ofcourse kinE fetching
	Float3 forces_prev[MAX_COMPOUND_PARTICLES];
	Float3 vels_prev[MAX_COMPOUND_PARTICLES]; // Get wierd change of outcome if i move this here??

	// Bonds
	SingleBond singlebonds[MAX_SINGLEBONDS_IN_COMPOUND];
	AngleBond anglebonds[MAX_ANGLEBONDS_IN_COMPOUND];
	DihedralBond dihedrals[MAX_DIHEDRALBONDS_IN_COMPOUND];
	ImproperDihedralBond impropers[MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND];

	CompoundInteractionBoundary interaction_boundary;
	int centerparticle_index = -1;			// Index of particle initially closest to CoM

	uint16_t bonded_compound_ids[max_bonded_compounds];	// *2-2because it should exclude itself from both sides
};








struct ParticleReference {
	ParticleReference() {}	// TODO: i dont think i need this.

#ifdef LIMAKERNELDEBUGMODE
	// Used by moleculebuilder only
	ParticleReference(int compound_id, int local_id_compound, int gro_id, uint8_t compoundid_local_to_bridge) :
		compound_id(compound_id), local_id_compound(local_id_compound), 
		gro_id(gro_id), compoundid_local_to_bridge(compoundid_local_to_bridge)
	{}
	int gro_id = -1;
#else
	// Used by moleculebuilder only
	ParticleReference(int compound_id, int local_id_compound, uint8_t compoundid_local_to_bridge) :
		compound_id(compound_id), local_id_compound(local_id_compound),
		compoundid_local_to_bridge(compoundid_local_to_bridge) 
	{}
#endif

	int compound_id;	// global
	int local_id_compound;	// id of particle
	uint8_t compoundid_local_to_bridge = 255;

	//int global_id = -1; // For debug
};

struct CompoundBridge {
	CompoundBridge() {}	
	
	ParticleReference particle_refs[MAX_PARTICLES_IN_BRIDGE]{};
	uint8_t n_particles = 0;					

	uint8_t n_singlebonds = 0;
	SingleBond singlebonds[MAX_SINGLEBONDS_IN_BRIDGE];

	uint8_t n_anglebonds = 0;
	AngleBond anglebonds[MAX_ANGLEBONDS_IN_BRIDGE];

	uint8_t n_dihedrals = 0;
	DihedralBond dihedrals[MAX_DIHEDRALBONDS_IN_BRIDGE];

	uint8_t n_improperdihedrals = 0;
	ImproperDihedralBond impropers[MAX_IMPROPERDIHEDRALBONDS_IN_BRIDGE];

	static_assert(MAX_COMPOUNDS < UINT16_MAX, "CompoundBridge cannot handle such large compound ids");
	uint8_t n_compounds;

	uint16_t compound_ids[MAX_COMPOUNDS_IN_BRIDGE];


	// -------------- Device functions ------------- //
	__device__ void loadMeta(CompoundBridge* bridge) {
		n_particles = bridge->n_particles;
		n_singlebonds = bridge->n_singlebonds;
		n_anglebonds = bridge->n_anglebonds;
		n_dihedrals = bridge->n_dihedrals;
		n_improperdihedrals = bridge->n_improperdihedrals;
		n_compounds = bridge->n_compounds;

		for (int i = 0; i < 4; i++)
			compound_ids[i] = bridge->compound_ids[i];
		
		//compound_id_left = bridge->compound_id_left;
		//compound_id_right = bridge->compound_id_right;
	}
	__device__ void loadData(CompoundBridge* bridge) {
		if (threadIdx.x < n_particles) {
			//atom_types[threadIdx.x] = bridge->atom_types[threadIdx.x];
			particle_refs[threadIdx.x] = bridge->particle_refs[threadIdx.x];
		}
		
		for (int i = 0; (i * blockDim.x) < n_singlebonds; i++) {
			int index = i * blockDim.x + threadIdx.x;
			if (index < n_singlebonds)
				singlebonds[index] = bridge->singlebonds[index];
		}
		for (int i = 0; (i * blockDim.x) < n_anglebonds; i++) {
			int index = i * blockDim.x + threadIdx.x;
			if (index < n_anglebonds)
				anglebonds[index] = bridge->anglebonds[index];
		}
		for (int i = 0; (i * blockDim.x) < n_dihedrals; i++) {
			int index = i * blockDim.x + threadIdx.x;
			if (index < n_dihedrals)
				dihedrals[index] = bridge->dihedrals[index];
		}
		if (threadIdx.x < n_improperdihedrals) {
			impropers[threadIdx.x] = bridge->impropers[threadIdx.x];
		}
	}
};

struct CompoundBridgeBundleCompact {
	CompoundBridgeBundleCompact() {}
	CompoundBridgeBundleCompact(const std::vector<CompoundBridge>& bridges) {
		for (int i = 0; i < bridges.size(); i++) {
			compound_bridges[i] = bridges[i];
		}
		n_bridges = static_cast<int>(bridges.size());
	}


	CompoundBridge compound_bridges[MAX_COMPOUNDBRIDGES];
	int n_bridges = 0;
};






struct ForceField_NB {
	struct ParticleParameters {	//Nonbonded
		float mass = -1;		//[kg/mol]	or 
		float sigma = -1;		// []
		float epsilon = -1;		// J/mol [kg*nm^2 / s^2]
	};

	ParticleParameters particle_parameters[MAX_ATOM_TYPES];
};









#pragma warning (pop)
