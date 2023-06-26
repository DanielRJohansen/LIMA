#pragma once
#include <iostream>

#include "Constants.cuh"
#include "LimaTypes.cuh"
#include <vector>
#include <array>

// Hyper-fast objects for kernel, so we turn all safety off here!
#pragma warning (push)
#pragma warning (disable : 26495)







//--------------------------- THE FOLLOWING IS FOR HANDLING INTRAMOLECULAR FORCES ---------------------------//

struct NBAtomtype {
	NBAtomtype(){}
	NBAtomtype(float m, float s, float e) : mass(m), sigma(s), epsilon(e) {}
	float mass = 0.f;			// kg / mol
	float sigma = 0.f;		// nm
	float epsilon = 0.f;		// J/mol
};


struct SingleBond {	// IDS and indexes are used interchangeably here!
	SingleBond(){}
	SingleBond(std::array<int, 2> ids); // Used when loading topology only
	SingleBond(int id1, int id2, float b0, float kb);
	SingleBond(std::array<uint32_t, 2> ids, float b0, float kb);

	SingleBond(uint32_t particleindex_a, uint32_t particleindex_b) {
		atom_indexes[0] = particleindex_a;
		atom_indexes[1] = particleindex_b;
	}
	SingleBond(float ref_dist, float kb, uint32_t particleindex_a, uint32_t particleindex_b) :	// TODO this is duplicate, remove!
		b0(ref_dist), kb(kb) {
		atom_indexes[0] = particleindex_a;
		atom_indexes[1] = particleindex_b;
	}

	float b0 = 0.f;
	float kb = 0.f;
	uint32_t atom_indexes[2] = {0,0};	// Relative to the compund - NOT ABSOLUTE INDEX. Used in global table with compunds start-index
	const static int n_atoms = 2;
	//bool invertLJ = false;		// When sharing multiple bonded connections
};

struct AngleBond {
	AngleBond() {}
	AngleBond(std::array<int, 3> ids); // Used when loading topology only
	AngleBond(int id1, int id2, int id3, float theta_0, float k_theta);	//todo: remove
	AngleBond(std::array<uint32_t, 3> ids, float theta_0, float k_theta);

	float theta_0 = 0.f;
	float k_theta = 0.f;
	uint32_t atom_indexes[3] = {0,0,0}; // i,j,k angle between i and k
	const static int n_atoms = 3;
};

struct DihedralBond {
	const static int n_atoms = 4;
	DihedralBond() {}
	DihedralBond(std::array<int, 4> ids); // Used when loading topology only
	DihedralBond(int id1, int id2, int id3, int id4, float phi_0, float k_phi, int n);
	DihedralBond(std::array<uint32_t, 4> ids, float phi0, float kphi, int n);

	float phi_0 = 0.f;
	float k_phi = 0.f;
	uint8_t n = 0;		// n parameter, how many energy equilibriums does the dihedral have // OPTIMIZE: maybe float makes more sense, to avoid conversion in kernels?
	uint32_t atom_indexes[4] = {0,0,0,0};
};

struct ImproperDihedralBond {
	ImproperDihedralBond() {}
	ImproperDihedralBond(std::array<uint32_t, 4> ids); // Used when loading topology only
	ImproperDihedralBond(std::array<uint32_t, 4> ids, float psi_0, float k_psi);

	//TODO: make const?
	float psi_0 = 0.f;
	float k_psi = 0.f;

	uint32_t atom_indexes[4] = { 0,0,0,0 };
	const static int n_atoms = 4;
};





// ------------------------------------------------- COMPOUNDS ------------------------------------------------- //
class NeighborList {
public:
	enum NEIGHBOR_TYPE {COMPOUND, SOLVENT};


	__host__ bool addId(uint16_t new_id, NEIGHBOR_TYPE nt);
	__host__ bool removeId(uint16_t neighbor_id, NEIGHBOR_TYPE nt);

	__host__ void addGridnode(uint16_t gridnode_id);
	__host__ void removeGridnode(uint16_t gridnode_id);



	__device__ void loadMeta(NeighborList* nl_ptr) {	// Called from thread 0
		n_compound_neighbors = nl_ptr->n_compound_neighbors;
		n_solvent_neighbors = nl_ptr->n_solvent_neighbors;
		n_gridnodes = nl_ptr->n_gridnodes;
	}
	__device__ void loadData(NeighborList* nl_ptr) {
		if (threadIdx.x < n_compound_neighbors)			// DANGER Breaks when threads < mAX_COMPOUND_Ns
			neighborcompound_ids[threadIdx.x] = nl_ptr->neighborcompound_ids[threadIdx.x];
		for (int i = threadIdx.x;  i < n_solvent_neighbors; i += blockDim.x) {// Same as THREADS_PER_COMPOUNDBLOCK
			neighborsolvent_ids[i] = nl_ptr->neighborsolvent_ids[i];
			i += blockDim.x;	
		}

		if (threadIdx.x < n_gridnodes) {
			gridnode_ids[threadIdx.x] = nl_ptr->gridnode_ids[threadIdx.x];
		}
	}


	uint16_t neighborcompound_ids[NEIGHBORLIST_MAX_COMPOUNDS];
	uint16_t n_compound_neighbors = 0;
	uint16_t neighborsolvent_ids[NEIGHBORLIST_MAX_SOLVENTS];
	uint16_t n_solvent_neighbors = 0;

	static const int max_gridnodes = 48;
	uint16_t gridnode_ids[max_gridnodes];
	uint8_t n_gridnodes = 0;

	int associated_id = -1;

};






struct CompoundCoords {
	__device__ void loadData(CompoundCoords& coords) {
		if (threadIdx.x == 0) { origo = coords.origo; };
		rel_positions[threadIdx.x] = coords.rel_positions[threadIdx.x];
	}

	NodeIndex origo{};								// [nm]
	Coord rel_positions[MAX_COMPOUND_PARTICLES]{};	// [lm]
	static constexpr int firststep_prev = STEPS_PER_LOGTRANSFER - 1;


	__host__ void static copyInitialCoordConfiguration(CompoundCoords* coords,
		CompoundCoords* coords_prev, CompoundCoords* coordarray_circular_queue);
	//__host__ Float3 getAbsolutePositionLM(int particle_id);
};

const int MAX_SINGLEBONDS_IN_COMPOUND = 64;
const int MAX_ANGLEBONDS_IN_COMPOUND = 128;
const int MAX_DIHEDRALBONDS_IN_COMPOUND = 256;
const int MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND = 16;
struct CompoundState {							// Maybe delete this soon?
	__device__ void setMeta(int n_p) {
		n_particles = n_p;
	}
	/*__device__ void loadData(CompoundState* state) {
		if (threadIdx.x < n_particles)
			positions[threadIdx.x] = state->positions[threadIdx.x];
	}*/
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

	//__host__ SolventCoord static createFromPositionNM(const Float3& solvent_pos);
	//__host__ __device__ Float3 getAbsolutePositionLM() const {
	//	return ((origo * NANO_TO_LIMA).toFloat3() + rel_position.toFloat3());
	//}
};






// blocks are notcentered 
struct SolventBlock {
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
	static constexpr float block_len = 1; // [nm]	// TODO: Remove this, there should only be node_len of BoxGrid
	static_assert((static_cast<int>(BOX_LEN_NM) % static_cast<int>(block_len)) == 0, "Illegal box dimension");

	NodeIndex origo{};
	uint16_t n_solvents = 0;
	Coord rel_pos[MAX_SOLVENTS_IN_BLOCK];	// Pos rel to lower left forward side of block, or floor() of pos
	uint32_t ids[MAX_SOLVENTS_IN_BLOCK];
};


static constexpr int BOXGRID_NODE_LEN_pico = 1200;
static constexpr float BOXGRID_NODE_LEN = static_cast<float>(BOXGRID_NODE_LEN_pico * PICO_TO_LIMA);	// [lm]
static constexpr int32_t BOXGRID_NODE_LEN_i = BOXGRID_NODE_LEN_pico * PICO_TO_LIMA;
static const int BOXGRID_N_NODES = static_cast<int>(BOX_LEN / BOXGRID_NODE_LEN);
template <typename NodeType>
class BoxGrid {
public:
	//static const int blocks_per_dim = static_cast<int>(BOX_LEN_NM) / SolventBlock::block_len;
	static const int blocks_total = BOXGRID_N_NODES * BOXGRID_N_NODES * BOXGRID_N_NODES;
	NodeType blocks[blocks_total];
	static const int first_step_prev = STEPS_PER_SOLVENTBLOCKTRANSFER - 1;


	// This function assumes the user has used PBC
	__host__ NodeType* getBlockPtr(const NodeIndex& index3d) {
		if (index3d.x >= BOXGRID_N_NODES || index3d.y >= BOXGRID_N_NODES || index3d.z >= BOXGRID_N_NODES
		|| index3d.x < 0 || index3d.y < 0 || index3d.z < 0) {
			printf("BAD 3D index\n"); exit(1); }
		return getBlockPtr(BoxGrid::get1dIndex(index3d));
	}

	// This function assumes the user has used PBC
	__device__ __host__ NodeType* getBlockPtr(const int index1d) {
		return &blocks[index1d];
	}

	__device__ __host__ static int get1dIndex(const NodeIndex& index3d) {
		static const int bpd = BOXGRID_N_NODES;
		return index3d.x + index3d.y * bpd + index3d.z * bpd * bpd;
	}
	__device__ static NodeIndex get3dIndex(int index1d) {
		static const int bpd = BOXGRID_N_NODES;
		auto z = index1d / (bpd * bpd);
		index1d -= z * bpd * bpd;
		auto y = index1d / bpd;
		index1d -= y * bpd;
		auto x = index1d;
		return NodeIndex{ x, y, z };
	}	
};

class SolventBlockGrid : public BoxGrid<SolventBlock> {
public:
	__host__ bool addSolventToGrid(const SolventCoord& coord, uint32_t solvent_id) {
		// TODO: Implement safety feature checking and failing if PBC is not met!
		return getBlockPtr(coord.origo)->addSolvent(coord.rel_position, solvent_id);
	}
};

template <int size>
struct SolventTransferqueue {
	Coord rel_positions[size];
	Coord rel_positions_prev[size];	// Adjust rel positions to the new block's origo BEFORE putting it here!
	uint32_t ids[size];
	int n_elements = 0;

	// Do NOT call on queue residing in global memory
	__device__ bool addElement(const Coord& pos, const Coord& pos_prev, uint32_t id) {
		if (n_elements >= size) { 
			printf("\nTried to add too many solvents in outgoing transferqueue\n"); 
			return false;
		}
		rel_positions[n_elements] = pos;
		rel_positions_prev[n_elements] = pos_prev;
		ids[n_elements] = id;
		n_elements++;
		return true;
	}

	// Insert relative to thread calling.
	__device__ void fastInsert(const Coord& relpos, const Coord& relpos_prev, const int id) {
		//rel_positions[threadIdx.x]		= relpos - transfer_dir * static_cast<int32_t>(NANO_TO_LIMA);
		//rel_positions_prev[threadIdx.x] = relpos_prev - transfer_dir * static_cast<int32_t>(NANO_TO_LIMA);
		rel_positions[threadIdx.x] = relpos;
		rel_positions_prev[threadIdx.x] = relpos_prev;
		ids[threadIdx.x] = id;
	}

	//__device__ void setUsedSize(int usedsize) {
	//	n_elements = usedsize;
	//}
};

struct SolventBlockTransfermodule {
	// Only use directly (full plane contact) adjecent blocks
	static const int n_queues = 6;			// or, adjecent_solvent_blocks
	static const int max_queue_size = 64;	// Maybe this is a bit dangerous

	// I NEED TO FIGURE OUT PREV_POS FOR NON-REMAINING SOLVENTS!!!!
	// Each queue will be owned solely by 1 adjecent solventblock
	SolventTransferqueue<max_queue_size> transfer_queues[n_queues];

	Coord remain_relpos_prev[MAX_SOLVENTS_IN_BLOCK];
	int n_remain = 0;

	/// <param name="transfer_direction">Relative to the originating block</param>
	__device__ static int getQueueIndex(const NodeIndex& transfer_direction) {
		// Fucking magic yo...
		// Only works if transfer_direction.len() == 1
		// First op leaves a space at index 3:
		//{-z, -y, -x, _ x, y, z}
		//const int tmp_index = transfer_direction.dot(NodeIndex{ 1, 2, 3 });
		int tmp_index = transfer_direction.x * 1 + transfer_direction.y * 2 + transfer_direction.z * 3;
		// Shift positive values left.
		tmp_index = tmp_index > 0 ? tmp_index + 2 : tmp_index + 3;
		return tmp_index;
	}

	// Maybe too convoluted...
	// I THINK THIS IS WRONG:, it should be using temp_index no?
	//static NodeIndex getDirectionFromQueueIndex(const int index) {
	//	const int tmp_index = index > 2 ? index - 2 : index - 3;
	//	NodeIndex direction{};
	//	direction.x += index * (abs(index) == 1);
	//	direction.y += index/2 * (abs(index) == 2);
	//	direction.x += index/3 * (abs(index) == 3);
	//	return direction;
	//}


	//Coord* getQueuePtr(const Coord& transfer_direction) {
	//	// Fucking magic yo...
	//	auto index = transfer_direction.dot(Coord{ 1, 2, 3 }) + 2;
	//}
};

using STransferQueue = SolventTransferqueue< SolventBlockTransfermodule::max_queue_size>;
using SRemainQueue = SolventTransferqueue<MAX_SOLVENTS_IN_BLOCK>;

namespace SolventBlockHelpers {
	//bool insertSolventcoordInGrid(SolventBlockGrid& grid, const SolventCoord& coord, uint32_t solvent_id);
	bool copyInitialConfiguration(const SolventBlockGrid& grid, const SolventBlockGrid& grid_prev,
		SolventBlockGrid* grid_circular_queue);

	void createSolventblockGrid(SolventBlockGrid** solventblockgrid_circularqueue);
	void setupBlockMetaOnHost(SolventBlockGrid* grid, SolventBlockGrid* grid_prev);	// Is this used?
	//__host__ void createSolventblockTransfermodules(SolventBlockTransfermodule** transfermodule_array);

	__device__ __host__ bool static isTransferStep(int step) {
		return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == SOLVENTBLOCK_TRANSFERSTEP;
	}
	__device__ bool static isFirstStepAfterTransfer(int step) {
		return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == 0;
	}

	// Only use this if if the index3d is inside the box!


	//__device__ static Float3 extractAbsolutePositionLM(const SolventBlock& block) {
	//	return (block.origo * static_cast<int32_t>(NANO_TO_LIMA) + block.rel_pos[threadIdx.x]).toFloat3();
	//}

}












namespace CoordArrayQueueHelpers {
	/*__host__ static void copyInitialCoordConfiguration(CompoundCoords* coords,
		CompoundCoords* coords_prev, CompoundCoords* coordarray_circular_queue) */

	__host__ __device__ static CompoundCoords* getCoordarrayRef(CompoundCoords* coordarray_circular_queue,
		uint32_t step, uint32_t compound_index) {
		const int index0_of_currentstep_coordarray = (step % STEPS_PER_LOGTRANSFER) * MAX_COMPOUNDS;
		return &coordarray_circular_queue[index0_of_currentstep_coordarray + compound_index];
	}


	__host__ __device__ static SolventBlockGrid* getSolventblockGridPtr(SolventBlockGrid* solventblockgrid_circular_queue,
		const int step) {
			const int index0_of_currentstep_blockarray = (step % STEPS_PER_SOLVENTBLOCKTRANSFER);
			return &solventblockgrid_circular_queue[index0_of_currentstep_blockarray];
	}

	__device__ static SolventBlock* getSolventBlockPtr(SolventBlockGrid* solventblockgrid_circular_queue,
		const int step, const int solventblock_id) {
		//const int index0_of_currentstep_blockarray = (step % STEPS_PER_SOLVENTBLOCKTRANSFER);
		//return solventblockgrid_circular_queue[index0_of_currentstep_blockarray].getBlockPtr(solventblock_id);
		return getSolventblockGridPtr(solventblockgrid_circular_queue, step)->getBlockPtr(solventblock_id);
	}

	// WTF is this??
	//__device__ static SolventCoord* getSolventcoordPtr(SolventCoord* solventcoordarray_circular_queue,
	//	const int step, const int solvent_id) {
	//	const int index0_of_currentstep_coordarray = (step % STEPS_PER_LOGTRANSFER) * MAX_SOLVENTS;
	//	return &solventcoordarray_circular_queue[index0_of_currentstep_coordarray + solvent_id];
	//}
}




struct Compound {
	__host__ __device__  Compound() {}

	uint8_t n_particles = 0;			
	Float3 forces[MAX_COMPOUND_PARTICLES];					// Carries forces from bridge_kernels // TODO: find another mechanism
	uint8_t atom_types[MAX_COMPOUND_PARTICLES];
	uint8_t atom_color_types[MAX_COMPOUND_PARTICLES];	// For drawing pretty spheres :)	//TODO Move somewhere else

	//bool is_in_bridge[MAX_COMPOUND_PARTICLES];	// TODO: implement this
	float potE_interim[MAX_COMPOUND_PARTICLES];
	//LJ_Ignores lj_ignore_list[MAX_COMPOUND_PARTICLES];

#ifdef LIMAKERNELDEBUGMODE
	uint32_t particle_global_ids[MAX_COMPOUND_PARTICLES];
#endif


	Float3 center_of_mass = Float3(0, 0, 0);

	// These n's can be made uint8
	uint16_t n_singlebonds = 0;
	SingleBond singlebonds[MAX_SINGLEBONDS_IN_COMPOUND];

	uint16_t n_anglebonds = 0;
	AngleBond anglebonds[MAX_ANGLEBONDS_IN_COMPOUND];

	uint16_t n_dihedrals = 0;
	DihedralBond dihedrals[MAX_DIHEDRALBONDS_IN_COMPOUND];

	uint8_t n_improperdihedrals = 0;
	ImproperDihedralBond impropers[MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND];

	int key_particle_index = -1;			// Index of particle initially closest to CoM
	float radius = 0;	// [nm] All particles in compound are PROBABLY within this radius 


	__device__ void loadMeta(Compound* compound) {
		n_particles = compound->n_particles;
		n_singlebonds = compound->n_singlebonds;
		n_anglebonds = compound->n_anglebonds;
		n_dihedrals = compound->n_dihedrals;
		n_improperdihedrals = compound->n_improperdihedrals;
	}

	__device__ void loadData(Compound* compound) {
		if (threadIdx.x < n_particles) {
			//prev_positions[threadIdx.x] = compound->prev_positions[threadIdx.x];
			atom_types[threadIdx.x] = compound->atom_types[threadIdx.x];
			//lj_ignore_list[threadIdx.x] = compound->lj_ignore_list[threadIdx.x];
			forces[threadIdx.x] = compound->forces[threadIdx.x];
			compound->forces[threadIdx.x] = Float3(0.f);

			potE_interim[threadIdx.x] = compound->potE_interim[threadIdx.x];
			//#ifdef LIMA_DEBUGMODE
			particle_global_ids[threadIdx.x] = compound->particle_global_ids[threadIdx.x];
			//#endif
		}
		else {
			//prev_positions[threadIdx.x] = Float3(-1.f);
			atom_types[threadIdx.x] = 0;
			//lj_ignore_list[threadIdx.x] = compound->lj_ignore_list[threadIdx.x];
			forces[threadIdx.x] = Float3(0.f);
			potE_interim[threadIdx.x] = 0.f;
			particle_global_ids[threadIdx.x] = 0;
		}

		for (int i = 0; (i * blockDim.x) < n_singlebonds; i++) {
			int index = i * blockDim.x + threadIdx.x;
			if (index < n_singlebonds)
				singlebonds[index] = compound->singlebonds[index];
		}
		for (int i = 0; (i * blockDim.x) < n_anglebonds; i++) {
			int index = i * blockDim.x + threadIdx.x;
			if (index < n_anglebonds)
				anglebonds[index] = compound->anglebonds[index];
		}
		for (int i = 0; (i * blockDim.x) < n_dihedrals; i++) {
			int index = i * blockDim.x + threadIdx.x;
			if (index < n_dihedrals)
				dihedrals[index] = compound->dihedrals[index];
		}
		if (threadIdx.x < n_improperdihedrals) {
			impropers[threadIdx.x] = compound->impropers[threadIdx.x];
		}
	}
};


using BondedParticlesLUT = FixedSizeMatrix<bool, MAX_COMPOUND_PARTICLES>;
using BondedParticlesLUTManager = FixedSizeMatrix<BondedParticlesLUT, 100>;









struct ParticleReference {
	ParticleReference() {}

#ifdef LIMAKERNELDEBUGMODE
	// Used by moleculebuilder only
	ParticleReference(int compound_id, int local_id_compound, int gro_id) :
		compound_id(compound_id), local_id_compound(local_id_compound), gro_id(gro_id)
	{}
	int gro_id = -1;
#else
	// Used by moleculebuilder only
	ParticleReference(int compound_id, int local_id_compound) : 
		compound_id(compound_id), local_id_compound(local_id_compound) {}
#endif

	int compound_id = -1;
	int local_id_compound = -1;

	//int global_id = -1; // For debug
};

struct CompoundBridge {
	CompoundBridge() {}	
	
	ParticleReference particle_refs[MAX_PARTICLES_IN_BRIDGE]{};
	uint8_t atom_types[MAX_PARTICLES_IN_BRIDGE]{};
	uint8_t n_particles = 0;					

	uint8_t n_singlebonds = 0;
	SingleBond singlebonds[MAX_SINGLEBONDS_IN_BRIDGE];

	uint8_t n_anglebonds = 0;
	AngleBond anglebonds[MAX_ANGLEBONDS_IN_BRIDGE];

	uint8_t n_dihedrals = 0;
	DihedralBond dihedrals[MAX_DIHEDRALBONDS_IN_BRIDGE];

	uint8_t n_improperdihedrals = 0;
	ImproperDihedralBond impropers[MAX_IMPROPERDIHEDRALBONDS_IN_BRIDGE];

	uint16_t compound_id_left{};
	uint16_t compound_id_right{};

	// -------------- Device functions ------------- //
	__device__ void loadMeta(CompoundBridge* bridge) {
		n_particles = bridge->n_particles;
		n_singlebonds = bridge->n_singlebonds;
		n_anglebonds = bridge->n_anglebonds;
		n_dihedrals = bridge->n_dihedrals;
		n_improperdihedrals = bridge->n_improperdihedrals;

		compound_id_left = bridge->compound_id_left;
		compound_id_right = bridge->compound_id_right;
	}
	__device__ void loadData(CompoundBridge* bridge) {
		if (threadIdx.x < n_particles) {
			atom_types[threadIdx.x] = bridge->atom_types[threadIdx.x];
			particle_refs[threadIdx.x] = bridge->particle_refs[threadIdx.x];
		}
		
		for (int i = 0; (i * blockDim.x) < n_singlebonds; i++) {	// TODO: Remove for loop, make static check that nparticles is larger than max bonds!
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
			compound_bridges[i] = bridges[i];//CompoundBridge{ bridges[i], false };
		}
		n_bridges = static_cast<int>(bridges.size());
	}


	CompoundBridge compound_bridges[COMPOUNDBRIDGES_IN_BUNDLE];
	int n_bridges = 0;
};






struct ForceField_NB {
	struct ParticleParameters {	//Nonbonded
		float mass = -1;		//[kg/mol]	or 
		float sigma = -1;
		float epsilon = -1;		// J/mol [kg*nm^2 / s^2]
	};

	ParticleParameters particle_parameters[MAX_ATOM_TYPES];
};






struct CompoundGridNode {
	__host__ void addNearbyCompound(int16_t compound_id);
	__host__ void addAssociatedCompound(int16_t compound_id);



	// Contains only compounds that are CLOSEST to this specific node
	static const int max_associated_compounds = 4;
	int16_t associated_ids[max_associated_compounds]{};
	int16_t n_associated_compounds = 0;

	// Compounds that are near this specific node
	// A particle belonging to this node coord, can iterate through this list
	// to find all appropriate nearby compounds;
	static const int max_nearby_compounds = 64;
	int16_t nearby_compound_ids[max_nearby_compounds]{};	// MAX_COMPOUNDS HARD LIMIT
	int16_t n_nearby_compounds = 0;
};

// Class for signaling compound origo's and quickly searching nearby compounds using a coordinate on the grid
class CompoundGrid : public BoxGrid<CompoundGridNode> {
public:
	CompoundGrid(){}

	void reset() {

	}

	// Function called by compound kernels before nlist update
	__device__ void signalOrigo(Coord origo) {
		if (threadIdx.x == 0) {
			//compound_origos[blockIdx.x] = origo;
		}
	}

	// Returns pointer to the list of origos on device. Receive on host
	//__host__ Coord* getOrigosPtr() { return compound_origos; }

	// Returns ptr to the list of CompoundGridNodes on device. Transmit to from host
	__host__ CompoundGridNode* getGridnodesPtr() { return blocks; }	// DOes this even make sense? Just use getBlockPtr?

//private:
	// Each compound in kernel will transmit their origo. Will be transferring from device to host by nlist
	//Coord compound_origos[MAX_COMPOUNDS];
};







#pragma warning (pop)