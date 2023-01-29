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

struct ParticleRef {		 // Maybe call map instead?
	ParticleRef() {}
	ParticleRef(int g, int c, int l) : global_id(g), compound_id(c), local_id_compound(l) {}

	int global_id = -1;		// refers to the id in conf file!

	// For designated compound
	int compound_id = -1;
	int local_id_compound = -1;		// Particles id in compound

	// For designated compound_bridge
	int bridge_id = -1;
	int local_id_bridge = -1;

	__host__ inline bool operator == (const ParticleRef a) const { return (global_id == a.global_id); }
};

struct NBAtomtype {
	NBAtomtype(){}
	NBAtomtype(float m, float s, float e) : mass(m), sigma(s), epsilon(e) {}
	float mass = 0.f;			// kg / mol
	float sigma = 0.f;		// nm
	float epsilon = 0.f;		// J/mol
};


struct PairBond {	// IDS and indexes are used interchangeably here!
	PairBond(){}
	PairBond(int id1, int id2, float b0, float kb);

	PairBond(uint32_t particleindex_a, uint32_t particleindex_b) {
		atom_indexes[0] = particleindex_a;
		atom_indexes[1] = particleindex_b;
	}
	PairBond(float ref_dist, float kb, uint32_t particleindex_a, uint32_t particleindex_b) :
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
	AngleBond(int id1, int id2, int id3, float theta_0, float k_theta);

	float theta_0 = 0.f;
	float k_theta = 0.f;
	uint32_t atom_indexes[3] = {0,0,0}; // i,j,k angle between i and k
	const static int n_atoms = 3;
};

struct DihedralBond {
	DihedralBond() {}
	DihedralBond(int id1, int id2, int id3, int id4, float phi_0, float k_phi, int n);

	float phi_0 = 0.f;
	float k_phi = 0.f;
	uint8_t n = 0;								// WTF is this? xD
	uint32_t atom_indexes[4] = {0,0,0,0};
	const static int n_atoms = 4;
};

struct GenericBond {				// ONLY used during creation, never on device!
	enum BONDTYPES { SINGLE, ANGLE, DIHEDRAL, PAIR };

	GenericBond() {}
	GenericBond(ParticleRef* particle_refs, int n);

	bool spansTwoCompounds();
	bool allParticlesExist();


	int compound_ids[2] = { -1,-1 };		// Can either span 1 or 2 compounds. If more, then undefined behaviour

	void* bond_ref = nullptr;

	ParticleRef particles[4] = {};
	int n_particles = 0;
};




// ------------------------------------------------- COMPOUNDS ------------------------------------------------- //
class NeighborList {
public:
	enum NEIGHBOR_TYPE {COMPOUND, SOLVENT};


	__host__ bool addId(uint16_t new_id, NEIGHBOR_TYPE nt);
	__host__ bool removeId(uint16_t neighbor_id, NEIGHBOR_TYPE nt);

	__device__ void loadMeta(NeighborList* nl_ptr) {	// Called from thread 0
		n_compound_neighbors = nl_ptr->n_compound_neighbors;
		n_solvent_neighbors = nl_ptr->n_solvent_neighbors;
	}
	__device__ void loadData(NeighborList* nl_ptr) {
		if (threadIdx.x < n_compound_neighbors)			// DANGER Breaks when threads < mAX_COMPOUND_Ns
			neighborcompound_ids[threadIdx.x] = nl_ptr->neighborcompound_ids[threadIdx.x];
		for (int i = threadIdx.x;  i < n_solvent_neighbors; i += blockDim.x) {// Same as THREADS_PER_COMPOUNDBLOCK
			neighborsolvent_ids[i] = nl_ptr->neighborsolvent_ids[i];
			i += blockDim.x;	
		}
	}


	uint16_t neighborcompound_ids[NEIGHBORLIST_MAX_COMPOUNDS];
	uint16_t n_compound_neighbors = 0;
	uint16_t neighborsolvent_ids[NEIGHBORLIST_MAX_SOLVENTS];
	uint16_t n_solvent_neighbors = 0;
	int associated_id = -1;

};







const int MAX_PAIRBONDS = 128;
const int MAX_ANGLEBONDS = 256;
const int MAX_DIHEDRALS = 384;
struct CompoundState {							// Maybe delete this soon?
	__device__ void setMeta(int n_p) {
		n_particles = n_p;
	}
	__device__ void loadData(CompoundState* state) {
		if (threadIdx.x < n_particles)
			positions[threadIdx.x] = state->positions[threadIdx.x];
	}


	Float3 positions[MAX_COMPOUND_PARTICLES];
	uint8_t n_particles = 0;
};













struct CompoundCoords {
	__device__ void loadData(CompoundCoords& coords) {
		if (threadIdx.x == 0) { origo = coords.origo; };
		rel_positions[threadIdx.x] = coords.rel_positions[threadIdx.x];
	}
	Coord origo{};									// [nm]
	Coord rel_positions[MAX_COMPOUND_PARTICLES]{};	// [lm]

	__host__ void static copyInitialCoordConfiguration(CompoundCoords* coords,
		CompoundCoords* coords_prev, CompoundCoords* coordarray_circular_queue);
	__host__ Float3 getAbsolutePositionLM(int particle_id);
};

struct SolventCoord {
	Coord origo{};									// [nm]
	Coord rel_position{};							// [lm]

	__host__ void static copyInitialCoordConfiguration(SolventCoord* coords,
		SolventCoord* coords_prev, SolventCoord* coordarray_circular_queue);

	__host__ SolventCoord static createFromPositionNM(const Float3& solvent_pos);
	__host__ __device__ Float3 getAbsolutePositionLM() const {
		return ((origo * NANO_TO_LIMA).toFloat3() + rel_position.toFloat3());
	}
};

namespace CoordArrayQueueHelpers {
	/*__host__ static void copyInitialCoordConfiguration(CompoundCoords* coords, 
		CompoundCoords* coords_prev, CompoundCoords* coordarray_circular_queue) */

	__device__ static CompoundCoords* getCoordarrayPtr(CompoundCoords* coordarray_circular_queue, int step, int compound_index) {
		const int index0_of_currentstep_coordarray = (step % STEPS_PER_LOGTRANSFER) * MAX_COMPOUNDS;
		return &coordarray_circular_queue[index0_of_currentstep_coordarray + compound_index];
	}


	__device__ static SolventCoord* getSolventcoordPtr(SolventCoord* solventcoordarray_circular_queue, const int step, const int solvent_id) {
		const int index0_of_currentstep_coordarray = (step % STEPS_PER_LOGTRANSFER) * MAX_SOLVENTS;
		return &solventcoordarray_circular_queue[index0_of_currentstep_coordarray + solvent_id];
	}
}


//struct Solvent {
//	__host__ __device__ Solvent() {}
//	__host__ Solvent(Float3 pos, Float3 pos_tsub1) : pos(pos), pos_tsub1(pos_tsub1) {}
//
//	SolventCoord pos;
//	SolventCoord pos_tsub1;
//};

struct Compound {
	__host__ __device__  Compound() {}	// {}


	uint8_t n_particles = 0;					// MAX 256 particles!!!!0
	//Float3 prev_positions[MAX_COMPOUND_PARTICLES];;			// Should this really belong to the compound and not the box?
	Float3 forces[MAX_COMPOUND_PARTICLES];					// Carries forces from bridge_kernels
	uint8_t atom_types[MAX_COMPOUND_PARTICLES];
	uint8_t atom_color_types[MAX_COMPOUND_PARTICLES];	// For drawing pretty spheres :)

	//LJ_Ignores lj_ignore_list[MAX_COMPOUND_PARTICLES];

#ifdef LIMA_DEBUGMODE
	uint32_t particle_global_ids[MAX_COMPOUND_PARTICLES];
#endif


	Float3 center_of_mass = Float3(0, 0, 0);
	//double radius = 0;

	uint16_t n_singlebonds = 0;
	PairBond singlebonds[MAX_PAIRBONDS];

	uint16_t n_anglebonds = 0;
	AngleBond anglebonds[MAX_ANGLEBONDS];

	uint16_t n_dihedrals = 0;
	DihedralBond dihedrals[MAX_DIHEDRALS];

	int key_particle_index = 404;		// particle which started at center. Used for PBC, applyhyperpos, and neighborlist search.
	float confining_particle_sphere = 0;		// All particles in compound are PROBABLY within this radius !!!! I think it is nm now, please change to lm!!


	//---------------------------------------------------------------------------------//

	// Only call this if the compound has already been assigned particles & bonds
	/*
	__host__ bool intersects(Compound a) {
		return (a.center_of_mass - center_of_mass).len() < (a.radius + radius + max_LJ_dist);
	}*/


	__host__ bool hasRoomForRes() {					// TODO: Implement, that it checks n atoms in res
		return ((int)n_particles + MAX_ATOMS_IN_RESIDUE) <= MAX_COMPOUND_PARTICLES;
	}
	//---------------------------------------------------------------------------------//

//	__device__ void loadMeta(Compound* compound);

	__device__ void loadMeta(Compound* compound) {
		n_particles = compound->n_particles;
		n_singlebonds = compound->n_singlebonds;
		n_anglebonds = compound->n_anglebonds;
		n_dihedrals = compound->n_dihedrals;
	}

	//__device__ void loadData(Compound* compound);
	__device__ void loadData(Compound* compound) {
		if (threadIdx.x < n_particles) {
			//prev_positions[threadIdx.x] = compound->prev_positions[threadIdx.x];
			atom_types[threadIdx.x] = compound->atom_types[threadIdx.x];
			//lj_ignore_list[threadIdx.x] = compound->lj_ignore_list[threadIdx.x];
			forces[threadIdx.x] = compound->forces[threadIdx.x];
			compound->forces[threadIdx.x] = Float3(0.f);
			//#ifdef LIMA_DEBUGMODE
			particle_global_ids[threadIdx.x] = compound->particle_global_ids[threadIdx.x];
			//#endif
		}
		else {
			//prev_positions[threadIdx.x] = Float3(-1.f);
			atom_types[threadIdx.x] = 0;
			//lj_ignore_list[threadIdx.x] = compound->lj_ignore_list[threadIdx.x];
			forces[threadIdx.x] = Float3(0.f);
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
	}

};

// Helper class containing some compound information we will not bring with us to the device
// This means, all functions are host only!
struct Compound_Carrier : public Compound {
	Compound_Carrier(int id) : compound_id(id) {}

	int compound_id = -1;
	CompoundState state;	// [nm], Absolute values 
	CompoundState state_tsub1;		// [nm], Absolute values 


	void addParticle(int atomtype_id, Float3 pos);// Is this ever used?
	void addParticle(int atomtype_id, Float3 pos, int atomtype_color_id, int global_id);
	void calcParticleSphere();

	void init() { calcCOM(); }
	Float3 calcCOM();

};

using BondedParticlesLUT = FixedSizeMatrix<bool, MAX_COMPOUND_PARTICLES>;
using BondedParticlesLUTManager = FixedSizeMatrix<BondedParticlesLUT, 100>;


struct CompoundBridgeBundleCompact;
struct CompoundCollection {
	CompoundCollection();
	int n_compounds = 0;
	//Compound* compounds = nullptr;
	std::vector<Compound_Carrier> compounds;
	CompoundBridgeBundleCompact* compound_bridge_bundle = nullptr;
	//CompoundBridgeBundle compound_bridge_bundle;	// Special compound, for special kernel. For now we only need one
	uint32_t n_atoms_total = 0;

	Float3 calcCOM();

	BondedParticlesLUTManager* bonded_particles_lut_manager = nullptr;

	~CompoundCollection() {
		//printf("Deleting\n");		// Huh, this deletes too early. I better implement properly at some point.
		//delete[] compounds;
	}
};

struct CompoundBridge {
	CompoundBridge() {}
	CompoundBridge(uint16_t id_left, uint16_t id_right): compound_id_left(id_left), compound_id_right(id_right) {
	}

	uint16_t compound_id_left{};
	uint16_t compound_id_right{};
	ParticleRef particle_refs[MAX_PARTICLES_IN_BRIDGE];
	uint8_t atom_types[MAX_PARTICLES_IN_BRIDGE];
	int n_particles = 0;




	GenericBond bonds[MAX_DIHEDRALBONDS_IN_BRIDGE*4];
	int n_bonds = 0;
	

	PairBond singlebonds[MAX_SINGLEBONDS_IN_BRIDGE];
	uint16_t n_singlebonds = 0;
	AngleBond anglebonds[MAX_ANGLEBONDS_IN_BRIDGE];
	uint16_t n_anglebonds = 0;
	DihedralBond dihedrals[MAX_DIHEDRALBONDS_IN_BRIDGE];
	uint16_t n_dihedrals = 0;

	bool bondBelongsInBridge(GenericBond* bond) const;
	bool particleAlreadyStored(ParticleRef* p_ref);
	void addParticle(ParticleRef* particle_ref, CompoundCollection* molecule);

	void addBondParticles(GenericBond* bond, CompoundCollection* molecule);



	void addGenericBond(PairBond pb);
	void addGenericBond(AngleBond ab);
	void addGenericBond(DihedralBond db);

private:

	template <typename T>
	void localizeIDs(T* bond, int n) {
		// First reassign the global indexes of the bond with local indexes of the bridge
		for (int p = 0; p < n; p++) {						
			for (int i = 0; i < n_particles; i++) {
				if (bond->atom_indexes[p] == particle_refs[i].global_id) {
					bond->atom_indexes[p] = particle_refs[i].local_id_bridge;
					break;
				}
			}
		}
	}
};




struct CompoundBridgeBundle {
	//ParticleRef particles[MAX_COMPOUND_PARTICLES * 2];
	//int n_particles = 0;
	CompoundBridge compound_bridges[COMPOUNDBRIDGES_IN_BUNDLE];
	int n_bridges = 0;

	bool addBridge(uint16_t left_c_id, uint16_t right_c_id);

	CompoundBridge* getBelongingBridge(GenericBond* bond);

};










struct ParticleRefCompact {
	ParticleRefCompact() {}
	ParticleRefCompact(ParticleRef pref) : compound_id(pref.compound_id), local_id_compound(pref.local_id_compound),
	global_id(pref.global_id) {}

	int compound_id = -1;
	int local_id_compound = -1;

	int global_id = -1; // temp
};

struct CompoundBridgeCompact {
	CompoundBridgeCompact() {}
	CompoundBridgeCompact(CompoundBridge* bridge, bool verbose);
	
	
	ParticleRefCompact particle_refs[MAX_PARTICLES_IN_BRIDGE]{};
	uint8_t atom_types[MAX_PARTICLES_IN_BRIDGE]{};
	uint8_t n_particles = 0;					

	uint16_t n_singlebonds = 0;
	PairBond singlebonds[MAX_SINGLEBONDS_IN_BRIDGE];

	uint16_t n_anglebonds = 0;
	AngleBond anglebonds[MAX_ANGLEBONDS_IN_BRIDGE];

	uint16_t n_dihedrals = 0;
	DihedralBond dihedrals[MAX_DIHEDRALBONDS_IN_BRIDGE];

	uint16_t compound_id_left{};
	uint16_t compound_id_right{};

	// -------------- Device functions ------------- //
	__device__ void loadMeta(CompoundBridgeCompact* bridge) {
		n_particles = bridge->n_particles;
		n_singlebonds = bridge->n_singlebonds;
		n_anglebonds = bridge->n_anglebonds;
		n_dihedrals = bridge->n_dihedrals;
	}
	__device__ void loadData(CompoundBridgeCompact* bridge) {
		if (threadIdx.x < n_particles) {
			atom_types[threadIdx.x] = bridge->atom_types[threadIdx.x];
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
	}
};

struct CompoundBridgeBundleCompact {
	CompoundBridgeBundleCompact() {}
	CompoundBridgeBundleCompact(CompoundBridgeBundle* bundle, bool verbose=false);
	CompoundBridgeCompact compound_bridges[COMPOUNDBRIDGES_IN_BUNDLE];
	int n_bridges = 0;
};



struct SlidingActionGroupModule {
	uint8_t pairs [64][64] = {};
};





struct ForceField_NB {
	struct ParticleParameters {	//Nonbonded
		float mass = -1;		//[kg/mol]	or 
		float sigma = -1;
		float epsilon = -1;		// J/mol [kg*nm^2 / s^2]
	};

	ParticleParameters particle_parameters[MAX_ATOM_TYPES];
};


#pragma warning (pop)





