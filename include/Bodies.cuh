#pragma once
#include <iostream>

#include "Constants.cuh"
#include "LimaTypes.cuh"
//#include "QuantomTypes.cuh"






// Hyper-fast objects for kernel, so we turn all safety off here!
#pragma warning (push)
#pragma warning (disable : 26495)




struct Solvent {
	__host__ __device__ Solvent() {}
	__host__ Solvent(Float3 pos, Float3 pos_tsub1) : pos(pos), pos_tsub1(pos_tsub1) {}

	Float3 pos;
	Float3 pos_tsub1;
};






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
	PairBond(int id1, int id2, float b0, float kb) : b0(b0), kb(kb) {
		// This is only for loading the forcefield, so the ID's refers to id's given in .conf file!
		atom_indexes[0] = id1;
		atom_indexes[1] = id2;
	}

	PairBond( uint32_t particleindex_a, uint32_t particleindex_b) {
		atom_indexes[0] = particleindex_a;
		atom_indexes[1] = particleindex_b;
	}
	PairBond(float ref_dist, float kb, uint32_t particleindex_a, uint32_t particleindex_b, bool invertLJ) :
		//reference_dist(ref_dist) {
		b0(ref_dist), kb(kb), invertLJ(invertLJ) {
		atom_indexes[0] = particleindex_a;
		atom_indexes[1] = particleindex_b;
	}
	
	//uint32_t bond_index;	
	//double reference_dist;
	float b0 = 0.f;
	float kb = 0.f;
	uint32_t atom_indexes[2] = {0,0};	// Relative to the compund - NOT ABSOLUTE INDEX. Used in global table with compunds start-index
	const static int n_atoms = 2;
	bool invertLJ = false;		// When sharing multiple bonded connections
};

struct AngleBond {
	AngleBond() {}
	AngleBond(int id1, int id2, int id3, float theta_0, float k_theta, bool invertLJ) : 
		theta_0(theta_0), k_theta(k_theta), invertLJ(invertLJ) {
		atom_indexes[0] = id1;
		atom_indexes[1] = id2;
		atom_indexes[2] = id3;
	}

	float theta_0 = 0.f;
	float k_theta = 0.f;
	uint32_t atom_indexes[3] = {}; // i,j,k angle between i and k
	const static int n_atoms = 3;
	bool invertLJ = false;		// When sharing multiple bonded connections
};

struct DihedralBond {
	DihedralBond() {}
	//DihedralBond(float phi_0, float k_phi) : phi_0(phi_0), k_phi(k_phi) {}
	DihedralBond(int id1, int id2, int id3, int id4, float phi_0, float k_phi, int n, bool invertLJ) : 
		phi_0(phi_0), k_phi(k_phi), n(n), invertLJ(invertLJ) {
		atom_indexes[0] = id1;
		atom_indexes[1] = id2;
		atom_indexes[2] = id3;
		atom_indexes[3] = id4;
	}
	float phi_0 = 0.f;
	float k_phi = 0.f;
	uint8_t n = 0;
	uint32_t atom_indexes[4] = {0,0,0,0};
	const static int n_atoms = 4;
	bool invertLJ = false;		// When sharing multiple bonded connections
};

struct GenericBond {				// ONLY used during creation, never on device!
	enum BONDTYPES { SINGLE, ANGLE, DIHEDRAL, PAIR };

	GenericBond() {}
	GenericBond(ParticleRef* particle_refs, int n) {
		for (int i = 0; i < n; i++) {
			particles[n_particles++] = particle_refs[i];


			if (compound_ids[0] == -1) { 
				compound_ids[0] = particles[i].compound_id; 
			}
			else if (compound_ids[0] != particles[i].compound_id) {
				compound_ids[1] = particles[i].compound_id;
			}
		}
		if (compound_ids[0] > compound_ids[1] && compound_ids[1] != -1) { std::swap(compound_ids[0], compound_ids[1]); }
	}

	bool spansTwoCompounds() {	// Simply check if any id has been put into index 1
		return (compound_ids[1] != -1);
	}
	bool allParticlesExist() {	// Returns false if any particle in bond does not exist, for example when ignoring hydrogens
		for (int i = 0; i < n_particles; i++) {
			if (particles[i].local_id_compound == -1)
				return false;
		}
		return true;
	}


	/*
	template<typename T>
	void assignBond(T* bond) {
		(T*) bond_ref = bond;
		for (int i = 0; i < bond->n_atoms; i++) {



		}
		//compound_ids[0] = compound_ids[0] == -1 ? bond->
	}
	*/


	int compound_ids[2] = { -1,-1 };		// Can either span 1 or 2 compounds. If more, then undefined behaviour

	void* bond_ref = nullptr;

	ParticleRef particles[4] = {};
	int n_particles = 0;
};



struct LJ_Ignores {	// Each particle is associated with 1 of these.
	LJ_Ignores() {
		for (int i = 0; i < max_ids; i++)
			global_ids[i] = UINT32_MAX;
	}

	static const int max_ids = 32;
	uint8_t local_ids[max_ids] = {};
	uint8_t compound_ids[max_ids] = {};

	uint32_t global_ids[max_ids] = {};
	uint8_t n_ignores = 0;


	uint16_t doublyconnected_id = 0;
	
	__host__ void addIgnoreTarget(int target_global_id) {
		global_ids[n_ignores] = target_global_id;
		n_ignores++;
	}

	__host__ void addIgnoreTarget(uint8_t target_id, uint8_t compound_id) {
		if (n_ignores == max_ids) {
			printf("Failed to add ignore target!\n");
			exit(0);
		}
		for (int i = 0; i < n_ignores; i++) {		// If target particle is already ignored, do nothing
			if (local_ids[i] == target_id && compound_ids[i] == compound_id)
				return;
		}

		local_ids[n_ignores] = target_id;
		compound_ids[n_ignores] = compound_id;
		n_ignores++;
	}
	
	__device__ bool ignore(uint8_t query_local_id, uint8_t query_compound_id) {
		for (int i = 0; i < n_ignores; i++) {
			if (local_ids[i] == query_local_id && compound_ids[i] == query_compound_id)
				return true;
		}
		return false;
	}

	__host__ bool checkAlreadyConnected(uint16_t global_id) {
		for (int i = 0; i < n_ignores; i++) {
			if (global_ids[i] == global_id)
				return true;
		}
		return false;
	}
	

	__host__ void assignDoublyConnectedID(uint16_t id) {
		doublyconnected_id = id;
	}

	__device__ bool doublyConnected(uint16_t dc_id) {
		return (dc_id == doublyconnected_id && doublyconnected_id != 0);
	}


	// For system 3
	__host__  void init() {
		for (int i = 0; i < max_ids; i++)
			global_ids[i] = UINT32_MAX;
	}
	__host__ bool isAlreadyConnected(uint32_t global_id) {
		for (int i = 0; i < n_ignores; i++) {
			if (global_ids[i] == global_id)
				return true;
		}
		return false;
	}

	__host__ void addConnection(uint32_t global_id) {
		for (int i = 0; i < max_ids; i++) {
			if (global_ids[i] == global_id) {
				//printf("MC found!\n");
				return;
			}
				
			if (global_ids[i] == UINT32_MAX) {
				global_ids[i] = global_id;
				n_ignores++;
				return;
			}
		}
		printf("Out of bounds force LJ ignores!\n");
		exit(1);
	}
};

// ------------------------------------------------- COMPOUNDS ------------------------------------------------- //


class NeighborList {
public:
	enum NEIGHBOR_TYPE {COMPOUND, SOLVENT};


	__host__ bool addId(uint16_t new_id, NEIGHBOR_TYPE nt) {
		switch (nt)
		{
		case NeighborList::COMPOUND:
			if (n_compound_neighbors < NEIGHBORLIST_MAX_COMPOUNDS) {
				neighborcompound_ids[n_compound_neighbors++] = new_id;
				return true;
			}
			printf("\nFailed to insert compound neighbor id %d!\n", new_id);
			exit(1);
			break;

		case NeighborList::SOLVENT:
			if (n_solvent_neighbors < NEIGHBORLIST_MAX_SOLVENTS) {
				neighborsolvent_ids[n_solvent_neighbors++] = new_id;
				//if (associated_id==0 || new_id == 0)
					//printf("%d added id %d\n", associated_id, new_id);
				return true;
			}

			printf("\nFailed to insert solvent neighbor id %d of %d\n", new_id, n_solvent_neighbors);
			exit(1);
			break;

		default:
			break;
		}
		exit(1);
		return false;
	}
	__host__ bool removeId(uint16_t neighbor_id, NEIGHBOR_TYPE nt) {

		switch (nt) {
		case NeighborList::SOLVENT:
			for (int i = 0; i < n_solvent_neighbors; i++) {
				if (neighborsolvent_ids[i] == neighbor_id) {
					neighborsolvent_ids[i] = neighborsolvent_ids[n_solvent_neighbors - 1];
					n_solvent_neighbors--;
					return true;
				}
			}
			printf("\n%d Failed to remove neighbor solvent ID: %d of %d total IDs. max+1: %d\n", associated_id, neighbor_id, n_solvent_neighbors, neighborcompound_ids[n_solvent_neighbors]);
			for (int i = 0; i < n_solvent_neighbors; i++) {
				//printf("%d\n", neighborsolvent_ids[i]);
			}
			break;
		case NeighborList::COMPOUND:
			for (int i = 0; i < n_compound_neighbors; i++) {
				if (neighborcompound_ids[i] == neighbor_id) {
					neighborcompound_ids[i] = neighborcompound_ids[n_compound_neighbors - 1];
					n_compound_neighbors--;
					return true;
				}
			}
			printf("\nFailed to remove neighbor compound %d of %d\n", neighbor_id, n_compound_neighbors);
			break;
		default:
			printf("Faulty neighbortype\n");
			break;
		}
		
		return false;
	}

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


struct Compound {
	__host__ __device__  Compound() {}	// {}


	uint8_t n_particles = 0;					// MAX 256 particles!!!!0
	Float3 prev_positions[MAX_COMPOUND_PARTICLES];;			// Should this really belong to the compound and not the box?
	Float3 forces[MAX_COMPOUND_PARTICLES];					// Carries forces from bridge_kernels
	uint8_t atom_types[MAX_COMPOUND_PARTICLES];
	uint8_t atom_color_types[MAX_COMPOUND_PARTICLES];	// For drawing pretty spheres :)

	LJ_Ignores lj_ignore_list[MAX_COMPOUND_PARTICLES];

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
	float confining_particle_sphere = 0;		// All particles in compound are PROBABLY within this radius


	//---------------------------------------------------------------------------------//

	
	__host__ void init() {	// Only call this if the compound has already been assigned particles & bonds
		center_of_mass = calcCOM();
		//printf("")
		//radius = singlebonds[0].reference_dist * n_particles * 0.5f;
		//center_of_mass.print('C');
		//printf("Radius %f\n", radius);
	}

	__host__ Float3 calcCOM() {
		Float3 com;
		for (int i = 0; i < n_particles; i++) {
			com += (prev_positions[i] * (1.f / (float)n_particles));
			//com += (particles[i].pos_tsub1 * (1.f / (float) n_particles));
		}
			
		return com;
	}
	/*
	__host__ bool intersects(Compound a) {
		return (a.center_of_mass - center_of_mass).len() < (a.radius + radius + max_LJ_dist);
	}*/

	


	__host__ void addParticle(int atomtype_id, Float3 pos) {
		if (n_particles == MAX_COMPOUND_PARTICLES) {
			printf("ERROR: Cannot add particle to compound!\n");
			exit(1);
		}

		atom_types[n_particles] = atomtype_id;
		prev_positions[n_particles] = pos;
		n_particles++;
	}
	__host__ void addParticle(int atomtype_id, Float3 pos, int atomtype_color_id, int global_id) {
		if (n_particles == MAX_COMPOUND_PARTICLES) {
			printf("ERROR: Cannot add particle to compound!\n");
			exit(1);
		}

		atom_types[n_particles] = atomtype_id;
		prev_positions[n_particles] = pos;
		atom_color_types[n_particles] = atomtype_color_id;
#ifdef LIMA_DEBUGMODE
		particle_global_ids[n_particles] = global_id;
		//printf("%d global id\n", global_id);
#endif
		n_particles++;
	}
	__host__ bool hasRoomForRes() {					// TODO: Implement, that it checks n atoms in res
		return ((int)n_particles + MAX_ATOMS_IN_RESIDUE) <= MAX_COMPOUND_PARTICLES;
	}
	__host__ void calcParticleSphere() {
		Float3 com = calcCOM();// calcCOM(compound);

		//float furthest = LONG_MIN;
		float furthest = FLT_MIN;
		//float closest = LONG_MAX;
		float closest = FLT_MAX;
		int closest_index = 0;

		for (int i = 0; i < n_particles; i++) {
			float dist = (prev_positions[i] - com).len();
			closest_index = dist < closest ? i : closest_index;
			closest = std::min(closest, dist);
			furthest = std::max(furthest, dist);
		}

		key_particle_index = closest_index;
		confining_particle_sphere = furthest;
	}
	//---------------------------------------------------------------------------------//

	__device__ void loadMeta(Compound* compound) {
		n_particles = compound->n_particles;
		n_singlebonds = compound->n_singlebonds;
		n_anglebonds = compound->n_anglebonds;
		n_dihedrals = compound->n_dihedrals;
	}
	__device__ void loadData(Compound* compound) {
		if (threadIdx.x < n_particles) {
			prev_positions[threadIdx.x] = compound->prev_positions[threadIdx.x];
			atom_types[threadIdx.x] = compound->atom_types[threadIdx.x];
			lj_ignore_list[threadIdx.x] = compound->lj_ignore_list[threadIdx.x];
			forces[threadIdx.x] = compound->forces[threadIdx.x];
			compound->forces[threadIdx.x] = Float3(0.f);
//#ifdef LIMA_DEBUGMODE
			particle_global_ids[threadIdx.x] = compound->particle_global_ids[threadIdx.x]; 
//#endif
		}
		else {
			prev_positions[threadIdx.x] = Float3(-1.f);
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


	//CompoundState* compound_state_ptr;
	//CompoundNeighborList* compound_neighborlist_ptr;


};



struct CompoundBridgeBundleCompact;
struct Molecule {
	Molecule();
	int n_compounds = 0;
	Compound* compounds;
	CompoundBridgeBundleCompact* compound_bridge_bundle;
	//CompoundBridgeBundle compound_bridge_bundle;	// Special compound, for special kernel. For now we only need one
	uint32_t n_atoms_total = 0;

	Float3 calcCOM();





	~Molecule() {
		//printf("Deleting\n");		// Huh, this deletes too early. I better implement properly at some point.
		//delete[] compounds;
	}
};

struct CompoundBridge {
	CompoundBridge() {}
	CompoundBridge(int id_left, int id_right): compound_id_left(id_left), compound_id_right(id_right) {
	}

	int compound_id_left, compound_id_right;
	ParticleRef particle_refs[MAX_PARTICLES_IN_BRIDGE];
	uint8_t atom_types[MAX_PARTICLES_IN_BRIDGE];
	int n_particles = 0;




	GenericBond bonds[MAX_DIHEDRALBONDS_IN_BRIDGE*4];
	int n_bonds = 0;
	

	PairBond singlebonds[MAX_SINGLEBONDS_IN_BRIDGE];
	int n_singlebonds = 0;
	AngleBond anglebonds[MAX_ANGLEBONDS_IN_BRIDGE];
	int n_anglebonds = 0;
	DihedralBond dihedrals[MAX_DIHEDRALBONDS_IN_BRIDGE];
	int n_dihedrals = 0;

	bool bondBelongsInBridge(GenericBond* bond) {
		return (compound_id_left == bond->compound_ids[0] && compound_id_right == bond->compound_ids[1]);
	}
	bool particleAlreadyStored(ParticleRef* p_ref) {
		for (int i = 0; i < n_particles; i++) {
			if (particle_refs[i] == *p_ref) {
				return true;
			}
		}
		return false;
	}
	void addParticle(ParticleRef* particle_ref, Molecule* molecule) {
		if (n_particles == MAX_PARTICLES_IN_BRIDGE) {
			printf("Too many particles in bridge\n");
			exit(0);
		}
		particle_ref->bridge_id = 0;
		particle_ref->local_id_bridge = n_particles;
		atom_types[n_particles] = molecule->compounds[particle_ref->compound_id].atom_types[particle_ref->local_id_compound];
		particle_refs[n_particles] = *particle_ref;
		n_particles++;
		//printf("Adding particle with global id: %d\n", particle_ref->global_id);
	}

	void addBondParticles(GenericBond* bond, Molecule* molecule) {
		for (int p = 0; p < bond->n_particles; p++) {
			if (!particleAlreadyStored(&bond->particles[p])) {
				addParticle(&bond->particles[p], molecule);
			}				
		}
	}

	template <typename T>
	void localizeIDs(T* bond, int n) {
		for (int p = 0; p < n; p++) {						// First reassign the global indexes of the bond with local indexes of the bridge
			for (int i = 0; i < n_particles; i++) {
				if (bond->atom_indexes[p] == particle_refs[i].global_id) {
					bond->atom_indexes[p] = particle_refs[i].local_id_bridge;
					break;
				}
			}
		}
	}
	/*void addSinglebond(PairBond pb) {
		localizeIDs(&pb, 2);
		singlebonds[n_singlebonds++] = pb;
		//printf("Singlebond added %d %d\n", singlebonds[n_singlebonds - 1].atom_indexes[0], singlebonds[n_singlebonds - 1].atom_indexes[1]);
	}
	void addAnglebond(AngleBond ab) {
		localizeIDs(&ab, 3);
		anglebonds[n_anglebonds++] = ab;
		//printf("Anglebond added %d %d %d\n", anglebonds[n_anglebonds - 1].atom_indexes[0], anglebonds[n_anglebonds - 1].atom_indexes[1], anglebonds[n_anglebonds - 1].atom_indexes[2]);
	}*/
	void addGenericBond(PairBond pb) {
		if (n_singlebonds == MAX_SINGLEBONDS_IN_BRIDGE) {
			printf("Cannot add bond to bridge\n");
			exit(0);
		}
		localizeIDs(&pb, 2);
		singlebonds[n_singlebonds++] = pb;
	}
	void addGenericBond(AngleBond ab) {
		if (n_anglebonds == MAX_ANGLEBONDS_IN_BRIDGE) {
			printf("Cannot add angle to bridge\n");
			exit(0);
		}
		localizeIDs(&ab, 3);
		anglebonds[n_anglebonds++] = ab;
	}
	void addGenericBond(DihedralBond db) {
		if (n_dihedrals == MAX_DIHEDRALBONDS_IN_BRIDGE) {
			printf("Cannot add dihedral to bridge\n");
			exit(0);
		}
		localizeIDs(&db, 4);
		dihedrals[n_dihedrals++] = db;
	}

};




struct CompoundBridgeBundle {
	//ParticleRef particles[MAX_COMPOUND_PARTICLES * 2];
	//int n_particles = 0;
	CompoundBridge compound_bridges[COMPOUNDBRIDGES_IN_BUNDLE];
	int n_bridges = 0;

	bool addBridge(int left_c_id, int right_c_id) {
		if (left_c_id > right_c_id) { std::swap(left_c_id, right_c_id); }
		if (n_bridges == COMPOUNDBRIDGES_IN_BUNDLE)
			return false;
		compound_bridges[n_bridges++] = CompoundBridge(left_c_id, right_c_id);
		return true;
	}

	CompoundBridge* getBelongingBridge(GenericBond* bond) {
		for (int i = 0; i < n_bridges; i++) {
			CompoundBridge* bridge = &compound_bridges[i];
			if (bridge->bondBelongsInBridge(bond)) {
				return bridge;
			}
		}
	}

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
	CompoundBridgeCompact(CompoundBridge* bridge) {
		n_particles = bridge->n_particles;
		
		for (int i = 0; i < n_particles; i++) {
			particle_refs[i] = ParticleRefCompact(bridge->particle_refs[i]);
			atom_types[i] = bridge->atom_types[i];
		}
		n_singlebonds = bridge->n_singlebonds;
		for (int i = 0; i < n_singlebonds; i++) {
			singlebonds[i] = bridge->singlebonds[i];
		}
		n_anglebonds = bridge->n_anglebonds;
		for (int i = 0; i < n_anglebonds; i++) {
			anglebonds[i] = bridge->anglebonds[i];
		}
		n_dihedrals = bridge->n_dihedrals;
		for (int i = 0; i < n_dihedrals; i++) {
			dihedrals[i] = bridge->dihedrals[i];
		}
		printf("Loading bridge with %d particles %d bonds %d angles %d dihedrals\n", n_particles, n_singlebonds, n_anglebonds, n_dihedrals);
	}
	
	
	ParticleRefCompact particle_refs[MAX_PARTICLES_IN_BRIDGE];
	uint8_t atom_types[MAX_PARTICLES_IN_BRIDGE];
	uint8_t n_particles = 0;					

	uint16_t n_singlebonds = 0;
	PairBond singlebonds[MAX_SINGLEBONDS_IN_BRIDGE];

	uint16_t n_anglebonds = 0;
	AngleBond anglebonds[MAX_ANGLEBONDS_IN_BRIDGE];

	uint16_t n_dihedrals = 0;
	DihedralBond dihedrals[MAX_DIHEDRALBONDS_IN_BRIDGE];


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
	CompoundBridgeBundleCompact(CompoundBridgeBundle* bundle);
	CompoundBridgeCompact compound_bridges[COMPOUNDBRIDGES_IN_BUNDLE];
	int n_bridges = 0;
};


#pragma warning (pop)





