#pragma once
//#include <iostream>

#include "Constants.h"
#include "LimaTypes.cuh"

#include <memory>
#include <vector>
#include <array>
#include <cuda_fp16.h>
#include <assert.h>
// Hyper-fast objects for kernel, so we turn all safety off here!
#pragma warning (push)
#pragma warning (disable : 26495)








// ------------------------------------------------- BondTypes ------------------------------------------------- //


struct SingleBond {
	struct Parameters {
		float b0 = 0.f;	// [lm]
		float kb = 0.f;	// [J/(mol*lm^2)] // V(bond) = 1/2 * kb * (r - b0)^2

		bool operator==(const Parameters&) const = default;
		bool HasZeroParam() const { return kb == 0.f; }
	};

	SingleBond(){}
	SingleBond(std::array<uint8_t, 2> ids, const Parameters&);

	Parameters params;
	uint8_t atom_indexes[2] = {0,0};	// Relative to the compund
	const static int nAtoms = 2;
};


struct AngleUreyBradleyBond {
	struct Parameters {
		float theta0 = 0.f;	// [rad]
		float kTheta = 0.f;	// [J/mol/rad^2]
		float ub0 = 0.f;	// [lm]
		float kUB = 0.f;	// [J/mol/lm^2]

		bool operator==(const Parameters&) const = default;
		bool HasZeroParam() const { return kTheta == 0.f; }
	};

	AngleUreyBradleyBond() {}
	AngleUreyBradleyBond(std::array<uint8_t, 3> ids, const Parameters&);

	Parameters params;
	uint8_t atom_indexes[3] = {0,0,0}; // i,j,k angle between i and k
	const static int nAtoms = 3;
};

struct DihedralBond {
	struct Parameters {
		float phi_0;		// [rad]
		float k_phi;		// [J/mol/rad^2]
		float n;			// [multiplicity] n parameter, how many energy equilibriums does the dihedral have // OPTIMIZE: maybe float makes more sense, to avoid conversion in kernels?

		bool operator==(const Parameters&) const = default;
		bool HasZeroParam() const { return k_phi == 0.f; }
	};
	const static int nAtoms = 4;
	DihedralBond() {}
	DihedralBond(std::array<uint8_t, 4> ids, const Parameters&);
	
	Parameters params;
	uint8_t atom_indexes[4] = {0,0,0,0};
};

struct ImproperDihedralBond {
	struct Parameters {
		float psi_0 = 0.f;	// [rad]
		float k_psi = 0.f;	// [J/mol/rad^2]

		bool operator==(const Parameters&) const = default;
		bool HasZeroParam() const { return k_psi == 0.f; }
	};

	ImproperDihedralBond() {}
	ImproperDihedralBond(std::array<uint8_t, 4> ids, const Parameters&);

	Parameters params;

	uint8_t atom_indexes[4] = { 0,0,0,0 };
	const static int nAtoms = 4;
};


// ------------------------------------------------- COMPOUNDS ------------------------------------------------- //







struct CompoundCoords {
	__device__ void loadData(const CompoundCoords& coords) {
		if (threadIdx.x == 0) { origo = coords.origo; };
		rel_positions[threadIdx.x] = coords.rel_positions[threadIdx.x];
	}
	
	NodeIndex origo{};								// [nm]
	Coord rel_positions[MAX_COMPOUND_PARTICLES];	// [lm]
};


//namespace CompoundcoordsCircularQueueUtils {
//	static const int queueLen = 3;
//
//	__host__ __device__ inline CompoundCoords* getCoordarrayRef(CompoundCoords* queue, uint32_t step, uint32_t compound_index) {
//		const int index0_of_currentstep_coordarray = (step % queueLen) * MAX_COMPOUNDS;
//		return &queue[index0_of_currentstep_coordarray + compound_index];
//	}
//};











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



















// Instead of having a single key_particle and an single radius, we now have multiple
struct CompoundInteractionBoundary {
	static const int k = 2;

	float radii[k];	// [nm]
	int key_particle_indices[k];
};

struct alignas(4) CompoundCompact {
	__host__ __device__ CompoundCompact() {}

	alignas(4) uint8_t atom_types[MAX_COMPOUND_PARTICLES];
	int n_particles = 0;

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

	__device__ void loadMeta(const CompoundCompact* const compound) {
		n_particles = compound->n_particles;
		n_singlebonds = compound->n_singlebonds;
		n_anglebonds = compound->n_anglebonds;
		n_dihedrals = compound->n_dihedrals;
		n_improperdihedrals = compound->n_improperdihedrals;
		n_bonded_compounds = compound->n_bonded_compounds;
	}

	__device__ void loadData(const CompoundCompact* const compound) {
		if (threadIdx.x < n_particles) {
			atom_types[threadIdx.x] = compound->atom_types[threadIdx.x];

			#ifdef LIMAKERNELDEBUGMODE
			particle_global_ids[threadIdx.x] = compound->particle_global_ids[threadIdx.x];
			#endif
		}
	}
};

const int MAX_SINGLEBONDS_IN_COMPOUND = MAX_COMPOUND_PARTICLES + 4;	// Due to AA such a TRP, thhere might be more bonds than atoms
const int MAX_ANGLEBONDS_IN_COMPOUND = 128;
const int MAX_DIHEDRALBONDS_IN_COMPOUND = 128 + 64 + 64 + 32;
const int MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND = 32;


struct CompoundInterimState {
	// Interims from the bridgekernel to compoundkernel
	float potE_interim[MAX_COMPOUND_PARTICLES];
	Float3 forces_interim[MAX_COMPOUND_PARTICLES];	// [GN/mol]

	// Used specifically for Velocity Verlet stormer, and ofcourse kinE fetching
	Float3 forces_prev[MAX_COMPOUND_PARTICLES];
	Float3 vels_prev[MAX_COMPOUND_PARTICLES]; // Get wierd change of outcome if i move this here??

	Coord coords[MAX_COMPOUND_PARTICLES];
};

// Rather large unique structures in global memory, that can be partly loaded when needed
struct Compound : public CompoundCompact {
	__host__ __device__ Compound() {}

	// Bonds
	SingleBond singlebonds[MAX_SINGLEBONDS_IN_COMPOUND];
	AngleUreyBradleyBond anglebonds[MAX_ANGLEBONDS_IN_COMPOUND];
	DihedralBond dihedrals[MAX_DIHEDRALBONDS_IN_COMPOUND];
	ImproperDihedralBond impropers[MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND];

	CompoundInteractionBoundary interaction_boundary;
	int centerparticle_index = -1;			// Index of particle initially closest to CoM

	uint16_t bonded_compound_ids[max_bonded_compounds];	// *2-2because it should exclude itself from both sides
	half atom_charges[MAX_COMPOUND_PARTICLES];	// [C/mol] - prolly move next to atomtypes to improve locality

	// For drawing pretty spheres :)
	char atomLetters[MAX_COMPOUND_PARTICLES];

	//bool is_in_bridge[MAX_COMPOUND_PARTICLES];	// TODO: implement this?

	int absoluteIndexOfFirstParticle = 0;
};








struct ParticleReference {
	ParticleReference() {}	// TODO: i dont think i need this.


	// Used by moleculebuilder only
	ParticleReference(int compound_id, int local_id_compound, uint8_t compoundid_local_to_bridge) :
		compound_id(compound_id), local_id_compound(local_id_compound),
		compoundid_local_to_bridge(compoundid_local_to_bridge) 
	{}

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
	AngleUreyBradleyBond anglebonds[MAX_ANGLEBONDS_IN_BRIDGE];

	uint8_t n_dihedrals = 0;
	DihedralBond dihedrals[MAX_DIHEDRALBONDS_IN_BRIDGE];

	uint8_t n_improperdihedrals = 0;
	ImproperDihedralBond impropers[MAX_IMPROPERDIHEDRALBONDS_IN_BRIDGE];

	static_assert(MAX_COMPOUNDS < UINT16_MAX, "CompoundBridge cannot handle such large compound ids");
	uint8_t n_compounds;
	uint16_t compound_ids[MAX_COMPOUNDS_IN_BRIDGE] = { UINT16_MAX,UINT16_MAX,UINT16_MAX,UINT16_MAX };


	// -------------- Device functions ------------- //
	__device__ void loadMeta(const CompoundBridge* const bridge) {
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
	__device__ void loadData(const CompoundBridge* const bridge) {
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




struct ForceField_NB {
	static const int MAX_TYPES = 64;

	struct ParticleParameters {	//Nonbonded
		float mass = -1;		//[kg/mol]	or 
		float sigma = -1;		// [lm]
		float epsilon = -1;		// [J/mol]
	};

	ParticleParameters particle_parameters[MAX_TYPES];
};


class UniformElectricField {
	Float3 field;	// [mV/limameter]

	public:
		UniformElectricField() {}
		/// <summary></summary>
		/// <param name="direction"></param>
		/// <param name="magnitude">[V/nm]</param>
		__host__ UniformElectricField(Float3 direction, float magnitude) 
			: field(direction.norm() * magnitude * KILO * LIMA_TO_NANO) {
			assert(direction.len() != 0.f);
		}
	
	/// <summary></summary>
	/// <param name="charge">[kC/mol]</param>
	/// <returns>[gigaN/mol]</returns>
	__device__ Float3 GetForce(float charge) const {
		return field * charge;
	}
};


#pragma warning (pop)
