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

namespace Bondtypes {
	struct SingleBond {
		struct Parameters {
			float b0 = 0.f;	// [nm]
			float kb = 0.f;	// [J/(mol*nm^2)] // V(bond) = 1/2 * kb * (r - b0)^2
			
			bool operator==(const Parameters&) const = default;
			bool HasZeroParam() const { return kb == 0.f; }
								
			/// <param name="b0">[nm]</param>
			/// <param name="kB">[kJ/mol/nm^2]</param>
			static Parameters CreateFromCharmm(float b0, float kB);
		};

		SingleBond() {}
		SingleBond(std::array<uint8_t, 2> ids, const Parameters&);

		Parameters params;
		uint8_t atom_indexes[2] = { 0,0 };	// Relative to the compund
		const static int nAtoms = 2;
	};

	struct PairBond {
		struct Parameters {
			float sigma = 0.f;
			float epsilon = 0.f;

			bool operator==(const Parameters&) const = default;
			bool HasZeroParam() const { return epsilon == 0.f; };

			
			/// <param name="sigma">[nm]</param>
			/// <param name="epsilon">[kJ/mol/nm]</param>
			static Parameters CreateFromCharmm(float sigma, float epsilon);
		};

		PairBond() {};
		PairBond(std::array<uint8_t, 2> ids, const Parameters&);

		Parameters params;
		uint8_t atom_indexes[2] = { 0,0 };	// Relative to the compund
		const static int nAtoms = 2;
	};

	struct AngleUreyBradleyBond {
		struct Parameters {
			float theta0 = 0.f;	// [rad]
			float kTheta = 0.f;	// [J/mol/rad^2]
			float ub0 = 0.f;	// [nm]
			float kUB = 0.f;	// [J/mol/nm^2]

			bool operator==(const Parameters&) const = default;
			bool HasZeroParam() const { return kTheta == 0.f; }
			
			
			/// <param name="t0">[degrees]</param>
			/// <param name="kTheta">[kJ/mol/rad^2]</param>
			/// <param name="ub0">[nm]</param>
			/// <param name="kUb">[kJ/molnm]</param>
			static Parameters CreateFromCharmm(float t0, float kTheta, float ub0, float kUb);
		};

		AngleUreyBradleyBond() {}
		AngleUreyBradleyBond(std::array<uint8_t, 3> ids, const Parameters&);

		Parameters params;
		uint8_t atom_indexes[3] = { 0,0,0 }; // i,j,k angle between i and k
		const static int nAtoms = 3;
	};

	struct DihedralBond {
		struct Parameters {
			float phi_0;		// [rad]
			float k_phi;		// [J/mol/rad^2]
			float n;			// [multiplicity] n parameter, how many energy equilibriums does the dihedral have // OPTIMIZE: maybe float makes more sense, to avoid conversion in kernels?

			bool operator==(const Parameters&) const = default;
			bool HasZeroParam() const { return k_phi == 0.f; }

			/// <param name="phi_0">[degress]</param>
			/// <param name="k_phi">[kJ/mol/rad^2]</param>
			static Parameters CreateFromCharmm(float phi0, float kPhi, int n);
		};
		const static int nAtoms = 4;
		DihedralBond() {}
		DihedralBond(std::array<uint8_t, 4> ids, const Parameters&);

		Parameters params;
		uint8_t atom_indexes[4] = { 0,0,0,0 };
	};

	struct ImproperDihedralBond {
		struct Parameters {
			float psi_0 = 0.f;	// [rad]
			float k_psi = 0.f;	// [J/mol/rad^2]

			bool operator==(const Parameters&) const = default;
			bool HasZeroParam() const { return k_psi == 0.f; }

			/// <param name="psi_0">[degrees]</param>
			/// <param name="k_psi">[kJ/mol/rad^2]</param>
			static Parameters CreateFromCharmm(float psi0, float kPsi);
		};

		ImproperDihedralBond() {}
		ImproperDihedralBond(std::array<uint8_t, 4> ids, const Parameters&);

		Parameters params;

		uint8_t atom_indexes[4] = { 0,0,0,0 };
		const static int nAtoms = 4;
	};
}
using namespace Bondtypes;
// ------------------------------------------------- COMPOUNDS ------------------------------------------------- //







struct CompoundCoords {
	__device__ void loadData(const CompoundCoords& coords) {
		if (threadIdx.x == 0) { origo = coords.origo; };
		rel_positions[threadIdx.x] = coords.rel_positions[threadIdx.x];
	}
	
	NodeIndex origo{};								// [nm]
	Coord rel_positions[MAX_COMPOUND_PARTICLES];	// [nm]
};




// struct with data that only the solvent itself needs
struct TinyMolState {
	Float3 vel_prev{};
	Float3 force_prev{};
	int tinymolTypeIndex = -1; // wrong place to have this
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


	// Use this to quickly lookup wheter a bondedparticleslut exists with another compound
	static const int max_bonded_compounds = MAX_COMPOUNDS_IN_BRIDGE * 2 - 2;
	int n_bonded_compounds = 0;

	__device__ void loadMeta(const CompoundCompact* const compound) {
		n_particles = compound->n_particles;
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


struct CompoundInterimState {
	// Interims from the bridgekernel to compoundkernel
	//float potE_interim[MAX_COMPOUND_PARTICLES];
	//Float3 forces_interim[MAX_COMPOUND_PARTICLES];	// [GN/mol]
	//__host__ Float3 sumForce(int particleIndex) const;
	//__host__ float sumPotentialenergy(int particleIndex) const;

	//ForceEnergy forceEnergyFarneighborShortrange[MAX_COMPOUND_PARTICLES];
	/*ForceEnergy forceEnergyImmediateneighborShortrange[MAX_COMPOUND_PARTICLES];
	ForceEnergy forceEnergyBonds[MAX_COMPOUND_PARTICLES];
	ForceEnergy forceEnergyBridge[MAX_COMPOUND_PARTICLES];*/

	// Used specifically for Velocity Verlet stormer, and ofcourse kinE fetching
	Float3 forces_prev[MAX_COMPOUND_PARTICLES];
	Float3 vels_prev[MAX_COMPOUND_PARTICLES]; // Get wierd change of outcome if i move this here??

	Coord coords[MAX_COMPOUND_PARTICLES];
};



struct BondgroupRef { // A particles ref to its position in a bondgroup
	int bondgroupId;
	int localIndexInBondgroup;

	bool operator<(const BondgroupRef& other) const {
		if (bondgroupId != other.bondgroupId)
			return bondgroupId < other.bondgroupId;
		return localIndexInBondgroup < other.localIndexInBondgroup;
	}
};

// Rather large unique structures in global memory, that can be partly loaded when needed
struct Compound : public CompoundCompact {
	__host__ __device__ Compound() {}

	CompoundInteractionBoundary interaction_boundary;
	int centerparticle_index = -1;			// Index of particle initially closest to CoM

	uint16_t bonded_compound_ids[max_bonded_compounds];	// *2-2because it should exclude itself from both sides
    float atom_charges[MAX_COMPOUND_PARTICLES];	// [C/mol] - prolly move next to atomtypes to improve locality

	// For drawing pretty spheres :)
	char atomLetters[MAX_COMPOUND_PARTICLES];

	float atomMasses[MAX_COMPOUND_PARTICLES];	// [kg/mol]

	int absoluteIndexOfFirstParticle = 0;

	struct BondgroupRefManager {
		static const int maxBondgroupApperances = 4;
		int nBondgroupApperances = 0;
		BondgroupRef bondgroupApperances[maxBondgroupApperances];
	} bondgroupReferences[MAX_COMPOUND_PARTICLES];
};

struct BondGroup {
	struct ParticleRef {
		int compoundId=0; // TODO: make uint16_t?
		int localIdInCompound=0; // TODO: make uint16_t?
	};

	static const int maxParticles = 64;
	static const int maxSinglebonds = 128;
	static const int maxAnglebonds = 128 + 64;
	static const int maxDihedralbonds = 256 + 64;
	static const int maxImproperdihedralbonds = 32;

	ParticleRef particles[maxParticles];
	SingleBond singlebonds[maxSinglebonds];
	AngleUreyBradleyBond anglebonds[maxAnglebonds];
	DihedralBond dihedralbonds[maxDihedralbonds];
	ImproperDihedralBond improperdihedralbonds[maxImproperdihedralbonds];

	int nParticles = 0;
	int nSinglebonds = 0;
	int nAnglebonds = 0;
	int nDihedralbonds = 0;
	int nImproperdihedralbonds = 0;
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

// Precomputed values for pairs of atomtypes
struct NonbondedInteractionParams {
	float sigma;
    float epsilon;
    //float chargeProduct;
};


struct ForceField_NB {
	static const int MAX_TYPES = 64;

	struct ParticleParameters {	//Nonbonded
		//float mass = -1;		//[kg/mol]	or 
		float sigma = -1;		// [nm]
		float epsilon = -1;		// [J/mol/nm]
	};

	ParticleParameters particle_parameters[MAX_TYPES];
};

struct ForcefieldTinymol {
	static const int MAX_TYPES = 16;

	// Can make mass and epsilon half
	struct TinyMolType {
		float sigma = -1;		// [nm]
		float epsilon = -1;		// [J/mol/nm]
		float mass = -1;		// [kg/mol]
		float charge = -1;		// [kC/mol]
	};

	TinyMolType types[MAX_TYPES];
};


class UniformElectricField {
	Float3 field;	// [mV/nm]

	public:
		UniformElectricField() {}
		/// <summary></summary>
		/// <param name="direction"></param>
		/// <param name="magnitude">[V/nm]</param>
		__host__ UniformElectricField(Float3 direction, float magnitude) 
			: field(direction.norm() * magnitude * KILO) {
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
