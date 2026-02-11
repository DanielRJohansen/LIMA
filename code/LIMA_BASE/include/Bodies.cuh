#pragma once

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

		constexpr SingleBond() {}
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
			static Parameters CreateFromCharmm(float t0, float kTheta, float ub0, float kUb, int func);
		};

		constexpr AngleUreyBradleyBond() {}
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
struct TinyMolParticleState {
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
	constexpr CompoundCompact() {}

	alignas(4) uint8_t atom_types[MAX_COMPOUND_PARTICLES];
	int n_particles = 0;

#if LIMAKERNELDEBUGMODE == 1
	uint32_t particle_global_ids[MAX_COMPOUND_PARTICLES];
#endif


	// Use this to quickly lookup wheter a bondedparticleslut exists with another compound
	static const int max_bonded_compounds = 4 * 2 - 2;
	int n_bonded_compounds = 0;

	__device__ void loadMeta(const CompoundCompact* const compound) {
		n_particles = compound->n_particles;
		n_bonded_compounds = compound->n_bonded_compounds;
	}

	__device__ void loadData(const CompoundCompact* const compound) {
		if (threadIdx.x < n_particles) {
			atom_types[threadIdx.x] = compound->atom_types[threadIdx.x];

			#if LIMAKERNELDEBUGMODE == 1
			particle_global_ids[threadIdx.x] = compound->particle_global_ids[threadIdx.x];
			#endif
		}
	}
};


struct CompoundInterimState {
	// Used specifically for Velocity Verlet stormer, and ofcourse kinE fetching
	Float3 forces_prev[MAX_COMPOUND_PARTICLES]; // [J/mol]
	Float3 vels_prev[MAX_COMPOUND_PARTICLES];

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
	static const int maxPairbonds = maxDihedralbonds;
	static const int maxImproperdihedralbonds = 32;

	ParticleRef particles[maxParticles];
	SingleBond singlebonds[maxSinglebonds];
	PairBond pairbonds[maxPairbonds];
	AngleUreyBradleyBond anglebonds[maxAnglebonds];
	DihedralBond dihedralbonds[maxDihedralbonds];
	ImproperDihedralBond improperdihedralbonds[maxImproperdihedralbonds];

	int nParticles = 0;
	int nSinglebonds = 0;
	int nPairbonds = 0;
	int nAnglebonds = 0;
	int nDihedralbonds = 0;
	int nImproperdihedralbonds = 0;
};

// TODO: OPTIM: THese should actually be cached in constant memory and accessed with a single id,
// because most tinymols are identical, just with different positions
struct BondgroupTinymol {
	
	static const int maxParticles = 4;
	static const int maxSinglebonds = 4;
	static const int maxAnglebonds = 4;

	
	//uint8_t particleIndicesRelativeToTinymol[maxParticles];
	// All indices are relative to the tinymol, so add the tinymols indexOfFirstInSolventlblock when accessing particle pos
	SingleBond singlebonds[maxSinglebonds];
	AngleUreyBradleyBond anglebonds[maxAnglebonds];
	int nParticles = 0;
	int nSinglebonds = 0;
	int nAnglebonds = 0;
};



struct ParticleReference {
	// Used by moleculebuilder only
	constexpr ParticleReference(int compound_id, int local_id_compound, uint8_t compoundid_local_to_bridge) :
		compound_id(compound_id), local_id_compound(local_id_compound),
		compoundid_local_to_bridge(compoundid_local_to_bridge) 
	{}

	int compound_id;	// global
	int local_id_compound;	// id of particle
	uint8_t compoundid_local_to_bridge = 255;

	//int global_id = -1; // For debug
};

struct NBParams {
	float sigmaHalf = -1;		// [nm]
	float epsilonSqrt = -1;		// [J/mol/nm]
	float charge = NAN;		// [kC/mol]
};

// Precomputed values for pairs of atomtypes
struct NonbondedInteractionParams {
	float sigma;
    float epsilon;
    float chargeProduct;
};


struct ForceField_NB {
	static const int MAX_TYPES = 64;

	struct ParticleParameters {	//Nonbonded
		//float mass = -1;		//[kg/mol]	or 
		// Values at a format for efficient parameter computation
		float sigmaHalf = -1;		// [nm]
		float epsilonSqrt = -1;		// [J/mol/nm]
	};

	ParticleParameters particle_parameters[MAX_TYPES];
};

struct ForcefieldTinymol {
    static const int MAX_TYPES = 16; // TODO OPTIM change to 4

	// Can make mass and epsilon half
	struct TinyMolType {
		float sigmaHalf = -1;		// [nm]
		float epsilonSqrt = -1;		// [J/mol/nm] // TODO: OPTIM: Should be 0 so the same logic handles missing data aswell as particles that doesnt interact with LJ
		float mass = -1;		// [kg/mol]
		float charge = -1;		// [kC/mol]
	};

	TinyMolType types[MAX_TYPES];
};

struct PData {
	Float3 position;
	NBParams params;
	constexpr bool Valid() const { return params.epsilonSqrt != -1.f; }
};

//struct PrecomputedSolventForcefield {
//	NonbondedInteractionParams ljParams[3]; // [O-O, O-H, H-H]
//	float chargeProducts[3]; // [O-O, O-H, H-H]
//};

struct PersistentCluster {
	static const int nParticles = 4;
	PData pqd[nParticles];
};
struct PersistentClusterMeta {
	int particleIdsGlobal[PersistentCluster::nParticles];
};

//struct PersistentCluster {
//	ParticleQuickData pqd[4];
//};

struct SuperCluster {
	static const int nPclusters = 4;
	static const int nParticles = PersistentCluster::nParticles * nPclusters;


	//Float3 positions[nParticles];
	PData pData[nParticles];

#if LIMAKERNELDEBUGMODE == 1
	Float3 center;
#endif

	/*float x[nParticles];
	float y[nParticles];
	float z[nParticles];*/

};

struct SCResult {
	ForceEnergy fe[SuperCluster::nParticles];

	__host__ bool operator!=(const SCResult& other) const {
		for (int i = 0; i < SuperCluster::nParticles; i++) {
			if (fe[i].force != other.fe[i].force ||
				fe[i].potE != other.fe[i].potE)
				return true;
		}
		return false;
	}
};


class BoolMatrix16x16 {
	uint16_t data[16]; // rowmajor

public:
	BoolMatrix16x16(){
		memset(data, 0, sizeof(data));
	}
	constexpr static bool Get(const uint16_t& row, int col) {
		return (row >> col) & 1;
	}
	template <typename T> constexpr static bool Get(const T& row, int col) = delete;

	constexpr uint16_t GetRow(int row) const {
		return data[row];
	}
	// TODO: Optim this with a SetRow
	constexpr void Set(int row, int col, bool val) {
		unsigned bit = 1u << col;
		unsigned mask = -static_cast<unsigned>(val);  // 0xFFFFFFFF if val==1, else 0
		unsigned old = data[row];

		data[row] = (old & ~bit) | (mask & bit);
		/*if (val)
			data[row] |= (1 << col);
		else
			data[row] &= ~(1 << col);		*/
	}

	__host__ void Print() const {
		for (int r = 0; r < 16; r++) {
			for (int c = 0; c < 16; c++) {
				printf("%d ", Get(data[r], c) ? 1 : 0);
			}
			printf("\n");
		}
		printf("\n");
	}
};
class NoMat {};// Needed as a nonlocal variant of the one above.


struct SuperClusterMeta {
	// Set by clustering kernel
	int pclusterIds[SuperCluster::nPclusters];
	
	//Float3 meanPos;


	// For debugging, find a way to remove in release automatically
	//int particlesIds[SuperCluster::nParticles];
	std::array<int, SuperCluster::nParticles> particlesIds;

	// Set by taskbuilder kernel
	int resultsStartIndex;
	int nResults;

	__host__ bool operator != (const SuperClusterMeta& other) const {
		if (resultsStartIndex != other.resultsStartIndex ||
			nResults != other.nResults)
			return true;
		for (int i = 0; i < SuperCluster::nPclusters; i++) {
			if (pclusterIds[i] != other.pclusterIds[i])
				return true;
		}
		return false;
	}
};

struct ScScTask {
	int scIds[2];
	int resultIndices[2];
	int nointeractionMatrixIndex = -1;

	__host__ constexpr bool operator!=(const ScScTask& other) const {
		for (int i = 0; i < 2; i++) {
			if (scIds[i] != other.scIds[i])
				return true;
			if (resultIndices[i] != other.resultIndices[i])
				return true;
		}
		return false;
	}
};


struct ParticleToCompoundOrSolventMapping {
	int compoundId = -1;
	int particleId = -1; // Relative to compound if in compound, otherwise solventid
	ParticleToCompoundOrSolventMapping() {}
	ParticleToCompoundOrSolventMapping(int solventId) {
		particleId = solventId;
	}
	ParticleToCompoundOrSolventMapping(int cid, int pid) {
		compoundId = cid;
		particleId = pid;
	}

	constexpr bool IsSolvent() const { return compoundId == -1; }
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
