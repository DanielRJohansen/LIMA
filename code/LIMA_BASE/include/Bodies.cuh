#pragma once

#include "Constants.h"
#include "LimaTypes.cuh"
#include <set>
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
		uint8_t idInBondgroup[2] = { 0,0 };	// Relative to the bondgroup
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











// ------------------------------------------------- Etc ------------------------------------------------- //



struct BondgroupRef { // A particles ref to its position in a bondgroup
	int bondgroupId;
	int localIndexInBondgroup;

	bool operator<(const BondgroupRef& other) const {
		if (bondgroupId != other.bondgroupId)
			return bondgroupId < other.bondgroupId;
		return localIndexInBondgroup < other.localIndexInBondgroup;
	}
};

struct BondGroup {
	struct ParticleRef {
		int pcid;
		int pid; // local to pcluster
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

	static_assert(maxParticles < UINT8_MAX, "bonds can't index their particles!");
};

struct NBParams {
	float sigmaHalf = -1;		// [nm]
	float epsilonSqrt = -1;		// [J/mol/nm]
	float charge = 0;		// [kC/mol]

	__host__ bool operator==(const NBParams& other) const = default;
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




// ------------------------------------------------- CLUSTERS ------------------------------------------------- //


struct PData {
	Float3 position;
	NBParams params;
	constexpr bool Valid() const { return params.epsilonSqrt != -1.f; }

	__host__ bool operator!=(const PData& other) const {
		return position != other.position || params != other.params;
	}
};

struct BondgroupRefManager {
	static const int maxBondgroupApperances = 4;
	int nBondgroupApperances = 0;
	BondgroupRef bondgroupApperances[maxBondgroupApperances];
	__host__ void Add(const BondgroupRef& bgRef) {
		if (nBondgroupApperances >= maxBondgroupApperances)
			throw std::runtime_error("Too many bondgroup apperances for a particle, increase maxBondgroupApperances or check your clustering");
		bondgroupApperances[nBondgroupApperances++] = bgRef;
	}
};

struct PersistentCluster {
	static const int maxParticles = 4;
	PData pqd[maxParticles];

	__host__ bool operator!=(const PersistentCluster& other) const {
		for (int i = 0; i < maxParticles; i++) {
			if (pqd[i] != other.pqd[i])
				return true;
		}
		return false;
	}
};
struct PersistentClusterMeta {
	int particleIdsGlobal[PersistentCluster::maxParticles]={ -1, -1, -1, -1 };
	float mass[PersistentCluster::maxParticles] = { 0,0,0,0 };		// [kg/mol]

	char atomLetter[PersistentCluster::maxParticles]; // For rendering
	bool isSolvent = false;
	int nParticles = 0;

	// I do not like this setup...
	BondgroupRefManager bondgroupReferences[PersistentCluster::maxParticles];
};

struct PersistentclusterInterimState {
	// Used specifically for Velocity Verlet stormer, and ofcourse kinE fetching
	Float3 forces_prev[PersistentCluster::maxParticles]; // [J/mol]
	Float3 vels_prev[PersistentCluster::maxParticles];
	//Coord coords[PersistentCluster::nParticles];
};

template <int size>
class StaticSet {	
	int data[size]; // is sorted
	static const int noVal = INT_MIN;
public:
	constexpr bool Contains(int value) const {
		for (int i = 0; i < size; i++) {
			if (data[i] == value)
				return true;
			if (data[i] > value || data[i] == noVal)
				return false;
		}
		return false;
	}

	static std::vector<StaticSet> Create(const std::vector<std::set<int>>& sets) {
		std::vector<StaticSet> result(sets.size());
		for (int i= 0; i < sets.size(); i++) {
			int j = 0;
			for (int val : sets[i]) {
				if (j >= size)
					throw std::runtime_error("Too many values in set, increase size or check your clustering");
				result[i].data[j++] = val;
			}
			for (; j < size; j++)
				result[i].data[j] = noVal;
		}
		return result;
	}
};

using ParticlesBondedToParticle = StaticSet<32>;
using PclustersBondedToPcluster = StaticSet<32>;

//struct PersistentCluster {
//	ParticleQuickData pqd[4];
//};




class BoolMatrix16x16 {
	uint16_t data[16]; // rowmajor

public:
	constexpr BoolMatrix16x16() {}
	constexpr void Clear() {
		for (int i = 0; i < 16; i++)
			data[i] = 0;
	}

	constexpr static bool Get(const uint16_t& row, int col) {
		return (row >> col) & 1;
	}
	template <typename T> constexpr static bool Get(const T& row, int col) = delete;

	constexpr uint16_t GetRow(int row) const {
		return data[row];
	}

	constexpr uint16_t SetRow(int row, uint16_t val) {
		return data[row] = val;
	}

	constexpr uint16_t GetColumn(int col) const {
		uint16_t out = 0;
		for (int row = 0; row < 16; ++row)
			out |= static_cast<uint16_t>(((data[row] >> col) & 1u) << row);
		return out;
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

	constexpr static void SetValueInRow(int col, uint16_t& rowData) {
		rowData |= (1 << col);
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


struct SuperCluster {
	//static const int maxPclusters = 4;
	static const int maxParticles = 16;

	//Float3 positions[nParticles];
	PData pData[maxParticles];

	__host__ bool operator!= (const SuperCluster& other) const {
		for (int i = 0; i < maxParticles; i++) {
			if (pData[i] != other.pData[i])
				return true;
		}
		return false;
	}
};

struct SuperClusterMeta {
	// Set by clustering kernel
	int _pclusterIds[SuperCluster::maxParticles];
	int indexInPcluster[SuperCluster::maxParticles];
	int globalParticleIds[SuperCluster::maxParticles];	
	int uniquePclusterIds[SuperCluster::maxParticles];
	int nUniquePcIds = 0;
	int nParticles;

	// Set by taskbuilder kernel
	int resultsStartIndex; // TODO: Is int always safe here??
	int nResults;
	

	//__host__ bool operator != (const SuperClusterMeta& other) const {
	//	if (nUniquePcIds != other.nUniquePcIds || nParticles != other.nParticles || )
	//		return true;


	//	if (resultsStartIndex != other.resultsStartIndex ||
	//		nResults != other.nResults)
	//		return true;
	//	for (int i = 0; i < SuperCluster::nPclusters; i++) {
	//		if (pclusterIds[i] != other.pclusterIds[i])
	//			return true;
	//	}
	//	return false;
	//}
};

struct SCResult {
	ForceEnergy fe[SuperCluster::maxParticles];

	__host__ bool operator!=(const SCResult& other) const {
		for (int i = 0; i < SuperCluster::maxParticles; i++) {
			if (fe[i].force != other.fe[i].force ||
				fe[i].potE != other.fe[i].potE)
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
	constexpr Float3 GetForce(float charge) const {
		return field * charge;
	}
};


#pragma warning (pop)
