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



struct BondgroupRef { // A particle's reference to its bonded-force result
	int bondgroupId;
	int indexInForceEnergiesBondgroups;

	bool operator<(const BondgroupRef& other) const {
		if (bondgroupId != other.bondgroupId)
			return bondgroupId < other.bondgroupId;
		return indexInForceEnergiesBondgroups < other.indexInForceEnergiesBondgroups;
	}
};

struct BondGroup {
	struct ParticleRef {
		int pcid;
		int pid; // local to pcluster
	};
	int indexOfFirstParticle = 0;
	int nParticles = 0;
	int indexOfFirstSinglebond = 0;
	int nSinglebonds = 0;
	int indexOfFirstPairbond = 0;
	int nPairbonds = 0;
	int indexOfFirstAnglebond = 0;
	int nAnglebonds = 0;
	int indexOfFirstDihedralbond = 0;
	int nDihedralbonds = 0;
	int indexOfFirstImproperdihedralbond = 0;
	int nImproperdihedralbonds = 0;
};

// Host-side SoA storage. BondGroup remains a compact device-friendly range descriptor.
struct BondGroups {
	std::vector<BondGroup> groups;
	std::vector<BondGroup::ParticleRef> particles;
	std::vector<SingleBond> singlebonds;
	std::vector<PairBond> pairbonds;
	std::vector<AngleUreyBradleyBond> anglebonds;
	std::vector<DihedralBond> dihedralbonds;
	std::vector<ImproperDihedralBond> improperdihedralbonds;

	size_t size() const { return groups.size(); }
	bool empty() const { return groups.empty(); }
	void clear() {
		groups.clear(); particles.clear(); singlebonds.clear(); pairbonds.clear(); anglebonds.clear(); dihedralbonds.clear(); improperdihedralbonds.clear();
	}
};

struct BondGroupsDevice {
	const BondGroup* groups;
	const BondGroup::ParticleRef* particles;
	const SingleBond* singlebonds;
	const PairBond* pairbonds;
	const AngleUreyBradleyBond* anglebonds;
	const DihedralBond* dihedralbonds;
	const ImproperDihedralBond* improperdihedralbonds;
};

struct NBParams {
	float sigmaHalf = -1;		// [nm]
	float epsilonSqrt = -1;		// [J/mol/nm]
	float charge = 0;		// [kC/mol]

	//__host__ bool operator==(const NBParams& other) const = default;
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
	NBParams params{};
	constexpr bool Valid() const { return params.epsilonSqrt != -1.f; }

	/*__host__ bool operator!=(const PData& other) const {
		return position != other.position || params != other.params;
	}*/
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

	//__host__ bool operator!=(const PersistentCluster& other) const {
	//	for (int i = 0; i < maxParticles; i++) {
	//		if (pqd[i] != other.pqd[i])
	//			return true;
	//	}
	//	return false;
	//}
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
	constexpr StaticSet() {
		for (int& value : data)
			value = noVal;
	}

	void AddOffset(int offset) {
		for (int& value : data)
			if (value != noVal) value += offset;
	}

	// The index-th value, or a value IsValue rejects past the last one
	constexpr int Get(int index) const { return data[index]; }
	static constexpr bool IsValue(int value) { return value != noVal; }

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
			result[i] = Create(sets[i]);
		}
		return result;
	}

	static StaticSet Create(const std::set<int>& values, int offset = 0) {
		StaticSet result;
		int index = 0;
		for (const int value : values) {
			if (index >= size)
				throw std::runtime_error("Too many values in set, increase size or check your clustering");
			result.data[index++] = value + offset;
		}
		return result;
	}

	// values must be sorted and unique
	static StaticSet CreateFromSorted(const std::vector<int>& values) {
		if (values.size() > size)
			throw std::runtime_error("Too many values in set, increase size or check your clustering");
		StaticSet result;
		for (int i = 0; i < values.size(); i++)
			result.data[i] = values[i];
		return result;
	}
};

using ParticlesBondedToParticle = StaticSet<32>;
using PclustersBondedToPcluster = StaticSet<32>;

//struct PersistentCluster {
//	ParticleQuickData pqd[4];
//};






struct SuperCluster {
	//static const int maxPclusters = 4;
	static const int maxParticles = 16;

	// AoS per particle, so the nonbonded kernels load any 4 consecutive particles in 3 sectors
	float4 posCharge[maxParticles];			// x, y, z [nm], charge [kC/mol]
	float2 sigmaEpsilon[maxParticles];		// sigmaHalf [nm], epsilonSqrt [J/mol/nm]. epsilonSqrt is -1 for padding

	__host__ __device__ void SetPdata(const PData& pdata, int index) {
		posCharge[index] = float4{ pdata.position.x, pdata.position.y, pdata.position.z, pdata.params.charge };
		sigmaEpsilon[index] = float2{ pdata.params.sigmaHalf, pdata.params.epsilonSqrt };
	}
	__device__ void LoadPdata(PData& pdata, int index) const {
		const float4 pq = posCharge[index];
		const float2 se = sigmaEpsilon[index];
		pdata.position = Float3(pq.x, pq.y, pq.z);
		pdata.params.sigmaHalf = se.x;
		pdata.params.epsilonSqrt = se.y;
		pdata.params.charge = pq.w;
	}
	__host__ __device__ Float3 Position(int index) const {
		return Float3(posCharge[index].x, posCharge[index].y, posCharge[index].z);
	}
	__device__ void SetPosition(int index, const Float3& position) {
		posCharge[index].x = position.x;
		posCharge[index].y = position.y;
		posCharge[index].z = position.z;
	}
	__host__ __device__ float SigmaHalf(int index) const { return sigmaEpsilon[index].x; }
	__host__ __device__ float EpsilonSqrt(int index) const { return sigmaEpsilon[index].y; }
	__host__ __device__ float Charge(int index) const { return posCharge[index].w; }
	__host__ __device__ bool Valid(int index) const { return sigmaEpsilon[index].y != -1.f; }
	//PData pData[maxParticles];

	//__host__ bool operator!= (const SuperCluster& other) const {
	//	for (int i = 0; i < maxParticles; i++) {
	//		if (pData[i] != other.pData[i])
	//			return true;
	//	}
	//	return false;
	//}
};

struct alignas(16) SuperClusterMeta {	
	// Set by clustering kernel
	int _pclusterIds[SuperCluster::maxParticles];
	int indexInPcluster[SuperCluster::maxParticles];
	int globalParticleIds[SuperCluster::maxParticles];	
	int uniquePclusterIds[SuperCluster::maxParticles];
	int nUniquePcIds = 0;
	int16_t nParticles;
	int16_t simulationId = 0;

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

// SoA layout so each component load is fully coalesced, and so potE (only computed on logging steps) occupies
// its own sectors which are never touched on non-logging steps
struct SCResult {
	float fx[SuperCluster::maxParticles];
	float fy[SuperCluster::maxParticles];
	float fz[SuperCluster::maxParticles];
	float potE[SuperCluster::maxParticles]; // Only written/read when withPotE

	template <bool withPotE>
	__device__ void Store(int pid, const ForceEnergy& fe) {
		fx[pid] = fe.force.x;
		fy[pid] = fe.force.y;
		fz[pid] = fe.force.z;
		if constexpr (withPotE)
			potE[pid] = fe.potE;
	}

	template <bool withPotE>
	__device__ ForceEnergy Load(int pid) const {
		if constexpr (withPotE)
			return ForceEnergy{ Float3{ fx[pid], fy[pid], fz[pid] }, potE[pid] };
		else
			return ForceEnergy{ Float3{ fx[pid], fy[pid], fz[pid] }, 0.f };
	}

	__host__ bool operator!=(const SCResult& other) const {
		for (int i = 0; i < SuperCluster::maxParticles; i++) {
			if (fx[i] != other.fx[i] || fy[i] != other.fy[i] || fz[i] != other.fz[i] || potE[i] != other.potE[i])
				return true;
		}
		return false;
	}
};

// Nonbonded work, see NbNonlocalKernel. Superclusters are split in quarters of 4 particles, and pairs are only computed
// for the 4x4 blocks of quarters that had a pair within the list radius when the tasks were built.
// An entry is quarter jQuarter of supercluster jScId, with the blocks it forms with the quarters of the task's own supercluster
struct QuarterEntry {
	int jScId;
	uint16_t noInteractions[4];	// Per own quarter: bit iLocal*4+jLocal set if the pair must not be computed (bonded, padding, or computed in the other order)
	uint8_t jQuarter;
	uint8_t ownQuarterMask;		// Bit i set if own quarter i has a pair within the list radius
	uint16_t _unused = 0;
};
static_assert(sizeof(QuarterEntry) == 16);

// The entries of one supercluster
struct QuarterEntryTask {
	int firstEntry = 0;
	int nEntries = 0;
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
