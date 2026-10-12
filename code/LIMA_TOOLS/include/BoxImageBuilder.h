#pragma once

#include <string.h>
#include <fstream>
#include <vector>
#include <array>
#include <memory>

#include "Bodies.cuh"
#include "Constants.h"
#include "Utilities.h"
#include "BoundaryConditionPublic.h"
#include "EngineCore.h"
#include "MoleculeGraph.h"
#include "MDFiles.h"

#include "Forcefield.h"
#include <future>
#include "tuple"
#include <set>
struct BoxImage;





// --------------------------------- Bond Factories --------------------------------- //


struct ParticleFactory {
	ParticleFactory(const TopologyFile::AtomsEntry& topologyAtom, const Float3& pos, int indexInGrofile, int activeLJParamIndex) :
		topologyAtom(topologyAtom), position(pos), indexInGrofile(indexInGrofile), activeLJParamIndex(activeLJParamIndex) {}

	const TopologyFile::AtomsEntry& topologyAtom;
	const Float3 position{};
	int indexInGrofile = -1; // 0-indexed

	const int activeLJParamIndex = -1;
};

template <int n_Atoms, typename ParamsType>
struct BondFactory {
	static const int nAtoms = n_Atoms;
	BondFactory() = default;
	BondFactory(const std::array<int, nAtoms>& ids, const ParamsType& parameters)
		: params(parameters), global_atom_indexes(ids) {}

	ParamsType params{};
	std::array<int, nAtoms> global_atom_indexes{};
};
using SingleBondFactory = BondFactory<2, SingleBond::Parameters>;
using PairBondFactory = BondFactory<2, PairBond::Parameters>;
using AngleBondFactory = BondFactory<3, AngleUreyBradleyBond::Parameters>;
using DihedralBondFactory = BondFactory<4, DihedralBond::Parameters>;
using ImproperDihedralBondFactory = BondFactory<4, ImproperDihedralBond::Parameters>;

struct ParticleToCompoundMapping {
	int compoundId = -1;
	int localIdInCompound = -1;
};
struct ParticleToBridgeMapping {
	int bridgeId = -1;
	int localIdInBridge = -1;
};

struct ParticleToPclusterMapping {
	int pcid;
	int pid; // pc local
};
using ParticleToPclusterMap = std::vector<ParticleToPclusterMapping>;

struct PersistentClusterFactory {
	std::vector<PersistentCluster> pClusters;
	std::vector<PersistentClusterMeta> pClusterMetas;
	ParticleToPclusterMap particleToPclusterMap;
	std::vector<ParticlesBondedToParticle> particleBondedToParticle;
	std::vector<PclustersBondedToPcluster> pclusterBondedToPcluster;
};

namespace LIMA_MOLECULEBUILD {
	class SuperTopology {

	public:
		struct MoleculeInstance {
			std::shared_ptr<const TopologyFile::Moleculetype> type;
			int particleOffset = 0;
			int nParticles = 0;

			// Index of the first bond of this molecule in the respective vectors. The bonds of a molecule are contiguous
			int firstSinglebond = 0;
			int firstPairbond = 0;
			int firstAnglebond = 0;
			int firstDihedralbond = 0;
			int firstImproperdihedralbond = 0;
		};

		SuperTopology(const TopologyFile::System& system, const GroFile& grofile, LIMAForcefield& forcefield);


		void VerifyBondsAreStable(const Float3& boxlen_nm, BoundaryConditionSelect bc_select, bool energyMinimizationMode) const;


		//Temporary, untill how i know how to deal with bonds in tinymols
		void RemoveBondsFromTinymol(const std::vector<ParticleToCompoundMapping>&);

		std::vector<ParticleFactory> particles;
		std::vector<SingleBondFactory> singlebonds;
		std::vector<PairBondFactory> pairbonds;
		std::vector<AngleBondFactory> anglebonds;
		std::vector<DihedralBondFactory> dihedralbonds;
		std::vector<ImproperDihedralBondFactory> improperdihedralbonds;
		std::vector<MoleculeInstance> moleculeInstances;
	};





	std::unique_ptr<BoxImage> buildMolecules(
		const GroFile& gro_file,
		const TopologyFile& top_file,
		const SimParams& simparams
	);
}




// Groups bonds into bondgroups of at most 64 particles. A bondgroup never spans multiple molecules
class BondGroupFactory {
	static constexpr int maxParticlesPerBondgroup = 64;
	BondGroups bondgroups;
	std::vector<int> particleGlobalIds;	// The global id of each particle in bondgroups.particles

public:
	BondGroupFactory(const LIMA_MOLECULEBUILD::SuperTopology& topology);

	// Warning: unfinished bondgroups, run the function below before using
	void AddPclusterRefs(const ParticleToPclusterMap& particleToPclusterMap);

	// Adds a reference to each particle's pcluster, for each bondgroup the particle is in, in order of bondgroup
	void AddBondgroupRefsToPclusters(const ParticleToPclusterMap& particleToPclusterMap, std::vector<PersistentClusterMeta>& pclusterMetas) const;
	BondGroups GetBondgroups();
};



// A translation unit between Gro file representation, and LIMA Box representation
struct BoxImage {	

	GroFile grofile;

	LIMA_MOLECULEBUILD::SuperTopology topology; // This is only used for debugging purposes

	BondGroups bondgroups;

	// Clusters
	std::vector<PersistentCluster> persistentClusters;
	std::vector<PersistentClusterMeta> persistentClustersMetadata;
	std::vector<ParticlesBondedToParticle> particleBondedToParticle;
	std::vector<PclustersBondedToPcluster> pclusterBondedToPcluster;
	std::vector<std::tuple<int, int>> gpidToPcidAndPid;
	//std::vector<ParticleToCompoundOrSolventMapping> particleToCompoundOrSolventMapping;
	int totalParticles = 0;
};
 
