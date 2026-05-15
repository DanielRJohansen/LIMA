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
	BondFactory(const std::array<int, nAtoms>& ids, const ParamsType& parameters)
		: params(parameters), global_atom_indexes(ids) {}

	ParamsType params;
	std::array<int, nAtoms> global_atom_indexes;
	std::string sourceLine;
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
using ParticleToCompoundMap = std::vector<ParticleToCompoundMapping>;
using ParticleToBridgeMap = std::vector<std::optional<ParticleToBridgeMapping>>;

struct ParticleToPclusterMapping {
	int pcid;
	int pid; // pc local
};
using ParticleToPclusterMap = std::vector<ParticleToPclusterMapping>;

struct PersistentClusterFactory {
	std::vector<PersistentCluster> pClusters;
	std::vector<PersistentClusterMeta> pClusterMetas;
	ParticleToPclusterMap particleToPclusterMap;
	std::vector<std::set<int>>particleBondedToParticle;
	std::vector<std::set<int>> pclusterBondedToPcluster;
};

namespace LIMA_MOLECULEBUILD {
	class SuperTopology {

		template <typename BondType, typename BondtypeFactory, typename BondTypeTopologyfile>
		void LoadBondsIntoTopology(const std::vector<BondTypeTopologyfile>& bondsInTopfile, 
			int atomIdOffset, LIMAForcefield& forcefield, std::vector<BondtypeFactory>& topology);

	public:
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
	};





	std::unique_ptr<BoxImage> buildMolecules(
		const GroFile& gro_file,
		const TopologyFile& top_file,
		VerbosityLevel vl,
		std::unique_ptr<LimaLogger>,
		bool ignore_hydrogens,
		const SimParams& simparams
	);
}




class BondGroupFactory {

	std::vector<BondGroup> bondgroups;
	std::vector<std::array<int, BondGroup::maxParticles>> particleGlobalIds;

	int FindLocalParticleId(int bgIndex, const int globalId) const;
	void AddBondParticles(int bgIndex, std::span<const int> globalIds, std::span<const uint8_t> localIds);

	template <int n>
	std::array<uint8_t, n> GetLocalIds(const std::array<int, n>& globalIds) const;

	
	bool AddBond(int bondgroupIndex, const SingleBondFactory&);
	bool AddBond(int bondgroupIndex, const PairBondFactory&);
	bool AddBond(int bondgroupIndex, const AngleBondFactory&);
	bool AddBond(int bondgroupIndex, const DihedralBondFactory&);
	bool AddBond(int bondgroupIndex, const ImproperDihedralBondFactory&);

	// Add bonds from a specific type to the bond group
	void AddBondsFromMap(int bgIndex, const auto& bondMap, auto& availableBondIds, const auto& bonds) {
		//TimeIt timer("addbondsfrommap");
		for (const int bondId : bondMap) {
			if (availableBondIds.contains(bondId)) {
				if (AddBond(bgIndex, bonds[bondId]))
					availableBondIds.erase(bondId);
			}
		}
	};

public:
	BondGroupFactory(const LIMA_MOLECULEBUILD::SuperTopology& topology);
	
	// Returns <nNewParticles, localParticleIds>, where localParticleIds may not be assigned yet..
	template <int n>
	std::tuple<int, std::array<uint8_t, n>> TryAssignLocalIds(int bgIndex, const std::array<int, n>& particleIds) const;



	// Warning: unfinished bondgroups, run the function below before using
	
	void AddPclusterRefs(const ParticleToPclusterMap& particleToPclusterMap);
	std::vector<std::set<BondgroupRef>> MakeParticleToBondgroupsMap(int nParticlesTotal) const;
	std::vector<BondGroup> GetBondgroups();

	//static std::vector<BondGroup> FinishBondgroups(const std::vector<BondGroupFactory>&);
};



// A translation unit between Gro file representation, and LIMA Box representation
struct BoxImage {	

	GroFile grofile;

	const ForceField_NB forcefield;


	LIMA_MOLECULEBUILD::SuperTopology topology; // This is only used for debugging purposes

	std::shared_ptr<LimaMoleculeGraph::MoleculeGraph> systemGraph;

	const std::vector<BondGroup> bondgroups;

	// Clusters
	std::vector<PersistentCluster> persistentClusters;
	std::vector<PersistentClusterMeta> persistentClustersMetadata;
	std::vector<std::set<int>> particleBondedToParticle;
	std::vector<std::set<int>> pclusterBondedToPcluster;
	std::vector<std::tuple<int, int>> gpidToPcidAndPid;
	//std::vector<ParticleToCompoundOrSolventMapping> particleToCompoundOrSolventMapping;
	int totalParticles = 0;
};
 
