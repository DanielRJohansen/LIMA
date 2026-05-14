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




class BondGroupFactory : public BondGroup {

	int FindLocalParticleId(const int globalId) const;
	void AddBondParticles(const ParticleToPclusterMap&, std::span<const int> globalIds);
	template <int n>
	std::array<uint8_t, n> GetLocalIds(const std::array<int, n>& globalIds) const;
public:
	BondGroupFactory() {}

	bool HasSpaceForParticlesInBond(const std::span<const int>& particleIds) const;

	//void AddParticles(const std::span<const uint32_t>& particleIds);

	void AddBond(const ParticleToPclusterMap&, const SingleBondFactory&, const PersistentCluster* pClusters = nullptr);
	void AddBond(const ParticleToPclusterMap&, const PairBondFactory&, const PersistentCluster* pClusters = nullptr);
	void AddBond(const ParticleToPclusterMap&, const AngleBondFactory&, const PersistentCluster* pClusters = nullptr);
	void AddBond(const ParticleToPclusterMap&, const DihedralBondFactory&, const PersistentCluster* pClusters = nullptr);
	void AddBond(const ParticleToPclusterMap&, const ImproperDihedralBondFactory&, const PersistentCluster* pClusters = nullptr);
	
	std::array<int, maxParticles> particleGlobalIds;
	//std::unordered_map<int, uint8_t> particleGlobalToLocalId;



	static std::vector<BondGroupFactory> MakeBondgroups(const LIMA_MOLECULEBUILD::SuperTopology&,
		const ParticleToPclusterMap&, const PersistentCluster* pClusters);

	static std::vector<std::set<BondgroupRef>> MakeParticleToBondgroupsMap(
		const std::vector<BondGroupFactory>&, int nParticlesTotal);

	static std::vector<BondGroup> FinishBondgroups(const std::vector<BondGroupFactory>&);
};



// A translation unit between Gro file representation, and LIMA Box representation
struct BoxImage {	

	GroFile grofile;

	const ForceField_NB forcefield;


	LIMA_MOLECULEBUILD::SuperTopology topology; // This is only used for debugging purposes

	std::shared_ptr<LimaMoleculeGraph::MoleculeGraph> systemGraph;

	const std::vector<NonbondedInteractionParams> nonbondedInteractionParams;

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
 
