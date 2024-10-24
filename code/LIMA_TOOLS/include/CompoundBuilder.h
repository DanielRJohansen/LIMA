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

#include "MDFiles.h"

struct BoxImage;





// --------------------------------- Bond Factories --------------------------------- //
template <int n_Atoms, typename ParamsType>
struct BondFactory {
	static const int nAtoms = n_Atoms;
	BondFactory(const std::array<uint32_t, nAtoms>& ids, const ParamsType& parameters) 
		: params(parameters), global_atom_indexes(ids) {}

	ParamsType params;
	std::array<uint32_t, nAtoms> global_atom_indexes;
	std::string sourceLine;
};

struct ParticleFactory {	
	const int activeLjtypeParameterIndex = -1;
	const TopologyFile::AtomsEntry* topAtom = nullptr;
	//const GroRecord* groAtom = nullptr;
	const Float3 position;
	const int uniqueResId;
};


using SingleBondFactory = BondFactory<2, SingleBond::Parameters>;
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
	void* aaa;
};

using ParticleToCompoundMap = std::vector<ParticleToCompoundMapping>;
using ParticleToBridgeMap = std::vector<std::optional<ParticleToBridgeMapping>>;

namespace LIMA_MOLECULEBUILD {
	struct Topology {
		std::vector<ParticleFactory> particles;

		std::vector<SingleBondFactory> singlebonds;
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


//struct ParticleInfo {
//
//
//	ParticleInfo() = default;
//	ParticleInfo(const GroRecord* groAtom, const TopologyFile::AtomsEntry* topAtom, int activeLjTypeParameterIndex, int uniqueResId)
//		: groAtom(groAtom), topAtom(topAtom), activeLjtypeParameterIndex(activeLjTypeParameterIndex), uniqueResId(uniqueResId), inUse(true) {}
//	const GroRecord* groAtom = nullptr;
//	const TopologyFile::AtomsEntry* topAtom = nullptr;
//	int activeLjtypeParameterIndex = -1;
//	int uniqueResId = -1;
//
//	// Only available once the compounds have been created
//	int compoundId = -1;
//	int localIdInCompound = -1;
//
//	// Only available once the bridges have been created
//	int bridgeId = -1;
//	int localIdInBridge = -1;
//
//	bool inUse = false;
//
//	int globalId = -1; // used for debugging
//	
//	std::string sourceLine;
//};

class CompoundFactory : public Compound, public CompoundInterimState {
public:
	CompoundFactory() {}
	CompoundFactory(const int id) : 
		id(id), local_atomtype_to_LATID_map(MAX_ATOM_TYPES, -1)
	{
		memset(forces_interim, 0, sizeof(forces_interim));
		memset(potE_interim, 0, sizeof(potE_interim));
	}

	void addParticle(const Float3& position, int atomtype_id, char atomLetter,
		int global_id, float boxlen_nm, BoundaryConditionSelect bc, float charge);


	void AddBond(const std::vector<ParticleToCompoundMapping>&, const SingleBondFactory&);
	void AddBond(const std::vector<ParticleToCompoundMapping>&, const AngleBondFactory&);
	void AddBond(const std::vector<ParticleToCompoundMapping>&, const DihedralBondFactory&);
	void AddBond(const std::vector<ParticleToCompoundMapping>&, const ImproperDihedralBondFactory&);

	bool hasRoomForRes(int n_particles_in_res) const {					// TODO: Implement, that it checks n atoms in res
		return ((int)n_particles + n_particles_in_res) <= MAX_COMPOUND_PARTICLES;
	}

	void addIdOfBondedCompound(int id);

	int id = -1;	// unique lima id

	Float3 positions[MAX_COMPOUND_PARTICLES];	// Extern positions [nm]
	int global_ids[MAX_COMPOUND_PARTICLES]{};		// For debug

	std::vector<int> local_atomtype_to_LATID_map;

	//std::unordered_set<int> connectedCompoundIds;
};


class BridgeFactory : public CompoundBridge {
public:
	BridgeFactory(int bridge_id, const std::vector<int>& _compound_ids) : bridge_id(bridge_id)
	{ 
		if (_compound_ids.size() > MAX_COMPOUNDS_IN_BRIDGE) {	// TODO: Move to .cpp and use format here
			throw std::runtime_error("Cannot add more compounds to a single bridge");
		}
		for (int i = 0; i < _compound_ids.size(); i++) {
			compound_ids[i] = _compound_ids[i];
		}
		n_compounds = _compound_ids.size();
	}

	void AddBond(const ParticleToCompoundMap&, ParticleToBridgeMap&, const SingleBondFactory&);
	void AddBond(const ParticleToCompoundMap&, ParticleToBridgeMap&, const AngleBondFactory&);
	void AddBond(const ParticleToCompoundMap&, ParticleToBridgeMap&, const DihedralBondFactory&);
	void AddBond(const ParticleToCompoundMap&, ParticleToBridgeMap&, const ImproperDihedralBondFactory&);

	const int bridge_id;
private:
//	std::vector<int> globalIdsOfParticlesInBridge;

	// Integrates the particle if it is not already, and returns its index relative to this bridge
	//uint8_t getBridgelocalIdOfParticle(ParticleInfo& particle_info);

	// Add particle to bridge and augment its particle info with the references to its position in this bridge
	// Only called from addBond, when a bond contains a particle that does NOT already exist in bridge
	void addParticle(const ParticleToCompoundMapping& p2cMapping, std::optional<ParticleToBridgeMapping>& p2bMapping);

	template <typename BondFactory_type>
	std::array<uint8_t, BondFactory_type::nAtoms> ConvertGlobalIdsToBridgelocalIds(
		const std::vector<ParticleToCompoundMapping>& p2cMap, 
		std::vector<std::optional<ParticleToBridgeMapping>>& p2bMap, 
		const BondFactory_type& bond);
};

// A translation unit between Gro file representation, and LIMA Box representation
struct BoxImage {
	const std::vector<CompoundFactory> compounds;
	const int32_t total_compound_particles;

	std::vector<BondedParticlesLUT> bpLutCollection;

	std::unique_ptr<CompoundBridgeBundleCompact> bridgebundle;

	const std::vector<Float3> solvent_positions;

	//const std::vector<ParticleInfo> particleinfos;

	GroFile grofile;

	const ForceField_NB forcefield;

	LIMA_MOLECULEBUILD::Topology topology; // This is only used for debugging purposes
};
 