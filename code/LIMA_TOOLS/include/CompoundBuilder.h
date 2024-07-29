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

namespace LIMA_MOLECULEBUILD {

	std::unique_ptr<BoxImage> buildMolecules(
		//Forcefield* ff,					// TODO: Can be removed when we dont need to do the stupid color lookup anymore
		const std::string& molecule_dir,	// We need access to .lff files aswell
		const GroFile& gro_file,
		const TopologyFile& top_file,
		VerbosityLevel vl,
		std::unique_ptr<LimaLogger>,
		bool ignore_hydrogens,
		const SimParams& simparams
	);
}









// --------------------------------- Bond Factories --------------------------------- //

struct SingleBondFactory {
	static const int n_atoms = 2;
	SingleBondFactory(const std::array<uint32_t, 2>& ids, const SingleBond::Parameters& parameters);
	SingleBond::Parameters params;
	uint32_t global_atom_indexes[2] = { 0,0 };
	//std::array<uint32_t, n_atoms> global_atom_indexes;
};

struct AngleBondFactory {
	static const int n_atoms = 3;
	AngleBondFactory(std::array<uint32_t, 3> ids, const AngleBond::Parameters& bondparameters);
	AngleBond::Parameters params;
	uint32_t global_atom_indexes[3] = { 0,0, 0 };
};

struct DihedralBondFactory {
	static const int n_atoms = 4;
	DihedralBondFactory(std::array<uint32_t, 4> ids, const DihedralBond::Parameters& bondparameters);
	DihedralBond::Parameters params;
	uint32_t global_atom_indexes[4] = { 0,0, 0, 0 };
};

struct ImproperDihedralBondFactory {
	static const int n_atoms = 4;
	ImproperDihedralBondFactory(std::array<uint32_t, 4> ids, const ImproperDihedralBond::Parameters& bondparameters);
	ImproperDihedralBond::Parameters params;
	uint32_t global_atom_indexes[4] = { 0,0, 0, 0 };
};






















struct ParticleInfo {
	ParticleInfo() = default;
	ParticleInfo(const GroRecord* groAtom, const TopologyFile::AtomsEntry* topAtom, int activeLjTypeParameterIndex, int uniqueResId)
		: groAtom(groAtom), topAtom(topAtom), activeLjtypeParameterIndex(activeLjTypeParameterIndex), uniqueResId(uniqueResId) {}
	const GroRecord* groAtom = nullptr;
	const TopologyFile::AtomsEntry* topAtom = nullptr;
	int activeLjtypeParameterIndex = -1;
	int uniqueResId = -1;

	// Only available once the compounds have been created
	int compoundId = -1;
	uint8_t localIdInCompound = -1;

	// Only available once the bridges have been created
	int bridgeId = -1;
	uint8_t localIdInBridge = -1;
};

class CompoundFactory : public Compound {
public:
	CompoundFactory() {}
	CompoundFactory(const int id) : 
		id(id), local_atomtype_to_LATID_map(MAX_ATOM_TYPES, -1)
	{
		for (int i = 0; i < MAX_COMPOUND_PARTICLES; i++) {
			potE_interim[i] = 0.f;
		}
	}

	void addParticle(const Float3& position, int atomtype_id, char atomLetter, 
		int global_id, float boxlen_nm, BoundaryConditionSelect bc, float charge) {
		if (n_particles >= MAX_COMPOUND_PARTICLES) {
			throw std::runtime_error("Failed to add particle to compound");
		}

		// Hyperposition each compound relative to particle 0, so we can find key_particle and radius later
		Float3 hyperpos = position;
		if (n_particles > 0) {
			//LIMAPOSITIONSYSTEM::applyHyperposNM<BoundaryCondition>(&positions[0], &hyperpos);
			BoundaryConditionPublic::applyHyperposNM(positions[0], hyperpos, boxlen_nm, bc);
		}

		// Variables only present in factory
		positions[n_particles] = hyperpos;
		global_ids[n_particles] = global_id;

		// Variables present in Compound
		atom_types[n_particles] = atomtype_id;
		atomLetters[n_particles] = atomLetter;

		atom_charges[n_particles] = charge;

		n_particles++;
	}
	int id = -1;	// unique lima id


	void AddBond(const std::vector<ParticleInfo>&, const SingleBondFactory&);
	void AddBond(const std::vector<ParticleInfo>&, const AngleBondFactory&);
	void AddBond(const std::vector<ParticleInfo>&, const DihedralBondFactory&);
	void AddBond(const std::vector<ParticleInfo>&, const ImproperDihedralBondFactory&);

	bool hasRoomForRes(int n_particles_in_res) const {					// TODO: Implement, that it checks n atoms in res
		return ((int)n_particles + n_particles_in_res) <= MAX_COMPOUND_PARTICLES;
	}

	void addIdOfBondedCompound(int id);

	Float3 positions[MAX_COMPOUND_PARTICLES];	// Extern positions [nm]
	int global_ids[MAX_COMPOUND_PARTICLES]{};		// For debug


	std::vector<int> local_atomtype_to_LATID_map;



	std::vector<int> indicesOfBondconnectedCompounds;
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

	void AddBond(std::vector<ParticleInfo>&, const SingleBondFactory&);
	void AddBond(std::vector<ParticleInfo>&, const AngleBondFactory&);
	void AddBond(std::vector<ParticleInfo>&, const DihedralBondFactory&);
	void AddBond(std::vector<ParticleInfo>&, const ImproperDihedralBondFactory&);

	bool containsCompound(int compound_id) const {
		for (int i = 0; i < n_compounds; i++) {
			if (compound_ids[i] == compound_id)
				return true;
		}
		return false;
	}


	const int bridge_id;
private:
	// Integrates the particle if it is not already, and returns its index relative to this bridge
	uint8_t getBridgelocalIdOfParticle(ParticleInfo& particle_info);

	// Add particle to bridge and augment its particle info with the references to its position in this bridge
	// Only called from addBond, when a bond contains a particle that does NOT already exist in bridge
	void addParticle(ParticleInfo&);
};

// A translation unit between Gro file representation, and LIMA Box representation
struct BoxImage {
	const std::vector<CompoundFactory> compounds;
	const int32_t total_compound_particles;

	std::unique_ptr<BondedParticlesLUTManager> bp_lut_manager;

	std::unique_ptr<CompoundBridgeBundleCompact> bridgebundle;

	const std::vector<Float3> solvent_positions;

	const float box_size;

	//const ParticleInfoTable particleinfotable;
	const std::vector<ParticleInfo> particleinfos;

	GroFile grofile;

	const ForceField_NB forcefield;
};
