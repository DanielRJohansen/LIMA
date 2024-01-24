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
		const ParsedGroFile& gro_file,
		const ParsedTopologyFile& top_file,
		VerbosityLevel vl,
		std::unique_ptr<LimaLogger>,
		bool ignore_hydrogens,
		BoundaryConditionSelect bc_select
	);
}






struct Residue {
	// Ballpark max size. Dont use for hard constraints, this is only for loosely reserving memory
	static const int max_size_soft = 24;

	Residue(int gro_id, int global_id, const std::string& name, int chain_id) :
		gro_id(gro_id), global_id(global_id), name(name), chain_id(chain_id) {
		//atoms.reserve(max_size_soft);
		atoms_globalid.reserve(max_size_soft);
	}

	const int gro_id;					// NOT UNIQUE. Used ONLY to spot when a new residue occurs in a gro file
	const int global_id;				// Unique id given by LIMA
	const std::string name;				// 3 letter residue name
	const int chain_id;					// unique, given by lima


	std::vector<int> atoms_globalid;	// GID of atoms in residue
	//std::vector<int> bondedresidue_ids;	// GIDs of residue with which this shares a singlebond
};

struct ParticleInfo {
	// All member variables are in the order they get assigned

	// Added when the entry is created, as the nonbonded.lff file is read.
	const int global_id;
	const int gro_id;
	const int chain_id;
	const int atomtype_id;
	const std::string atomname;
	const int residue_groid;
	const int unique_res_id;
	// Added when determining which residues of the chains are connected
	std::vector<int> singlebonds_indices;	// indices of singlebonds of which this atom is in

	// Added when creating compounds
	int compound_index = -1;
	uint8_t local_id_compound = 255;

	// Added when allocated to a bridge
	int bridge_id = -1;
	uint8_t local_id_bridge = 255;
};
using ParticleInfoTable = std::vector<ParticleInfo>;



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

	void addParticle(const Float3& position, int atomtype_id, int atomtype_color_id, int global_id, float boxlen_nm, BoundaryConditionSelect bc) {
		if (n_particles >= MAX_COMPOUND_PARTICLES) {
			throw std::runtime_error("Failed to add particle to compound");
		}

		// Hyperposition each compound relative to particle 0, so we can find key_particle and radius later
		Float3 hyperpos = position;
		if (n_particles > 0) {
			//LIMAPOSITIONSYSTEM::applyHyperposNM<BoundaryCondition>(&positions[0], &hyperpos);
			BoundaryConditionPublic::applyHyperposNM(&positions[0], &hyperpos, boxlen_nm, bc);
		}

		// Variables only present in factory
		positions[n_particles] = hyperpos;
		global_ids[n_particles] = global_id;

		// Variables present in Compound
		atom_types[n_particles] = atomtype_id;
		atom_color_types[n_particles] = atomtype_color_id;	// wtf is this

		n_particles++;
	}
	int id = -1;	// unique lima id


	void addBond(const ParticleInfoTable&, const SingleBondFactory&);
	void addBond(const ParticleInfoTable&, const AngleBondFactory&);
	void addBond(const ParticleInfoTable&, const DihedralBondFactory&);
	void addBond(const ParticleInfoTable&, const ImproperDihedralBondFactory&);

	bool hasRoomForRes(int n_particles_in_res) const {					// TODO: Implement, that it checks n atoms in res
		return ((int)n_particles + n_particles_in_res) <= MAX_COMPOUND_PARTICLES;
	}

	void addIdOfBondedCompound(int id);

	Float3 positions[MAX_COMPOUND_PARTICLES];	// Extern positions [nm]
	int global_ids[MAX_COMPOUND_PARTICLES]{};		// For debug


	std::vector<int> local_atomtype_to_LATID_map;
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


	void addBond(ParticleInfoTable&, const SingleBondFactory&);
	void addBond(ParticleInfoTable&, const AngleBondFactory&);
	void addBond(ParticleInfoTable&, const DihedralBondFactory&);
	void addBond(ParticleInfoTable&, const ImproperDihedralBondFactory&);

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

	const ParticleInfoTable particleinfotable;

	ParsedGroFile grofile;
};
