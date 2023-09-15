#pragma once



#include <string.h>
#include <fstream>
#include <vector>
#include <array>
#include <memory>

#include "Bodies.cuh"
#include "Constants.h"
#include "Forcefield.cuh"
#include "Utilities.h"


struct CompoundCollection;

namespace LIMA_MOLECULEBUILD {
	// This is the only function to be called from outside :)
	CompoundCollection buildMolecules(
		Forcefield* ff,					// Can be removed when we dont need to do the stupid color lookup anymore
		const std::string& work_dir, 
		VerbosityLevel vl,
		const string& gro_path, 
		const string& topol_path, 
		std::unique_ptr<LimaLogger> logger,
		bool ignore_hydrogens = true
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
	// Added when determining which residues of the chains are connected
	std::vector<int> singlebonds_indices;	// indices of singlebonds of which this atom is in

	// Added when creating compounds
	int compound_index = -1;
	int local_id_compound = -1;

	// Added when allocated to a bridge
	int bridge_id = -1;
	int local_id_bridge = -1;
};
using ParticleInfoTable = std::vector<ParticleInfo>;



class CompoundFactory : public Compound {
public:
	CompoundFactory() {}
	CompoundFactory(const int id) : id(id) {
		for (int i = 0; i < MAX_COMPOUND_PARTICLES; i++) {
			potE_interim[i] = 0.f;
		}
	}

	void addParticle(const Float3& position, int atomtype_id, int atomtype_color_id, int gro_id);
	int id = -1;	// unique lima id


	void addBond(const ParticleInfoTable&, const SingleBond&);
	void addBond(const ParticleInfoTable&, const AngleBond&);
	void addBond(const ParticleInfoTable&, const DihedralBond&);
	void addBond(const ParticleInfoTable&, const ImproperDihedralBond&);

	bool hasRoomForRes(int n_particles_in_res) const {					// TODO: Implement, that it checks n atoms in res
		return ((int)n_particles + n_particles_in_res) <= MAX_COMPOUND_PARTICLES;
	}

	Float3 positions[MAX_COMPOUND_PARTICLES];	// Extern positions [nm]
	int global_ids[MAX_COMPOUND_PARTICLES]{};		// For debug
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


	void addBond(ParticleInfoTable&, const SingleBond&);
	void addBond(ParticleInfoTable&, const AngleBond&);
	void addBond(ParticleInfoTable&, const DihedralBond&);
	void addBond(ParticleInfoTable&, const ImproperDihedralBond&);

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
	uint32_t getBridgelocalIdOfParticle(ParticleInfo& particle_info);

	// Add particle to bridge and augment its particle info with the references to its position in this bridge
	// Only called from addBond, when a bond contains a particle that does NOT already exist in bridge
	void addParticle(ParticleInfo&);
};


struct CompoundCollection {
	const std::vector<CompoundFactory> compounds;
	const int32_t total_compound_particles;

	std::unique_ptr<BondedParticlesLUTManager> bp_lut_manager;

	std::unique_ptr<CompoundBridgeBundleCompact> bridgebundle;

	const std::vector<Float3> solvent_positions;
};
