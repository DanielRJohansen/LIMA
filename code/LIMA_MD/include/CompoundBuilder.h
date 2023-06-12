#pragma once


#include "LIMA_BASE/include/Bodies.cuh"
#include <string.h>
#include <fstream>
#include <vector>
#include <array>

#include "LIMA_BASE/include/Constants.cuh"
#include "LIMA_ENGINE/include/Forcefield.cuh"
#include "LIMA_BASE/include/Utilities.h"




struct CompoundCollection;

namespace LIMA_MOLECULEBUILD {
	// This is the only function to be called from outside :)
	CompoundCollection buildMolecules(
		Forcefield* ff, 
		const std::string& work_dir, 
		VerbosityLevel vl,
		const string& gro_path, 
		const string& topol_path, 
		bool ignore_hydrogens=true);
}


















enum TopologyMode { INACTIVE, ATOMS, BOND, ANGLE, DIHEDRAL };

struct GroRecord {
	int residue_number{};
	std::string residue_name{};
	std::string atom_name{};
	int atom_number{};
	Float3 position{};
	Float3 velocity{};
};


struct AtomEntry {	// Vague defined class, temporarily holding information
	AtomEntry(const int id, const Float3& position, const std::string name) :
		id(id), position(position), name(name) {}

	const int id;				// Unique 1-indexed id, given by gro (I hope its unique..)
	const Float3 position;
	const std::string name;
	//const Float3 velocity;	// to be added later...
};


struct Residue {
	// Ballpark max size. Dont use for hard constraints, this is only for loosely reserving memory
	static const int max_size_soft = 24;

	Residue(int gro_id, int id, const std::string& name) :
		gro_id(gro_id), id(id), name(name) {
		atoms.reserve(max_size_soft);
	}

	const int gro_id;					// NOT UNIQUE. Used ONLY to spot when a new residue occurs in a gro file
	const int id;						// Unique id given by LIMA
	const std::string name;				// 3 letter residue name

	std::vector<AtomEntry> atoms;
	std::vector<int> bondedresidue_ids;	// Lima Ids of residue with which this shares a singlebond
};


struct ParticleInfo {
	bool isUsed = false;						// false if discarded at load, or simply never found
	int gro_id = -1;

	std::vector<int> singlebonds_indices;		// lima indices of the bond in moleculebuilder
	std::vector<int> anglebonds_indices;		// lima indices of the bond in moleculebuilder
	std::vector<int> dihedralbonds_indices;		// lima indices of the bond in moleculebuilder

	// First available when added to compound
	int compound_index = -1;
	int local_id_compound = -1;	// index in compound

	// First available (if) when added to bridge
	int local_id_bridge = -1;	// index in bridge
};


class CompoundFactory : public Compound {
public:
	CompoundFactory() {}
	CompoundFactory(const int id) : id(id) {}

	void addParticle(const Float3& position, int atomtype_id, int atomtype_color_id, int gro_id);
	int id = -1;	// unique lima id

	void addSingleBond(const std::array<ParticleInfo, 2>& particle_info, const SingleBond& bondtype);
	void addAngleBond(const std::array<ParticleInfo, 3>& particle_info, const AngleBond& bondtype);
	void addDihedralBond(const std::array<ParticleInfo, 4>& particle_info, const DihedralBond& bondtype);

	bool hasRoomForRes(int n_particles_in_res) const {					// TODO: Implement, that it checks n atoms in res
		return ((int)n_particles + n_particles_in_res) <= MAX_COMPOUND_PARTICLES;
	}

	Float3 positions[MAX_COMPOUND_PARTICLES];	// Extern positions [nm]
	int gro_ids[MAX_COMPOUND_PARTICLES]{};		// For debug
};


class BridgeFactory : public CompoundBridge {
public:
	BridgeFactory(const std::array<int, 2>& compound_ids) {
		compound_id_left = compound_ids[0];
		compound_id_right = compound_ids[1];
	}

	// Augments the particle_info with local_id_bridge if necessary
	void addSingleBond(std::array<ParticleInfo*, 2> particle_info, const SingleBond& bondtype);
	// Augments the particle_info with local_id_bridge if necessary
	void addAngleBond(std::array<ParticleInfo*, 3> particle_info, const AngleBond& bondtype);
	// Augments the particle_info with local_id_bridge if necessary
	void addDihedralBond(std::array<ParticleInfo*, 4> particle_info, const DihedralBond& bondtype);


private:
	// Integrates the particle if it is not already, and returns its index relative to this bridge
	int getBridgelocalIdOfParticle(ParticleInfo& particle_info);
};


struct CompoundCollection {
	const std::vector<CompoundFactory> compounds;
	const int32_t total_compound_particles;

	std::unique_ptr<BondedParticlesLUTManager> bp_lut_manager;

	std::unique_ptr<CompoundBridgeBundleCompact> bridgebundle;

	const std::vector<Float3> solvent_positions;
};