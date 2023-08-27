#pragma once


#include "LIMA_BASE/include/Bodies.cuh"
#include <string.h>
#include <fstream>
#include <vector>
#include <array>

#include "LIMA_BASE/include/Constants.cuh"
#include "LIMA_BASE/include/Forcefield.cuh"
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
		std::unique_ptr<LimaLogger> logger,
		bool ignore_hydrogens = true
		);
}


















enum TopologyMode { INACTIVE, ATOMS, BOND, ANGLE, DIHEDRAL };

struct GroRecord {
	int residue_number{};
	std::string residue_name{};
	std::string atom_name{};
	int gro_id{};
	Float3 position{};
	Float3 velocity{};
};


struct AtomEntry {	// Vague defined class, temporarily holding information
	AtomEntry(const int global_id, const int gro_id, const Float3& position, const std::string name) :
		global_id(global_id), gro_id(gro_id), position(position), name(name) {}

	const int global_id;			// Unique, given by LIMA
	const int gro_id;				// Not unique but cyclic, given by gromacs
	const Float3 position;			// nm
	const std::string name;
	//const Float3 velocity;	// to be added later...
};


struct Residue {
	// Ballpark max size. Dont use for hard constraints, this is only for loosely reserving memory
	static const int max_size_soft = 24;

	Residue(int gro_id, int global_id, const std::string& name, int chain_id) :
		gro_id(gro_id), global_id(global_id), name(name), chain_id(chain_id) {
		atoms.reserve(max_size_soft);
	}

	const int gro_id;					// NOT UNIQUE. Used ONLY to spot when a new residue occurs in a gro file
	const int global_id;				// Unique id given by LIMA
	const std::string name;				// 3 letter residue name
	const int chain_id;					// unique, given by lima

	std::vector<AtomEntry> atoms;
	std::vector<int> bondedresidue_ids;	// Lima Ids of residue with which this shares a singlebond
};

struct ProteinChain {
	struct Topology {
		std::vector<SingleBond> singlebonds;
		std::vector<AngleBond> anglebonds;
		std::vector<DihedralBond> dihedralbonds;
		std::vector<ImproperDihedralBond> improperdihedralbonds;
	};

	std::vector<Residue> residues;
	Topology topology;
};

struct ParticleInfo {
	// Given once created from .gro file
	bool isUsed = false;						// false if discarded at load, or simply never found
	int gro_id = -1;
	int chain_id = -1;


	std::vector<int> singlebonds_indices;		// lima indices of the bond in moleculebuilder
	//std::vector<int> anglebonds_indices;		// lima indices of the bond in moleculebuilder
	std::vector<int> dihedralbonds_indices;		// lima indices of the bond in moleculebuilder

	// First available when added to compound
	int compound_index = -1;
	int local_id_compound = -1;	// index in compound

	// First available (if) when added to bridge
	int bridge_id = -1;
	int local_id_bridge = -1;	// index in bridge
	
};


class ParticleInfoLUT {
	const int max_chains = 1000;	// TODO: This is a shitty temporary solution 
	int current_size = 0;
	std::vector<ParticleInfo> particle_info;

	// chains x gro_id
	std::vector<std::vector<int>> id_map;


public:
	ParticleInfoLUT(int size) {
		particle_info.resize(size);

		id_map.resize(max_chains);
		for (auto& e : id_map) {
			e.resize(size);
			for (auto& mapping : e) {
				mapping = -1;
			}
		}
	}

	// If the numbering in the .gro file is funky we can increase the capacity, but that is very expensive
	void increaseCapacity(int minsize) {
		particle_info.resize(minsize * 2);

		for (auto& e : id_map) {
			e.resize(minsize);
		}
	}

	ParticleInfo& get(int global_id) {
		assert(global_id < particle_info.size());
		return particle_info[global_id];
	}
	const ParticleInfo& getConst(int global_id) const {
		assert(global_id < particle_info.size());
		return particle_info[global_id];
	}

	// This fucking bullshit is necessary since topology uses the non-unique gro-ids to determine the bonds
	ParticleInfo& get(int chain_id, int gro_id);

	int size() { return particle_info.size(); }
};


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

	template <typename Bondtype>
	void addBond(const ParticleInfoLUT&, const Bondtype&);

	template <> void addBond(const ParticleInfoLUT&, const SingleBond&);
	template <> void addBond(const ParticleInfoLUT&, const AngleBond&);
	template <> void addBond(const ParticleInfoLUT&, const DihedralBond&);
	template <> void addBond(const ParticleInfoLUT&, const ImproperDihedralBond&);

	bool hasRoomForRes(int n_particles_in_res) const {					// TODO: Implement, that it checks n atoms in res
		return ((int)n_particles + n_particles_in_res) <= MAX_COMPOUND_PARTICLES;
	}

	Float3 positions[MAX_COMPOUND_PARTICLES];	// Extern positions [nm]
	int gro_ids[MAX_COMPOUND_PARTICLES]{};		// For debug
};


class BridgeFactory : public CompoundBridge {
public:
	BridgeFactory(int bridge_id, const std::array<int, 2>& compound_ids) :bridge_id(bridge_id){
		compound_id_left = compound_ids[0];
		compound_id_right = compound_ids[1];
	}

	// Augments the particle_info with local_id_bridge if necessary
	template <typename Bondtype>
	void addBond(ParticleInfoLUT&, const Bondtype&);

	template <> void addBond(ParticleInfoLUT&, const SingleBond&);
	template <> void addBond(ParticleInfoLUT&, const AngleBond&);
	template <> void addBond(ParticleInfoLUT&, const DihedralBond&);
	template <> void addBond(ParticleInfoLUT&, const ImproperDihedralBond&);

private:
	// Integrates the particle if it is not already, and returns its index relative to this bridge
	uint32_t getBridgelocalIdOfParticle(ParticleInfo& particle_info);
	const int bridge_id;
};


struct CompoundCollection {
	const std::vector<CompoundFactory> compounds;
	const int32_t total_compound_particles;

	std::unique_ptr<BondedParticlesLUTManager> bp_lut_manager;

	std::unique_ptr<CompoundBridgeBundleCompact> bridgebundle;

	const std::vector<Float3> solvent_positions;
};