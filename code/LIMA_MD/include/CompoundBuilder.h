#pragma once


#include "Bodies.cuh"
#include <string.h>
#include <fstream>
#include <vector>
#include <array>

#include "Constants.cuh"
#include "Forcefield.cuh"
#include "Utilities.h"

using namespace std;


enum LineType { atom, singlebond, anglebond, torsionbond };

struct ParsedLine {
	LineType linetype;
	int pdb_index;
	string id;
	string monomer;
	Float3 position;
	string atom;
};


struct Topology {

	void addAtomsdataEntry(uint32_t nr, string type, uint32_t res_nr) { atoms_data.push_back({ nr, type, res_nr }); }
	void addBondsdataEntry(uint32_t ai, uint32_t aj, int funct) { bonds_data.push_back({ ai, aj, funct }); }


	// Data types templated directly on GMX default .top files sections
	struct atoms_data_entry {
		uint32_t nr{};
		string type{};
		uint32_t res_nr{};
	};

	struct bonds_data_entry {
		uint32_t ai{};
		uint32_t aj{};
		int funct{};
	};

	vector<atoms_data_entry> atoms_data;
	vector<bonds_data_entry> bonds_data;
};

enum TopologyMode { INACTIVE, ATOMS, BOND, ANGLE, DIHEDRAL };

// The compound-builder must be unique to a single conf.gro-topol.top file pair!
using namespace std;
class CompoundBuilder
{
	struct IDMap {				// Delete for particle refs!
		IDMap() {}	
		IDMap(int global_id, int compound_id, int local_id) : global_id(global_id), compound_id(compound_id), local_id(local_id) {}
		int global_id = -1, compound_id = -1, local_id = -1;	// local_id is local to the compound
	};
public:
	CompoundBuilder() {}
	CompoundBuilder(Forcefield* ff, VerbosityLevel vl = SILENT);
	CompoundCollection buildCompoundCollection(string gro_path, string itp_path, uint32_t max_residue_id=UINT32_MAX, uint32_t min_residue_id=0, bool ignore_hydrogens=true);

	vector<Float3> getSolventPositions(string gro_path);


private:
	Forcefield* forcefield = nullptr;
	//IDMap* particle_id_maps;
	ParticleRef* particle_id_maps = nullptr;
	CompoundBridgeBundle* compound_bridge_bundle = nullptr;

	uint16_t unique_doublyconnected_id = 1;

	VerbosityLevel verbosity_level = SILENT;
	//uint32_t** bonded_interactions_list;	// Contains each particles list of (larger) ids of particles with which it shares a bonded interaction
	//LJ_Ignores* bonded_interactions_list = nullptr;

	// Count to keep track of multiple topol sections named "dihedral". This can be done smarter..
	int dihedral_sections_count = 0;


	struct Record_ATOM;
	void loadParticles(CompoundCollection* molecule, vector<Record_ATOM>* pdb_data, uint32_t max_monomer_id = INT32_MAX, uint32_t min_residue_id=0, bool ignore_protons =false);
	void loadTopology(CompoundCollection* molecule, vector<vector<string>>* itp_data);
	void calcParticleSphere(Compound* compound);


	bool setMode(vector<string>& entry, TopologyMode& current_mode);
	void loadMaps(ParticleRef* maps, vector<string>* record, int n);
	void addGeneric(CompoundCollection* molecule, vector<string>* record, TopologyMode mode);
	void addBond(CompoundCollection* molecule, ParticleRef* maps, vector<string>* record);
	void addAngle(CompoundCollection* molecule, ParticleRef* maps, vector<string>* record);
	void addDihedral(CompoundCollection* molecule, ParticleRef* maps, vector<string>* record);
	void distributeLJIgnores(CompoundCollection* molecule, ParticleRef* maps, int n);
	//bool checkIfFirstBondedInteraction(CompoundCollection* molecule, ParticleRef* maps, int n);



	// Topology related stuff
	void assignMoleculeIDs();
	bool areResiduesBonded(uint32_t res1, uint32_t res2) const;
	uint32_t getFirstAtomindexOfRes(vector<Record_ATOM>& atom_data, uint32_t res_id);
	vector<uint32_t> makeResidueIdToAtomindexMap();
	vector<vector<uint32_t>> makeBondedatomsLUT(const string& top_path );

	vector<uint32_t> residueId_to_firstatomindex;	// GMX 1-indexing
	vector<vector<uint32_t>> bonded_atoms;		// GMX 1-indexing

	//Topology topology;
	vector<Record_ATOM> atom_data;

	ParsedLine parseLine(int line_index);
	ParsedLine parseAtom(string line);
	ParsedLine parseConnection(string line);

	vector<vector<string>> parseTOP(string path);
	Topology parseTop1(string path);
	vector<Record_ATOM> parsePDB(string path);
	vector<Record_ATOM> parseGRO(string path);


	struct ResidueComboId;

	ResidueComboId parseResidueID(string s);
	inline bool isAsciiNumber(char c) { return (c > 47 && c < 58); }
	void countElements(CompoundCollection* molecule);
	vector<string> splitAtomnameFromId(vector<string> words);







	struct ResidueComboId {	// Combined name and id.
		ResidueComboId(){}
		ResidueComboId(int id, string name) : id(id), name(name) { 
			valid = true; 
		}
		bool valid = false;
		int id = 0;
		string name;
	};





	struct Record_ATOM {	//	https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
		Record_ATOM(int a, string b, char c, string d, char e, int f, char g, Float3 h) {
			atom_serial_number = a;
			atom_name = b;
			alternate_location_indicator = c;
			residue_name = d;
			chain_identifier = e;
			residue_seq_number = f;
			code_for_insertions_of_residues = g;
			coordinate = h;
		}

		int atom_serial_number = 0;
		string atom_name = "";
		char alternate_location_indicator = 0;
		string residue_name = "";
		char chain_identifier = 0;
		uint32_t residue_seq_number = 0;
		char code_for_insertions_of_residues = 0;
		Float3 coordinate;						// [nm]
		int moleculeID = -1;	// Given by LIMA
	};




};















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
		gro_id(gro_id), id(id), name(name) { atoms.reserve(max_size_soft); }

	const int gro_id;					// NOT UNIQUE. Used ONLY to spot when a new residue occurs in a gro file
	const int id;						// Unique id given by LIMA
	const std::string name;				// 3 letter residue name

	std::vector<AtomEntry> atoms;
	std::vector<int> bondedresidue_ids;	// Lima Ids of residue with which this shares a singlebond
};

struct ParticleInfo {
	bool isUsed = false;						// false if discarded at load, or simply never found

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

	void addParticle(const Float3& position, int atomtype_id, int atomtype_color_id);
	int id = -1;	// unique lima id

	void addSingleBond(std::array<ParticleInfo, 2> particle_info, const SingleBond& bondtype);
	void addAngleBond(std::array<ParticleInfo, 3> particle_info, const AngleBond& bondtype);
	void addDihedralBond(std::array<ParticleInfo, 4> particle_info, const DihedralBond& bondtype);


	Float3 positions[MAX_COMPOUND_PARTICLES];	// Extern positions [nm]
};

class BridgeFactory : public CompoundBridgeCompact {
public:
	BridgeFactory( const std::array<int, 2> compound_ids) {
		compound_id_left = compound_ids[0];
		compound_id_right = compound_ids[1];
	}

	void addSingleBond(std::array<ParticleInfo, 2> particle_info, const SingleBond& bondtype);
	void addAngleBond(std::array<ParticleInfo, 3> particle_info, const AngleBond& bondtype);
	void addDihedralBond(std::array<ParticleInfo, 4> particle_info, const DihedralBond& bondtype);


private:
	// Integrates the particle if it is not already, and returns its index relative to this bridge
	int getBridgelocalIdOfParticle(ParticleInfo& particle_info);
};

struct CompoundCollection2 {
	const std::vector<CompoundFactory>& compounds;
	const int64_t total_compound_particles;

	BondedParticlesLUTManager* bp_lut_manager;

	const CompoundBridgeBundleCompact bridgebundle;


	const std::vector<Float3>& solvent_positions;
};


class MoleculeBuilder {
public:
	MoleculeBuilder(Forcefield* ff, const std::string& work_dir="", VerbosityLevel vl = SILENT);

	CompoundCollection2 buildMolecules(const string& gro_path, const string& topol_path, bool ignore_hydrogens = true);

private:
	// Members, seen in the order they are filled in
	std::vector<Residue> residues;
	std::vector<Float3> solvent_positions;
	int64_t n_particles_in_residues = 0;

	std::vector<std::array<int, 2>> singlebonds;	// gro id's of particles
	std::vector<std::array<int, 3>> anglebonds;		// gro id's of particles
	std::vector<std::array<int, 4>> dihedralbonds;	// gro id's of particles
	std::vector<ParticleInfo> particle_info;		// Uses gro 1-indexed 

	std::vector<CompoundFactory> compounds;
	std::vector<BridgeFactory> compound_bridges;

	BondedParticlesLUTManager* bp_lut_manager = nullptr;






	LimaLogger logger;
	const VerbosityLevel verbosity_level;
	const Forcefield* forcefield;






	// ------------------------------------ HELPER FUNCTIONS ------------------------------------ //

	// Only works for fixed positions gro files!! Meaning max, 99.999 particles
	void loadResiduesAndSolvents(const std::string gro_path);

	// Loads all bonds, and sets references to these for each particle - HEAVY!
	void loadTopology(const std::string& topol_path);

	/// <summary>
	/// Goes through all the residues, and fills the bondedresidue_ids vector. 
	/// The particle_bonds_lut is used to determine whether particles of two different 
	/// share a bond.
	/// </summary>
	void matchBondedResidues();

	// TODO: Put correct references in the bridges too!
	void createCompoundsAndBridges();

	void distributeBondsToCompoundsAndBridges();

	template <int n_ids>
	bool spansTwoCompounds(std::array<int, n_ids> bond_ids);
};