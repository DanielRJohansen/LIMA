#pragma once


#include "Bodies.cuh"
#include <string.h>
#include <fstream>
#include <vector>
#include <array>

#include "Constants.cuh"
#include "Forcefield.cuh"


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


	enum TopologyMode { INACTIVE, ATOMS, BOND, ANGLE, DIHEDRAL };
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





