// Superstructure to the Filehandler. 
// Provides special functionality for MD files such as .gro .itp .top .pdb and more
#pragma once

#include <filesystem>
#include "Bodies.cuh"
#include "Filehandling.h"


struct GroRecord {
	int residue_number{};
	std::string residue_name{};
	std::string atom_name{};
	int gro_id{};
	Float3 position{};
	Float3 velocity{};
};

struct ParsedGroFile {
	std::string title;
	int n_atoms{ -1 };
	std::vector<GroRecord>atoms;
	Float3 box_size{};
};

struct ParsedTopologyFile {	//.top or .itp
	struct MoleculetypeEntry {
		std::string name;
		int nrexcl;
	};

	// Variable names are from .itp file
	struct AtomsEntry {
		int nr;	// Not guarantee'd to be unique, atleast not with multiple files!
		std::string type;
		int resnr;
		std::string residue;
		std::string atom;
		int cgnr;
		float charge;
		float mass;
	};


	struct SingleBond {	// IDS and indexes are used interchangeably here!
		int atom_indexes[2] = { 0,0 };
		int funct;
	};
	struct Pair {
		int atom_indexes[2] = { 0,0 };
		int funct;
	};
	struct AngleBond {
		int atom_indexes[3] = { 0,0,0 }; // i,j,k angle between i and k
		int funct;
	};
	struct DihedralBond {
		int atom_indexes[4] = { 0,0,0,0 };
		int funct;
	};
	struct ImproperDihedralBond {
		int atom_indexes[4] = { 0,0,0,0 };
		int funct;
	};

	template <typename EntryType>
	struct Section {
		std::string title;
		std::string legend;
		std::vector<EntryType> entries;
	};


	std::string title;
	int n_atoms{ -1 };

	//std::vector<MoleculetypeEntry>moleculetypes;
	Section<MoleculetypeEntry> moleculetypes{ "hey", "ho" };

	Section<AtomsEntry> atoms{ "[ atoms ]",
		";   nr       type  resnr residue  atom   cgnr       charge     mass" };

	Section<SingleBond> singlebonds{ "[ bonds ]",
		";  ai    aj funct            c0            c1            c2            c3" };
	Section<Pair> pairs{ "[ pairs ]",
	";  ai    aj funct            c0            c1            c2            c3" };
	Section<AngleBond> anglebonds{ "[ angles ]",
	";  ai    aj    ak funct            c0            c1            c2            c3" };
	Section<DihedralBond> dihedralbonds{ "[ dihedrals ]",
	";  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5" };
	Section<ImproperDihedralBond> improperdihedralbonds{ "[ dihedrals ]",
	";  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5" };

	Float3 box_size{};
};


namespace MDFiles {
	namespace fs = std::filesystem;

	ParsedGroFile loadGroFile(const fs::path& path);
	ParsedTopologyFile loadTopologyFile(const fs::path& path);
}