#pragma once

#include "LIMA_BASE/include/Filehandling.h"
#include "LIMA_BASE/include/Bodies.cuh"

#include <string>
#include <vector>
#include <map>


#define ATOMTYPE_SOL 0

using std::string;


class Forcefield {
	struct Topology {
		std::vector<SingleBond> singlebonds;
		std::vector<AngleBond> anglebonds;
		std::vector<DihedralBond> dihedralbonds;
		std::vector<ImproperDihedralBond> improperdihedralbonds;
	};

public:
	Forcefield(VerbosityLevel vl);


	void loadForcefield(std::string molecule_dir);
	int getAtomtypeID(int global_id) const;

	const SingleBond& getSinglebondtype(int bond_index, std::array<int, 2> gro_ids) const;
	const AngleBond& getAnglebondtype(int bond_index, std::array<int, 3> gro_ids) const;
	const DihedralBond& getDihedralbondtype(int bond_index, std::array<int, 4> gro_ids) const;
	const ImproperDihedralBond& getImproperdihedralbondtype(int bond_index, std::array<int, 4> gro_ids) const;
	

	ForceField_NB getNBForcefield() const {
		return forcefield;
	}


	int atomTypeToIndex(const char& atom) const {
		if (atom == 'C')
			return 1;
		if (atom == 'O')
			return 2;
		if (atom == 'N')
			return 3;
		if (atom == 'H')
			return 4;
		if (atom == 'P')
			return 5;
		if (atom == 'S')
			return 6;
		printf("Unable to find atom %c\n", atom);
		exit(1);
	}

	bool forcefield_loaded = false;

private:
	ForceField_NB forcefield;
	std::map<int, int> groIdToAtomtypeMap;

	Topology topology;


	VerbosityLevel vl = SILENT;


	enum STATE { INACTIVE, FF_NONBONDED, NB_ATOMTYPES, BONDS, ANGLES, DIHEDRALS, IMPROPERDIHEDRALS };
	static STATE setState(std::string s, STATE current_state) {
		if (s == "ff_nonbonded")
			return FF_NONBONDED;
		if (s == "atoms")
			return NB_ATOMTYPES;
		if (s == "bonds")
			return BONDS;
		if (s == "angles")
			return ANGLES;
		if (s == "dihedrals")
			return DIHEDRALS;
		if (s == "improperdihedrals")
			return IMPROPERDIHEDRALS;
		return current_state;
	}

	bool newParseTitle(std::vector<std::string> row) {
		return (row[0][0] == '#');
	}

	std::vector<NBAtomtype> loadAtomTypes(const SimpleParsedFile& nonbonded_parsed);
	std::map<int, int> loadAtomTypeMap(const SimpleParsedFile& nonbonded_parsed);

	Topology loadTopology(const SimpleParsedFile& bonded_parsed);

	void loadAtomypesIntoForcefield(const std::vector<NBAtomtype>& atomtypes);

};






