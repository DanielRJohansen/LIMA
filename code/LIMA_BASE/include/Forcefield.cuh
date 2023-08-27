#pragma once

#include "Filehandling.h"
#include "Bodies.cuh"

#include <string>
#include <vector>
#include <array>
#include <map>
//#include <span>

constexpr int ATOMTYPE_SOLVENT = 0;

using std::string;


class Forcefield {
public:
	struct Topology {
		std::vector<SingleBond> singlebonds;
		std::vector<AngleBond> anglebonds;
		std::vector<DihedralBond> dihedralbonds;
		std::vector<ImproperDihedralBond> improperdihedralbonds;
	};

	Forcefield(VerbosityLevel vl);


	void loadForcefield(std::string molecule_dir);
	int getAtomtypeID(int global_id) const;


	const Topology& getTopology() const { return topology; }



	ForceField_NB getNBForcefield() const {
		return forcefield_nb;
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
	ForceField_NB forcefield_nb;
	std::map<int, int> globaldToAtomtypeMap;

	Topology topology;


	VerbosityLevel vl = SILENT;

	std::vector<NBAtomtype> loadAtomTypes(const SimpleParsedFile& nonbonded_parsed);
	std::map<int, int> loadAtomTypeMap(const SimpleParsedFile& nonbonded_parsed);

	Topology loadTopology(const SimpleParsedFile& bonded_parsed);

	void loadAtomypesIntoForcefield(const std::vector<NBAtomtype>& atomtypes);

};






