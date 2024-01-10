#pragma once

#include "Filehandling.h"
#include "Bodies.cuh"

#include <string>
#include <vector>
#include <array>
#include <map>

constexpr int ATOMTYPE_SOLVENT = 0;

class Forcefield {
public:

	Forcefield(VerbosityLevel vl, const std::string& molecule_dir);

	ForceField_NB getNBForcefield() const {
		return forcefield_nb;
	}
	const ForceField_NB& getNBForcefieldRef() const {
		return forcefield_nb;
	}

	static int atomTypeToIndex(const char& atom) {
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
		throw std::runtime_error("Forcefield failed");
	}

	bool forcefield_loaded = false;

private:
	ForceField_NB forcefield_nb;


	VerbosityLevel vl = SILENT;

	std::vector<NBAtomtype> loadAtomTypes(const SimpleParsedFile& nonbonded_parsed);
	std::map<int, int> loadAtomTypeMap(const SimpleParsedFile& nonbonded_parsed);

	void loadAtomypesIntoForcefield(const std::vector<NBAtomtype>& atomtypes);
};