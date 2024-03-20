#pragma once

#include "Filehandling.h"
#include "Bodies.cuh"

#include <string>
#include <map>

constexpr int ATOMTYPE_SOLVENT = 0;

class Forcefield {
public:

	Forcefield(VerbosityLevel vl, const std::string& molecule_dir);

	ForceField_NB constexpr getNBForcefield() const {
		return forcefield_nb;
	}
	const ForceField_NB& getNBForcefieldRef() const {
		return forcefield_nb;
	}

	bool forcefield_loaded = false;

private:
	ForceField_NB forcefield_nb;


	VerbosityLevel vl = SILENT;

	std::vector<NBAtomtype> loadAtomTypes(const SimpleParsedFile& nonbonded_parsed);
	std::map<int, int> loadAtomTypeMap(const SimpleParsedFile& nonbonded_parsed);

	void loadAtomypesIntoForcefield(const std::vector<NBAtomtype>& atomtypes);
};