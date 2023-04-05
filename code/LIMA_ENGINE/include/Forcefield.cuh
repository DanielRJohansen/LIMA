#pragma once

#include "Bodies.cuh"
//#include "Simulation.cuh"
#include "Filehandling.h"
#include <string>
#include <vector>

using std::string;
using std::vector;




#define ATOMTYPE_SOL 0







class Forcefield {
public:
	Forcefield(VerbosityLevel vl);


	void loadForcefield(std::string molecule_dir);
	int getAtomtypeID(int global_id) const;
	SingleBond* getBondType(std::array<int, 2> ids);
	AngleBond* getAngleType(std::array<int, 3> ids);
	DihedralBond* getDihedralType(std::array<int, 4> ids);



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

private:
	ForceField_NB forcefield;


	int* nb_atomtype_ids = nullptr;
	int n_atoms = 0;

	NBAtomtype* nb_atomtypes = nullptr;
	int n_nb_atomtypes = 0;

	std::vector<SingleBond> topol_bonds;
	std::vector<AngleBond> topol_angles;
	std::vector<DihedralBond> topol_dihedrals;

	VerbosityLevel vl = SILENT;


	enum STATE { INACTIVE, FF_NONBONDED, NB_ATOMTYPES, BONDS, ANGLES, DIHEDRALS };
	static STATE setState(string s, STATE current_state) {
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
		return current_state;
	}

	bool newParseTitle(vector<string> row) {
		return (row[0][0] == '#');
	}

	//StringMap parseNBAtomtypeMaps(vector<vector<string>> forcefield_rows) {}

	

	
	NBAtomtype* parseAtomTypes(vector<vector<string>> summary_rows);

	int* parseAtomTypeIDs(vector<vector<string>> forcefield_rows);

	std::vector<SingleBond> parseBonds(vector<vector<string>> forcefield_rows);
	std::vector<AngleBond> parseAngles(vector<vector<string>> forcefield_rows);
	std::vector<DihedralBond> parseDihedrals(vector<vector<string>> forcefield_rows);
	
	void loadAtomypesIntoForcefield();

};






