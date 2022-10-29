#pragma once

#include "Bodies.cuh"
//#include "Simulation.cuh"
#include "Filehandling.h"
#include <string>
#include <vector>

using std::string;


const int MAX_ATOM_TYPES = 32;


#define ATOMTYPE_SOL 0


struct ParticleParameters {	//Nonbonded
	float mass = -1;		//[kg/mol]	or 
	float sigma = -1;
	float epsilon = -1;		// J/mol [kg*nm^2 / s^2]
};


struct ForceField {
	ParticleParameters particle_parameters[MAX_ATOM_TYPES];
};



class ForceFieldMaker {
public:
	ForceFieldMaker() {}


	void buildForcefield();
	int getAtomtypeID(int global_id);
	PairBond* getBondType(int id1, int id2);
	AngleBond* getAngleType(int id1, int id2, int id3);
	DihedralBond* getDihedralType(int id1, int id2, int id3, int id4);



	ForceField getNBForcefield() {
		return forcefield;
	}


	int atomTypeToIndex(char atom) {
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
	ForceField forcefield;


	int* nb_atomtype_ids = nullptr;
	int n_atoms = 0;

	NBAtomtype* nb_atomtypes = nullptr;
	int n_nb_atomtypes = 0;

	PairBond* topol_bonds = nullptr;
	int n_topol_bonds = 0;
	AngleBond* topol_angles = nullptr;
	int n_topol_angles = 0;
	DihedralBond* topol_dihedrals = nullptr;
	int n_topol_dihedrals = 0; 

#ifdef __linux__
	string ff_dir = "../../Simulation/Forcefield/";
#else
	//string ff_dir = "C:\\Users\\Daniel\\git_repo\\Quantom\\";
	string sim_dir = "C:\\PROJECTS\\Quantom\\Simulation\\";
	string ff_dir = "C:\\PROJECTS\\Quantom\\Simulation\\Forcefield\\";
#endif

	//enum STATE { INACTIVE, NB_ATOMTYPE_MAPPINGS, FF_NONBONDED, NB_ATOMTYPES, PAIRTYPES };
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

	PairBond* parseBonds(vector<vector<string>> forcefield_rows);
	AngleBond* parseAngles(vector<vector<string>> forcefield_rows);
	DihedralBond* parseDihedrals(vector<vector<string>> forcefield_rows);
	
	void loadAtomypesIntoForcefield();

};






