#pragma once

#include "LIMA_ENGINE/include/Forcefield.cuh"
#include "LIMA_BASE/include/Printer.h"
#include "LIMA_ENGINE/include/EngineUtils.cuh"

using namespace LIMA_Print;

const int min_reserve_size = 10000;	// This should NOT be permanent...

Forcefield::Forcefield(VerbosityLevel vl) : vl(vl) {};


void Forcefield::loadForcefield(string molecule_dir) {
	if (vl >= CRITICAL_INFO) { printH2("Building forcefield"); }

	vector<vector<string>> summary_rows = Filehandler::readFile(molecule_dir + "/LIMA_ffnonbonded_filtered.txt");
	vector<vector<string>> forcefield_rows = Filehandler::readFile(molecule_dir + "/LIMA_ffbonded_filtered.txt");


	nb_atomtypes = parseAtomTypes(summary_rows);					// 1 entry per type in compressed forcefield
	loadAtomypesIntoForcefield();


	groIdToAtomtypeMap = parseAtomTypeIDs(forcefield_rows);	// 1 entry per atom in conf

	topol_bonds = parseBonds(forcefield_rows);
	topol_angles = parseAngles(forcefield_rows);
	topol_dihedrals = parseDihedrals(forcefield_rows);
	topol_improperdihedrals = parseImproperDihedrals(forcefield_rows);

	forcefield_loaded = true;

	if (vl >= CRITICAL_INFO) {
		printf("Nonbonded parameters size: %llu bytes\n", sizeof(ForceField_NB));
		printH2("Finished building forcefield");
	}
}


int Forcefield::getAtomtypeID(int gro_id) const {
	if (groIdToAtomtypeMap.count(gro_id) == 0 || gro_id == 0) {	// 0 is an error, as atoms are 1-indexed
		printf("Attempting to fetch atomtype of non-loaded atom with global_id %d\n", gro_id);
		exit(0);
	}
	return groIdToAtomtypeMap.find(gro_id)->second;
}

template <int array_length>
bool isMatch(const uint32_t* topolbonds, const std::array<int, array_length> query_ids) {
	for (int i = 0; i < array_length; i++) {
		if (topolbonds[i] != query_ids[i]) { return false; }
	}
	return true;
}

const SingleBond* Forcefield::getBondType(std::array<int, 2> ids) const {
	for (auto& singlebond : topol_bonds) {
		if (isMatch(singlebond.atom_indexes, ids)) {
			return &singlebond;
		}
	}
	printf("Bond not found with ids %d %d\n", ids[0], ids[1]);
	exit(0);
}

const AngleBond* Forcefield::getAngleType(std::array<int, 3> ids) const {
	for (auto& anglebond : topol_angles) {
		if (isMatch(anglebond.atom_indexes, ids)) {
			return &anglebond;
		}
	}

	printf("Angle not found with ids %d %d %d\n", ids[0], ids[1], ids[2]);
	exit(0);
}

const DihedralBond* Forcefield::getDihedralType(std::array<int, 4> ids) const {
	for (auto& dihedralbond : topol_dihedrals) {
		if (isMatch(dihedralbond.atom_indexes, ids)) {
			return &dihedralbond;
		}
	}

	printf("Dihedral not found with ids %d %d %d %d\n", ids[0], ids[1], ids[2], ids[3]);
	exit(0);
}

std::vector<NBAtomtype> Forcefield::parseAtomTypes(vector<vector<string>> summary_rows) {
	nb_atomtypes.reserve(10000);
	STATE current_state = INACTIVE;

	for (vector<string> row : summary_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		

		if (current_state == FF_NONBONDED) {
			//for (string e : row)
				//cout << e << '\t';
			//printf("\n");
			// Row is type, id, weight [g], sigma [nm], epsilon [J/mol]

			//atomtypes[ptr++] = NBAtomtype(stof(row[2]), stof(row[3]), stof(row[4]));
			nb_atomtypes.push_back(NBAtomtype{ stof(row[2]), stof(row[3]), stof(row[4]) });
		}			
	}
	n_nb_atomtypes = nb_atomtypes.size();
	if (vl >= V1) { printf("%d NB_Atomtypes loaded\n", n_nb_atomtypes); }
	return nb_atomtypes;
}

std::map<int, int> Forcefield::parseAtomTypeIDs(vector<vector<string>> forcefield_rows) {	// returns the nonbonded atomtype
	std::map<int, int> groidToType;

	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == NB_ATOMTYPES) {
			int gro_id = stoi(row[0]);
			int atomtype_id = stoi(row[1]);


			groidToType.insert(std::pair<int, int>(gro_id, atomtype_id));
		}
			
	}
	if (vl >= V1) { printf("%d NB_Atomtype_IDs loaded\n", groidToType.size()); }

	return groidToType;
}

std::vector<SingleBond> Forcefield::parseBonds(vector<vector<string>> forcefield_rows) {
	std::vector<SingleBond> singlebonds;
	singlebonds.reserve(min_reserve_size);
	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == BONDS) {
			singlebonds.push_back(SingleBond(
				stoi(row[0]), 
				stoi(row[1]), 
				stof(row[4]) * NANO_TO_LIMA,						// convert [nm] to [lm]*/
				stof(row[5]) / (NANO_TO_LIMA * NANO_TO_LIMA)		// convert [J/(mol * nm^2)] to [J/(mol * nm * lm) I dont know why, but one of the "nm" might be the direction unitvector?? i am confused...
			));		
		}
	}
	if (vl >= V1) { printf("%d bonds loaded\n", singlebonds.size()); }
	return singlebonds;
}

std::vector<AngleBond> Forcefield::parseAngles(vector<vector<string>> forcefield_rows) {
	std::vector<AngleBond> anglebonds;
	anglebonds.reserve(min_reserve_size);

	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == ANGLES) {
			anglebonds.push_back(AngleBond(stoi(row[0]), stoi(row[1]), stoi(row[2]), stof(row[6]) , stof(row[7])));		// Assumes radians here
		}

	}
	if (vl >= V1) { printf("%d angles loaded\n", anglebonds.size()); }
	return anglebonds;
}

std::vector<DihedralBond> Forcefield::parseDihedrals(vector<vector<string>> forcefield_rows) {
	std::vector<DihedralBond> dihedralbonds;
	dihedralbonds.reserve(min_reserve_size);

	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			
		//	if (has_been_enabled)	// To deal with the wierd dihedrals at the bottom of the topol.top
			//	break;
			continue;
		}

		if (current_state == DIHEDRALS) {
			dihedralbonds.push_back(DihedralBond(
				stoi(row[0]), 
				stoi(row[1]), 
				stoi(row[2]), 
				stoi(row[3]), 
				stof(row[8]),
				abs(stof(row[9])), // MIGHT HAVE TO DO AN ABS() ON K_PHI, SINCE IT IS NEGATIVE SOMETIMES??? WHAT THE FUCKKKKKKKKKK CHEMISTS?????!?!?!
				stoi(row[10])
			));			
			//has_been_enabled = true;
		}
	}
	if (vl >= V1) { printf("%d dihedrals loaded\n", dihedralbonds.size()); }
	return dihedralbonds;
}

std::vector<ImproperDihedralBond> Forcefield::parseImproperDihedrals(const vector<vector<string>>& forcefield_rows) {
	std::vector<ImproperDihedralBond> improperdihedrals;
	improperdihedrals.reserve(min_reserve_size);

	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == IMPROPERDIHEDRALS) {
			std::array<uint32_t,4> ids{ stoul(row[0]), stoi(row[1]), stoi(row[2]), stoi(row[3]) };
			float psi0 = stof(row[8]);
			float kpsi = stof(row[9]);
			improperdihedrals.push_back(ImproperDihedralBond(ids, psi0, kpsi));
		}
	}
	if (vl >= V1) { printf("%d improper dihedrals loaded\n", improperdihedrals.size()); }
	return improperdihedrals;
}

void Forcefield::loadAtomypesIntoForcefield() {
	static const float mass_min = 0.001f;	// [kg/mol]
	static const float sigma_min = 0.001f;
	static const float epsilon_min = 0.001f;

	for (int i = 0; i < n_nb_atomtypes; i++) {
		forcefield.particle_parameters[i].mass = nb_atomtypes[i].mass * 1e-3f;				// Convert g/mol to kg/mol
		forcefield.particle_parameters[i].sigma = nb_atomtypes[i].sigma * NANO_TO_LIMA;		// Convert from [nm] to [lm]
		forcefield.particle_parameters[i].epsilon = nb_atomtypes[i].epsilon;				// Interpreted as kg*lm^2/ls^2 

		bool illegal_parameter = (forcefield.particle_parameters[i].mass < mass_min) || (forcefield.particle_parameters[i].sigma < sigma_min) || (forcefield.particle_parameters[i].epsilon < epsilon_min);

		if ((vl >= V2) || illegal_parameter) { printf("Mass %f Sigma %f Epsilon %f\n", nb_atomtypes[i].mass, nb_atomtypes[i].sigma, nb_atomtypes[i].epsilon); }
	}

}
