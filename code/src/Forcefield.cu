#pragma once


#include "Forcefield.cuh"





void ForceFieldMaker::buildForcefield() {
	printf("#################################### BUILDING FORCEFIELD ####################################\n");

//	vector<vector<string>> summary_rows = Filehandler::readFile("C:\\Users\\Daniel\\git_repo\\Quantom\\ForcefieldSummary.txt", INT_MAX, true);
//	vector<vector<string>> forcefield_rows = Filehandler::readFile("C:\\Users\\Daniel\\git_repo\\Quantom\\Forcefield.txt", INT_MAX, true);

	//vector<vector<string>> summary_rows = Filehandler::readFile(ff_dir + "ForcefieldSummary.txt", INT_MAX, true);
	//vector<vector<string>> forcefield_rows = Filehandler::readFile(ff_dir + "Forcefield.txt", INT_MAX, true);

	vector<vector<string>> summary_rows = Filehandler::readFile(ff_dir + "LIMA_ffnonbonded_filtered.txt", INT_MAX, true);
	vector<vector<string>> forcefield_rows = Filehandler::readFile(ff_dir + "LIMA_ffbonded_filtered.txt", INT_MAX, true);


	nb_atomtypes = parseAtomTypes(summary_rows);					// 1 entry per type in compressed forcefield
	loadAtomypesIntoForcefield();

	
	//printf("%f %f %f\n", forcefield.particle_parameters[1].mass, forcefield.particle_parameters[1].sigma, forcefield.particle_parameters[1].epsilon);
	//printf("%f %f %f\n", forcefield1.particle_parameters[1].mass, forcefield1.particle_parameters[1].sigma, forcefield1.particle_parameters[1].epsilon);

	nb_atomtype_ids = parseAtomTypeIDs(forcefield_rows);				// 1 entry per atom in conf

	topol_bonds = parseBonds(forcefield_rows);
	topol_angles = parseAngles(forcefield_rows);
	topol_dihedrals = parseDihedrals(forcefield_rows);

	printf("Nonbonded parameters size: %llu bytes\n", sizeof(ForceField));

	printf("\n\n############################# FINISHED BUILDING FORCEFIELD #############################\n\n\n");

	//exit(0);
}


int ForceFieldMaker::getAtomtypeID(int global_id) {
	if (global_id > n_atoms || global_id == 0) {	// 0 is an error, as atoms are 1-indexed
		printf("Attempting to fetch atomtype of non-loaded atom with global_id %d\n", global_id);
		exit(0);
	}
	return nb_atomtype_ids[global_id];
}

PairBond* ForceFieldMaker::getBondType(int id1, int id2) {
	for (int i = 0; i < n_topol_bonds; i++) {
		if (topol_bonds[i].atom_indexes[0] == id1 && topol_bonds[i].atom_indexes[1] == id2) {
			return &topol_bonds[i];
		}
	}
	printf("Bond not found with ids %d %d\n", id1, id2);
	exit(0);
}

AngleBond* ForceFieldMaker::getAngleType(int id1, int id2, int id3) {
	for (int i = 0; i < n_topol_angles; i++) {
		if (topol_angles[i].atom_indexes[0] == id1 && topol_angles[i].atom_indexes[1] == id2 && topol_angles[i].atom_indexes[2] == id3) {
			return &topol_angles[i];
		}
	}
	printf("Angle not found with ids %d %d %d\n", id1, id2, id3);
	exit(0);
}

DihedralBond* ForceFieldMaker::getDihedralType(int id1, int id2, int id3, int id4) {
	for (int i = 0; i < n_topol_dihedrals; i++) {
		if (topol_dihedrals[i].atom_indexes[0] == id1 && topol_dihedrals[i].atom_indexes[1] == id2 && topol_dihedrals[i].atom_indexes[2] == id3 && topol_dihedrals[i].atom_indexes[3] == id4) {
			return &topol_dihedrals[i];
		}
	}
	printf("Dihedral not found with ids %d %d %d %d\n", id1, id2, id3, id4);
	exit(0);
}












NBAtomtype* ForceFieldMaker::parseAtomTypes(vector<vector<string>> summary_rows) {
	NBAtomtype* atomtypes = new NBAtomtype[10000];
	int ptr = 0;
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
			atomtypes[ptr++] = NBAtomtype(stof(row[2]), stof(row[3]), stof(row[4]));
		}			
	}
	n_nb_atomtypes = ptr;
	printf("%d NB_Atomtypes loaded\n", ptr);
	return atomtypes;
}

int* ForceFieldMaker::parseAtomTypeIDs(vector<vector<string>> forcefield_rows) {	// returns the nonbonded atomtype
	int* atomtype_ids = new int[10000];
	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == NB_ATOMTYPES) {
			atomtype_ids[stoi(row[0])] = stoi(row[1]);
			n_atoms++;
		}
			
	}
	printf("%d NB_Atomtype_IDs loaded\n", n_atoms);
	return atomtype_ids;
}

PairBond* ForceFieldMaker::parseBonds(vector<vector<string>> forcefield_rows) {
	PairBond* bonds = new PairBond[10000];
	int ptr = 0;
	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == BONDS) {
			bonds[ptr++] = PairBond(stoi(row[0]), stoi(row[1]), stof(row[4]), stof(row[5]));
		}
	}
	n_topol_bonds = ptr;
	printf("%d bonds loaded\n", ptr);
	return bonds;
}

AngleBond* ForceFieldMaker::parseAngles(vector<vector<string>> forcefield_rows) {
	AngleBond* angles = new AngleBond[10000];
	int ptr = 0;
	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == ANGLES) {
			angles[ptr++] = AngleBond(stoi(row[0]), stoi(row[1]), stoi(row[2]), stof(row[6]) , stof(row[7]));		// Assumes radians here
		}

	}
	n_topol_angles = ptr;
	printf("%d angles loaded\n", ptr);
	return angles;
}

DihedralBond* ForceFieldMaker::parseDihedrals(vector<vector<string>> forcefield_rows) {
	DihedralBond* dihedrals = new DihedralBond[10000];
	int ptr = 0;
	STATE current_state = INACTIVE;
	//bool has_been_enabled = false;


	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			
		//	if (has_been_enabled)	// To deal with the wierd dihedrals at the bottom of the topol.top
			//	break;
			continue;
		}

		if (current_state == DIHEDRALS) {
			dihedrals[ptr++] = DihedralBond(stoi(row[0]), stoi(row[1]), stoi(row[2]), stoi(row[3]), stof(row[8]), abs(stof(row[9])), stoi(row[10]));			// MIGHT HAVE TO DO AN ABS() ON K_PHI, SINCE IT IS NEGATIVE SOMETIMES??? WHAT THE FUCKKKKKKKKKK CHEMISTS?????!?!?!
			//has_been_enabled = true;
		}
	}
	n_topol_dihedrals = ptr;
	printf("%d dihedrals loaded\n", ptr);
	return dihedrals;
}






void ForceFieldMaker::loadAtomypesIntoForcefield() {
	for (int i = 0; i < n_nb_atomtypes; i++) {
		forcefield.particle_parameters[i].mass = nb_atomtypes[i].mass * 1e-3;		// Convert g/mol to kg/mol
		forcefield.particle_parameters[i].sigma = nb_atomtypes[i].sigma;
		forcefield.particle_parameters[i].epsilon = nb_atomtypes[i].epsilon;
		printf("Mass %f Sigma %f Epsilon %f\n", nb_atomtypes[i].mass, nb_atomtypes[i].sigma, nb_atomtypes[i].epsilon);
	}
}
