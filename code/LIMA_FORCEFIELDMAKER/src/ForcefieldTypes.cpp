#include "LIMA_FORCEFIELDMAKER/include/ForcefieldTypes.h"




using std::cout, std::endl;


vector<NB_Atomtype> NB_Atomtype::filterUnusedTypes(vector<NB_Atomtype> forcefield, vector<string> active_types, Map* map, LimaLogger& logger, bool print_mappings) {
	vector<NB_Atomtype> filtered_list;

	filtered_list.push_back(forcefield[0]);				// Solvent is not present in topology, so we force it to be added here!
	for (string type : active_types) {
		bool found = false;
		string alias = type;



	TRY_AGAIN:


		if (typeIsPresent(&forcefield, alias)) {		// If type from conf exists, in forcefield, we follow happy path
			NB_Atomtype record = findRecord(&forcefield, alias);

			if (!typeIsPresent(&filtered_list, alias)) {			// If type is not present in filtered list, add it with a simulation-specific type-id.
				record.atnum_local = static_cast<int>(filtered_list.size());
				filtered_list.push_back(record);
			}
			if (type != alias) {									// If type != alias, add it to map
				map->addMap(type, alias);
			}
		}
		else {												// Type from conf does not exist, so we change the alias and try again
			if (alias.length() > 1) {
				alias.pop_back();
				goto TRY_AGAIN;
			}
			else {
				printf("Failed");
				exit(404);
			}
		}

	}
	
	logger.print(std::format("Filtered ff down to {} entries ({} bytes)\n", 
		filtered_list.size(), sizeof(float) * 3 * filtered_list.size()));

	if (print_mappings) {
		logger.print(std::format("Aliases ({}): \t", map->n_mappings));
		for (int i = 0; i < map->n_mappings; i++) {
			logger.print(std::format("\t{} -> {}\n", map->mappings[i].left, map->mappings[i].right));
		}
		logger.print("\n");
	}

	return filtered_list;
}


















vector<Atom> Atom::parseTopolAtoms(vector<vector<string>>& rows) {
	STATE current_state = INACTIVE;

	vector<Atom> records;

	for (vector<string> row : rows) {

		if (row.size() == 3) {
			current_state = setState(row[1], current_state);
			continue;
		}



		switch (current_state)
		{
		case INACTIVE:
			break;
		case ATOMS:
			records.push_back(Atom(stoi(row[0]), row[1], row[4]));
			break;
		default:
			break;
		}

		if (current_state == FINISHED)
			break;
	}
	printf("%lld atoms found in topology file\n", records.size());
	return records;
}















vector<Bondtype> Bondtype::parseFFBondtypes(vector<vector<string>> rows) {
	FTHelpers::STATE current_state = FTHelpers::INACTIVE;

	vector<Bondtype> records;

	for (vector<string> row : rows) {
		if (row.size() == 2) {
			current_state = FTHelpers::setState(row[1], current_state);
			continue;
		}


		switch (current_state) {
		case FTHelpers::FF_BONDTYPES:

			records.push_back(Bondtype(row[0], row[1], stof(row[2]), stof(row[3])));
			break;
		default:
			break;
		}
	}
	printf("%lld single bonds read from file\n", records.size());
	return records;
}

vector<Bondtype> Bondtype::parseTopolBondtypes(vector<vector<string>> rows) {
	FTHelpers::STATE current_state = FTHelpers::INACTIVE;
	vector<Bondtype> records;

	for (vector<string> row : rows) {
		if (row.size() == 3) {
			if (row[0][0] == '[') {
				current_state = FTHelpers::setState(row[1], current_state);
				continue;
			}
		}



		switch (current_state) {
		case FTHelpers::FF_BONDTYPES:
			records.push_back(Bondtype(stoi(row[0]), stoi(row[1])));
			break;
		default:
			break;
		}
	}
	printf("%lld bonds found in topology file\n", records.size());
	return records;
}

void Bondtype::assignTypesFromAtomIDs(vector<Bondtype>* topol_bonds, vector<Atom> atoms) {
	for (int i = 0; i < topol_bonds->size(); i++) {
		Bondtype* bond = &topol_bonds->at(i);


		bond->type1 = atoms.at(bond->id1 - size_t{ 1 }).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
		bond->type2 = atoms.at(bond->id2 - size_t{ 1 }).atomtype_bond;
		bond->sort();
		//cout << bond->type1 << '\t' << bond->type2 << endl;;
	}
}



Bondtype* Bondtype::findBestMatchInForcefield(Bondtype* query_type, vector<Bondtype>* forcefield) {
	if (forcefield->size() == 0) { throw std::exception("No angletypes in forcefield!"); }
	float best_likeness = 0;
	Bondtype* best_bond = &((*forcefield)[0]);
	for (Bondtype& bond : *forcefield) {
		float likeness = FTHelpers::calcLikeness(query_type->type1, bond.type1) * FTHelpers::calcLikeness(query_type->type2, bond.type2);
		if (likeness > best_likeness) {
			best_likeness = likeness;
			best_bond = &bond;
		}
	}
	if (best_likeness > 0.01f)
		return best_bond;

	cout << "Failed to match bond types.\n Closest match " << best_bond->type1 << "    " << best_bond->type2;
	printf("Likeness %f\n", best_likeness);
	printf("Topol ids: %d %d\n", query_type->id1, query_type->id2);
	cout << query_type->type1 << '\t' << query_type->type2 << endl;
	exit(0);
}





vector<Angletype> Angletype::parseFFAngletypes(vector<vector<string>> rows) {
	FTHelpers::STATE current_state = FTHelpers::INACTIVE;
	vector<Angletype> angletypes;

	for (vector<string> row : rows) {

		if (row.size() == 2) {
			current_state = FTHelpers::setState(row[1], current_state);
			continue;
		}



		switch (current_state) {
		case FTHelpers::FF_ANGLETYPES:
			angletypes.push_back(Angletype(row[0], row[1], row[2], stof(row[3]), stof(row[4])));
			break;
		default:
			break;
		}


		if (current_state == FTHelpers::FF_DIHEDRALTYPES)
			break;
	}
	printf("%lld angletypes in forcefield\n", angletypes.size());
	return angletypes;
}

vector<Angletype> Angletype::parseTopolAngletypes(vector<vector<string>> rows) {
	FTHelpers::STATE current_state = FTHelpers::INACTIVE;
	vector<Angletype> records;

	for (vector<string> row : rows) {
		if (row.size() == 3) {
			if (row[0][0] == '[') {
				current_state = FTHelpers::setState(row[1], current_state);
				continue;
			}
		}



		switch (current_state)
		{
		case FTHelpers::FF_ANGLETYPES:
			records.push_back(Angletype(stoi(row[0]), stoi(row[1]), stoi(row[2])));
			break;
		default:
			break;
		}


		if (current_state == FTHelpers::FF_DIHEDRALTYPES)
			break;
	}
	printf("%lld angles found in topology file\n", records.size());
	return records;
}

void Angletype::assignTypesFromAtomIDs(vector<Angletype>* topol_angles, vector<Atom> atoms) {
	for (int i = 0; i < topol_angles->size(); i++) {
		Angletype* angle = &topol_angles->at(i);

		angle->type1 = atoms.at(angle->id1 - size_t{1}).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
		angle->type2 = atoms.at(angle->id2 - size_t{1}).atomtype_bond;
		angle->type3 = atoms.at(angle->id3 - size_t{1}).atomtype_bond;
		angle->sort();
		//cout << bond->type1 << '\t' << bond->type2 << endl;;
	}
}



Angletype* Angletype::findBestMatchInForcefield(Angletype* query_type, vector<Angletype>* FF_angletypes) {
	if (FF_angletypes->size() == 0) { throw std::exception("No angletypes in forcefield!"); }
	float best_likeness = 0;
	Angletype* best_angle = &((*FF_angletypes)[0]);
	for (Angletype& angle : *FF_angletypes) {
		float likeness = FTHelpers::calcLikeness(query_type->type1, angle.type1) * FTHelpers::calcLikeness(query_type->type2, angle.type2) * FTHelpers::calcLikeness(query_type->type3, angle.type3);
		if (likeness > best_likeness) {
			best_likeness = likeness;
			best_angle = &angle;
		}
		if (likeness == 1)
			break;
	}
	if (best_likeness > 0.01f)
		return best_angle;

	cout << "\n\n\nFailed to match angle types.\n Closest match " << best_angle->type1 << "    " << best_angle->type2 << "    " << best_angle->type3 << endl;
	printf("Likeness %f\n", best_likeness);
	printf("Topol ids: %d %d %d\n", query_type->id1, query_type->id2, query_type->id3);
	cout << query_type->type1 << '\t' << query_type->type2 << '\t' << query_type->type3 << endl;
	exit(0);
}






















vector<Dihedraltype> Dihedraltype::parseFFDihedraltypes(vector<vector<string>> rows) {
	FTHelpers::STATE current_state = FTHelpers::INACTIVE;
	vector<Dihedraltype> dihedraltypes;

	for (vector<string> row : rows) {

		if (row.size() == 2) {
			current_state = FTHelpers::setState(row[1], current_state);
			continue;
		}

		switch (current_state) {
		case FTHelpers::FF_DIHEDRALTYPES:
			dihedraltypes.push_back(Dihedraltype(row[0], row[1], row[2], row[3], stof(row[4]), stof(row[5]), stoi(row[6])));
			break;
		default:
			break;
		}

	}
	printf("%zd dihedraltypes in forcefield\n", dihedraltypes.size());
	return dihedraltypes;
}

vector<Dihedraltype> Dihedraltype::parseTopolDihedraltypes(vector<vector<string>> rows) {
	FTHelpers::STATE current_state = FTHelpers::INACTIVE;
	vector<Dihedraltype> records;

	bool entered_dihedrals_first_time = false;	// The top contains maybe improper dihedrals? ANyways these mess up, as they are simply named as dihedrals

	for (vector<string> row : rows) {
		if (row.size() == 3) {
			if (row[0][0] == '[') {
				current_state = FTHelpers::setState(row[1], current_state);
				if (entered_dihedrals_first_time)
					break;
				continue;
			}
		}



		switch (current_state) {
		case FTHelpers::FF_DIHEDRALTYPES:
			if (row.size() != 5)
				continue;
			records.push_back(Dihedraltype(stoi(row[0]), stoi(row[1]), stoi(row[2]), stoi(row[3])));
			entered_dihedrals_first_time = true;
			break;
		default:
			break;
		}
	}
	printf("%zd dihedrals found in topology file\n", records.size());
	return records;
}

void Dihedraltype::assignTypesFromAtomIDs(vector<Dihedraltype>* topol_dihedrals, vector<Atom> atoms) {
	for (int i = 0; i < topol_dihedrals->size(); i++) {
		Dihedraltype* dihedral = &topol_dihedrals->at(i);

		//printf("Accessing atoms %d %d %d %d\n", dihedral->id1, dihedral->id2, dihedral->id3, dihedral->id4);
		dihedral->type1 = atoms.at(static_cast<size_t>(dihedral->id1) - 1).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
		dihedral->type2 = atoms.at(static_cast<size_t>(dihedral->id2) - 1).atomtype_bond;
		dihedral->type3 = atoms.at(static_cast<size_t>(dihedral->id3) - 1).atomtype_bond;
		dihedral->type4 = atoms.at(static_cast<size_t>(dihedral->id4) - 1).atomtype_bond;
		dihedral->sort();
		//cout << bond->type1 << '\t' << bond->type2 << endl;;
	}
}

Dihedraltype* Dihedraltype::findBestMatchInForcefield(Dihedraltype* query_type, vector<Dihedraltype>* forcefield) {
	if (forcefield->size() == 0) { throw std::exception("No elements in forcefield"); }
	float best_likeness = 0;
	Dihedraltype* best_match = &((*forcefield)[0]);
	for (Dihedraltype& dihedral : *forcefield) {
		float likeness = FTHelpers::calcLikeness(query_type->type1, dihedral.type1)
			* FTHelpers::calcLikeness(query_type->type2, dihedral.type2)
			* FTHelpers::calcLikeness(query_type->type3, dihedral.type3)
			* FTHelpers::calcLikeness(query_type->type4, dihedral.type4);
		if (likeness > best_likeness) {
			best_likeness = likeness;
			best_match = &dihedral;
		}
		if (likeness == 1)
			break;
	}
	if (best_likeness > 0.001f)
		return best_match;

	cout << "\n\n\nFailed to match angle types.\nClosest match: " << best_match->type1 << "    " << best_match->type2 << "    " << best_match->type3 << "    " << best_match->type4 << endl;
	printf("Likeness %f\n", best_likeness);
	printf("Topol ids: %d %d %d %d\n", query_type->id1, query_type->id2, query_type->id3, query_type->id4);
	cout << query_type->type1 << '\t' << query_type->type2 << '\t' << query_type->type3 << '\t' << query_type->type4 << endl;
	exit(0);
}
