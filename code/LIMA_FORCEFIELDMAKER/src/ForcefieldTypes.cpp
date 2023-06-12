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
				throw std::exception("No alias found");
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
















Atom::STATE Atom::setState(string s, STATE current_state) {
	if (s == "atoms")
		return ATOMS;
	if (s == "bonds")
		return FINISHED;
	return INACTIVE;
}

vector<Atom> Atom::parseTopolAtoms(vector<vector<string>>& rows, bool verbose) {
	STATE current_state = INACTIVE;

	vector<Atom> records;

	AtomMap atommap;	// I guess maps are not ideal, since we always get increasing gro_ids, meaning it is always worst case..


	for (vector<string> row : rows) {

		if (row.size() == 3) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == ATOMS) {
			const int gro_id = stoi(row[0]);
			const string bonded_type = row[1];
			const string nbonded_type = row[4];

			records.push_back(Atom(gro_id, bonded_type, nbonded_type));

			atommap.insert(std::pair(gro_id, Atom(gro_id, bonded_type, nbonded_type)));
		}

		if (current_state == FINISHED)
			break;
	}

	if (verbose) {
		printf("%lld atoms found in topology file\n", records.size());
	}
	
	return records;
}

bool Atom::assignAtomtypeID(Atom& atom, vector<NB_Atomtype>& forcefield, const string& alias) {
	for (NB_Atomtype force_parameters : forcefield) {
		if (force_parameters.type == alias) {
			atom.atomtype_id = force_parameters.atnum_local;
			return true;
		}
	}

	return false;
}

void Atom::assignAtomtypeIDs(vector<Atom>* atoms, vector<NB_Atomtype>* forcefield, Map* map) {
	for (int i = 0; i < atoms->size(); i++) {
		Atom* atom = &atoms->at(i);
		string alias = map->mapRight(atom->atomtype);

		bool success = assignAtomtypeID(*atom, *forcefield, alias);
	}
}




void Singlebondtype::assignTypesFromAtomIDs(vector<Singlebondtype>* topol_bonds, vector<Atom> atoms) {
	for (int i = 0; i < topol_bonds->size(); i++) {
		Singlebondtype* bond = &topol_bonds->at(i);


		bond->bonded_typenames[0] = atoms.at(bond->gro_ids[0] - size_t{ 1 }).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
		bond->bonded_typenames[1] = atoms.at(bond->gro_ids[1] - size_t{ 1 }).atomtype_bond;
		bond->sort();
		//cout << bond->type1 << '\t' << bond->type2 << endl;;
	}
}



Singlebondtype* Singlebondtype::findBestMatchInForcefield(Singlebondtype* query_type, vector<Singlebondtype>* forcefield) {
	if (forcefield->size() == 0) { throw std::exception("No angletypes in forcefield!"); }
	float best_likeness = 0;
	Singlebondtype* best_bond = &((*forcefield)[0]);
	for (Singlebondtype& bond : *forcefield) {
		float likeness = FTHelpers::calcLikeness(query_type->bonded_typenames[0], bond.bonded_typenames[0]) * FTHelpers::calcLikeness(query_type->bonded_typenames[1], bond.bonded_typenames[1]);
		if (likeness > best_likeness) {
			best_likeness = likeness;
			best_bond = &bond;
		}
	}
	if (best_likeness > 0.01f)
		return best_bond;

	cout << "Failed to match bond types.\n Closest match " << best_bond->bonded_typenames[0] << "    " << best_bond->bonded_typenames[1];
	printf("Likeness %f\n", best_likeness);
	printf("Topol ids: %d %d\n", query_type->gro_ids[0], query_type->gro_ids[1]);
	cout << query_type->bonded_typenames[0] << '\t' << query_type->bonded_typenames[1] << endl;
	exit(0);
}


void Anglebondtype::assignTypesFromAtomIDs(vector<Anglebondtype>* topol_angles, vector<Atom> atoms) {
	for (int i = 0; i < topol_angles->size(); i++) {
		Anglebondtype* angle = &topol_angles->at(i);

		angle->bonded_typenames[0] = atoms.at(angle->gro_ids[0] - size_t{ 1 }).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
		angle->bonded_typenames[1] = atoms.at(angle->gro_ids[1] - size_t{ 1 }).atomtype_bond;
		angle->bonded_typenames[2] = atoms.at(angle->gro_ids[2] - size_t{ 1 }).atomtype_bond;
		angle->sort();
		//cout << bond->type1 << '\t' << bond->type2 << endl;;
	}
}



Anglebondtype* Anglebondtype::findBestMatchInForcefield(Anglebondtype* query_type, vector<Anglebondtype>* FF_angletypes) {
	if (FF_angletypes->size() == 0) { throw std::exception("No angletypes in forcefield!"); }
	float best_likeness = 0;
	Anglebondtype* best_angle = &((*FF_angletypes)[0]);
	for (Anglebondtype& angle : *FF_angletypes) {
		float likeness = FTHelpers::calcLikeness(query_type->bonded_typenames[0], angle.bonded_typenames[0]) 
			* FTHelpers::calcLikeness(query_type->bonded_typenames[1], angle.bonded_typenames[1]) 
			* FTHelpers::calcLikeness(query_type->bonded_typenames[2], angle.bonded_typenames[2]);

		if (likeness > best_likeness) {
			best_likeness = likeness;
			best_angle = &angle;
		}
		if (likeness == 1)
			break;
	}
	if (best_likeness > 0.01f)
		return best_angle;

	cout << "\n\n\nFailed to match angle types.\n Closest match " << best_angle->bonded_typenames[0] << "    " << best_angle->bonded_typenames[1] << "    " << best_angle->bonded_typenames[2] << endl;
	printf("Likeness %f\n", best_likeness);
	printf("Topol ids: %d %d %d\n", query_type->gro_ids[0], query_type->gro_ids[1], query_type->gro_ids[2]);
	cout << query_type->bonded_typenames[0] << '\t' << query_type->bonded_typenames[1] << '\t' << query_type->bonded_typenames[2] << endl;
	exit(0);
}






















vector<Dihedralbondtype> Dihedralbondtype::parseFFDihedraltypes(vector<vector<string>> rows, bool verbose) {
	FTHelpers::STATE current_state = FTHelpers::INACTIVE;
	vector<Dihedralbondtype> dihedraltypes;

	for (vector<string> row : rows) {

		if (row.size() == 2) {
			current_state = FTHelpers::setState(row[1], current_state);
			continue;
		}

		switch (current_state) {
		case FTHelpers::FF_DIHEDRALTYPES:
			dihedraltypes.push_back(Dihedralbondtype(row[0], row[1], row[2], row[3], stof(row[4]), stof(row[5]), stoi(row[6])));
			break;
		default:
			break;
		}

	}
	if (verbose) {
		printf("%zd dihedraltypes in forcefield\n", dihedraltypes.size());
	}
	
	return dihedraltypes;
}

vector<Dihedralbondtype> Dihedralbondtype::parseTopolDihedraltypes(vector<vector<string>> rows, bool verbose) {
	FTHelpers::STATE current_state = FTHelpers::INACTIVE;
	vector<Dihedralbondtype> records;

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
			records.push_back(Dihedralbondtype(stoi(row[0]), stoi(row[1]), stoi(row[2]), stoi(row[3])));
			entered_dihedrals_first_time = true;
			break;
		default:
			break;
		}
	}

	if (verbose) {
		printf("%zd dihedrals found in topology file\n", records.size());
	}
	
	return records;
}

void Dihedralbondtype::assignTypesFromAtomIDs(vector<Dihedralbondtype>* topol_dihedrals, vector<Atom> atoms) {
	for (int i = 0; i < topol_dihedrals->size(); i++) {
		Dihedralbondtype* dihedral = &topol_dihedrals->at(i);

		//printf("Accessing atoms %d %d %d %d\n", dihedral->id1, dihedral->id2, dihedral->id3, dihedral->id4);
		dihedral->type1 = atoms.at(static_cast<size_t>(dihedral->id1) - 1).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
		dihedral->type2 = atoms.at(static_cast<size_t>(dihedral->id2) - 1).atomtype_bond;
		dihedral->type3 = atoms.at(static_cast<size_t>(dihedral->id3) - 1).atomtype_bond;
		dihedral->type4 = atoms.at(static_cast<size_t>(dihedral->id4) - 1).atomtype_bond;
		dihedral->sort();
		//cout << bond->type1 << '\t' << bond->type2 << endl;;
	}
}

Dihedralbondtype* Dihedralbondtype::findBestMatchInForcefield(Dihedralbondtype* query_type, vector<Dihedralbondtype>* forcefield) {
	if (forcefield->size() == 0) { throw std::exception("No elements in forcefield"); }
	float best_likeness = 0;
	Dihedralbondtype* best_match = &((*forcefield)[0]);
	for (Dihedralbondtype& dihedral : *forcefield) {
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
