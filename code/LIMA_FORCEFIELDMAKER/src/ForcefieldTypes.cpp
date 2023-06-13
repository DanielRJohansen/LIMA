#include "LIMA_FORCEFIELDMAKER/include/ForcefieldTypes.h"




using std::cout, std::endl;


vector<NB_Atomtype> NB_Atomtype::filterUnusedTypes(const vector<NB_Atomtype>& forcefield, const vector<string>& active_types, Map& map, LimaLogger& logger, bool print_mappings) {
	vector<NB_Atomtype> filtered_list;

	filtered_list.push_back(forcefield[0]);				// Solvent is not present in topology, so we force it to be added here!
	for (const string& type : active_types) {
		string alias = type;

		while (alias.length() > 0) {
			if (typeIsPresent(forcefield, alias)) {		// If type from conf exists, in forcefield, we follow happy path
				NB_Atomtype record = *findRecord(forcefield, alias);

				if (!typeIsPresent(filtered_list, alias)) {			// If type is not present in filtered list, add it with a simulation-specific type-id.
					record.atnum_local = static_cast<int>(filtered_list.size());
					filtered_list.push_back(record);
				}
				if (type != alias) {									// If type != alias, add it to map
					map.addMap(type, alias);
				}
				break;
			}
			else {												// Type from conf does not exist, so we change the alias and try again
				alias.pop_back();
				
				if (alias.length() == 0) {
					throw std::exception("No alias found");
				}
			}
		}

	}
	
	logger.print(std::format("Filtered ff down to {} entries ({} bytes)\n", 
		filtered_list.size(), sizeof(float) * 3 * filtered_list.size()));

	if (print_mappings) {
		logger.print(std::format("Aliases ({}): \t", map.n_mappings));
		for (int i = 0; i < map.n_mappings; i++) {
			logger.print(std::format("\t{} -> {}\n", map.mappings[i].left, map.mappings[i].right));
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

	AtomTable atomtable;	// I guess maps are not ideal, since we always get increasing gro_ids, meaning it is always worst case..


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

			atomtable.insert(std::pair(gro_id, Atom(gro_id, bonded_type, nbonded_type)));
		}

		if (current_state == FINISHED)
			break;
	}

	if (verbose) {
		printf("%lld atoms found in topology file\n", records.size());
	}
	
	return records;
}

void Atom::assignAtomtypeIDs(vector<Atom>& atoms, const vector<NB_Atomtype>& forcefield, const Map& map) {
	for (Atom& atom : atoms) {
		const string alias = map.mapRight(atom.atomtype);

		for (const NB_Atomtype& force_parameters : forcefield) {
			if (force_parameters.type == alias) {
				atom.atomtype_id = force_parameters.atnum_local;
				break;
			}
		}
	}
}





void Singlebondtype::sort() {
	if (!FTHelpers::isSorted(&bonded_typenames[0], &bonded_typenames[1])) {
		swap(bonded_typenames[0], bonded_typenames[1]);
		//std::swap(id1, id2);	// TODO: Should this not be like this???
	}
}

void Anglebondtype::sort() {
	// The middle type must always be in the middle, so we only sort the outer 2 types
	// No need to sort if our edges are the same
	if (bonded_typenames[0] != bonded_typenames[2]) {
		if (!FTHelpers::isSorted(&bonded_typenames[0], &bonded_typenames[2])) {
			swap(bonded_typenames[0], bonded_typenames[2]);
		}
	}
}

void Dihedralbondtype::sort() {
	if (bonded_typenames[0] != bonded_typenames[3]) {
		if (!FTHelpers::isSorted(&bonded_typenames[0], &bonded_typenames[3])) {
			flip();
		}
	}
	else {			// In case the outer two is identical, we check the inner two.
		if (!FTHelpers::isSorted(&bonded_typenames[1], &bonded_typenames[2])) {
			flip();
		}
	}
}
