#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

///////////////////////////////// READ HERE FIRST /////////////////////////////////
// ffbonded.itp and ffnonbonden.itp has different atomtypes. 
// The purpose if this bit of code is to compress all forcefield parameters so they can fit on a GPU more easily
//

using namespace std;

const float water_mass = 15.999000 + 2 * 1.008000;
const float water_sigma = 0.302905564168 + 2 * 0.040001352445;
const float water_epsilon = (0.50208f + 2.f * 0.19246f) *1000.f;		// Convert kJ->J


struct FTHelpers {
	static bool charIsNumber(char c) {
		return ((int)c > 47 && (int)c < 58);
	}
	static bool charIsNumberAbove1(char c) {
		return ((int)c > 49 && (int)c < 58);
	}
	static float calcLikeness(const string& a,const string& b) {
		float likeness = 0;
		float point_scale = 1.f / max(a.length(), b.length());

		for (int i = 0; i < min(a.length(), b.length()); i++) {
			if (a[i] == b[i])
				likeness += point_scale;
			else
				break;
		}
		return likeness;
	}

	static bool isSorted(string* leftmost, string* rightmost) {
		int ptr = 0;
		while (leftmost->length() > ptr && rightmost->length() > ptr) {
			if ((int)(*leftmost)[ptr] < (int)(*rightmost)[ptr]) {
				return true;
			}
			else if ((int)(*leftmost)[ptr] > (int)(*rightmost)[ptr]) {
				return false;
			}
			ptr++;
		}
		if (leftmost->length() > rightmost->length()) {
			return false;
		}
		return true;
	}

	static vector<string> parseConf(vector<vector<string>> rows) {		// KNOWN ERROR HERE. SOMETIMES SPACES IN FILES ARE TABS, AND IT DOESN'T RECOGNISE A ROW AS A PROPER ENTRY!!
		vector<string> atom_types;

		for (vector<string> row : rows) {
			if (row.size() != 6) {
				continue;
			}

			bool type_already_found = false;
			for (string atom_type : atom_types) {
				if (atom_type == row[1]) {
					type_already_found = true;
					break;
				}
			}
			if (!type_already_found)
				atom_types.push_back(row[1]);
		}

		printf("%d atom types in conf\n", atom_types.size());

		return atom_types;
	}


	enum STATE { INACTIVE, FF_NONBONDED, FF_BONDTYPES, FF_ANGLETYPES, FF_DIHEDRALTYPES, FF_PAIRTYPES };
	static STATE setState(string s, STATE current_state) {
		if (s == "ff_nonbonded")
			return FF_NONBONDED;
		if (s == "ff_bondtypes" || s == "bonds")
			return FF_BONDTYPES;
		if (s == "ff_angletypes" || s == "angles")
			return FF_ANGLETYPES;
		if (s == "ff_dihedraltypes" || s == "dihedrals")
			return FF_DIHEDRALTYPES;
		if (s == "pairs")
			return FF_PAIRTYPES;

		return current_state;
	}
};


struct Map {
	struct Mapping {
		Mapping () {}
		Mapping (string l, string r) : left(l), right(r) {}
		string left, right;
	};



	Map() {
		mappings = new Mapping[10000];
	}
	~Map() {
		delete[] mappings;
	}

	Mapping* mappings;
	int n_mappings = 0; 

	bool mapExists(string l) {
		for (int i = 0; i < n_mappings; i++)
			if (l == mappings[i].left)
				return true;
		return false;
	}
	void addMap(string l, string r) {
		if (!mapExists(l))
			mappings[n_mappings++] = Mapping(l, r);
	}
	string mapRight(const string& query) {
		for (int i = 0; i < n_mappings; i++) {
			if (mappings[i].left == query) {
				return mappings[i].right;
			}

		}
		return query;
		//cout << "No map for query" << query << endl;
		//exit(0);
	}
};


struct NB_Atomtype {
	NB_Atomtype() {}
	NB_Atomtype(string t) : type(t){}		// This is for loading form the conf file, for easy comparisons
	NB_Atomtype(string t, float mass) : type(t), mass(mass) {}		// This if for forcefield merging
	NB_Atomtype(string t, float mass, float sigma, float epsilon) : type(t),  mass(mass), sigma(sigma), epsilon(epsilon) {}		// For LIMA_ffnonbonded
	NB_Atomtype(string t, int atnum, float mass, float sigma, float epsilon) : type(t), atnum(atnum), mass(mass), sigma(sigma), epsilon(epsilon) {}		// For random dudes ff_nonbonded
	// Official parameters
	string type = "";
	int atnum = -1;					// atnum given by input file (CHARMM)
	int atnum_local = 0;			// atnum specific to simulation
	float mass = -1;		// [g/mol]
	float sigma = -1;	// [nm]
	float epsilon = -1;	// J/mol

	// LIMA parameters
	bool is_present_in_simulation = false;
	int simulation_specific_id;





	//bool operator==(const NB_Atomtype a) { return (type == a.type); }
	//static bool sameType(NB_Atomtype a, NB_Atomtype b) { return (a.type == b.type); }


	// Parser functions

	static vector<NB_Atomtype> parseNonbonded(vector<vector<string>> rows) {		// KNOWN ERROR HERE. SOMETIMES SPACES IN FILES ARE TABS, AND IT DOESN'T RECOGNISE A ROW AS A PROPER ENTRY!!
		FTHelpers::STATE current_state = FTHelpers::INACTIVE;

		vector<NB_Atomtype> records;

		records.push_back(NB_Atomtype("WATER", 0, water_mass, water_sigma, water_epsilon));		// Solvent type always first!

		
		for (vector<string> row : rows) {

			if (row.size() == 2) {
				current_state =FTHelpers::setState(row[1], current_state);
				continue;
			}


			switch (current_state) {
			case FTHelpers::FF_NONBONDED:

				records.push_back(NB_Atomtype(
					row[0], 
					stof(row[1]), 
					stof(row[2]), 
					stof(row[3])
				));
				break;
			default:
				break;
			}
		}
		printf("%d atom types read from file\n", records.size());
		return records;
	}
	static NB_Atomtype findRecord(vector<NB_Atomtype>* records, string type) {
		for (NB_Atomtype record : *records) {
			if (record.type == type) {
				return record;
			}
		}
		return NB_Atomtype();
	}
	static bool typeIsPresent(vector<NB_Atomtype>* records, string type) {
		return (findRecord(records, type).type != "");
	}
	static vector<NB_Atomtype> filterUnusedTypes(vector<NB_Atomtype> forcefield, vector<string> active_types, Map* map) {
		vector<NB_Atomtype> filtered_list;

		filtered_list.push_back(forcefield[0]);				// Solvent is not present in topology, so we force it to be added here!
		//string* mappings = new string[10000];
		for (string type : active_types) {
			bool found = false;
			string alias = type;



			TRY_AGAIN:


			if (typeIsPresent(&forcefield, alias)) {		// If type from conf exists, in forcefield, we follow happy path
				NB_Atomtype record = findRecord(&forcefield, alias);

				if (!typeIsPresent(&filtered_list, alias)) {			// If type is not present in filtered list, add it with a simulation-specific type-id.
					record.atnum_local = filtered_list.size();
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

			/*
			NB_Atomtype record = findRecord(&forcefield, alias);

			if (record.atnum != -1) {
				if (!typeIsPresent(&filtered_list, alias)) {
					record.atnum_local = filtered_list.size();
					filtered_list.push_back(record);
				}
				if (type != alias) {
					map->addMap(type, alias);
				}					
			}
			else {
				if (alias.length() > 0) {
					alias.pop_back();
					goto TRY_AGAIN;
				}
				else {
					printf("Failed");
					exit(404);
				}
			}
			*/
		}
		printf("Filtered ff down to %d entries (%d bytes)\n", filtered_list.size(), sizeof(float) * 3 * filtered_list.size());
		printf("Aliases (%d): \n", map->n_mappings);
		for (int i = 0; i < map->n_mappings; i++) {
			//cout << map->mappings[i].left << "    " << map->mappings[i].right << endl;
		}
		printf("\n\n\n");
		return filtered_list;
	}
};


// This is for bonded atoms!!!!!!!!!!!
struct Atom {
	Atom() {}
	Atom(int id, string type_b, string type_nb) : id(id), atomtype_bond(type_b), atomtype(type_nb) {
		//convertToZeroindexed();
	}
	int id;										// Come from topol.top file
	string atomtype;
	string atomtype_bond;
	int atomtype_id = -1;				// Asigned later
	//float charge;



	enum STATE { INACTIVE, ATOMS, FINISHED };
	static STATE setState(string s, STATE current_state) {
		if (s == "atoms")
			return ATOMS;
		if (s == "bonds")
			return FINISHED;
		return current_state;
	}

	static vector<Atom> parseTopolAtoms(vector<vector<string>> rows) {
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
		printf("%d atoms found in topology file\n", records.size());
		return records;
	}


	static bool assignAtomtypeID(Atom& atom, vector<NB_Atomtype>& forcefield, const string& alias) {
		for (NB_Atomtype force_parameters : forcefield) {
			if (force_parameters.type == alias) {
				atom.atomtype_id = force_parameters.atnum_local;
				return true;
			}			
		}
		//printf("Alias not found!\n");	// TODO: FIX. irght now i dont care tho
		//exit(1);
		return false;
	}

	static void assignAtomtypeIDs(vector<Atom>* atoms, vector<NB_Atomtype>* forcefield, Map* map) {
		for (int i = 0; i < atoms->size(); i++) {
			Atom* atom = &((*atoms).at(i));
			string alias = map->mapRight(atom->atomtype);

			bool success = assignAtomtypeID(*atom, *forcefield, alias);
		}
	}
};



struct Bondtype {
	Bondtype() {}
	Bondtype(string t1, string t2) : type1(t1), type2(t2) {
		sort();
	}
	Bondtype(string t1, string t2, float b0, float kb) : type1(t1), type2(t2), b0(b0), kb(kb) {
		sort();
	}
	Bondtype(int id1, int id2) : id1(id1), id2(id2) {
		//convertToZeroIndexed();
	}
	
	string type1, type2;
	int id1, id2;			// bonds from .top only has these values! 
	float b0;
	float kb;

	
	void sort() {
		int ptr = 0; 
		//cout << "Sorting " << type1 << "    " << type2 << endl;
		while (type1.length() > ptr && type2.length() > ptr) {

			if ((int)type1[ptr] < (int)type2[ptr]) {
				return;
			}
			else if ((int)type1[ptr] == (int)type2[ptr]) {
				// Do nothing, move ptr to the right
			}
			else {
				//printf("Swapping %d %d   ", (int) type1[ptr], (int) type2[ptr]);
				//cout << type1[ptr] << "   " << type2[ptr] << endl;
				//cout << type1 << "    " << type2 << endl << endl;
				swap(type1, type2);
				return;
			}
			ptr++;
		}
		if (type1.length() > type2.length()) {
			swap(type1, type2);
		}
		//printf("\n");
	}


	static vector<Bondtype> parseFFBondtypes(vector<vector<string>> rows) {
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
		printf("%d single bonds read from file\n", records.size());
		return records;
	}

	static vector<Bondtype> parseTopolBondtypes(vector<vector<string>> rows) {
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
		printf("%d bonds found in topology file\n", records.size());
		return records;	
	}

	static void assignTypesFromAtomIDs(vector<Bondtype>* topol_bonds, vector<Atom> atoms) {
		for (int i = 0; i < topol_bonds->size(); i++) {
			Bondtype* bond = &topol_bonds->at(i);


			bond->type1 = atoms.at(bond->id1 - 1).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
			bond->type2 = atoms.at(bond->id2 - 1).atomtype_bond;
			bond->sort();
			//cout << bond->type1 << '\t' << bond->type2 << endl;;
		}
	}




	static Bondtype getBondFromTypes(Bondtype* query_type, vector<Bondtype>* FF_bondtypes) { 		
		float best_likeness = 0;
		Bondtype best_bond;
		for (Bondtype bond : *FF_bondtypes) {
			float likeness = FTHelpers::calcLikeness(query_type->type1, bond.type1) * FTHelpers::calcLikeness(query_type->type2, bond.type2);
			if (likeness > best_likeness) {
				best_likeness = likeness;
				best_bond = bond;
			}
		}
		if (best_likeness > 0.01f)
			return best_bond;

		cout << "Failed to match bond types.\n Closest match " << best_bond.type1 << "    " << best_bond.type2;
		printf("Likeness %f\n", best_likeness);		
		printf("Topol ids: %d %d\n", query_type->id1, query_type->id2);
		cout << query_type->type1 << '\t' << query_type->type2 << endl;
		exit(0);
	}

	static void assignFFParametersFromBondtypes(vector<Bondtype>* topol_bonds, vector<Bondtype>* FF_bondtypes) {
		for (int i = 0; i < topol_bonds->size(); i++) {
			if (!(i%10))
			printf("\rAssigning FF parameters to bond %06d of %06d", i, topol_bonds->size());
			Bondtype* bond = &topol_bonds->at(i);

			Bondtype appropriateForcefield = getBondFromTypes(bond, FF_bondtypes);	// This might not return the correct one, as it tries to fit the atomtypes_bond to atomtypes_bond known in the CHARMM forcefield
			
			bond->kb = appropriateForcefield.kb;
			bond->b0 = appropriateForcefield.b0;
		}
		printf("\n");
	}
};




struct Angletype {
	Angletype() {}
	Angletype(string t1, string t2, string t3) : type1(t1), type2(t2), type3(t3) {
		sort();
	}
	Angletype(string t1, string t2, string t3, float t0, float kt) : type1(t1), type2(t2), type3(t3), theta0(t0), ktheta(kt) {
		sort();
	}
	Angletype(int id1, int id2, int id3) : id1(id1), id2(id2), id3(id3) {
	}

	string type1, type2, type3;			// left, middle, right
	int id1, id2, id3;			// bonds from .top only has these values! 
	float theta0;
	float ktheta;


	void sort(string* leftmost, string* rightmost) {
		int ptr = 0;
		while (leftmost->length() > ptr && rightmost->length() > ptr) {
			if ((int)(*leftmost)[ptr] < (int)(*rightmost)[ptr]) {
				return;
			}
			else if ((int)(*leftmost)[ptr] > (int)(*rightmost)[ptr]) {
				swap(*leftmost, *rightmost);
				return;
			}
			ptr++;
		}
		if (leftmost->length() > rightmost->length()) {
			swap(*leftmost, *rightmost);
		}
	}
	void sort() {
		sort(&type1, &type3);		// We can ONLY sort these two, as type2 MUSTTT be in the middle!!
		//sort(&type1, &type2);
		//sort(&type2, &type3);
	}




	static vector<Angletype> parseFFAngletypes(vector<vector<string>> rows) {
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
		printf("%d angletypes in forcefield\n", angletypes.size());
		return angletypes;
	}

	static vector<Angletype> parseTopolAngletypes(vector<vector<string>> rows) {
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
		printf("%d angles found in topology file\n", records.size());
		return records;
	}

	static void assignTypesFromAtomIDs(vector<Angletype>* topol_angles, vector<Atom> atoms) {
		for (int i = 0; i < topol_angles->size(); i++) {
			Angletype* angle = &topol_angles->at(i);

			angle->type1 = atoms.at(angle->id1 - 1).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
			angle->type2 = atoms.at(angle->id2 - 1).atomtype_bond;
			angle->type3 = atoms.at(angle->id3 - 1).atomtype_bond;
			angle->sort();
			//cout << bond->type1 << '\t' << bond->type2 << endl;;
		}
	}



	static Angletype getAngleFromTypes(Angletype* query_type, vector<Angletype>* FF_angletypes) {
		float best_likeness = 0;
		Angletype best_angle;
		for (Angletype angle : *FF_angletypes) {
			float likeness = FTHelpers::calcLikeness(query_type->type1, angle.type1) * FTHelpers::calcLikeness(query_type->type2, angle.type2) * FTHelpers::calcLikeness(query_type->type3, angle.type3);
			if (likeness > best_likeness) {
				best_likeness = likeness;
				best_angle = angle;
			}
			if (likeness == 1)
				break;
		}
		if (best_likeness > 0.01f)
			return best_angle;

		cout << "\n\n\nFailed to match angle types.\n Closest match " << best_angle.type1 << "    " << best_angle.type2 << "    " << best_angle.type3 << endl;
		printf("Likeness %f\n", best_likeness);
		printf("Topol ids: %d %d %d\n", query_type->id1, query_type->id2, query_type->id3);
		cout << query_type->type1 << '\t' << query_type->type2 << '\t' << query_type->type3 << endl;
		exit(0);
	}
	
	static void assignFFParametersFromAngletypes(vector<Angletype>* topol_angles, vector<Angletype>* FF_angletypes) {
		for (int i = 0; i < topol_angles->size(); i++) {
			printf("\rAssigning FF parameters to angle %06d of %06d", i, topol_angles->size());

			Angletype* angle = &topol_angles->at(i);

			Angletype appropriateForcefield = getAngleFromTypes(angle, FF_angletypes);	// This might not return the correct one, as it tries to fit the atomtypes_bond to atomtypes_bond known in the CHARMM forcefield

			angle->theta0 = appropriateForcefield.theta0;
			angle->ktheta = appropriateForcefield.ktheta;
		}
		printf("\n");

	}
	
};




struct Dihedraltype {
	Dihedraltype() {}
	Dihedraltype(const string& t1, const string& t2, const string& t3, const string& t4) : type1(t1), type2(t2), type3(t3), type4(t4) {
		sort();
	}
	Dihedraltype(const string& t1, const string& t2, const string& t3, const string& t4, float phi0, float kphi, int n) : type1(t1), type2(t2), type3(t3), type4(t4), phi0(phi0), kphi(kphi), n(n) {
		sort();
	}
	Dihedraltype(int id1, int id2, int id3, int id4) : id1(id1), id2(id2), id3(id3), id4(id4) {
	}

	string type1, type2, type3, type4;			// left, lm, rm, right
	int id1, id2, id3, id4;			// bonds from .top only has these values! 
	float phi0;
	float kphi;
	int n;


	void flip() {
		swap(type1, type4);
		swap(type2, type3);
	}
	void sort() {		
		if (type1 != type4) {
			if (!FTHelpers::isSorted(&type1, &type4)) {
				flip();
			}
		}
		else {			// In case the outer two is identical, we check the inner two.
			if (!FTHelpers::isSorted(&type2, &type3)) {
				flip();
			}
		}
	}



	static vector<Dihedraltype> parseFFDihedraltypes(vector<vector<string>> rows) {
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
		printf("%d dihedraltypes in forcefield\n", dihedraltypes.size());
		return dihedraltypes;
	}

	static vector<Dihedraltype> parseTopolDihedraltypes(vector<vector<string>> rows) {
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
		printf("%d dihedrals found in topology file\n", records.size());
		return records;
	}

	static void assignTypesFromAtomIDs(vector<Dihedraltype>* topol_dihedrals, vector<Atom> atoms) {
		//printf("GOt here!\n");
		for (int i = 0; i < topol_dihedrals->size(); i++) {
			Dihedraltype* dihedral = &topol_dihedrals->at(i);

			//printf("Accessing atoms %d %d %d %d\n", dihedral->id1, dihedral->id2, dihedral->id3, dihedral->id4);
			dihedral->type1 = atoms.at(dihedral->id1 - 1).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
			dihedral->type2 = atoms.at(dihedral->id2 - 1).atomtype_bond;
			dihedral->type3 = atoms.at(dihedral->id3 - 1).atomtype_bond;
			dihedral->type4 = atoms.at(dihedral->id4 - 1).atomtype_bond;
			dihedral->sort();
			//cout << bond->type1 << '\t' << bond->type2 << endl;;
		}
	}



	static Dihedraltype getDihedralFromTypes(Dihedraltype* query_type, vector<Dihedraltype>* forcefield) {
		float best_likeness = 0;
		Dihedraltype best_match;
		for (Dihedraltype dihedral : *forcefield) {
			float likeness = FTHelpers::calcLikeness(query_type->type1, dihedral.type1) * FTHelpers::calcLikeness(query_type->type2, dihedral.type2) * FTHelpers::calcLikeness(query_type->type3, dihedral.type3) * FTHelpers::calcLikeness(query_type->type4, dihedral.type4);
			if (likeness > best_likeness) {
				best_likeness = likeness;
				best_match = dihedral;
			}
			if (likeness == 1)
				break;
		}
		if (best_likeness > 0.001f)
			return best_match;

		cout << "\n\n\nFailed to match angle types.\nClosest match: " << best_match.type1 << "    " << best_match.type2 << "    " << best_match.type3 << "    " << best_match.type4 << endl;
		printf("Likeness %f\n", best_likeness);
		printf("Topol ids: %d %d %d %d\n", query_type->id1, query_type->id2, query_type->id3, query_type->id4);
		cout << query_type->type1 << '\t' << query_type->type2 << '\t' << query_type->type3 << '\t' << query_type->type4 << endl;
		exit(0);
	}

	static void assignFFParametersFromDihedraltypes(vector<Dihedraltype>* topol_dihedrals, vector<Dihedraltype>* forcefield) {
		for (int i = 0; i < topol_dihedrals->size(); i++) {
			printf("\rAssigning FF parameters to dihedral %06d of %06d", i, topol_dihedrals->size());

			Dihedraltype* dihedral = &topol_dihedrals->at(i);

			Dihedraltype appropriateForcefield = getDihedralFromTypes(dihedral, forcefield);	// This might not return the correct one, as it tries to fit the atomtypes_bond to atomtypes_bond known in the CHARMM forcefield

			dihedral->phi0 = appropriateForcefield.phi0;
			dihedral->kphi = appropriateForcefield.kphi;
			dihedral->n = appropriateForcefield.n;
		}
		printf("\n");
	}
};