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

//using namespace std;

using std::vector, std::string;

const float water_mass = 15.999000f + 2.f * 1.008000f;
const float water_sigma = 0.302905564168f + 2.f * 0.040001352445f;
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
		float point_scale = 1.f / std::max(a.length(), b.length());

		for (int i = 0; i < std::min(a.length(), b.length()); i++) {
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

		printf("%zd atom types in conf\n", atom_types.size());

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

	static string makeBondTag(std::vector<string>atom_ids) {
		string out = "";
		for (const auto& atom_id : atom_ids) {
			out = out + atom_id + "-";
		}
		out.pop_back();	// remove the last '-'
		return out;
	}
};


struct Map {
	struct Mapping {
		Mapping () {}
		Mapping (string l, string r) : left(l), right(r) {}
		string left{}, right{};
	};



	Map() {
		mappings = new Mapping[10000];
	}
	~Map() {
		delete[] mappings;
	}

	Mapping* mappings;
	int n_mappings = 0; 

	bool mapExists(std::string l) {
		for (int i = 0; i < n_mappings; i++)
			if (l == mappings[i].left)
				return true;
		return false;
	}
	void addMap(std::string l, std::string r) {
		if (!mapExists(l))
			mappings[n_mappings++] = Mapping(l, r);
	}
	std::string mapRight(const std::string& query) {
		for (int i = 0; i < n_mappings; i++) {
			if (mappings[i].left == query) {
				return mappings[i].right;
			}
		}
		return query;
	}
};


struct NB_Atomtype {
	NB_Atomtype() {}
	NB_Atomtype(std::string t) : type(t){}		// This is for loading form the conf file, for easy comparisons
	NB_Atomtype(std::string t, float mass) : type(t), mass(mass) {}		// This if for forcefield merging
	NB_Atomtype(std::string t, float mass, float sigma, float epsilon) : type(t),  mass(mass), sigma(sigma), epsilon(epsilon) {}		// For LIMA_ffnonbonded
	NB_Atomtype(std::string t, int atnum, float mass, float sigma, float epsilon) : type(t), atnum(atnum), mass(mass), sigma(sigma), epsilon(epsilon) {}		// For random dudes ff_nonbonded
	// Official parameters
	std::string type = "";
	int atnum = -1;					// atnum given by input file (CHARMM)
	int atnum_local = 0;			// atnum specific to simulation
	float mass = -1;		// [g/mol]
	float sigma = -1;	// [nm]
	float epsilon = -1;	// J/mol

	// LIMA parameters
	bool is_present_in_simulation = false;
	int simulation_specific_id = -1;





	//bool operator==(const NB_Atomtype a) { return (type == a.type); }
	//static bool sameType(NB_Atomtype a, NB_Atomtype b) { return (a.type == b.type); }


	// Parser functions

	static std::vector<NB_Atomtype> parseNonbonded(std::vector<std::vector<std::string>> rows) {		// KNOWN ERROR HERE. SOMETIMES SPACES IN FILES ARE TABS, AND IT DOESN'T RECOGNISE A ROW AS A PROPER ENTRY!!
		FTHelpers::STATE current_state = FTHelpers::INACTIVE;

		std::vector<NB_Atomtype> records;

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
		printf("%lld atom types read from file\n", records.size());
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
	static vector<NB_Atomtype> filterUnusedTypes(vector<NB_Atomtype> forcefield, vector<string> active_types, Map* map);
};


// This is for bonded atoms!!!!!!!!!!!
struct Atom {
	Atom() {}
	Atom(int id, string type_b, string type_nb) : id(id), atomtype_bond(type_b), atomtype(type_nb) {
		//convertToZeroindexed();
	}
	int id=-1;										// Come from topol.top file
	string atomtype{};
	string atomtype_bond{};
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

	static vector<Atom> parseTopolAtoms(vector<vector<string>>& rows);


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
	
	string type1{}, type2{};
	int id1{}, id2{};			// bonds from .top only has these values! 
	float b0{};
	float kb{};

	void sort() {
		if (!FTHelpers::isSorted(&type1, &type2)) {
			swap(type1, type2);
		}
	}
	
	
	//void sort() {
	//	int ptr = 0; 
	//	//cout << "Sorting " << type1 << "    " << type2 << endl;
	//	while (type1.length() > ptr && type2.length() > ptr) {

	//		if ((int)type1[ptr] < (int)type2[ptr]) {
	//			return;
	//		}
	//		else if ((int)type1[ptr] == (int)type2[ptr]) {
	//			// Do nothing, move ptr to the right
	//		}
	//		else {
	//			//printf("Swapping %d %d   ", (int) type1[ptr], (int) type2[ptr]);
	//			//cout << type1[ptr] << "   " << type2[ptr] << endl;
	//			//cout << type1 << "    " << type2 << endl << endl;
	//			swap(type1, type2);
	//			return;
	//		}
	//		ptr++;
	//	}
	//	if (type1.length() > type2.length()) {
	//		swap(type1, type2);
	//	}
	//	//printf("\n");
	//}


	static vector<Bondtype> parseFFBondtypes(vector<vector<string>> rows);

	static vector<Bondtype> parseTopolBondtypes(vector<vector<string>> rows);

	static void assignTypesFromAtomIDs(vector<Bondtype>* topol_bonds, vector<Atom> atoms);


	static Bondtype getBondFromTypes(Bondtype* query_type, vector<Bondtype>* FF_bondtypes);

	static void assignFFParametersFromBondtypes(vector<Bondtype>* topol_bonds, vector<Bondtype>* FF_bondtypes);
};




struct Angletype {
	Angletype() {}
	Angletype(string t1, string t2, string t3) : type1(t1), type2(t2), type3(t3) {
		sort();
	}
	Angletype(string t1, string t2, string t3, float t0, float kt) : type1(t1), type2(t2), type3(t3), theta0(t0), ktheta(kt) {
		sort();
	}
	Angletype(int id1, int id2, int id3) : id1(id1), id2(id2), id3(id3) {}

	string type1{}, type2{}, type3{};			// left, middle, right
	int id1{}, id2{}, id3{};			// bonds from .top only has these values! 
	float theta0{};
	float ktheta{};

	void assignForceVariables(const Angletype& a) {
		theta0 = a.theta0;
		ktheta = a.ktheta;
	}


	/*void sort(string* leftmost, string* rightmost) {
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
	}*/

	void sort() {
		// The middle type must always be in the middle, so we only sort the outer 2 types
		// No need to sort if our edges are the same
		if (type1 != type3) {	
			if (!FTHelpers::isSorted(&type1, &type3)) {
				swap(type1, type3);
			}
		}
	}


	static vector<Angletype> parseFFAngletypes(vector<vector<string>> rows);

	static vector<Angletype> parseTopolAngletypes(vector<vector<string>> rows);

	static void assignTypesFromAtomIDs(vector<Angletype>* topol_angles, vector<Atom> atoms);



	static Angletype* getAngleFromTypes(Angletype* query_type, vector<Angletype>* FF_angletypes);
	
	static void assignFFParametersFromAngletypes(vector<Angletype>* topol_angles, vector<Angletype>* forcefield);
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

	void assignForceVariables(const Dihedraltype& a) {
		phi0 = a.phi0;
		kphi = a.kphi;
		n = a.n;
	}

	string type1, type2, type3, type4;			// left, lm, rm, right
	int id1{}, id2{}, id3{}, id4{};			// bonds from .top only has these values! 
	float phi0{};
	float kphi{};
	int n{};


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



	static vector<Dihedraltype> parseFFDihedraltypes(vector<vector<string>> rows);

	static vector<Dihedraltype> parseTopolDihedraltypes(vector<vector<string>> rows);

	static void assignTypesFromAtomIDs(vector<Dihedraltype>* topol_dihedrals, vector<Atom> atoms);

	static Dihedraltype* getDihedralFromTypes(Dihedraltype* query_type, vector<Dihedraltype>* forcefield);

	static void assignFFParametersFromDihedraltypes(vector<Dihedraltype>* topol_dihedrals, vector<Dihedraltype>* forcefield);
};