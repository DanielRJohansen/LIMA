#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <format>
#include <functional>
#include "LIMA_BASE/include/Utilities.h"

///////////////////////////////// READ HERE FIRST /////////////////////////////////
// ffbonded.itp and ffnonbonden.itp has different atomtypes. 
// The purpose if this bit of code is to compress all forcefield parameters so they can fit on a GPU more easily
//





const float water_mass = 15.999000f + 2.f * 1.008000f;
const float water_sigma = 0.302905564168f + 2.f * 0.040001352445f;
const float water_epsilon = (0.50208f + 2.f * 0.19246f) *1000.f;		// Convert kJ->J

#include <map>

#include <span>

namespace FTHelpers {
	using std::vector, std::string;

	static bool charIsNumber(char c) {
		return ((int)c > 47 && (int)c < 58);
	}
	static bool charIsNumberAbove1(char c) {
		return ((int)c > 49 && (int)c < 58);
	}
	static float _calcLikeness(const string& query_type, const string& forcefield_type) {

		// Edgecase: perfect match
		if (query_type == forcefield_type) { return 1.f; }

		// Edgecase: wildcard
		if (forcefield_type == "X") {
			return 0.9;
		}

		float likeness = 0;
		float point_scale = 1.f / std::max(query_type.length(), forcefield_type.length());

		for (int i = 0; i < std::min(query_type.length(), forcefield_type.length()); i++) {
			if (query_type[i] == forcefield_type[i])
				likeness += point_scale;
			else
				break;
		}
		return likeness;
	}
		
			
	template <class DerivedType>
	static float calcLikeness(DerivedType query_type, const DerivedType &forcefield_type) {
		float likeness_unflipped = 1.f;
		for (int i = 0; i < DerivedType::n_atoms; i++) {
			likeness_unflipped *= FTHelpers::_calcLikeness(query_type.bonded_typenames[i], forcefield_type.bonded_typenames[i]);
		}

		query_type.flip();

		float likeness_flipped = 1.f;
		for (int i = 0; i < DerivedType::n_atoms; i++) {
			likeness_flipped *= FTHelpers::_calcLikeness(query_type.bonded_typenames[i], forcefield_type.bonded_typenames[i]);
		}

		return std::max(likeness_flipped, likeness_unflipped);
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


	static string makeBondTag(std::span<string>atom_ids) {
		string out = "";
		for (const auto& atom_id : atom_ids) {
			out = out + atom_id + "-";
		}
		out.pop_back();	// remove the last '-'
		return out;
	}

	template <typename GenericBondType>
	static void assignForceVariablesFromForcefield(vector<GenericBondType>& topol_bonds, const vector<GenericBondType>& forcefield) {
		std::unordered_map<string, const GenericBondType> forcefieldMap;

		for (int i = 0; i < topol_bonds.size(); i++) {
			//std::cout << std::format("\rAssigning FF parameters to {} {} of {}", 
			//	GenericBondType::getBondtype(), i+1, topol_bonds->size());

			GenericBondType* bond = &topol_bonds.at(i);

			// Try to see if we have already searched for a dihedral with an identical order of identical atom types
			string tag = FTHelpers::makeBondTag(bond->bonded_typenames);
			auto cachedFF = forcefieldMap.find(tag);

			if (cachedFF != forcefieldMap.end()) {
				const GenericBondType& appropriateForcefieldType = cachedFF->second;
				bond->assignForceVariables(appropriateForcefieldType);
			}
			else {
				const GenericBondType appropriateForcefieldType = GenericBondType::findBestMatchInForcefield(*bond, forcefield);
				forcefieldMap.insert({ tag, appropriateForcefieldType });
				bond->assignForceVariables(appropriateForcefieldType);
			}
		}
	}
};



// Use composition, so we both have the type which is simulation agnostic and the atom which is sim specific, 
// and has a reference or id of a type. The same should be done for each bond!
struct NB_Atomtype {
	NB_Atomtype(std::string t) : type(t){}		// This is for loading form the conf file, for easy comparisons
	NB_Atomtype(std::string t, float mass) : type(t), mass(mass) {}		// This if for forcefield merging
	NB_Atomtype(std::string t, float mass, float sigma, float epsilon) : type(t),  mass(mass), sigma(sigma), epsilon(epsilon) {}		// For LIMA_ffnonbonded
	NB_Atomtype(std::string t, int atnum, float mass, float sigma, float epsilon) : type(t), mass(mass), sigma(sigma), epsilon(epsilon) {}		// For random dudes ff_nonbonded
	NB_Atomtype(int gro_id, const std::string& atomtype) : gro_id(gro_id), type(atomtype) {}
	// Official parameters
	std::string type = "";
	//int atnum = -1;					// atnum given by input file (CHARMM)
	int atnum_local = 0;			// atnum specific to simulation	//TODO: rename to simulation id?
	float mass = -1;				// [g/mol]
	float sigma = -1;				// [nm]
	float epsilon = -1;				// J/mol

	int gro_id = -1;

	// LIMA parameters
	bool is_present_in_simulation = false;
	int simulation_specific_id = -1;


	static const NB_Atomtype* findRecord(const std::vector<NB_Atomtype>& records, const std::string& type) {
		for (const NB_Atomtype& record : records) {
			if (record.type == type) {
				return &record;
			}
		}
		return nullptr;
	}
	static bool typeIsPresent(const std::vector<NB_Atomtype>& records, std::string type) {
		return (findRecord(records, type) != nullptr);
	}
	
	/*static std::vector<NB_Atomtype> filterUnusedTypes(const std::vector<NB_Atomtype>& forcefield, 
		const std::vector<std::string>& active_types, Map& map, LimaLogger& logger, bool print_mappings);*/
};

// This is for bonded atoms!!!!!!!!!!!
struct Atom {
	Atom(int id, std::string atomtype, std::string atomname) : gro_id(id), atomname(atomname), atomtype(atomtype) {}
	const int gro_id=-1;										// Come from topol.top file
	const std::string atomtype{};	
	const std::string atomname{};	// I dunno what this is for
	int atomtype_id = -1;				// Asigned later
	//float charge;
};

using AtomTable = std::map<int, Atom>;













template <int n_atoms>	// n atoms in bond
struct BondtypeBase {
	BondtypeBase(const std::array<std::string, n_atoms>& typenames) : bonded_typenames(typenames) {}
	BondtypeBase(const std::array<int, n_atoms>& ids) : gro_ids(ids) {}	// TODO: remove this
	BondtypeBase(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: bonded_typenames(typenames), gro_ids(ids) {}

	//virtual void sort() = 0;
	virtual void flip() = 0;

	template <class DerivedType>
	static void assignTypesFromAtomIDs(std::vector<DerivedType>& topol_bonds, const AtomTable& atoms) {
		for (auto& bond : topol_bonds) {
			for (int i = 0; i < n_atoms; i++) {
				if (!atoms.contains(bond.gro_ids[i])) { throw "Atomtable does not contain bondedtype"; }
				bond.bonded_typenames[i] = atoms.find(bond.gro_ids[i])->second.atomtype_bond;
			}
			bond.sort();
		}
	}

	template <class DerivedType>
	static const DerivedType findBestMatchInForcefield(const DerivedType& query_type, const std::vector<DerivedType>& forcefield) {
		if (forcefield.size() == 0) { throw std::exception("No angletypes in forcefield!"); }

		//query_type.sort();

		float best_likeness = 0;
		DerivedType best_bond = forcefield.at(0);
		for (const DerivedType& ff_bondtype : forcefield) {
			/*float likeness = 1.f;
			for (int i = 0; i < n_atoms; i++) {
				likeness *= FTHelpers::calcLikeness(query_type.bonded_typenames[i], ff_bondtype.bonded_typenames[i]);
			}*/
			const float likeness = FTHelpers::calcLikeness(query_type, ff_bondtype);


			if (likeness > best_likeness) {
				best_likeness = likeness;
				best_bond = ff_bondtype;
			}
		}
		if (best_likeness > 0.01f)
			return best_bond;

		std::cout << "Failed to match bond types.\n Closest match ";
		for (auto& name : best_bond.bonded_typenames) {
			std::cout << name << " ";
		}
		// << best_bond.bonded_typenames[0] << "    " << best_bond.bonded_typenames[1];	//TODO: make this generic
		printf("\nLikeness %f\n", best_likeness);
		printf("Query typenames: ");
		for (auto& name : query_type.bonded_typenames) {
			std::cout << name << " ";
		}
		printf("\nQuery gro_ids: ");
		for (auto& id : query_type.gro_ids) {
			std::cout << std::to_string(id) << " ";
		}
		//std::cout << query_type.bonded_typenames[0] << '\t' << query_type.bonded_typenames[1] << std::endl;
		exit(0);
	}


	std::array<std::string, n_atoms> bonded_typenames;
	std::array<int, n_atoms> gro_ids;
};

struct Singlebondtype : public BondtypeBase<2>{
	static const int n_atoms = 2;
	Singlebondtype(const std::array<std::string, n_atoms>& typenames, float b0, float kb) : BondtypeBase(typenames), b0(b0), kb(kb) {
	}
	Singlebondtype(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: BondtypeBase(ids, typenames) {}

	float b0{};
	float kb{};

	void assignForceVariables(const Singlebondtype& a) {
		b0 = a.b0;
		kb = a.kb;
	}
	void flip() {
		swap(bonded_typenames[0], bonded_typenames[1]);
	}
	//void sort();
};




struct Anglebondtype : public BondtypeBase<3> {
	static const int n_atoms = 3;
	Anglebondtype(const std::array<std::string, n_atoms>& typenames, float t0, float kt)
		: BondtypeBase(typenames), theta0(t0), ktheta(kt) 
	{}
	//Anglebondtype(const std::array<int, n_atoms>& ids) : BondtypeBase(ids) {}
	Anglebondtype(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: BondtypeBase(ids, typenames) {}

	float theta0{};	// [rad]
	float ktheta{};

	void assignForceVariables(const Anglebondtype& a) {
		theta0 = a.theta0;
		ktheta = a.ktheta;
	}
	void flip() {
		swap(bonded_typenames[0], bonded_typenames[2]);
	}
	//void sort();
};

struct Dihedralbondtype : public BondtypeBase<4> {
	static const int n_atoms = 4;
	Dihedralbondtype(const std::array<std::string, n_atoms>& typenames, float phi0, float kphi, int n) 
		: BondtypeBase(typenames), phi0(phi0), kphi(kphi), n(n) {}
	Dihedralbondtype(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: BondtypeBase(ids, typenames) {}

	float phi0{};
	float kphi{};
	int n{};

	void assignForceVariables(const Dihedralbondtype& a) {
		phi0 = a.phi0;
		kphi = a.kphi;
		n = a.n;
	}

	void flip() {
		std::swap(bonded_typenames[0], bonded_typenames[3]);
		std::swap(bonded_typenames[1], bonded_typenames[2]);
	}
	//void sort();
};

struct Improperdihedralbondtype : public BondtypeBase<4> {
	static const int n_atoms = 4;
	// i j k l - https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html
	
	Improperdihedralbondtype(const std::array<std::string, n_atoms>& typenames, float psi0, float kpsi)
		: BondtypeBase(typenames), psi0(psi0), kpsi(kpsi){
	}
	Improperdihedralbondtype(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: BondtypeBase(ids, typenames) {}

	float psi0{};
	float kpsi{};

	void assignForceVariables(const Improperdihedralbondtype& a) {
		psi0 = a.psi0;
		kpsi = a.kpsi;
	}

	void flip() {
		//std::swap(bonded_typenames[0], bonded_typenames[3]);
		//std::swap(bonded_typenames[1], bonded_typenames[2]);
	}
	//void sort() {
	//	// Not sure what to do here yet
	//}
};