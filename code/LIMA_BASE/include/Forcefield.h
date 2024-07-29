#pragma once

#include <string>

#include "Utilities.h"
#include "LimaTypes.cuh"
#include "MDFiles.h"
#include "ForcefieldTypes.h"





namespace ForcefieldHelpers{
	using std::string;

	float _calcLikeness(const string& query_type, const string& forcefield_type);


	template <class GenericBondType>
	static float calcLikeness(std::array<std::string, GenericBondType::nAtoms> query, const GenericBondType& forcefield_type) {
		float likeness = 1.f;
		for (int i = 0; i < GenericBondType::nAtoms; i++) {
			likeness *= _calcLikeness(query[i], forcefield_type.bonded_typenames[i]);
		}

		return likeness;
	}



	// TODONOW QueryTyoe should only be the atomnames, we dont have the other info
	template <typename GenericBondType>
	static int findBestMatchInForcefield_Index(const std::array<std::string, GenericBondType::nAtoms>& query, 
		const std::vector<GenericBondType>& forcefield) {
		if (forcefield.size() == 0) { throw std::runtime_error("No bonds in forcefield!"); }

		float best_likeness = 0;
		int bestBondIndex = 0;
		for (int i = 0; i < forcefield.size(); i++) {
			const GenericBondType& ff_bondtype = forcefield[i];
			const float likeness = calcLikeness(query, ff_bondtype);

			if (likeness > best_likeness) {
				best_likeness = likeness;
				bestBondIndex = i;
			}
		}

		if (best_likeness > 0.001f) {
			return bestBondIndex;
		}

		std::cout << "Failed to match bond types.\n Closest match ";
		for (const auto& name : forcefield[bestBondIndex].bonded_typenames) {
			std::cout << name << " ";
		}
		if constexpr (std::is_same_v<GenericBondType, DihedralbondType>) {
			std::cout << "Dihedral type\n";
		}
		else {
			std::cout << "Improper type\n";
		}
		printf("\nLikeness %f\n", best_likeness);
		printf("Query typenames: ");
		for (auto& name : query) {
			std::cout << name << " ";
		}
		//printf("\nQuery gro_ids: ");
		//for (auto& id : query_type.global_ids) {
		//	std::cout << std::to_string(id) << " ";
		//}
		throw std::runtime_error("\nfindBestMatchInForcefield failed");
	}

};

template <typename GenericBondType>
class ParameterDatabase {
public:
	const GenericBondType& get(const std::array<std::string, GenericBondType::nAtoms>& query) {
		const std::string key = MakeKey(query);
		if (fastLookup.count(key) != 0)
			return parameters[fastLookup.find(key)->second];

		//GenericBondType temp(query);	// TEMP, supply a public sort()
		const int index = ForcefieldHelpers::findBestMatchInForcefield_Index(query, parameters);
		fastLookup.insert({ key, index });
		return parameters[index];
	}


	// Due to types such as this:
	// X	CT1	OH1	X	9	0.00	0.58576	3
	// We cannot sort and have a single definition of each type, as the names in the topol entry at the 'X'
	// places may be in either order and match here
	void insert(GenericBondType element) {
		__insert(element);
		FlipBondedtypeNames<GenericBondType>(element.bonded_typenames);
		__insert(element);
	}

private:
	void __insert(const GenericBondType& element) {
		if (fastLookup.count(MakeKey(element.bonded_typenames)) != 0)
			return;
		parameters.push_back(element);
		fastLookup.insert({ MakeKey(element.bonded_typenames), parameters.size() - 1 });
	}
	std::vector<GenericBondType> parameters;
	std::unordered_map<std::string, int> fastLookup;	// Map a query to an index in parameters

	std::string MakeKey(const std::array<std::string, GenericBondType::nAtoms>& arr) {
		std::string key(GenericBondType::nAtoms * 5, ' '); // 5 is the maximum length of a typename
		

		for (size_t i = 0; i < GenericBondType::nAtoms; ++i) {
			std::strcpy(&key[i * 5], arr[i].c_str());
		}
		return key;
	}
};


class LjParameterDatabase {
	std::unordered_map<std::string, AtomType> atomTypes;

	std::vector<AtomType> activeAtomTypes;
	std::unordered_map<std::string, int> fastLookup;	// Map a query to an index in activeParameters

	bool finished = false;

	public:
		LjParameterDatabase();
		void insert(AtomType element);
		int GetActiveIndex(const std::string& query);
		std::vector<AtomType>& GetActiveParameters() {
			finished = true;			
			return activeAtomTypes;
		}
};

class LIMAForcefield {
public:
	LIMAForcefield();
	LIMAForcefield(const LIMAForcefield&) = delete;


	int GetActiveLjParameterIndex(const std::string& query) {
		return ljParameters.GetActiveIndex(query);
	}
	ForceField_NB GetActiveLjParameters() {
		ForceField_NB forcefield{};
		const auto& activeParameters = ljParameters.GetActiveParameters();
		for (int i = 0; i < activeParameters.size(); i++)
			forcefield.particle_parameters[i] = activeParameters[i].parameters;
		return forcefield;
	}


	template<typename BondParamType>
	const BondParamType& GetBondParameters(const auto& query) {
		if constexpr (std::is_same<BondParamType, SingleBond::Parameters>::value) {
			return singlebondParameters.get(query).params;
		}
		else if constexpr (std::is_same<BondParamType, AngleBond::Parameters>::value ) {
			return anglebondParameters.get(query).params;
		}
		else if constexpr (std::is_same<BondParamType, DihedralBond::Parameters>::value ) {
			return dihedralbondParameters.get(query).params;
		}
		else if constexpr (std::is_same<BondParamType, ImproperDihedralBond::Parameters>::value ) {
			return improperdihedralbondParameters.get(query).params;
		}
		else {
			throw std::runtime_error("Unsupported bond type");
			//static_assert(false, "Unsupported bond type");
			//static BondParamType dummy; // Fallback to handle unsupported cases&
			//return dummy;
		}
	}


private:
	LjParameterDatabase ljParameters;

	ParameterDatabase<SinglebondType> singlebondParameters;
	ParameterDatabase<AnglebondType> anglebondParameters;
	ParameterDatabase<DihedralbondType> dihedralbondParameters;
	ParameterDatabase<ImproperDihedralbondType> improperdihedralbondParameters;

	void LoadFileIntoForcefield(const fs::path& path);
};