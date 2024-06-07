#pragma once

#include <string>

#include "Utilities.h"
#include "LimaTypes.cuh"
#include "MDFiles.h"
#include "ForcefieldTypes.h"





namespace ForcefieldHelpers{
	using std::string;

	float _calcLikeness(const string& query_type, const string& forcefield_type);


	template <class DerivedType>
	static float calcLikeness(DerivedType query_type, const DerivedType& forcefield_type) {
		float likeness = 1.f;
		for (int i = 0; i < DerivedType::n_atoms; i++) {
			likeness *= _calcLikeness(query_type.bonded_typenames[i], forcefield_type.bonded_typenames[i]);
		}

		return likeness;
	}




	template <typename GenericBondType>
	static int findBestMatchInForcefield_Index(const GenericBondType& query_type, const std::vector<GenericBondType>& forcefield) {
		if (forcefield.size() == 0) { throw std::runtime_error("No bonds in forcefield!"); }

		float best_likeness = 0;
		int bestBondIndex = 0;
		for (int i = 0; i < forcefield.size(); i++) {
			const GenericBondType& ff_bondtype = forcefield[i];
			const float likeness = calcLikeness(query_type, ff_bondtype);

			if (likeness > best_likeness) {
				best_likeness = likeness;
				bestBondIndex = i;
			}
		}

		if (best_likeness > 0.05f) {
			return bestBondIndex;
		}

		std::cout << "Failed to match bond types.\n Closest match ";
		for (const auto& name : forcefield[bestBondIndex].bonded_typenames) {
			std::cout << name << " ";
		}
		if constexpr (std::is_same_v<GenericBondType, Dihedralbondtype>) {
			std::cout << "Dihedral type\n";
		}
		else {
			std::cout << "Improper type\n";
		}
		printf("\nLikeness %f\n", best_likeness);
		printf("Query typenames: ");
		for (auto& name : query_type.bonded_typenames) {
			std::cout << name << " ";
		}
		printf("\nQuery gro_ids: ");
		for (auto& id : query_type.global_ids) {
			std::cout << std::to_string(id) << " ";
		}
		throw std::runtime_error("\nfindBestMatchInForcefield failed");
	}

};

template <typename GenericBondType>
class ParameterDatabase {
public:
	const GenericBondType& get(const std::array<std::string, GenericBondType::n_atoms>& query) {
		const std::string key = MakeKey(query);
		if (fastLookup.count(key) != 0)
			return parameters[fastLookup.find(key)->second];

		GenericBondType temp(query);	// TEMP, supply a public sort()
		const int index = ForcefieldHelpers::findBestMatchInForcefield_Index(temp, parameters);
		fastLookup.insert({ key, index });
		return parameters[index];
	}

	void insert(const GenericBondType& element) {
		if (fastLookup.count(MakeKey(element.bonded_typenames)) != 0)
			return;
		parameters.push_back(element);
		fastLookup.insert({ MakeKey(element.bonded_typenames), parameters.size()-1});
	}

private:
	std::vector<GenericBondType> parameters;
	std::unordered_map<std::string, int> fastLookup;	// Map a query to an index in parameters

	std::string MakeKey(const std::array<std::string, GenericBondType::n_atoms>& arr) {
		std::string key(GenericBondType::n_atoms * 5, ' '); // 5 is the maximum length of a typename
		

		for (size_t i = 0; i < GenericBondType::n_atoms; ++i) {
			std::strcpy(&key[i * 5], arr[i].c_str());
		}
		return key;
	}
};


class LjParameterDatabase {
	ParameterDatabase<LJParameter> ljParameters;

	std::vector<LJParameter> activeLjParameters;
	std::unordered_map<std::string, int> fastLookup;	// Map a query to an index in activeParameters

	bool finished = false;

	public:
		LjParameterDatabase();
		void insert(LJParameter element);
		int GetActiveIndex(const std::string& query);
		std::vector<LJParameter>& GetActiveParameters() {
			finished = true;			
			return activeLjParameters;
		}
};

// TODO move to the other file, delete this file
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
		if constexpr (std::is_same<BondParamType, Singlebondtype>::value) {
			return singlebondParameters.get(query);
		}
		else if constexpr (std::is_same<BondParamType, Anglebondtype>::value ) {
			return anglebondParameters.get(query);
		}
		else if constexpr (std::is_same<BondParamType, Dihedralbondtype>::value ) {
			return dihedralbondParameters.get(query);
		}
		else if constexpr (std::is_same<BondParamType, Improperdihedralbondtype>::value ) {
			return improperdihedralbondParameters.get(query);
		}
		else {
			static BondParamType dummy; // Fallback to handle unsupported cases
			return dummy;
		}
	}


private:
	/*ParameterDatabase<LJParameter> ljParameters;
	std::vector<LJParameter> activeLjParameters;*/
	LjParameterDatabase ljParameters;

	ParameterDatabase<Singlebondtype> singlebondParameters;
	ParameterDatabase<Anglebondtype> anglebondParameters;
	ParameterDatabase<Dihedralbondtype> dihedralbondParameters;
	ParameterDatabase<Improperdihedralbondtype> improperdihedralbondParameters;

	void loadFileIntoForcefield(const SimpleParsedFile& parsedfile);
};