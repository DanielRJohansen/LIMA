#pragma once

#include <string>

#include "Utilities.h"
#include "LimaTypes.cuh"
#include "MDFiles.h"
#include "ForcefieldTypes.h"

struct NB_Atomtype;

namespace LimaForcefieldBuilder {

	/// <summary>
	/// Create ffbonded.lff and ffnonbonded.lff files.
	/// </summary>
	/// <param name="molecule_dir">Dir where conf and topol are, with all .itp include files</param>
	/// <param name="output_dir">Where .lff files will be created</param>
	/// <param name="conf_name">name of main .gro file in molecule_dir</param>
	/// <param name="topol_name">name of main .top/.itp file in molecule_dir</param>
	/// <param name="envmode"></param>
	void buildForcefield(const std::string& molecule_dir, const std::string& output_dir,
		const ParsedTopologyFile& topol_file, EnvMode envmode);
}



template <typename GenericBondType>
class ParameterDatabase {
public:
	const GenericBondType& get(const std::array<std::string, GenericBondType::n_atoms>& query) {
		GenericBondType temp(query);	// TEMP, supply a public sort()
		const std::string key = MakeKey(temp.bonded_typenames);
		if (fastLookup.count(key) != 0)
			return parameters[fastLookup.find(key)->second];

		int index = FTHelpers::findBestMatchInForcefield_Index(temp, parameters);
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
		std::ostringstream key;
		for (size_t i = 0; i < GenericBondType::n_atoms; ++i) {
			key << (i ? "|" : "") << arr[i];  // Inline conditional adds separator for all but the first element
		}
		return key.str();
	}
};


class LjParameterDatabase {
	ParameterDatabase<LJParameter> ljParameters;

	std::vector<LJParameter> activeLjParameters;
	std::unordered_map<std::string, int> fastLookup;	// Map a query to an index in activeParameters

	bool finished = false;

	public:
		void insert(LJParameter element) {
			if (finished)
				throw std::runtime_error("Cannot insert after finishing");
			ljParameters.insert(element);
		}
		int GetActiveIndex(const std::string& query) {
			if (finished)
				throw std::runtime_error("Cannot query for more active types after finishing");
			if (fastLookup.count(query) != 0)
				return fastLookup.find(query)->second;

			const LJParameter& parameter = ljParameters.get({ query });
			// First check of the parameter is already in the active list (since it might be 
			// lexicographically similar but not identical, we cant expect the fastlookup above to always catch)
			for (int i = 0; i < activeLjParameters.size(); i++) {
				if (parameter.bonded_typenames == activeLjParameters[i].bonded_typenames) {
					fastLookup.insert({ query, i });
					return i;
				}
			}

			activeLjParameters.push_back(parameter);
			fastLookup.insert({ query, activeLjParameters.size() - 1 });
			return activeLjParameters.size() - 1;
		}
		std::vector<LJParameter>& GetActiveParameters() {
			finished = true;			
			return activeLjParameters;
		}
};

class LIMAForcefield {
public:
	LIMAForcefield();
	LIMAForcefield(const LIMAForcefield&) = delete;

	/*const Singlebondtype& GetSinglebondParameters(const std::array<std::string, 2>& query) {
		return singlebondParameters.get(query);
	}*/

	int GetActiveLjParameterIndex(const std::string& query) {
		return ljParameters.GetActiveIndex(query);
	}
	std::vector<LJParameter>& GetActiveLjParameters() {
		return ljParameters.GetActiveParameters();	// After this, the Get() above should no longer be legal
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