#pragma once

#include <string>

#include "Utilities.h"
#include "LimaTypes.cuh"
#include "MDFiles.h"
#include <span>



/////
// These are the full representation of what information a AtomType or Bondtype contains.
// Instances of such bonds are defined more compactly in bodies.cuh
/////

struct AtomType {
	std::string name{};
	int atNum{};
	ForceField_NB::ParticleParameters parameters{};
	float charge{}; // [kilo C/mol]
	char ptype{};
};

struct SinglebondType {
	static const int nAtoms = 2;
	std::array<std::string, 2> bonded_typenames; // i,j
	int func{};
	SingleBond::Parameters params;
};

struct AnglebondType {
	static const int nAtoms = 3;
	std::array<std::string, 3> bonded_typenames; // i,j,k
	int func{};
	AngleBond::Parameters params;
	float ub0;	// No idea what this is
	float cub;	// No idea what this is
};

struct DihedralbondType {
	static const int nAtoms = 4;
	std::array<std::string, 4> bonded_typenames; // ijkl
	int func{};
	DihedralBond::Parameters params;
};

struct ImproperDihedralbondType {
	static const int nAtoms = 4;
	std::array<std::string, 4> bonded_typenames; // i j k l - https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html
	int func{};
	ImproperDihedralBond::Parameters params;
};

template <typename bond>
void FlipBondedtypeNames(std::array<std::string, bond::nAtoms>& bondedTypeNames) {
	if constexpr (std::is_same<bond, SinglebondType>::value) {
		swap(bondedTypeNames[0], bondedTypeNames[1]);
	}
	else if constexpr (std::is_same<bond, AnglebondType>::value) {
		swap(bondedTypeNames[0], bondedTypeNames[2]);
	}
	else if constexpr (std::is_same<bond, DihedralbondType>::value) {

		std::swap(bondedTypeNames[0], bondedTypeNames[3]);
		std::swap(bondedTypeNames[1], bondedTypeNames[2]);

	}
	else if constexpr (std::is_same<bond, ImproperDihedralbondType>::value) {
		std::swap(bondedTypeNames[0], bondedTypeNames[3]);
	}
	else {
		throw std::runtime_error("Illegal bond type");
	}
}

///-----------------------------------------------------------------------------------///



template <typename GenericBondType>
class ParameterDatabase {
public:
	const GenericBondType& get(const std::array<std::string, GenericBondType::nAtoms>& query) {
		const std::string key = MakeKey(query);
		if (fastLookup.count(key) != 0)
			return parameters[fastLookup.find(key)->second];

		const int index = findBestMatchInForcefield_Index(query);
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

	int findBestMatchInForcefield_Index(const std::array<std::string, GenericBondType::nAtoms>& query);
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