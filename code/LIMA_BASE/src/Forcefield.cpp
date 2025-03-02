#include <vector>

#include "Forcefield.h"
#include "MDFiles.h"

#include <format>
#include <tuple>

using std::string;

class AtomtypeDatabase {
	std::unordered_map<std::string, AtomType> atomTypes;

	std::shared_ptr<std::vector<AtomType>> activeAtomTypes;
	std::unordered_map<std::string, int> fastLookup;	// Map a query to an index in activeParameters

	bool finished = false;

public:
	AtomtypeDatabase() {
		activeAtomTypes = std::make_shared<std::vector<AtomType>>();
	}
	void insert(AtomType element);
	int GetActiveIndex(const std::string& query);
	std::vector<AtomType>& GetActiveParameters() {
		finished = true;
		return *activeAtomTypes;
	}

	AtomType& GetAtomType(const std::string& query) {
		return atomTypes.at(query);
	}

	// Debug only
	std::unordered_map<std::string, AtomType>& _getAll() {
		return atomTypes;
	}
};


void AtomtypeDatabase::insert(AtomType element) {
	if (finished)
		throw std::runtime_error("Cannot insert after finishing");

	// Only add a type the first time we see it. TODO: We should also do this for all bonded types
	if (!atomTypes.contains(element.name))
		atomTypes.insert({ element.name, element });
}
int AtomtypeDatabase::GetActiveIndex(const std::string& query) {
	if (finished)
		throw std::runtime_error("Cannot query for more active types after finishing");
	if (fastLookup.count(query) != 0)
		return fastLookup.find(query)->second;
	if (!atomTypes.contains(query))
		throw std::runtime_error(std::format("Failed to find atomtype [{}]", query));

	const AtomType& parameter = atomTypes.at({ query });
	// First check of the parameter is already in the active list (since it might be 
	// lexicographically similar but not identical, we cant expect the fastlookup above to always catch)
	for (int i = 0; i < activeAtomTypes->size(); i++) {
		if (parameter.name == (*activeAtomTypes)[i].name) {
			fastLookup.insert({ query, i });
			return i;
		}
	}
	if (activeAtomTypes->size() == ForceField_NB::MAX_TYPES)
		throw std::runtime_error("Cant handle so many unique atomtypes");
	activeAtomTypes->push_back(parameter);
	fastLookup.insert({ query, activeAtomTypes->size() - 1 });
	return activeAtomTypes->size() - 1;
}












template <typename GenericBondType>
class ParameterDatabase {
public:
	const std::vector<typename GenericBondType::Parameters>& get(const std::array<std::string, GenericBondType::nAtoms>& query) {
		const std::string key = MakeKey(query);

		if (fastLookup.count(key) == 0) {
			const std::vector<int> indices = findBestMatchesInForcefield_Index(query);

			// If we found no matches, return an empty list
			if (indices.empty())
				return emptyParams;

			std::vector<typename GenericBondType::Parameters> bondparams;
			for (const auto& index : indices) {
				bondparams.push_back(parameters[index].params);
			}

			fastLookup.insert({ key, bondparams });
		}

		return fastLookup.find(key)->second;
	}

	// Due to types such as this:
	// X	CT1	OH1	X	9	0.00	0.58576	3
	// We cannot sort and have a single definition of each type, as the names in the topol entry at the 'X'
	// places may be in either order and match here
	void insert(const GenericBondType& element) {
		__insert(element);
		GenericBondType flippedBond = element;
		FlipBondedtypeNames<GenericBondType>(flippedBond.bonded_typenames);
		if (flippedBond.bonded_typenames != element.bonded_typenames)
			__insert(flippedBond);
	}

private:

	template <typename bond>
	void FlipBondedtypeNames(std::array<std::string, bond::nAtoms>& bondedTypeNames) {
		if constexpr (std::is_same<bond, SinglebondType>::value) {
			swap(bondedTypeNames[0], bondedTypeNames[1]);
		}
		else if constexpr (std::is_same<bond, PairbondType>::value) {
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

	void __insert(const GenericBondType& element) {
		if (locked) {
			throw std::runtime_error("Cannot add to database after it has been locked");
		}
		for (const auto& param : parameters) {
			if (param.bonded_typenames == element.bonded_typenames && param.params == element.params) {
				// I think this will be a problem, either because we include multiple files for the forcefiel, either 
				// custom or there are duplicates in the ffnabonded.itp, or it may just be a problem with the forcefield itself
				//throw std::runtime_error("Duplicate bond type found in forcefield");


				// This does happen in Slipids, then duplicates denoted with  ;*new...... So lets just ignore the problem, and ignore the duplicates
				return;
			}
		}
		parameters.push_back(element);
	}

	const std::vector<typename GenericBondType::Parameters> emptyParams{}; // Return this if we find no matches
	std::vector<GenericBondType> parameters;
	std::unordered_map<std::string, std::vector<typename GenericBondType::Parameters>> fastLookup;	// Map a query to an index in parameters

	bool locked = false; // Once this is true, we may no longer add to the database

	std::string MakeKey(const std::array<std::string, GenericBondType::nAtoms>& arr) {
		std::string key(GenericBondType::nAtoms * 5, ' '); // 5 is the maximum length of a typename


		for (size_t i = 0; i < GenericBondType::nAtoms; ++i) {
			std::strcpy(&key[i * 5], arr[i].c_str());
		}
		return key;
	}


	// Returns {isMatch, wildcardCoint}, so 0 is perfect match, 3 in a dihedral, is a poor match
	constexpr std::tuple<bool, int> DetermineMatchDegree(const std::span<const std::string>& query, const std::span<const std::string>& typeInForcefield) {
		int wildcardCount = 0;
		const auto wildcard = "X";

		for (int i = 0; i < query.size(); i++) {
			if (query[i] == typeInForcefield[i]) {
				continue;
			}
			else if (typeInForcefield[i] == wildcard) {
				wildcardCount++;
				continue;
			}
			return { false, 0 };
		}
		return { true, wildcardCount };
	}

	std::vector<int> findBestMatchesInForcefield_Index(const std::array<std::string, GenericBondType::nAtoms>& query) {
		if (parameters.size() == 0) {
			throw std::runtime_error("No bonds in forcefield!");
		}

		//int bestBondIndex = 0;
		std::vector<int> bestBondIndices{};	// Can be multiple bonds, if they are equally good matches
		int lowestWildcardDegree = INT_MAX;


		for (size_t i = 0; i < parameters.size(); i++) {
			auto [isMatch, wildcardDegree] = DetermineMatchDegree(query, parameters[i].bonded_typenames);

			if (isMatch) {
				if (wildcardDegree < lowestWildcardDegree) {
					lowestWildcardDegree = wildcardDegree;
					bestBondIndices.clear();
					bestBondIndices.push_back(i);
				}
				else if (wildcardDegree == lowestWildcardDegree) {
					bestBondIndices.push_back(i);
				}
			}
		}
		return bestBondIndices;
	}
};
template class ParameterDatabase<SinglebondType>;
template class ParameterDatabase<PairbondType>;
template class ParameterDatabase<AnglebondType>;
template class ParameterDatabase<DihedralbondType>;
template class ParameterDatabase<ImproperDihedralbondType>;














LIMAForcefield::LIMAForcefield() {};

LIMAForcefield::LIMAForcefield(const GenericItpFile& file) {

	ljParameters = std::make_unique<AtomtypeDatabase>();
	tinymolTypes = std::make_unique<AtomtypeDatabase>();

	singlebondParameters = std::make_unique<ParameterDatabase<SinglebondType>>();
	pairbondParameters = std::make_unique<ParameterDatabase<PairbondType>>();
	anglebondParameters = std::make_unique<ParameterDatabase<AnglebondType>>();
	dihedralbondParameters = std::make_unique<ParameterDatabase<DihedralbondType>>();
	improperdihedralbondParameters = std::make_unique<ParameterDatabase<ImproperDihedralbondType>>();

	LoadFileIntoForcefield(file);

	// TEMP while we force solvents to be singleparticle
	if constexpr (!AllAtom) {
		if (tinymolTypes->_getAll().contains("OW"))
			tinymolTypes->_getAll().at("OW").mass += 2.f * tinymolTypes->_getAll().at("HW").mass;
	}
}

LIMAForcefield::~LIMAForcefield() {}

int LIMAForcefield::GetActiveLjParameterIndex(const std::string& query) {
	return ljParameters->GetActiveIndex(query);
}
ForceField_NB LIMAForcefield::GetActiveLjParameters() {
	ForceField_NB forcefieldNB{};
	const std::vector<AtomType>& activeParameters = ljParameters->GetActiveParameters();//  *activeLJParamtypes;
	if (activeParameters.size() > ForceField_NB::MAX_TYPES)
		throw std::runtime_error("Too many atom types");
	for (int i = 0; i < activeParameters.size(); i++)
		forcefieldNB.particle_parameters[i] = activeParameters[i].parameters;
	return forcefieldNB;
}

std::vector<NonbondedInteractionParams> LIMAForcefield::GetNonbondedInteractionParams() const {
    std::vector<NonbondedInteractionParams> nonbondedInteractionParams(ForceField_NB::MAX_TYPES* ForceField_NB::MAX_TYPES, NonbondedInteractionParams{0,0});

	const std::vector<AtomType>& activeParameters = ljParameters->GetActiveParameters();
	for (int i = 0; i < activeParameters.size(); i++) {
		for (int j = 0; j < activeParameters.size(); j++) {
			nonbondedInteractionParams[i * ForceField_NB::MAX_TYPES + j] = NonbondedInteractionParams{
				(activeParameters[i].parameters.sigmaHalf + activeParameters[j].parameters.sigmaHalf),
                activeParameters[i].parameters.epsilonSqrt * activeParameters[j].parameters.epsilonSqrt
                    //,activeParameters[i].charge* activeParameters[j].charge
			};
		}
	}
	return nonbondedInteractionParams;
}


int LIMAForcefield::GetActiveTinymoltypeIndex(const std::string& query) {
	return tinymolTypes->GetActiveIndex(query);
}

ForcefieldTinymol LIMAForcefield::GetTinymolTypes() {
	ForcefieldTinymol forcefieldTinymol{};
	const std::vector<AtomType>& activeParameters = tinymolTypes->GetActiveParameters();
	if (activeParameters.size() > ForcefieldTinymol::MAX_TYPES)
		throw std::runtime_error("Too many atom types");
	for (int i = 0; i < activeParameters.size(); i++) {
		const AtomType& at = activeParameters[i];
		forcefieldTinymol.types[i] = ForcefieldTinymol::TinyMolType{ at.parameters.sigmaHalf, at.parameters.epsilonSqrt, at.mass, at.charge };
	}
	return forcefieldTinymol;
}

void LIMAForcefield::LoadFileIntoForcefield(const GenericItpFile& file) 
{
	for (const auto& line : file.GetSection(TopologySection::includes)) {	
		throw std::runtime_error("Includes not supported here");
	}

	for (const auto& line : file.GetSection(TopologySection::atomtypes)) {
		std::istringstream iss(line);
		float sigma, epsilon;
		AtomType atomtype{};
		iss >> atomtype.name >> atomtype.atNum
			>> atomtype.mass				// [g]         // // we take the one from topology, not FF file
			>> atomtype.charge				// [e]
			>> atomtype.ptype
			>> sigma	// [nm]
			>> epsilon;	// [kJ/mol]

		atomtype.charge *= elementaryChargeToKiloCoulombPerMole;
		atomtype.mass /= static_cast<float>(KILO);
		atomtype.parameters.epsilonSqrt = std::sqrt(epsilon * KILO);
		atomtype.parameters.sigmaHalf = sigma * 0.5f;

		ljParameters->insert(atomtype);
		tinymolTypes->insert(atomtype);
	}
	for (const auto& line : file.GetSection(TopologySection::bondtypes)) {
		std::istringstream iss(line);

		SinglebondType bondtype{};
		float b0, kb;
		iss >> bondtype.bonded_typenames[0] >> bondtype.bonded_typenames[1] >> bondtype.func
			>> b0	// [nm]
			>> kb;	// [kJ/mol/nm^2]

		bondtype.params = SingleBond::Parameters::CreateFromCharmm(b0, kb);

		singlebondParameters->insert(bondtype);
	}
	for (const auto& line : file.GetSection(TopologySection::pairtypes)) {
		std::istringstream iss(line);

		PairbondType pairbondtype{};
		float sigma, epsilon;
		iss >> pairbondtype.bonded_typenames[0] >> pairbondtype.bonded_typenames[1] >> pairbondtype.func
			>> sigma	// [nm]
			>> epsilon;	// [kJ/mol]

		pairbondtype.params = PairBond::Parameters::CreateFromCharmm(sigma, epsilon);

		pairbondParameters->insert(pairbondtype);
	}
	for (const auto& line : file.GetSection(TopologySection::angletypes)) {
		std::istringstream iss(line);

		AnglebondType anglebondtype{};
		float t0, kT, ub0, kUB;
		iss >> anglebondtype.bonded_typenames[0] >> anglebondtype.bonded_typenames[1] >> anglebondtype.bonded_typenames[2] >> anglebondtype.func
			>> t0		// [degrees]
			>> kT		// [kJ/mol/rad^2]
			>> ub0		// [nm]
			>> kUB;		// [kJ/mol/nm^2]

		if (anglebondtype.func != 5) {// 5 is for UB, for all other we assume harmonic. TEMP LONG TODO This is a bandaid, we need to parse ALL forcefield types based on the func...
			ub0 = 0;
			kUB = 0;
		}

		anglebondtype.params = AngleUreyBradleyBond::Parameters::CreateFromCharmm(t0, kT, ub0, kUB);

		anglebondParameters->insert(anglebondtype);
	}
	for (const auto& line : file.GetSection(TopologySection::dihedraltypes)) {
		std::istringstream iss(line);

		DihedralbondType dihedralbondtype{};
		float phi0, kphi;
		int n;
		iss >> dihedralbondtype.bonded_typenames[0] >> dihedralbondtype.bonded_typenames[1] >> dihedralbondtype.bonded_typenames[2] >> dihedralbondtype.bonded_typenames[3]
			>> dihedralbondtype.func
			>> phi0	// [degrees]
			>> kphi	// [kJ/mol]
			>> n;

		dihedralbondtype.params = DihedralBond::Parameters::CreateFromCharmm(phi0, kphi, n);

		dihedralbondParameters->insert(dihedralbondtype);
	}
	for (const auto& line : file.GetSection(TopologySection::impropertypes)) {
		std::istringstream iss(line);

		ImproperDihedralbondType improperdihedralbondtype{};
		float psi0, kPsi;
		iss >> improperdihedralbondtype.bonded_typenames[0] >> improperdihedralbondtype.bonded_typenames[1] >> improperdihedralbondtype.bonded_typenames[2] >> improperdihedralbondtype.bonded_typenames[3]
			>> improperdihedralbondtype.func
			>> psi0		// [degrees]
			>> kPsi;	// [2 * kJ/mol]

		improperdihedralbondtype.params = ImproperDihedralBond::Parameters::CreateFromCharmm(psi0, kPsi);

		improperdihedralbondParameters->insert(improperdihedralbondtype);
	}
}

template<typename GenericBond>
const std::vector<typename GenericBond::Parameters>& LIMAForcefield::_GetBondParameters(const auto& query) {
	if constexpr (std::is_same<GenericBond, SingleBond>::value) {
		return singlebondParameters->get(query);
	}
	else if constexpr (std::is_same<GenericBond, PairBond>::value) 
	{
		// Pairbonds are special, as they are only defined for special interactions. In other cases we simply combine the LJ params for the types
		if (pairbondParameters->get(query).empty()) {
			
			const AtomType& left = ljParameters->GetAtomType(query[0]);
			const AtomType& right = ljParameters->GetAtomType(query[1]);

			// TODO: Would prefer to have this computation in a file specialized for it..
			const float sigma = left.parameters.sigmaHalf + right.parameters.sigmaHalf;
			const float epsilon = left.parameters.epsilonSqrt * right.parameters.epsilonSqrt;

			pairbondParameters->insert(PairbondType{
				.bonded_typenames = query,
				.func = 1,
				.params = PairBond::Parameters{sigma, epsilon}
				});
		}
		return pairbondParameters->get(query);
	}
	else if constexpr (std::is_same<GenericBond, AngleUreyBradleyBond>::value) {
		return anglebondParameters->get(query);
	}
	else if constexpr (std::is_same<GenericBond, DihedralBond>::value) {
		return dihedralbondParameters->get(query);
	}
	else if constexpr (std::is_same<GenericBond, ImproperDihedralBond>::value) {
		return improperdihedralbondParameters->get(query);
	}
	else {
		throw std::runtime_error("Unsupported bond type");
	}
}

template<typename GenericBond>
const std::vector<typename GenericBond::Parameters>& LIMAForcefield::GetBondParameters(const auto& query) {
	const auto& parameters = _GetBondParameters<GenericBond>(query);
	if (!parameters.empty())
		return parameters;

	throw std::runtime_error("Failed to find bond parameters");
}
template const std::vector<SingleBond::Parameters>& LIMAForcefield::GetBondParameters<SingleBond>(const std::array<std::string, SingleBond::nAtoms>&);
template const std::vector<PairBond::Parameters>& LIMAForcefield::GetBondParameters<PairBond>(const std::array<std::string, PairBond::nAtoms>&);
template const std::vector<AngleUreyBradleyBond::Parameters>& LIMAForcefield::GetBondParameters<AngleUreyBradleyBond>(const std::array<std::string, AngleUreyBradleyBond::nAtoms>&);
template const std::vector<DihedralBond::Parameters>& LIMAForcefield::GetBondParameters<DihedralBond>(const std::array<std::string, DihedralBond::nAtoms>&);
template const std::vector<ImproperDihedralBond::Parameters>& LIMAForcefield::GetBondParameters<ImproperDihedralBond>(const std::array<std::string, ImproperDihedralBond::nAtoms>&);


