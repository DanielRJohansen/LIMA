#include <vector>
#include "Format.h"

#include "Forcefield.h"
#include "MDFiles.h"

#include <format>
#include <tuple>
#include <algorithm>

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

	std::optional<AtomType> GetAtomType(const std::string& query) {
		if (!atomTypes.contains(query))
			return std::nullopt;			
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
		throw std::runtime_error(Lima::Format("Failed to find atomtype [{}]", query));

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
	fastLookup.insert({ query, static_cast<int>(activeAtomTypes->size()) - 1 });
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
		// Only entries with identical typenames can be duplicates, so index by typenames to keep loading linear
		std::vector<size_t>& sameNames = indicesByTypenames[JoinTypenames(element.bonded_typenames)];
		for (const size_t index : sameNames) {
			if (parameters[index].params == element.params) {
				// I think this will be a problem, either because we include multiple files for the forcefiel, either 
				// custom or there are duplicates in the ffnabonded.itp, or it may just be a problem with the forcefield itself
				//throw std::runtime_error("Duplicate bond type found in forcefield");


				// This does happen in Slipids, then duplicates denoted with  ;*new...... So lets just ignore the problem, and ignore the duplicates
				return;
			}
		}
		sameNames.push_back(parameters.size());
		parameters.push_back(element);
	}

	const std::vector<typename GenericBondType::Parameters> emptyParams{}; // Return this if we find no matches
	std::vector<GenericBondType> parameters;
	std::unordered_map<std::string, std::vector<typename GenericBondType::Parameters>> fastLookup;	// Map a query to an index in parameters
	std::unordered_map<std::string, std::vector<size_t>> indicesByTypenames;	// Indices in parameters, for duplicate detection

	static std::string JoinTypenames(const std::array<std::string, GenericBondType::nAtoms>& names) {
		std::string key;
		for (const auto& name : names) {
			key += name;
			key += ' ';
		}
		return key;
	}

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














namespace {
	// GROMACS' cmap spline setup (grompp setup_cmap, spline1d and interpolate1d), in degrees and kJ/mol like the itp
	void Spline1d(double dx, const double* y, int n, double* u, double* y2) {
		y2[0] = 0.0;
		u[0] = 0.0;
		for (int i = 1; i < n - 1; i++) {
			const double p = 0.5 * y2[i - 1] + 2.0;
			y2[i] = -0.5 / p;
			const double q = (y[i + 1] - 2.0 * y[i] + y[i - 1]) / dx;
			u[i] = (3.0 * q / dx - 0.5 * u[i - 1]) / p;
		}
		y2[n - 1] = 0.0;
		for (int i = n - 2; i >= 0; i--)
			y2[i] = y2[i] * y2[i + 1] + u[i];
	}

	void Interpolate1d(double xmin, double dx, const double* ya, const double* y2a, double x, double& y, double& y1) {
		const int ix = static_cast<int>((x - xmin) / dx);
		const double a = (xmin + (ix + 1) * dx - x) / dx;
		const double b = (x - xmin - ix * dx) / dx;
		y = a * ya[ix] + b * ya[ix + 1] + ((a * a * a - a) * y2a[ix] + (b * b * b - b) * y2a[ix + 1]) * (dx * dx) / 6.0;
		y1 = (ya[ix + 1] - ya[ix]) / dx - (3.0 * a * a - 1.0) / 6.0 * dx * y2a[ix] + (3.0 * b * b - 1.0) / 6.0 * dx * y2a[ix + 1];
	}

	// energies: n*n values [kJ/mol], phi-major, starting at -180 degrees
	CmapGrid MakeCmapGrid(int n, const std::vector<double>& energies) {
		const int n2 = 2 * n;
		std::vector<double> tmpU(n2), tmpU2(n2), tmpYY(n2), tmpY1(n2), tmpT2(n2 * n2), tmpGrid(n2 * n2);
		const double dx = 360.0 / n;
		const double xmin = -180.0 - dx * n / 2;

		// The grid extended by half a period on each side, so the splines are periodic within the original range
		for (int i = 0; i < n2; i++) {
			const int ii = (i + n - n / 2) % n;
			for (int j = 0; j < n2; j++) {
				const int jj = (j + n - n / 2) % n;
				tmpGrid[i * n2 + j] = energies[ii * n + jj];
			}
		}
		for (int i = 0; i < n2; i++)
			Spline1d(dx, &tmpGrid[n2 * i], n2, tmpU.data(), &tmpT2[n2 * i]);

		CmapGrid grid{ n, std::vector<float4>(n * n) };
		for (int i = n / 2; i < n + n / 2; i++) {
			const int ii = i - n / 2;
			const double phi = ii * dx - 180.0;
			for (int j = n / 2; j < n + n / 2; j++) {
				const int jj = j - n / 2;
				const double psi = jj * dx - 180.0;
				for (int k = 0; k < n2; k++)
					Interpolate1d(xmin, dx, &tmpGrid[n2 * k], &tmpT2[n2 * k], psi, tmpYY[k], tmpY1[k]);
				double v, v1, v2, v12;
				Spline1d(dx, tmpYY.data(), n2, tmpU.data(), tmpU2.data());
				Interpolate1d(xmin, dx, tmpYY.data(), tmpU2.data(), phi, v, v1);
				Spline1d(dx, tmpY1.data(), n2, tmpU.data(), tmpU2.data());
				Interpolate1d(xmin, dx, tmpY1.data(), tmpU2.data(), phi, v2, v12);

				// Derivatives per degree, scaled to a grid cell, and kJ -> J
				grid.points[ii * n + jj] = float4{ static_cast<float>(energies[ii * n + jj] * KILO), static_cast<float>(v1 * dx * KILO),
					static_cast<float>(v2 * dx * KILO), static_cast<float>(v12 * dx * dx * KILO) };
			}
		}
		return grid;
	}

	std::string JoinCmapTypenames(const std::array<std::string, CmapBond::nAtoms>& names) {
		std::string key;
		for (const auto& name : names) {
			key += name;
			key += ' ';
		}
		return key;
	}
}

LIMAForcefield::LIMAForcefield() {};

std::optional<int> LIMAForcefield::GetCmapGridIndex(const std::array<std::string, CmapBond::nAtoms>& atomtypes) const {
	const auto found = cmapGridIndices.find(JoinCmapTypenames(atomtypes));
	if (found == cmapGridIndices.end())
		return std::nullopt;
	return found->second;
}

LIMAForcefield::LIMAForcefield(const GenericItpFile& file) {

	ljParameters = std::make_unique<AtomtypeDatabase>();

	singlebondParameters = std::make_unique<ParameterDatabase<SinglebondType>>();
	pairbondParameters = std::make_unique<ParameterDatabase<PairbondType>>();
	anglebondParameters = std::make_unique<ParameterDatabase<AnglebondType>>();
	dihedralbondParameters = std::make_unique<ParameterDatabase<DihedralbondType>>();
	improperdihedralbondParameters = std::make_unique<ParameterDatabase<ImproperDihedralbondType>>();

	LoadFileIntoForcefield(file);
}

LIMAForcefield::~LIMAForcefield() {}

int LIMAForcefield::GetActiveLjParameterIndex(const std::string& query) {
	return ljParameters->GetActiveIndex(query);
}
NBParams LIMAForcefield::GetLjParameters(const std::string& query) const {
	auto params = ljParameters->GetAtomType(query);

	if (!params.has_value())
		throw std::runtime_error(Lima::Format("Failed to find atomtype [{}]", query));

	NBParams nbparams{};
	nbparams.sigmaHalf = params->parameters.sigmaHalf;
	nbparams.epsilonSqrt = params->parameters.epsilonSqrt;
	nbparams.charge = params->charge;
	return nbparams;
}

std::optional<AtomType> LIMAForcefield::GetAtomtype(const std::string& query) const {
	return ljParameters->GetAtomType(query);
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
                    ,activeParameters[i].charge* activeParameters[j].charge
			};
		}
	}
	return nonbondedInteractionParams;
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

		anglebondtype.params = AngleUreyBradleyBond::Parameters::CreateFromCharmm(t0, kT, ub0, kUB, anglebondtype.func);

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
	{
		// An entry spans several lines joined by '\': 5 atomtypes, funct, the grid size in phi and psi, then the energies [kJ/mol]
		std::string entries;
		for (const auto& line : file.GetSection(TopologySection::cmaptypes)) {
			std::string data = line.substr(0, line.find(';'));
			std::replace(data.begin(), data.end(), '\\', ' ');
			entries += data;
			entries += ' ';
		}
		std::istringstream iss(entries);
		std::array<std::string, CmapBond::nAtoms> names;
		int func = 0, nPhi = 0, nPsi = 0;
		while (iss >> names[0] >> names[1] >> names[2] >> names[3] >> names[4] >> func >> nPhi >> nPsi) {
			if (nPhi != nPsi || nPhi < 4)
				throw std::runtime_error(Lima::Format("Unsupported cmaptypes grid of {}x{} points", nPhi, nPsi));
			std::vector<double> energies(nPhi * nPsi);
			for (double& energy : energies)
				if (!(iss >> energy))
					throw std::runtime_error("Incomplete cmaptypes grid");
			// Like other types, the first definition wins
			if (cmapGridIndices.try_emplace(JoinCmapTypenames(names), static_cast<int>(cmapGrids.size())).second)
				cmapGrids.push_back(MakeCmapGrid(nPhi, energies));
		}
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
			
			const AtomType left = *ljParameters->GetAtomType(query[0]);
			const AtomType right = *ljParameters->GetAtomType(query[1]);

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


