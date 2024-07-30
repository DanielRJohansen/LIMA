#include <vector>

#include "Forcefield.h"
#include "MDFiles.h"



#include <filesystem>
#include <format>

using std::string;
namespace lfs = Filehandler;




const float water_mass = (15.999000f + 2.f * 1.008000f);	// [g]
//const float water_sigma = 1.7398 * rminToSigma * AngToNm;	// Value guessed from param19.inp: OH2      0.0000    -0.0758    1.7398 !ST2   water oxygen
const float water_sigma = 0.22f;	// Made up value. Works better than the one above. I guess i need to implement proper tip3 at some point?
const float water_epsilon = 0.1591f * kcalToJoule;



LjParameterDatabase::LjParameterDatabase() {
	// First add solvent
	insert(AtomType{"solvent", 0, ForceField_NB::ParticleParameters{water_mass * 1e-3f, water_sigma * NANO_TO_LIMA, water_epsilon}, 0.f, 'A'} ); // TODO: Stop doing this, read from the proper file
	GetActiveIndex("solvent");	// Make sure solvent always are active
	assert(GetActiveIndex("solvent") == ATOMTYPE_SOLVENT);
}



void LjParameterDatabase::insert(AtomType element) {
	if (finished)
		throw std::runtime_error("Cannot insert after finishing");

	// Only add a type the first time we see it. TODO: We should also do this for all bonded types
	if (!atomTypes.contains(element.name))
		atomTypes.insert({ element.name, element });
}
int LjParameterDatabase::GetActiveIndex(const std::string& query) {
	if (finished)
		throw std::runtime_error("Cannot query for more active types after finishing");
	if (fastLookup.count(query) != 0)
		return fastLookup.find(query)->second;
	if (!atomTypes.contains(query))
		throw std::runtime_error(std::format("Failed to find atomtype [{}]", query));


	const AtomType& parameter = atomTypes.at({ query });
	// First check of the parameter is already in the active list (since it might be 
	// lexicographically similar but not identical, we cant expect the fastlookup above to always catch)
	for (int i = 0; i < activeAtomTypes.size(); i++) {
		if (parameter.name == activeAtomTypes[i].name) {
			fastLookup.insert({ query, i });
			return i;
		}
	}

	activeAtomTypes.push_back(parameter);
	fastLookup.insert({ query, activeAtomTypes.size() - 1 });
	return activeAtomTypes.size() - 1;
}


LIMAForcefield::LIMAForcefield() {
	bool verbose = false;
	const char ignore_atomtype = IGNORE_HYDROGEN ? 'H' : '.';

#ifdef __linux__
	const fs::path ff_dir = "/opt/LIMA/resources/Forcefields/charm27.ff";
#else
	const fs::path ff_dir = "C:/Users/Daniel/git_repo/LIMA/resources/Forcefields/charmm27.ff";
#endif

	LoadFileIntoForcefield(ff_dir / "../lima_custom_forcefield.itp");
	LoadFileIntoForcefield(ff_dir / "ffnonbonded.itp");
	LoadFileIntoForcefield(ff_dir / "ffbonded.itp");
	/*LoadFileIntoForcefield(ff_dir / "ffnanonbonded.itp");
	LoadFileIntoForcefield(ff_dir / "ffnabonded.itp");
	LoadFileIntoForcefield(ff_dir / "ions.itp");*/

}













namespace ForcefieldHelpers {
	float _calcLikeness(const std::string& query_type, const std::string& forcefield_type) {

		// Edgecase: perfect match
		if (query_type == forcefield_type) { return 1.f; }

		// Edgecase: wildcard
		if (forcefield_type == "X") {
			return 0.9f;
		}

		float likeness = 0;
		float point_scale = 1.f / std::max(query_type.length(), forcefield_type.length());

		for (size_t i = 0; i < std::min(query_type.length(), forcefield_type.length()); i++) {
			if (query_type[i] == forcefield_type[i])
				likeness += point_scale;
			else
				break;
		}
		return likeness;
	}


	static float calcLikeness(std::span<const std::string> query,
		std::span<const std::string> typeInForcefield) {
		float likeness = 1.f;
		for (int i = 0; i < query.size(); i++) {
			likeness *= _calcLikeness(query[i], typeInForcefield[i]);
		}

		return likeness;
	}
};


void LIMAForcefield::LoadFileIntoForcefield(const fs::path& path) {

	GenericItpFile file{ path };


	if (file.GetSection(TopologySection::atomtypes)) {
		for (const auto& line : file.GetSection(TopologySection::atomtypes)->get().lines) {
			std::istringstream iss(line);

			AtomType atomtype{};
			iss >> atomtype.name >> atomtype.atNum
				>> atomtype.parameters.mass		// [g]
				>> atomtype.charge				// [e]
				>> atomtype.ptype
				>> atomtype.parameters.sigma	// [nm]
				>> atomtype.parameters.epsilon;	// [kJ/mol]

			atomtype.charge *= elementaryChargeToKiloCoulombPerMole;

			atomtype.parameters.mass /= static_cast<float>(KILO);
			atomtype.parameters.sigma *= NANO_TO_LIMA;
			atomtype.parameters.epsilon *= KILO;

			ljParameters.insert(atomtype);
		}
	}
	if (file.GetSection(TopologySection::bondtypes)) {
		for (const auto& line : file.GetSection(TopologySection::bondtypes)->get().lines) {
			std::istringstream iss(line);

			SinglebondType bondtype{};
			iss >> bondtype.bonded_typenames[0] >> bondtype.bonded_typenames[1] >> bondtype.func
				>> bondtype.params.b0	// [nm]
				>> bondtype.params.kb;	// [kJ/mol/nm^2]

			bondtype.params.b0 *= NANO_TO_LIMA;
			bondtype.params.kb *= 1. / NANO_TO_LIMA / NANO_TO_LIMA * KILO; // Convert to J/mol

			//SortBondedtypeNames<SingleBond>(bondtype.bonded_typenames);

			singlebondParameters.insert(bondtype);
		}
	}
	if (file.GetSection(TopologySection::angletypes)) {
		for (const auto& line : file.GetSection(TopologySection::angletypes)->get().lines) {
			std::istringstream iss(line);

			AnglebondType anglebondtype{};
			iss >> anglebondtype.bonded_typenames[0] >> anglebondtype.bonded_typenames[1] >> anglebondtype.bonded_typenames[2] >> anglebondtype.func
				>> anglebondtype.params.theta_0		// [degrees]
				>> anglebondtype.params.k_theta;	// [kJ/mol/rad^2]

			anglebondtype.params.theta_0 *= DEG_TO_RAD;
			anglebondtype.params.k_theta *= KILO; // Convert to J/mol/rad^2
			//SortBondedtypeNames<AngleBond>(anglebondtype.bonded_typenames);

			anglebondParameters.insert(anglebondtype);
		}
	}
	if (file.GetSection(TopologySection::dihedraltypes)) {
		for (const auto& line : file.GetSection(TopologySection::dihedraltypes)->get().lines) {
			std::istringstream iss(line);

			DihedralbondType dihedralbondtype{};
			float phi0;
			float kphi;
			float n;
			iss >> dihedralbondtype.bonded_typenames[0] >> dihedralbondtype.bonded_typenames[1] >> dihedralbondtype.bonded_typenames[2] >> dihedralbondtype.bonded_typenames[3]
				>> dihedralbondtype.func
				>> phi0	// [degrees]
				>> kphi	// [kJ/mol]
				>> n;

			dihedralbondtype.params.phi_0 = phi0 * DEG_TO_RAD;
			dihedralbondtype.params.k_phi = kphi * KILO / 2.f; // Convert to J/mol// Convert to J/mol TODO: Move the /2 from here to the force calculation to save an op there
			dihedralbondtype.params.n = n;

			dihedralbondParameters.insert(dihedralbondtype);
		}
	}
	if (file.GetSection(TopologySection::impropertypes)) {
		for (const auto& line : file.GetSection(TopologySection::impropertypes)->get().lines) {
			std::istringstream iss(line);

			ImproperDihedralbondType improperdihedralbondtype{};
			float psi0;	// [degrees]
			float kpsi;	// [kJ/mol]
			iss >> improperdihedralbondtype.bonded_typenames[0] >> improperdihedralbondtype.bonded_typenames[1] >> improperdihedralbondtype.bonded_typenames[2] >> improperdihedralbondtype.bonded_typenames[3]
				>> improperdihedralbondtype.func
				>> psi0		// [degrees]
				>> kpsi;	// [2 * kJ/mol]

			improperdihedralbondtype.params.psi_0 = psi0 * DEG_TO_RAD;
			improperdihedralbondtype.params.k_psi = kpsi * KILO; 
			//SortBondedtypeNames<ImproperDihedralBond>(improperdihedralbondtype.bonded_typenames);
			improperdihedralbondParameters.insert(improperdihedralbondtype);
		}
	}
}

template <typename GenericBondType>
int ParameterDatabase<GenericBondType>::findBestMatchInForcefield_Index(const std::array<std::string, GenericBondType::nAtoms>& query) {
	if (parameters.size() == 0) {
		throw std::runtime_error("No bonds in forcefield!");
	}

	float best_likeness = 0;
	int bestBondIndex = 0;
	for (size_t i = 0; i < parameters.size(); i++) {
		const GenericBondType& ff_bondtype = parameters[i];
		const float likeness = ForcefieldHelpers::calcLikeness(query, parameters[i].bonded_typenames);

		if (likeness > best_likeness) {
			best_likeness = likeness;
			bestBondIndex = static_cast<int>(i);
		}
	}

	if (best_likeness > 0.001f) {
		return bestBondIndex;
	}

	std::cout << "Failed to match bond types.\n Closest match ";
	for (const auto& name : parameters[bestBondIndex].bonded_typenames) {
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
	for (const auto& name : query) {
		std::cout << name << " ";
	}

	throw std::runtime_error("\nfindBestMatchInForcefield failed");
}
template class ParameterDatabase<SinglebondType>;
template class ParameterDatabase<AnglebondType>;
template class ParameterDatabase<DihedralbondType>;
template class ParameterDatabase<ImproperDihedralbondType>;