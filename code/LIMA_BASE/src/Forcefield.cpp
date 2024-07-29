#include <vector>

#include "Filehandling.h"
#include "ForcefieldTypes.h"
#include "Forcefield.h"
#include "MDFiles.h"



#include <filesystem>
#include <assert.h>
#include <format>
#include <unordered_set>

using std::string;
namespace lfs = Filehandler;




const float water_mass = (15.999000f + 2.f * 1.008000f);	// [g]
//const float water_sigma = 1.7398 * rminToSigma * AngToNm;	// Value guessed from param19.inp: OH2      0.0000    -0.0758    1.7398 !ST2   water oxygen
const float water_sigma = 0.22f;	// Made up value. Works better than the one above. I guess i need to implement proper tip3 at some point?
const float water_epsilon = 0.1591f * kcalToJoule;







std::vector<std::string> getFiles() {
	std::vector<std::string> files;

	// Some files are commented out because it is NOT clear whether they are using rmin or rmin/2
#ifdef __linux__
	const std::string ff_dir = "/opt/LIMA/resources/Forcefields/charmm36-mar2019.ff";
#else
	const std::string ff_dir = "C:/Users/Daniel/git_repo/LIMA/resources/Forcefields/charmm36-mar2019.ff";
#endif

	files.push_back(ff_dir + "/par_all35_ethers.prm");
	files.push_back(ff_dir + "/par_all36_carb.prm");
	files.push_back(ff_dir + "/par_all36_lipid.prm");
	files.push_back(ff_dir + "/par_all36_na.prm");	
	files.push_back(ff_dir + "/par_all36m_prot.prm");
	//files.push_back(ff_dir + "/par_all36_cgenff.prm");
	files.push_back(ff_dir + "/par_all22_prot.prm");
	files.push_back(ff_dir + "/../lima_custom_forcefield.prm");

	return files;
}




LjParameterDatabase::LjParameterDatabase() {
	// First add solvent
	insert(AtomType{"solvent", 0, ForceField_NB::ParticleParameters{water_mass * 1e-3f, water_sigma * NANO_TO_LIMA, water_epsilon}, 0.f, 'A'} );
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



//
//void LIMAForcefield::loadFileIntoForcefield(const SimpleParsedFile& parsedfile) {
//	std::unordered_map<std::string, float> atomnameToMassMap;
//
//	for (const SimpleParsedFile::Row& row : parsedfile.rows) {
//
//		if (row.section == "ATOMS") {
//			assert(row.words.size() >= 4);
//
//			const string& atomtype = row.words[2];
//			const float mass = stof(row.words[3]);		// Should this come from topol too?
//
//			atomnameToMassMap.insert(std::pair{ atomtype, mass });
//
//			//forcefield.nb_atoms.emplace_back(NB_Atomtype(atomtype, atnum, mass, sigma, epsilon));
//			//forcefield.atomToTypeMap.insert(std::pair(atomtype, NB_Atomtype(atomtype, atnum, mass, sigma, epsilon)));
//		}
//		if (row.section == "pairtypes") {
//			// TODO: Fill out this logic
//		}
//		else if (row.section == "BONDS") {
//			assert(row.words.size() >= 4);
//
//			const std::array<string, 2> bondedatom_typenames{ row.words[0], row.words[1] };
//			const float kb = stof(row.words[2]) * kcalToJoule / AngToNm / AngToNm * 2.f / (NANO_TO_LIMA* NANO_TO_LIMA);		// * 2 to go from V(bond) = Kb(b - b0)**2 -> V(bond) = 0.5*Kb(b - b0)**2
//			const float b0 = stof(row.words[3]) * AngToNm * NANO_TO_LIMA;	// A to nm
//
//			auto bond = SinglebondType{ bondedatom_typenames, -1,  {b0, kb} };
//
//			//auto bond = Singlebondtype(bondedatom_typenames, b0, kb);
//			SortBondedtypeNames<SingleBond>(bond.bonded_typenames);
//			singlebondParameters.insert(bond);
//		}
//		else if (row.section == "ANGLES") {
//			assert(row.words.size() >= 5);
//
//			const std::array<string, 3> angle_typenames{ row.words[0], row.words[1], row.words[2] };
//			const float ktheta = stof(row.words[3]) * kcalToJoule * 2.f;	// * 2 to convert Ktheta(Theta - Theta0)**2 -> 0.5*Ktheta(Theta - Theta0)**2
//			const float theta0 = stof(row.words[4]) * degreeToRad;
//
//			//auto bond = Anglebondtype(angle_typenames, theta0, ktheta);
//			auto bond = AnglebondType{ angle_typenames, -1, {theta0, ktheta}, 0.f, 0.f };
//
//			//SortBondedtypeNames<AngleBond>(bond.bonded_typenames);
//			anglebondParameters.insert(bond);
//		}
//		else if (row.section == "DIHEDRALS") {
//			assert(row.words.size() >= 7);
//
//			const std::array<string, 4> dihedral_typenames{ row.words[0], row.words[1], row.words[2], row.words[3] };
//			const float kphi = stof(row.words[4]) * kcalToJoule;
//			const int n = stoi(row.words[5]);
//			const float phi0 = stof(row.words[6]) * degreeToRad;
//
//			//auto bond = Dihedralbondtype(dihedral_typenames, phi0, kphi, n);
//			auto bond = DihedralbondType{ dihedral_typenames, -1, {phi0, kphi, n} };
//			//SortBondedtypeNames<DihedralBond>(bond.bonded_typenames);
//			dihedralbondParameters.insert(bond);
//		}
//		else if (row.section == "IMPROPER") {
//			assert(row.words.size() >= 7);
//
//			const std::array<string, 4> improper_dihedral_typenames{ row.words[0], row.words[1], row.words[2], row.words[3] };
//			const float kpsi = stof(row.words[4]) * kcalToJoule * 2.f;	// * 2 to convert Kpsi(psi - psi0)**2 -> 0.5*Kpsi(psi - psi0)**2
//			const float psi0 = stof(row.words[6]) * degreeToRad;
//			//auto bond = Improperdihedralbondtype(improper_dihedral_typenames, psi0, kpsi);
//			auto bond = ImproperDihedralbondType{ improper_dihedral_typenames, -1, {psi0, kpsi} };
//			//SortBondedtypeNames<ImproperDihedralBond>(bond.bonded_typenames);
//			improperdihedralbondParameters.insert(bond);
//		}
//		else if (row.section == "NONBONDED") {
//			assert(row.words.size() >= 4);
//
//			const string& atomTypeName = row.words[0];
//			const float epsilon = abs(stof(row.words[2])) * kcalToJoule;	// For some fucked reason the eps is *inconsistently* negative...
//			const float sigma = stof(row.words[3]) * 2.f * rminToSigma * AngToNm;	// rmin/2 [A] -> sigma [lm] 
//			const float sigmaOld = stof(row.words[3]) / rminToSigma * AngToNm;	// rmin/2 [A] -> sigma [lm]  // WRONG
//			if (!atomnameToMassMap.count(atomTypeName))
//				throw std::runtime_error(std::format("Failed to find atomtype [{}]", atomTypeName));
//			const float mass = atomnameToMassMap.at(atomTypeName) / static_cast<float>(KILO);		// [g/mol] -> [kg/mol]
//
//			//LJParameter entry{ { atomtype }, ForceField_NB::ParticleParameters{mass, sigma * NANO_TO_LIMA, epsilon} };
//			//ljParameters.insert(entry);
//
//			AtomType atomtype{ atomTypeName , 0, ForceField_NB::ParticleParameters{mass, sigma * NANO_TO_LIMA, epsilon} , 0.f, 'A' };
//			ljParameters.insert(atomtype);
//			
//
//			// Not yet used
//			//const float epsilon_1_4 = stof(row.words[6]) * 2;	 // rmin/2 -> sigma
//			//const float sigma_1_4 = stof(row.words[6]) * 2;	 // rmin/2 -> sigma			
//		}
//	}
//}

LIMAForcefield::LIMAForcefield() {
	bool verbose = false;
	const char ignore_atomtype = IGNORE_HYDROGEN ? 'H' : '.';
	// Check if filtered files already exists, if so return

	//std::vector<std::string> files = getFiles();

	//for (auto& file_path : files) {
	//	const SimpleParsedFile parsedfile = lfs::parsePrmFile(file_path, verbose);
	//	loadFileIntoForcefield(parsedfile);
	//}



#ifdef __linux__
	const fs::path ff_dir = "/opt/LIMA/resources/Forcefields/charm27.ff";
#else
	const fs::path ff_dir = "C:/Users/Daniel/git_repo/LIMA/resources/Forcefields/charmm27.ff";
#endif

	LoadFileIntoForcefield(ff_dir / "ffnonbonded.itp");
	LoadFileIntoForcefield(ff_dir / "ffbonded.itp");
	/*LoadFileIntoForcefield(ff_dir / "ffnanonbonded.itp");
	LoadFileIntoForcefield(ff_dir / "ffnabonded.itp");
	LoadFileIntoForcefield(ff_dir / "ions.itp");*/

}
















float ForcefieldHelpers::_calcLikeness(const string& query_type, const string& forcefield_type) {

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
			iss >> anglebondtype.bonded_typenames[0] >> anglebondtype.bonded_typenames[1] >> anglebondtype.bonded_typenames[2]
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
			dihedralbondtype.params.k_phi = kphi * KILO; // Convert to J/mol
			dihedralbondtype.params.n = n;
			//SortBondedtypeNames<DihedralBond>(dihedralbondtype.bonded_typenames);

			dihedralbondParameters.insert(dihedralbondtype);
		}
	}
	if (file.GetSection(TopologySection::impropertypes)) {
		for (const auto& line : file.GetSection(TopologySection::impropertypes)->get().lines) {
			std::istringstream iss(line);

			ImproperDihedralbondType improperdihedralbondtype{};
			float psi0;
			float kpsi;
			iss >> improperdihedralbondtype.bonded_typenames[0] >> improperdihedralbondtype.bonded_typenames[1] >> improperdihedralbondtype.bonded_typenames[2] >> improperdihedralbondtype.bonded_typenames[3]
				>> psi0		// [degrees]
				>> kpsi;	// [2 * kJ/mol]

			improperdihedralbondtype.params.psi_0 = psi0 * DEG_TO_RAD;
			improperdihedralbondtype.params.k_psi = kpsi * KILO / 2.f; // Convert to J/mol TODO: Move the /2 from here to the force calculation to save an op there
			//SortBondedtypeNames<ImproperDihedralBond>(improperdihedralbondtype.bonded_typenames);
			improperdihedralbondParameters.insert(improperdihedralbondtype);
		}
	}
}