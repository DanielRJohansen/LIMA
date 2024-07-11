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
	insert(LJParameter{ {"solvent"}, ForceField_NB::ParticleParameters{water_mass * 1e-3f, water_sigma * NANO_TO_LIMA, water_epsilon} });
	GetActiveIndex("solvent");	// Make sure solvent always are active
	assert(GetActiveIndex("solvent") == ATOMTYPE_SOLVENT);
}



void LjParameterDatabase::insert(LJParameter element) {
	if (finished)
		throw std::runtime_error("Cannot insert after finishing");
	ljParameters.insert(element);
}
int LjParameterDatabase::GetActiveIndex(const std::string& query) {
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




void LIMAForcefield::loadFileIntoForcefield(const SimpleParsedFile& parsedfile) {
	std::unordered_map<std::string, float> atomnameToMassMap;

	for (const SimpleParsedFile::Row& row : parsedfile.rows) {

		if (row.section == "ATOMS") {
			assert(row.words.size() >= 4);

			const string& atomtype = row.words[2];
			const float mass = stof(row.words[3]);		// Should this come from topol too?

			atomnameToMassMap.insert(std::pair{ atomtype, mass });

			//forcefield.nb_atoms.emplace_back(NB_Atomtype(atomtype, atnum, mass, sigma, epsilon));
			//forcefield.atomToTypeMap.insert(std::pair(atomtype, NB_Atomtype(atomtype, atnum, mass, sigma, epsilon)));
		}
		if (row.section == "pairtypes") {
			// TODO: Fill out this logic
		}
		else if (row.section == "BONDS") {
			assert(row.words.size() >= 4);

			const std::array<string, 2> bondedatom_typenames{ row.words[0], row.words[1] };
			const float kb = stof(row.words[2]) * kcalToJoule / AngToNm / AngToNm * 2.f;		// * 2 to go from V(bond) = Kb(b - b0)**2 -> V(bond) = 0.5*Kb(b - b0)**2
			const float b0 = stof(row.words[3]) * AngToNm;	// A to nm

			singlebondParameters.insert(Singlebondtype(bondedatom_typenames, b0, kb));
		}
		else if (row.section == "ANGLES") {
			assert(row.words.size() >= 5);

			const std::array<string, 3> angle_typenames{ row.words[0], row.words[1], row.words[2] };
			const float ktheta = stof(row.words[3]) * kcalToJoule * 2.f;	// * 2 to convert Ktheta(Theta - Theta0)**2 -> 0.5*Ktheta(Theta - Theta0)**2
			const float theta0 = stof(row.words[4]) * degreeToRad;

			anglebondParameters.insert(Anglebondtype(angle_typenames, theta0, ktheta));
		}
		else if (row.section == "DIHEDRALS") {
			assert(row.words.size() >= 7);

			const std::array<string, 4> dihedral_typenames{ row.words[0], row.words[1], row.words[2], row.words[3] };
			const float kphi = stof(row.words[4]) * kcalToJoule;
			const int n = stoi(row.words[5]);
			const float phi0 = stof(row.words[6]) * degreeToRad;

			dihedralbondParameters.insert(Dihedralbondtype(dihedral_typenames, phi0, kphi, n));
		}
		else if (row.section == "IMPROPER") {
			assert(row.words.size() >= 7);

			const std::array<string, 4> improper_dihedral_typenames{ row.words[0], row.words[1], row.words[2], row.words[3] };
			const float kpsi = stof(row.words[4]) * kcalToJoule * 2.f;	// * 2 to convert Kpsi(psi - psi0)**2 -> 0.5*Kpsi(psi - psi0)**2
			const float psi0 = stof(row.words[6]) * degreeToRad;

			improperdihedralbondParameters.insert(Improperdihedralbondtype(improper_dihedral_typenames, psi0, kpsi));
		}
		else if (row.section == "NONBONDED") {
			assert(row.words.size() >= 4);

			const string& atomtype = row.words[0];
			const float epsilon = abs(stof(row.words[2]) * kcalToJoule);	// For some fucked reason the eps is *inconsistently* negative...
			const float sigma = stof(row.words[3]) * rminToSigma * AngToNm;	// rmin/2 [A] -> sigma [lm] 
			if (!atomnameToMassMap.count(atomtype))
				throw std::runtime_error(std::format("Failed to find atomtype [{}]", atomtype));
			const float mass = atomnameToMassMap.at(atomtype) * 1e-3f;		// [g] -> [kg]

			LJParameter entry{ { atomtype }, ForceField_NB::ParticleParameters{mass, sigma * NANO_TO_LIMA, epsilon} };

			ljParameters.insert(entry);

			// Not yet used
			//const float epsilon_1_4 = stof(row.words[6]) * 2;	 // rmin/2 -> sigma
			//const float sigma_1_4 = stof(row.words[6]) * 2;	 // rmin/2 -> sigma			
		}
	}
}

LIMAForcefield::LIMAForcefield() {
	bool verbose = false;
	const char ignore_atomtype = IGNORE_HYDROGEN ? 'H' : '.';
	// Check if filtered files already exists, if so return

	std::vector<std::string> files = getFiles();

	for (auto& file_path : files) {
		const SimpleParsedFile parsedfile = lfs::parsePrmFile(file_path, verbose);
		loadFileIntoForcefield(parsedfile);
	}
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