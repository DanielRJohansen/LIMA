#include "Printer.h"
#include "Forcefield.cuh"
#include <filesystem>
//#include "LIMA_ENGINE/include/EngineUtils.cuh"


using namespace LIMA_Print;

Forcefield::Forcefield(VerbosityLevel vl, const std::string& molecule_dir) : vl(vl) {
	if (vl >= CRITICAL_INFO) { printH2("Building forcefield"); }


	std::string nonbonded_path = Filehandler::fileExists(Filehandler::pathJoin(molecule_dir, "custom_ffnonbonded.lff"))
		? Filehandler::pathJoin(molecule_dir, "custom_ffnonbonded.lff")
		: Filehandler::pathJoin(molecule_dir, "ffnonbonded.lff");


	const SimpleParsedFile nonbonded_parsed = Filehandler::parseLffFile(nonbonded_path, vl >= V1);

	// First load the nb forcefield
	auto atomtypes = loadAtomTypes(nonbonded_parsed);					// 1 entry per type in compressed forcefield
	loadAtomypesIntoForcefield(atomtypes);

	// Find mappings between the atoms in the simulation and the nb forcefield
	//globaldToAtomtypeMap = loadAtomTypeMap(nonbonded_parsed);	// 1 entry per atom in conf

	forcefield_loaded = true;

	if (vl >= CRITICAL_INFO) {
		printf("Nonbonded parameters size: %llu bytes\n", sizeof(ForceField_NB));
		printH2("Finished building forcefield");
	}
};

template <int n>
bool isMatch(const uint32_t* topolbonds, const std::array<int, n> query_ids) {
	for (int i = 0; i < query_ids.size(); i++) {
		if (topolbonds[i] != query_ids[i]) {
			return false;
		}
	}
	return true;
}

std::vector<NBAtomtype> Forcefield::loadAtomTypes(const SimpleParsedFile& parsedfile) {
	std::vector<NBAtomtype> atomtypes;
	atomtypes.reserve(200);

	for (auto& row : parsedfile.rows) {
		if (row.section == "atomtypes") {
			// Row is type, id, weight [g], sigma [nm], epsilon [J/mol]
			float mass = stof(row.words[2]);
			float sigma = stof(row.words[3]);
			float epsilon = stof(row.words[4]);
			atomtypes.emplace_back(NBAtomtype(mass, sigma, epsilon));
		}
	}

	assert(atomtypes.size() < MAX_ATOM_TYPES);
	if (vl >= V1) { printf("%zu NB_Atomtypes loaded\n", atomtypes.size()); }
	return atomtypes;
}

void Forcefield::loadAtomypesIntoForcefield(const std::vector<NBAtomtype>& atomtypes) {
	static const float mass_min = 0.001f;	// [kg/mol]
	static const float sigma_min = 0.001f;
	static const float epsilon_min = 0.001f;

	for (int i = 0; i < atomtypes.size(); i++) {
		forcefield_nb.particle_parameters[i].mass = atomtypes[i].mass * 1e-3f;				// Convert g/mol to kg/mol
		forcefield_nb.particle_parameters[i].sigma = atomtypes[i].sigma * NANO_TO_LIMA;		// Convert from [nm] to [lm]
		forcefield_nb.particle_parameters[i].epsilon = atomtypes[i].epsilon;				// Interpreted as kg*lm^2/ls^2 

		bool illegal_parameter = (forcefield_nb.particle_parameters[i].mass < mass_min) || (forcefield_nb.particle_parameters[i].sigma < sigma_min) || (forcefield_nb.particle_parameters[i].epsilon < epsilon_min);

		if ((vl >= V2) || illegal_parameter) { printf("Mass %f Sigma %f Epsilon %f\n", atomtypes[i].mass, atomtypes[i].sigma, atomtypes[i].epsilon); }
	}

}
