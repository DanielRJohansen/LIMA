#pragma once

#include "Printer.h"
#include "Forcefield.cuh"
//#include "LIMA_ENGINE/include/EngineUtils.cuh"


using namespace LIMA_Print;

const int min_reserve_size = 10000;	// This should NOT be permanent...

Forcefield::Forcefield(VerbosityLevel vl) : vl(vl) {};


void Forcefield::loadForcefield(string molecule_dir) {
	if (vl >= CRITICAL_INFO) { printH2("Building forcefield"); }

	SimpleParsedFile nonbonded_parsed = Filehandler::parseLffFile(Filehandler::pathJoin(molecule_dir, "ffnonbonded_filtered.lff"), vl>= V1);
	SimpleParsedFile bonded_parsed = Filehandler::parseLffFile(Filehandler::pathJoin(molecule_dir, "ffbonded_filtered.lff"), vl >= V1);

	// First load the nb forcefield
	auto atomtypes = loadAtomTypes(nonbonded_parsed);					// 1 entry per type in compressed forcefield
	loadAtomypesIntoForcefield(atomtypes);

	// Find mappings between the atoms in the simulation and the nb forcefield
	groIdToAtomtypeMap = loadAtomTypeMap(nonbonded_parsed);	// 1 entry per atom in conf

	// Load the topology with their included forcefield parameters
	topology = loadTopology(bonded_parsed);

	forcefield_loaded = true;

	if (vl >= CRITICAL_INFO) {
		printf("Nonbonded parameters size: %llu bytes\n", sizeof(ForceField_NB));
		printH2("Finished building forcefield");
	}
}


int Forcefield::getAtomtypeID(int gro_id) const {
	if (groIdToAtomtypeMap.count(gro_id) == 0 || gro_id == 0) {	// 0 is an error, as atoms are 1-indexed
		printf("Attempting to fetch atomtype of non-loaded atom with global_id %d\n", gro_id);
		exit(0);
	}
	return groIdToAtomtypeMap.find(gro_id)->second;
}

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

	if (vl >= V1) { printf("%d NB_Atomtypes loaded\n", atomtypes.size()); }
	return atomtypes;
}

std::map<int, int> Forcefield::loadAtomTypeMap(const SimpleParsedFile& parsedfile) {	// returns the nonbonded atomtype
	std::map<int, int> groidToType;

	for (auto& row : parsedfile.rows) {
		if (row.section == "atomtype_map") {
			const int gro_id = stoi(row.words[0]);
			const int atomtype_id = stoi(row.words[1]);
			groidToType.insert({ gro_id, atomtype_id });
		}	
	}
	if (vl >= V1) { printf("%d NB_Atomtype_IDs loaded\n", groidToType.size()); }

	return groidToType;
}


Forcefield::Topology Forcefield::loadTopology(const SimpleParsedFile& parsedfile)
{
	Topology topology;

	for (auto& row : parsedfile.rows) {
		if (row.section == "singlebonds") {
			assert(row.words.size() == 6);

			std::array<uint32_t, 2> gro_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 2; i++) {
				gro_ids[i] = stoi(row.words[i]);
			}
			const float b0 = stof(row.words[4]) * NANO_TO_LIMA;							// convert [nm] to [lm]
			// Units of kb is [J/mol/nm^2]. One nm is for the error in forcecalc, and one distance in integration
			const float kb = stof(row.words[5]) / (NANO_TO_PICO * NANO_TO_LIMA);		// convert [J/(mol * nm^2)] to [J/(mol *  * lm)
			//const float kb = stof(row.words[5]) / (NANO_TO_LIMA * NANO_TO_LIMA);

			topology.singlebonds.emplace_back(SingleBond{ gro_ids, b0, kb });
		}
		else if (row.section == "anglebonds") {
			assert(row.words.size() == 8);

			std::array<uint32_t, 3> gro_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 3; i++) {
				gro_ids[i] = stoi(row.words[i]);
			}
			const float theta0 = stof(row.words[6]);
			const float ktheta = stof(row.words[7]);
			topology.anglebonds.emplace_back(AngleBond{ gro_ids, theta0, ktheta });
		}
		else if (row.section == "dihedralbonds") {
			assert(row.words.size() == 11);

			std::array<uint32_t, 4> gro_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 4; i++) {
				gro_ids[i] = stoi(row.words[i]);
			}
			const float phi0 = stof(row.words[8]);
			const float kphi = stof(row.words[9]);
			const int multiplicity = stoi(row.words[10]);
			topology.dihedralbonds.emplace_back(DihedralBond{ gro_ids, phi0, kphi, multiplicity });
		}
		else if (row.section == "improperdihedralbonds") {
			assert(row.words.size() == 10);

			std::array<uint32_t, 4> gro_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 4; i++) {
				gro_ids[i] = stoi(row.words[i]);
			}
			const float psi0 = stof(row.words[8]);
			const float kpsi = stof(row.words[9]);
			topology.improperdihedralbonds.emplace_back(ImproperDihedralBond{ gro_ids, psi0, kpsi });
		}
	}

	return topology;
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
