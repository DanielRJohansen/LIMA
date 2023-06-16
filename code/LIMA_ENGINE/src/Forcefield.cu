#pragma once

#include "LIMA_ENGINE/include/Forcefield.cuh"
#include "LIMA_BASE/include/Printer.h"
#include "LIMA_ENGINE/include/EngineUtils.cuh"

using namespace LIMA_Print;

const int min_reserve_size = 10000;	// This should NOT be permanent...

Forcefield::Forcefield(VerbosityLevel vl) : vl(vl) {};


void Forcefield::loadForcefield(string molecule_dir) {
	if (vl >= CRITICAL_INFO) { printH2("Building forcefield"); }

	SimpleParsedFile nonbonded_parsed = Filehandler::parseLffFile(Filehandler::pathJoin(molecule_dir, "ffnonbonded_filtered.lff"));
	SimpleParsedFile bonded_parsed = Filehandler::parseLffFile(Filehandler::pathJoin(molecule_dir, "ffbonded_filtered.lff"));

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

template <int array_length>
bool isMatch(const uint32_t* topolbonds, const std::array<int, array_length> query_ids) {
	for (int i = 0; i < array_length; i++) {
		if (topolbonds[i] != query_ids[i]) { return false; }
	}
	return true;
}

const SingleBond& Forcefield::getSinglebondtype(int bond_index, std::array<int, 2> gro_ids) const
{
	const auto& bond = topology.singlebonds[bond_index];
	assert(isMatch(bond.atom_indexes, gro_ids));
	return bond;
}
const AngleBond& Forcefield::getAnglebondtype(int bond_index, std::array<int, 3> gro_ids) const
{
	const auto& bond = topology.anglebonds[bond_index];
	assert(isMatch(bond.atom_indexes, gro_ids));
	return bond;
}
const DihedralBond& Forcefield::getDihedralbondtype(int bond_index, std::array<int, 4> gro_ids) const
{
	const auto& bond = topology.dihedralbonds[bond_index];
	assert(isMatch(bond.atom_indexes, gro_ids));
	return bond;
}
const ImproperDihedralBond& Forcefield::getImproperdihedralbondtype(int bond_index, std::array<int, 4> gro_ids) const
{
	const auto& bond = topology.improperdihedralbonds[bond_index];
	assert(isMatch(bond.atom_indexes, gro_ids));
	return bond;
}
//
//const SingleBond* Forcefield::getBondType(std::array<int, 2> ids) const {
//	for (auto& singlebond : topol_bonds) {
//		if (isMatch(singlebond.atom_indexes, ids)) {
//			return &singlebond;
//		}
//	}
//	printf("Bond not found with ids %d %d\n", ids[0], ids[1]);
//	exit(0);
//}
//
//const AngleBond* Forcefield::getAngleType(std::array<int, 3> ids) const {
//	for (auto& anglebond : topol_angles) {
//		if (isMatch(anglebond.atom_indexes, ids)) {
//			return &anglebond;
//		}
//	}
//
//	printf("Angle not found with ids %d %d %d\n", ids[0], ids[1], ids[2]);
//	exit(0);
//}
//
//const DihedralBond* Forcefield::getDihedralType(std::array<int, 4> ids) const {
//	for (auto& dihedralbond : topol_dihedrals) {
//		if (isMatch(dihedralbond.atom_indexes, ids)) {
//			return &dihedralbond;
//		}
//	}
//
//	printf("Dihedral not found with ids %d %d %d %d\n", ids[0], ids[1], ids[2], ids[3]);
//	exit(0);
//}

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


	//for (vector<string> row : summary_rows) {
	//	if (newParseTitle(row)) {
	//		current_state = setState(row[1], current_state);
	//		continue;
	//	}

	//	

	//	if (current_state == FF_NONBONDED) {
	//		//for (string e : row)
	//			//cout << e << '\t';
	//		//printf("\n");

	//		//atomtypes[ptr++] = NBAtomtype(stof(row[2]), stof(row[3]), stof(row[4]));
	//		nb_atomtypes.push_back(NBAtomtype{ stof(row[2]), stof(row[3]), stof(row[4]) });
	//	}			
	//}
	//n_nb_atomtypes = nb_atomtypes.size();
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

		//if (current_state == NB_ATOMTYPES) {
		//	int gro_id = stoi(row[0]);
		//	int atomtype_id = stoi(row[1]);


		//	groidToType.insert(std::pair<int, int>(gro_id, atomtype_id));
		//}
			
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

			std::array<int, 2> gro_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 2; i++) {
				gro_ids[i] = stoi(row.words[i]);
			}
			const float b0 = stof(row.words[4]) * NANO_TO_LIMA;						// convert [nm] to [lm]*/
			const float kb = stof(row.words[5]) / (NANO_TO_LIMA * NANO_TO_LIMA);		// convert [J/(mol * nm^2)] to [J/(mol * nm * lm)

			topology.singlebonds.emplace_back(SingleBond{ gro_ids, b0, kb });
		}
		else if (row.section == "anglebonds") {
			assert(row.words.size() == 8);

			std::array<int, 3> gro_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 3; i++) {
				gro_ids[i] = stoi(row.words[i]);
			}
			const float theta0 = stof(row.words[6]);
			const float ktheta = stof(row.words[7]);
			topology.anglebonds.emplace_back(AngleBond{ gro_ids, theta0, ktheta });
		}
		else if (row.section == "dihedrals") {
			assert(row.words.size() == 11);

			std::array<uint32_t, 4> gro_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 4; i++) {
				gro_ids[i] = stoi(row.words[i]);
			}
			const float phi0 = stof(row.words[7]);
			const float kphi = stof(row.words[8]);
			const int multiplicity = stoi(row.words[9]);
			topology.dihedralbonds.emplace_back(DihedralBond{ gro_ids, phi0, kphi, multiplicity });
		}
		else if (row.section == "impropers") {
			assert(row.words.size() == 10);

			std::array<uint32_t, 4> gro_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 4; i++) {
				gro_ids[i] = stoi(row.words[i]);
			}
			const float psi0 = stof(row.words[7]);
			const float kpsi = stof(row.words[8]);
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
		forcefield.particle_parameters[i].mass = atomtypes[i].mass * 1e-3f;				// Convert g/mol to kg/mol
		forcefield.particle_parameters[i].sigma = atomtypes[i].sigma * NANO_TO_LIMA;		// Convert from [nm] to [lm]
		forcefield.particle_parameters[i].epsilon = atomtypes[i].epsilon;				// Interpreted as kg*lm^2/ls^2 

		bool illegal_parameter = (forcefield.particle_parameters[i].mass < mass_min) || (forcefield.particle_parameters[i].sigma < sigma_min) || (forcefield.particle_parameters[i].epsilon < epsilon_min);

		if ((vl >= V2) || illegal_parameter) { printf("Mass %f Sigma %f Epsilon %f\n", atomtypes[i].mass, atomtypes[i].sigma, atomtypes[i].epsilon); }
	}

}
