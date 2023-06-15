#include <vector>

#include "LIMA_FORCEFIELDMAKER/include/ForcefieldTypes.h"
#include "LIMA_FORCEFIELDMAKER/include/ForcefieldMaker.h"
#include "LIMA_FORCEFIELDMAKER/include/Filehandling.h"	// TODO: Delete file

#include "LIMA_BASE/include/Filehandling.h"


#include <filesystem>
#include <assert.h>
#include <format>


using std::string;
using FileRows = std::vector<std::vector<std::string>>;

struct Forcefield {
	// Contains only 1 entry for each type that exists
	//std::vector<NB_Atomtype> nb_atoms;
	std::unordered_map<std::string, NB_Atomtype> atomToTypeMap;
	std::vector<Singlebondtype> singlebonds;
	std::vector<Anglebondtype> anglebonds;
	std::vector<Dihedralbondtype> dihedralbonds;
	std::vector<Improperdihedralbondtype> improperdeihedralbonds;


};

// TODO: Change all of these so they don't have both the atomtypename and the gro_ids?
struct Topology {
	// Contains only 1 entry for each entry in the topology file
	//std::vector<NB_Atomtype> nb_atoms;
	AtomTable atomtable;
	std::vector<Singlebondtype> singlebonds;
	std::vector<Anglebondtype> anglebonds;
	std::vector<Dihedralbondtype> dihedralbonds;
	std::vector<Improperdihedralbondtype> improperdeihedralbonds;
};

struct AtomtypeMapping {
	const int gro_id;
	const int atomtype_id; // simulation specific
};

//#ifdef __linux__
//	string sim_path = "../../Simulation";
//	//string sim_path = "/home/lima/Desktop/LIMA/Simulation";
//#else
//	string sim_path = "C:\\PROJECTS\\QUANTOM\\Simulation";
//#endif

//string mol_path = FileHelpers::pathJoin(sim_path, "Molecule");
//string forcefield_path = FileHelpers::pathJoin(sim_path, "Forcefield");

ForcefieldMaker::ForcefieldMaker(const string& workdir, EnvMode envmode, const string& ff_dir, const string& conf_file, const string& topol_file) :
	molecule_dir(FileHelpers::pathJoin(workdir, "molecule")),
	forcefield_dir(ff_dir),
	logger(LimaLogger::LogMode::compact, envmode, "forcefieldmaker", workdir),
	m_verbose(envmode != Headless)
{
	ff_bonded_path = FileHelpers::pathJoin(ff_dir, "ffbonded.itp");
	ff_nonbonded_path = FileHelpers::pathJoin(ff_dir, "ffnonbonded.itp");
	assertPath(ff_bonded_path);
	assertPath(ff_nonbonded_path);

	conf_path = FileHelpers::pathJoin(molecule_dir, conf_file);
	topol_path = FileHelpers::pathJoin(molecule_dir, topol_file);
	assertPath(conf_path);
	assertPath(topol_path);
}

void loadFFnonbondedIntoForcefield(const SimpleParsedFile& parsedfile, Forcefield& forcefield) {
	for (const SimpleParsedFile::Row& row : parsedfile.rows) {

		if (row.section == "atomtypes") {
			if (row.words.size() < 5) {
				int a = 0;
			}
			assert(row.words.size() >= 5);

			const string& atomtype = row.words[0];
			const int atnum = stoi(row.words[1]);
			const float mass = stof(row.words[2]);		// Should this come from topol too?
			const float charge = stof(row.words[3]);	// Comes from topol
			const float sigma = stof(row.words[5]);
			const float epsilon = stof(row.words[6]);

			//forcefield.nb_atoms.emplace_back(NB_Atomtype(atomtype, atnum, mass, sigma, epsilon));
			forcefield.atomToTypeMap.insert(std::pair(atomtype, NB_Atomtype(atomtype, atnum, mass, sigma, epsilon)));
		}
		else if (row.section == "pairtypes") {
			// TODO: Fill out this logic
		}		
	}
}

void loadFFbondedIntoForcefield(const SimpleParsedFile& parsedfile, Forcefield& forcefield) {
	for (const SimpleParsedFile::Row& row : parsedfile.rows) {

		if (row.section == "bondtypes") {
			assert(row.words.size() >= 5);

			const std::array<string, 2> bondedatom_typenames{ row.words[0], row.words[1] };
			const float b0 = stof(row.words[3]);
			const float kb = stof(row.words[4]);

			forcefield.singlebonds.emplace_back(Singlebondtype(bondedatom_typenames, b0, kb));
		}
		else if (row.section == "angletypes") {
			assert(row.words.size() >= 6);

			const std::array<string, 3> angle_typenames{ row.words[0], row.words[1], row.words[2] };
			const float theta0 = stof(row.words[4]);
			const float ktheta = stof(row.words[5]);

			forcefield.anglebonds.emplace_back(Anglebondtype(angle_typenames, theta0, ktheta));
		}
		else if (row.section == "dihedraltypes") {
			assert(row.words.size() >= 7);

			const std::array<string, 4> dihedral_typenames{ row.words[0], row.words[1], row.words[2], row.words[3] };
			const float phi0 = stof(row.words[5]);
			const float kphi = stof(row.words[6]);
			const float n = stoi(row.words[7]);

			forcefield.dihedralbonds.emplace_back(Dihedralbondtype(dihedral_typenames, phi0, kphi, n));
		}
		else if (row.section == "improperdihedraltypes") {
			assert(row.words.size() >= 7);

			const std::array<string, 4> improper_dihedral_typenames{ row.words[0], row.words[1], row.words[2], row.words[3] };
			const float psi0 = stof(row.words[5]);
			const float kpsi = stof(row.words[6]);

			forcefield.improperdeihedralbonds.emplace_back(Improperdihedralbondtype(improper_dihedral_typenames, psi0, kpsi));
		}
	}
}

template <int n>
std::pair<std::array<int, n>, std::array<string, n>> getGroidsAndTypenames(const std::vector<string>& words, const Topology& topology) {
	std::array<int, n> gro_ids;// { stoi(row.words[0]), stoi(row.words[1]) };
	std::array<string, n> atomtypes;
	
	for (int i = 0; i < n; i++) {
		gro_ids[i] = stoi(words[i]);
		if (topology.atomtable.find(gro_ids[i]) == topology.atomtable.end()) {
			// TODO:
			// I guess this just means that there was no match in the declared atoms of the topology, and thus we can simply skip this bond?
			throw "Couldn't find appropriate atomtype from gro id";
		}
		atomtypes[i] = topology.atomtable.find(gro_ids[i])->second.atomtype;
	}

	return std::pair(gro_ids, atomtypes);
}

Topology loadTopology(const SimpleParsedFile& parsedfile) {
	Topology topology;

	for (const SimpleParsedFile::Row& row : parsedfile.rows) {

		if (row.section == "atoms") {
			assert(row.words.size() >= 8);

			const int gro_id = stoi(row.words[0]);
			const string atomtype = row.words[1];
			const string atomname = row.words[4];		// nb type??
			const float charge = stof(row.words[6]);	// not currently used
			const float mass = stof(row.words[7]);		// not currently used

			topology.atomtable.insert(std::pair(gro_id, Atom(gro_id, atomtype, atomname)));
		}
		else if (row.section == "bonds") {
			assert(row.words.size() >= 3);

			auto groidid_atomtype_pair = getGroidsAndTypenames<2>(row.words, topology);
			topology.singlebonds.emplace_back(Singlebondtype{ groidid_atomtype_pair.first, groidid_atomtype_pair.second});
		}
		else if (row.section == "angles") {
			assert(row.words.size() >= 4);

			auto groidid_atomtype_pair = getGroidsAndTypenames<3>(row.words, topology);
			topology.anglebonds.emplace_back(Anglebondtype{ groidid_atomtype_pair.first, groidid_atomtype_pair.second });
		}
		else if (row.section == "dihedrals") {
			assert(row.words.size() >= 5);

			auto groidid_atomtype_pair = getGroidsAndTypenames<4>(row.words, topology);
			topology.dihedralbonds.emplace_back(Dihedralbondtype{ groidid_atomtype_pair.first, groidid_atomtype_pair.second });
		}
		else if (row.section == "improperdihedrals") {
			assert(row.words.size() >= 5);

			auto groidid_atomtype_pair = getGroidsAndTypenames<4>(row.words, topology);
			topology.improperdeihedralbonds.emplace_back(Improperdihedralbondtype{ groidid_atomtype_pair.first, groidid_atomtype_pair.second });
		}
	}
	return topology;
}

// Adds entries of forcefield.atomtypes to a filtered list, if there is any reference to said atomtype in the topology
const std::vector<NB_Atomtype> filterAtomtypes(const Topology& topology, Forcefield& forcefield) {
	std::vector<NB_Atomtype> atomtypes_filtered;

	atomtypes_filtered.emplace_back(NB_Atomtype("WATER", 0, water_mass, water_sigma, water_epsilon));		// Solvent type always first!


	for (const auto& record : topology.atomtable) {
		const string& atomtype = record.second.atomtype;

		NB_Atomtype& atomtype_ff = forcefield.atomToTypeMap.find(atomtype)->second;

		if (!atomtype_ff.is_present_in_simulation) {
			atomtype_ff.is_present_in_simulation = true;
			atomtypes_filtered.push_back(atomtype_ff);
		}
	}
	return atomtypes_filtered;
}

void fillTBondParametersFromForcefield(const Forcefield& forcefield, Topology& topology) {
	FTHelpers::assignForceVariablesFromForcefield(topology.singlebonds, forcefield.singlebonds);
	FTHelpers::assignForceVariablesFromForcefield(topology.anglebonds, forcefield.anglebonds);
	FTHelpers::assignForceVariablesFromForcefield(topology.dihedralbonds, forcefield.dihedralbonds);
	FTHelpers::assignForceVariablesFromForcefield(topology.improperdeihedralbonds, forcefield.improperdeihedralbonds);
}

int findIndexOfAtomtype(const string& query_atomtype_name, const std::vector<NB_Atomtype>& atomtypes) {
	for (int i = 0; i < atomtypes.size(); i++) {
		if (query_atomtype_name == atomtypes[i].type) {
			return i;
		}
	}

	throw "Failed to find atomtype";
}

const std::vector<AtomtypeMapping> mapGroidsToSimulationspecificAtomtypeids(const Topology& topology, const std::vector<NB_Atomtype>& atomtypes_filtered) {
	std::vector<AtomtypeMapping> map;

	for (const auto& e : topology.atomtable) {
		const Atom& atom = e.second;

		int filted_atomtype_id = findIndexOfAtomtype(atom.atomtype, atomtypes_filtered);

		map.push_back(AtomtypeMapping{ atom.gro_id, filted_atomtype_id });
	}
	return map;
}

void printFFNonbonded(const string& path, const std::vector<NB_Atomtype>& filtered_atomtypes) {}
void printFFBonded(const string& path, const Topology& topology) {}

void ForcefieldMaker::prepSimulationForcefield() {
	// Check if filtered files already exists, if so return


	Forcefield forcefield;
	
	// Load the forcefields. 
	const SimpleParsedFile nb_parsedfile = Filehandler::parseItpFile(ff_nonbonded_path);
	loadFFnonbondedIntoForcefield(nb_parsedfile, forcefield);

	const SimpleParsedFile bonded_parsedfile = Filehandler::parseItpFile(ff_bonded_path);
	loadFFbondedIntoForcefield(bonded_parsedfile, forcefield);

	// Load the topology
	const SimpleParsedFile topology_parsedfile = Filehandler::parseTopFile(topol_path);
	Topology topology = loadTopology(topology_parsedfile);

	// Filter for the atomtypes used in this simulation and map to them
	const std::vector<NB_Atomtype> atomtypes_filtered = filterAtomtypes(topology, forcefield);
	const std::vector<AtomtypeMapping> atomtype_map = mapGroidsToSimulationspecificAtomtypeids(topology, atomtypes_filtered);

	// Match the topology with the forcefields.
	fillTBondParametersFromForcefield(forcefield, topology);

	printFFNonbonded(Filehandler::pathJoin(molecule_dir, "LIMA_ffnonbonded_filtered.txt"), atomtypes_filtered);
	printFFBonded(Filehandler::pathJoin(molecule_dir, "LIMA_ffbonded_filtered.txt"), topology);
	logger.finishSection("Prepare Forcefield has finished");
}







