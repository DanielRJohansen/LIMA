#include <vector>

#include "LIMA_BASE/include/Filehandling.h"

#include "LIMA_FORCEFIELDMAKER/include/ForcefieldTypes.h"
#include "LIMA_FORCEFIELDMAKER/include/ForcefieldMaker.h"



#include <filesystem>
#include <assert.h>
#include <format>

using std::string, std::cout, std::endl, std::to_string;
using FileRows = std::vector<std::vector<std::string>>;

struct BondedTypes {
	// Contains only 1 entry for each type that exists
	std::unordered_map<std::string, NB_Atomtype> atomToTypeMap;	// THis is a bad name
	std::vector<Singlebondtype> singlebonds;
	std::vector<Anglebondtype> anglebonds;
	std::vector<Dihedralbondtype> dihedralbonds;
	std::vector<Improperdihedralbondtype> improperdeihedralbonds;


};

// TODO: Change all of these so they don't have both the atomtypename and the gro_ids?
struct Topology {
	// Contains only 1 entry for each entry in the topology file
	AtomTable atomtable;
	std::vector<Singlebondtype> singlebonds{};
	std::vector<Anglebondtype> anglebonds{};
	std::vector<Dihedralbondtype> dihedralbonds{};
	std::vector<Improperdihedralbondtype> improperdeihedralbonds{};
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
	molecule_dir(Filehandler::pathJoin(workdir, "molecule")),
	forcefield_dir(ff_dir),
	logger(LimaLogger::LogMode::compact, envmode, "forcefieldmaker", workdir),
	m_verbose(envmode != Headless)
{
	ff_bonded_path = Filehandler::pathJoin(ff_dir, "ffbonded.itp");
	ff_nonbonded_path = Filehandler::pathJoin(ff_dir, "ffnonbonded.itp");
	Filehandler::assertPath(ff_bonded_path);
	Filehandler::assertPath(ff_nonbonded_path);

	conf_path = Filehandler::pathJoin(molecule_dir, conf_file);
	topol_path = Filehandler::pathJoin(molecule_dir, topol_file);
	Filehandler::assertPath(conf_path);
	Filehandler::assertPath(topol_path);
}

void loadFFnonbondedIntoForcefield(const SimpleParsedFile& parsedfile, BondedTypes& forcefield) {
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

void loadFFbondedIntoForcefield(const SimpleParsedFile& parsedfile, BondedTypes& forcefield) {
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
bool getGroidsAndTypenames(const std::vector<string>& words, const Topology& topology, std::array<int, n>& gro_ids, std::array<string, n>& atomtypes) {	
	for (int i = 0; i < n; i++) {
		gro_ids[i] = stoi(words[i]);
		if (topology.atomtable.find(gro_ids[i]) == topology.atomtable.end()) {
			return false;	// All atoms of bond is not available in simulation input files (topology atoms)
		}
		atomtypes[i] = topology.atomtable.find(gro_ids[i])->second.atomtype;
	}

	return true;
}

Topology loadTopology(const SimpleParsedFile& parsedfile, const char ignored_atom) {
	Topology topology;

	for (const SimpleParsedFile::Row& row : parsedfile.rows) {

		if (row.section == "atoms") {
			assert(row.words.size() >= 8);

			const int gro_id = stoi(row.words[0]);
			const string atomtype = row.words[1];
			const string atomname = row.words[4];		// nb type??
			const float charge = stof(row.words[6]);	// not currently used
			const float mass = stof(row.words[7]);		// not currently used

			if (atomtype[0] == ignored_atom) {	// Typically for ignoring hydrogen
				continue;
			}

			topology.atomtable.insert(std::pair(gro_id, Atom(gro_id, atomtype, atomname)));
		}
		else if (row.section == "bonds") {
			assert(row.words.size() >= 3);

			std::array<int, 2> gro_ids;
			std::array<string, 2> atomtypes;
			if (getGroidsAndTypenames(row.words, topology, gro_ids, atomtypes)) {
				topology.singlebonds.emplace_back(Singlebondtype{ gro_ids, atomtypes});
			}
			
		}
		else if (row.section == "angles") {
			assert(row.words.size() >= 4);

			std::array<int, 3> gro_ids;
			std::array<string, 3> atomtypes;
			if (getGroidsAndTypenames(row.words, topology, gro_ids, atomtypes)) {
				topology.anglebonds.emplace_back(Anglebondtype{ gro_ids, atomtypes });
			}
		}
		else if (row.section == "dihedrals") {
			assert(row.words.size() >= 5);

			std::array<int, 4> gro_ids;
			std::array<string, 4> atomtypes;
			if (getGroidsAndTypenames(row.words, topology, gro_ids, atomtypes)) {
				topology.dihedralbonds.emplace_back(Dihedralbondtype{ gro_ids, atomtypes });
			}
		}
		else if (row.section == "improperdihedrals") {
			assert(row.words.size() >= 5);

			std::array<int, 4> gro_ids;
			std::array<string, 4> atomtypes;
			if (getGroidsAndTypenames(row.words, topology, gro_ids, atomtypes)) {
				topology.improperdeihedralbonds.emplace_back(Improperdihedralbondtype{ gro_ids, atomtypes });
			}
		}
	}
	return topology;
}

// Adds entries of forcefield.atomtypes to a filtered list, if there is any reference to said atomtype in the topology
const std::vector<NB_Atomtype> filterAtomtypes(const Topology& topology, BondedTypes& forcefield) {
	std::vector<NB_Atomtype> atomtypes_filtered;

	atomtypes_filtered.emplace_back(NB_Atomtype("WATER", 0, water_mass, water_sigma, water_epsilon));		// Solvent type always first!


	for (const auto& record : topology.atomtable) {
		const string& atomtype = record.second.atomtype;

		NB_Atomtype& atomtype_ff = forcefield.atomToTypeMap.find(atomtype)->second;

		if (!atomtype_ff.is_present_in_simulation) {
			atomtype_ff.is_present_in_simulation = true;
			atomtype_ff.atnum_local = atomtypes_filtered.size();
			atomtypes_filtered.push_back(atomtype_ff);
		}
	}
	return atomtypes_filtered;
}

void fillTBondParametersFromForcefield(const BondedTypes& forcefield, Topology& topology) {
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


namespace FFPrintHelpers {
	static string titleH1(string text) {
		return "// ----------------################ " + text + " ################---------------- //\n\n";
	}
	static string titleH2(string text) {
		return "// ---------------- " + text + " ---------------- //\n";
	}
	static string titleH3(string text) {
		return "// " + text + " //\n";
	}
	static string parserTitle(string text) {
		return "# " + text + '\n';
	}
	static string endBlock() {
		return "\n\n\n\n";
	}
};


const char delimiter = ' ';

void printFFNonbonded(const string& path, const std::vector<AtomtypeMapping>& atomtype_map, const std::vector<NB_Atomtype>& filtered_atomtypes) {
	std::ofstream file(path, std::ofstream::out);
	if (!file.is_open()) {
		cout << "Failed to open file " << path << endl;
		exit(0);
	}

	file << FFPrintHelpers::titleH1("Forcefield Non-bonded");
	file << FFPrintHelpers::titleH2("GRO_id to simulation-specific atomtype map");
	file << FFPrintHelpers::titleH3("{gro_id \t atomtype}");
	file << FFPrintHelpers::parserTitle("atomtype_map");
	for (auto& mapping : atomtype_map) {
		file << to_string(mapping.gro_id) << delimiter << to_string(mapping.atomtype_id) << endl;
	}
	file << FFPrintHelpers::endBlock();

	file << FFPrintHelpers::titleH2("Non-bonded parameters");
	file << FFPrintHelpers::titleH3("{atom_type \t type_id \t mass [g/mol] \t sigma [nm] \t epsilon [J/mol]}");
	file << FFPrintHelpers::parserTitle("atomtypes");
	for (auto& atomtype : filtered_atomtypes) {
		file << atomtype.type << delimiter << to_string(atomtype.atnum_local) << delimiter << to_string(atomtype.mass) << delimiter << to_string(atomtype.sigma) << delimiter << to_string(atomtype.epsilon) << endl;
	}
	file << FFPrintHelpers::endBlock();

	file.close();
}

void printFFBonded(const string& path, const Topology& topology) {
	std::ofstream file(path, std::ofstream::out);
	if (!file.is_open()) {
		printf(("Failed to open file\n"));
		exit(0);
	}

	file << FFPrintHelpers::titleH1("Forcefield Bonded");
	file << FFPrintHelpers::titleH2("Singlebonds");
	file << FFPrintHelpers::titleH3("{ID_p1 \t ID_p2 \t Atomtype \t Atomtype \t b_0 [nm] \t k_b [J/(mol * nm^2)]}");
	file << FFPrintHelpers::parserTitle("singlebonds");
	for (auto& bond : topology.singlebonds) {
		file << to_string(bond.gro_ids[0]) << delimiter << to_string(bond.gro_ids[1]) << delimiter
			<< bond.bonded_typenames[0] << delimiter << bond.bonded_typenames[1] << delimiter
			<< to_string(bond.b0) << delimiter << to_string(bond.kb) << endl;
	}
	file << FFPrintHelpers::endBlock();



	file << FFPrintHelpers::titleH2("Anglebonds");
	file << FFPrintHelpers::titleH3("{Atom-IDs \t Atomtypes \t theta_0 [rad] \t k_theta [J/(mol * rad^2)}");
	file << FFPrintHelpers::parserTitle("anglebonds");
	for (auto& angle : topology.anglebonds) {
		file << to_string(angle.gro_ids[0]) << delimiter << to_string(angle.gro_ids[1]) << delimiter << to_string(angle.gro_ids[2]) << delimiter
			<< angle.bonded_typenames[0] << delimiter << angle.bonded_typenames[1] << delimiter << angle.bonded_typenames[2] << delimiter
			<< to_string(angle.theta0) << delimiter << to_string(angle.ktheta) << endl;
	}
	file << FFPrintHelpers::endBlock();



	file << FFPrintHelpers::titleH2("Dihedrals");
	file << FFPrintHelpers::titleH3("{Atom IDs \t Atomtypes \t phi_0 [rad] \t k_phi [J/(mol * rad^2)] \t n}");
	file << FFPrintHelpers::parserTitle("dihedralbonds");
	for (auto& dihedral : topology.dihedralbonds) {
		for (auto& groid : dihedral.gro_ids) {
			file << to_string(groid) << delimiter;
		}
		for (auto& type : dihedral.bonded_typenames) {
			file << type << delimiter;
		}
		file << to_string(dihedral.phi0) << delimiter << to_string(dihedral.kphi) << delimiter << to_string(dihedral.n) << endl;
	}
	file << FFPrintHelpers::endBlock();



	file << FFPrintHelpers::titleH2("ImproperDihedrals");
	file << FFPrintHelpers::titleH3("{Atom IDs \t Atomtypes \t psi_0 [rad] \t k_psi [J/(mol * rad^2)]}");
	file << FFPrintHelpers::parserTitle("improperdihedralbonds");
	for (auto improper : topology.improperdeihedralbonds) {
		for (auto& groid : improper.gro_ids) {
			file << to_string(groid) << delimiter;
		}
		for (auto& type : improper.bonded_typenames) {
			file << type << delimiter;
		}
		file << to_string(improper.psi0) << delimiter << to_string(improper.kpsi) << endl;
	}
	file << FFPrintHelpers::endBlock();

	file.close();
}

void ForcefieldMaker::prepSimulationForcefield(const char ignored_atomtype) {
	// Check if filtered files already exists, if so return
	BondedTypes forcefield;
	
	// Load the forcefields. 
	const SimpleParsedFile nb_parsedfile = Filehandler::parseItpFile(ff_nonbonded_path);
	loadFFnonbondedIntoForcefield(nb_parsedfile, forcefield);

	const SimpleParsedFile bonded_parsedfile = Filehandler::parseItpFile(ff_bonded_path);
	loadFFbondedIntoForcefield(bonded_parsedfile, forcefield);

	// Load the topology
	const SimpleParsedFile topology_parsedfile = Filehandler::parseTopFile(topol_path);
	Topology topology = loadTopology(topology_parsedfile, ignored_atomtype);

	// Filter for the atomtypes used in this simulation and map to them
	const std::vector<NB_Atomtype> atomtypes_filtered = filterAtomtypes(topology, forcefield);
	const std::vector<AtomtypeMapping> atomtype_map = mapGroidsToSimulationspecificAtomtypeids(topology, atomtypes_filtered);

	// Match the topology with the forcefields.
	fillTBondParametersFromForcefield(forcefield, topology);

	printFFNonbonded(Filehandler::pathJoin(molecule_dir, "ffnonbonded_filtered.lff"), atomtype_map, atomtypes_filtered);
	printFFBonded(Filehandler::pathJoin(molecule_dir, "ffbonded_filtered.lff"), topology);
	logger.finishSection("Prepare Forcefield has finished");
}