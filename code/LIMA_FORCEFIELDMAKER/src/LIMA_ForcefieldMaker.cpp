#include <vector>

#include "Filehandling.h"
#include "ForcefieldTypes.h"
#include "ForcefieldMaker.h"



#include <filesystem>
#include <assert.h>
#include <format>
#include <unordered_set>

using std::string, std::cout, std::endl, std::to_string;
using FileRows = std::vector<std::vector<std::string>>;
namespace lfs = Filehandler;

namespace ForcefieldMakerTypes {
	struct BondedTypes {
		// Contains only 1 entry for each type that exists
		std::unordered_map<std::string, NB_Atomtype> atomToTypeMap;	// THis is a bad name
		std::vector<Singlebondtype> singlebonds;
		std::vector<Anglebondtype> anglebonds;
		std::vector<Dihedralbondtype> dihedralbonds;
		std::vector<Improperdihedralbondtype> improperdeihedralbonds;

		// This structure is for storing the mass of an atomtype, before we get to the NBATOMTYPE section
		std::unordered_map<std::string, float> atomnameToMassMap;

	};

	// TODO: Change all of these so they don't have both the atomtypename and the gro_ids?
	struct Topology {
		// Contains only 1 entry for each entry in the topology file
		AtomInfoTable atominfotable;

		std::unordered_set<std::string> active_atomtypes;

		std::vector<Singlebondtype> singlebonds;
		std::vector<Anglebondtype> anglebonds;
		std::vector<Dihedralbondtype> dihedralbonds;
		std::vector<Improperdihedralbondtype> improperdeihedralbonds;
	};

	struct AtomtypeMapping {
		AtomtypeMapping(int global, int gro, int chain_id, int res_id, int atomtype_id, const std::string& name, int unique_residue_id) 
			: global_id(global), gro_id(gro), chain_id(chain_id), residue_id(res_id), atomtype_id(atomtype_id), atomname(name), unique_residue_id(unique_residue_id) {}
		const int global_id;	// Given by LIMA
		const int gro_id;		// not unique
		const int chain_id;		// Unique
		const int residue_id;	// Unique within chain (i fucking hope..)
		const int atomtype_id;	// simulation specific
		const std::string atomname;
		const int unique_residue_id;
	};
}

using namespace ForcefieldMakerTypes;


const float water_mass = 15.999000f + 2.f * 1.008000f;
//const float water_sigma = 1.7398 * rminToSigma * AngToNm;	// Value guessed from param19.inp: OH2      0.0000    -0.0758    1.7398 !ST2   water oxygen
const float water_sigma = 0.22f;	// Made up value. Works better than the one above. I guess i need to implement proper tip3 at some point?
const float water_epsilon = 0.1591f * kcalToJoule;

static const NB_Atomtype Water_atomtype{ "WATER", 0, water_mass, water_sigma, water_epsilon };


void loadFileIntoForcefield(const SimpleParsedFile& parsedfile, BondedTypes& forcefield) {
	for (const SimpleParsedFile::Row& row : parsedfile.rows) {

		if (row.section == "ATOMS") {
			assert(row.words.size() >= 4);

			const string& atomtype = row.words[2];
			const float mass = stof(row.words[3]);		// Should this come from topol too?
			forcefield.atomnameToMassMap.insert(std::pair{ atomtype, mass });

			//forcefield.nb_atoms.emplace_back(NB_Atomtype(atomtype, atnum, mass, sigma, epsilon));
			//forcefield.atomToTypeMap.insert(std::pair(atomtype, NB_Atomtype(atomtype, atnum, mass, sigma, epsilon)));
		}
		else if (row.section == "pairtypes") {
			// TODO: Fill out this logic
		}
		else if (row.section == "BONDS") {
			if(row.words.size() < 4)
				int a = 0;

			assert(row.words.size() >= 4);

			const std::array<string, 2> bondedatom_typenames{ row.words[0], row.words[1] };
			const float kb = stof(row.words[2]) * kcalToJoule / AngToNm / AngToNm * 2.f;		// * 2 to go from V(bond) = Kb(b - b0)**2 -> V(bond) = 0.5*Kb(b - b0)**2
			const float b0 = stof(row.words[3]) * AngToNm;	// A to nm

			forcefield.singlebonds.emplace_back(Singlebondtype(bondedatom_typenames, b0, kb));
		}
		else if (row.section == "ANGLES") {
			assert(row.words.size() >= 5);

			const std::array<string, 3> angle_typenames{ row.words[0], row.words[1], row.words[2] };
			const float ktheta = stof(row.words[3]) * kcalToJoule * 2.f;	// * 2 to convert Ktheta(Theta - Theta0)**2 -> 0.5*Ktheta(Theta - Theta0)**2
			const float theta0 = stof(row.words[4]) * degreeToRad;

			forcefield.anglebonds.emplace_back(Anglebondtype(angle_typenames, theta0, ktheta));
		}
		else if (row.section == "DIHEDRALS") {
			if (row.words.size() < 7)
				int a = 0;
			assert(row.words.size() >= 7);

			const std::array<string, 4> dihedral_typenames{ row.words[0], row.words[1], row.words[2], row.words[3] };
			const float kphi = stof(row.words[4]) * kcalToJoule;
			const int n = stoi(row.words[5]);
			const float phi0 = stof(row.words[6]) * degreeToRad;
			
			forcefield.dihedralbonds.emplace_back(Dihedralbondtype(dihedral_typenames, phi0, kphi, n));
		}
		else if (row.section == "IMPROPER") {
			if (row.words.size() < 7)
				int a = 0;
			assert(row.words.size() >= 7);

			const std::array<string, 4> improper_dihedral_typenames{ row.words[0], row.words[1], row.words[2], row.words[3] };
			const float kpsi = stof(row.words[4]) * kcalToJoule * 2.f;	// * 2 to convert Kpsi(psi - psi0)**2 -> 0.5*Kpsi(psi - psi0)**2
			const float psi0 = stof(row.words[6]) * degreeToRad;

			// TODO: Check that the entry doesn't already exists. Current we sort these bonds, which i am really not very comfortable with, 
			// as the order here should be quite relevant?
			forcefield.improperdeihedralbonds.emplace_back(Improperdihedralbondtype(improper_dihedral_typenames, psi0, kpsi));
		}
		else if (row.section == "NONBONDED") {
			if (row.words.size() < 4)
				int a = 0;
			assert(row.words.size() >= 4);

			const string& atomtype = row.words[0];

			
			const float epsilon = abs(stof(row.words[2]) * kcalToJoule);	// For some fucked reason the eps is *inconsistently* negative...
			const float sigma = stof(row.words[3]) * rminToSigma * AngToNm;	// rmin/2 [A] -> sigma [nm]

			const float mass = forcefield.atomnameToMassMap.find(atomtype)->second;

			forcefield.atomToTypeMap.insert(std::pair(atomtype, NB_Atomtype(atomtype, -1, mass, sigma, epsilon)));


			// Not yet used
			//const float epsilon_1_4 = stof(row.words[6]) * 2;	 // rmin/2 -> sigma
			//const float sigma_1_4 = stof(row.words[6]) * 2;	 // rmin/2 -> sigma			
		}
	}
}

template <int n>
bool getGlobalIDsAndTypenames(const int (&gro_ids)[n], const AtomInfoTable& atominfotable, const int chain_id, std::array<int, n>& global_ids, std::array<string, n>& atomtypes) {
	for (int i = 0; i < n; i++) {
		const int gro_id = gro_ids[i];

		if (!atominfotable.exists(chain_id, gro_id)) {
			return false;
		}
		const Atom& atom = atominfotable.getConstRef(chain_id, gro_id);
		global_ids[i] = atom.global_id;
		atomtypes[i] = atom.atomtype;
	}

	return true;
}

template <int n>
bool getGlobalIDsAndTypenames(const std::vector<string>& words, const AtomInfoTable& atominfotable, const int chain_id, std::array<int, n>& global_ids, std::array<string, n>& atomtypes) {
	for (int i = 0; i < n; i++) {
		const int gro_id = stoi(words[i]);

		if (!atominfotable.exists(chain_id, gro_id)) {
			return false;
		}
		const Atom& atom = atominfotable.getConstRef(chain_id, gro_id);
		global_ids[i] = atom.global_id;
		atomtypes[i] = atom.atomtype;
	}

	return true;
}

#include "MDFiles.h"

void includeFileInTopology(Topology& topology, const ParsedTopologyFile& topolfile, int& current_chain_id, int& current_unique_residue_cnt) {
	for (const auto& atom : topolfile.atoms.entries) {
		if (current_chain_id == -1) {
			current_chain_id++;	// Assume there are no include topol files, and we simply read the atoms into chain 0
		}
		const int gro_id = atom.nr;

		if (atom.section_name.has_value()) {
			current_unique_residue_cnt++;
		}

		// Special case: The first set of atoms we read are from a topol file that starts with a non-residue such as a lipid
		// TODO: WARNING: DANGER: SERIOUS: This means reading 2 consecutive lipids will make the two lipids into the same residue!!!
		// Solved with adding lipid_section manually to all lipid files
		// current_unique_residue_cnt = std::max(current_unique_residue_cnt, 0);
		if (current_unique_residue_cnt == -1) {
			throw std::runtime_error("Something is wrong, we should have a res_id by now.");
		}

		topology.atominfotable.insert(current_chain_id, gro_id, atom.type, atom.atom, atom.resnr, current_unique_residue_cnt);
		//topology.atominfotable.insert(current_chain_id, gro_id, atom.atomtype, atomname, residue_id, current_unique_residue_cnt);
		topology.active_atomtypes.insert(atom.type);
	}

	for (const auto& bond : topolfile.singlebonds.entries) {
		std::array<int, 2> global_ids;
		std::array<string, 2> atomtypes;
		if (getGlobalIDsAndTypenames<2>(bond.atom_indexes, topology.atominfotable, current_chain_id, global_ids, atomtypes)) {
			topology.singlebonds.emplace_back(Singlebondtype{ global_ids, atomtypes });
		}
	}
	for (const auto& bond : topolfile.anglebonds.entries) {
		std::array<int, 3> global_ids;
		std::array<string, 3> atomtypes;
		if (getGlobalIDsAndTypenames<3>(bond.atom_indexes, topology.atominfotable, current_chain_id, global_ids, atomtypes)) {
			topology.anglebonds.emplace_back(Anglebondtype{ global_ids, atomtypes });
		}
	}
	for (const auto& bond : topolfile.dihedralbonds.entries) {
		std::array<int, 4> global_ids;
		std::array<string, 4> atomtypes;
		if (getGlobalIDsAndTypenames<4>(bond.atom_indexes, topology.atominfotable, current_chain_id, global_ids, atomtypes)) {
			topology.dihedralbonds.emplace_back(Dihedralbondtype{ global_ids, atomtypes });
		}
	}
	for (const auto& bond : topolfile.improperdihedralbonds.entries) {
		std::array<int, 4> global_ids;
		std::array<string, 4> atomtypes;
		if (getGlobalIDsAndTypenames<4>(bond.atom_indexes, topology.atominfotable, current_chain_id, global_ids, atomtypes)) {
			topology.improperdeihedralbonds.emplace_back(Improperdihedralbondtype{ global_ids, atomtypes });
		}
	}
}



// call this funciton initially with current_unique_residue_cnt=-1
Topology loadTopology(const std::string& molecule_dir, const ParsedTopologyFile& topol_file, const char ignored_atom)
{
	Topology topology{};
	int current_chain_id = -1;
	int current_unique_residue_cnt = -1;
	
	//const std::unique_ptr<ParsedTopologyFile> topolfile = MDFiles::loadTopologyFile(topol_path);

	// First load the top-level topol
	includeFileInTopology(topology, topol_file, current_chain_id, current_unique_residue_cnt);

	// Now load all includes (.itp)
	for (const auto& molecule : topol_file.molecules.entries) {
		current_chain_id++;
		includeFileInTopology(topology, *molecule.include_file, current_chain_id, current_unique_residue_cnt);
	}

	return topology;
}




// Adds entries of forcefield.atomtypes to a filtered list, if there is any reference to said atomtype in the topology
const std::vector<NB_Atomtype> filterAtomtypes(const Topology& topology, BondedTypes& forcefield) {
	std::vector<NB_Atomtype> atomtypes_filtered;

	atomtypes_filtered.emplace_back(Water_atomtype);

	for (const string& atomtype: topology.active_atomtypes) {

		if (!forcefield.atomToTypeMap.contains(atomtype)) {
			throw std::runtime_error(std::format("Failed to find atomtype {}", atomtype).c_str());
		}
		NB_Atomtype& atomtype_ff = forcefield.atomToTypeMap.find(atomtype)->second;

		if (!atomtype_ff.is_present_in_simulation) {
			atomtype_ff.is_present_in_simulation = true;
			atomtype_ff.atnum_local = static_cast<int>(atomtypes_filtered.size());
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

	// TODO: Some bonds (atleast dihedrals, see par_all35_ethers.prm have coefficients of 0.0, meaning we should disregard them here, and never write them to the .lff file...
	// although, we still need to preserve them not calcing LJ? So we need a new section for in-a-pseudo-bond-dont-calc-lj
}

int findIndexOfAtomtype(const string& query_atomtype_name, const std::vector<NB_Atomtype>& atomtypes) {
	for (int i = 0; i < atomtypes.size(); i++) {
		if (query_atomtype_name == atomtypes[i].type) {
			return i;
		}
	}

	throw std::runtime_error("Failed to find atomtype");
}

const std::vector<AtomtypeMapping> mapGroidsToSimulationspecificAtomtypeids(const Topology& topology, const std::vector<NB_Atomtype>& atomtypes_filtered) {
	std::vector<AtomtypeMapping> map;

	for (const Atom& atom : topology.atominfotable.getAllAtoms()) {
		const int filted_atomtype_id = findIndexOfAtomtype(atom.atomtype, atomtypes_filtered);
		map.push_back(AtomtypeMapping{ atom.global_id, atom.gro_id, atom.chain_id, atom.res_id, filted_atomtype_id, atom.atomname, atom.unique_res_id });
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
		throw std::runtime_error(std::format("Failed to open file {}", path).c_str());
	}

	file << FFPrintHelpers::titleH1("Forcefield Non-bonded");


	file << FFPrintHelpers::titleH2("Non-bonded parameters");
	file << FFPrintHelpers::titleH3("{atom_type \t type_id \t mass [g/mol] \t sigma [nm] \t epsilon [J/mol]}");
	file << FFPrintHelpers::titleH3("Potential(r) = 4 * Epsilon * ((sigma/r)^12 - (sigma/r)^6)");
	file << FFPrintHelpers::parserTitle("atomtypes");
	for (auto& atomtype : filtered_atomtypes) {
		file << atomtype.type << delimiter << to_string(atomtype.atnum_local) << delimiter << to_string(atomtype.mass) << delimiter << to_string(atomtype.sigma) << delimiter << to_string(atomtype.epsilon) << endl;
	}
	file << FFPrintHelpers::endBlock();


	file << FFPrintHelpers::titleH2("GRO_id to simulation-specific atomtype map");
	file << FFPrintHelpers::titleH3("{global_id \t gro_id \t chain_id \t residue_id \t atomtype_id \t atomname \t unique_residue_id}");
	file << FFPrintHelpers::parserTitle("atomtype_map");
	for (auto& mapping : atomtype_map) {
		file << to_string(mapping.global_id) << delimiter 
			<< to_string(mapping.gro_id) << delimiter 
			<< to_string(mapping.chain_id) << delimiter 
			<< to_string(mapping.residue_id) << delimiter 
			<< to_string(mapping.atomtype_id) << delimiter 
			<< mapping.atomname << delimiter
			<< mapping.unique_residue_id
			<< endl;
	}
	file << FFPrintHelpers::endBlock();


	file.close();
}

void printFFBonded(const string& path, const Topology& topology) {
	std::ofstream file(path, std::ofstream::out);
	if (!file.is_open()) {
		throw std::runtime_error("Failed to open file\n");
	}

	file << FFPrintHelpers::titleH1("Forcefield Bonded");
	{
		std::stringstream buffer;  // Create a stringstream to build the data

		file << FFPrintHelpers::titleH2("Singlebonds");
		file << FFPrintHelpers::titleH3("{IDs (global & unique) \t Atomtypes \t b_0 [nm] \t k_b [J/(mol * nm^2)]}");
		file << FFPrintHelpers::titleH3("Potential(r) = 0.5 * k_b * (r-b_0)^2");
		file << FFPrintHelpers::parserTitle("singlebonds");

		for (auto& bond : topology.singlebonds) {
			for (auto& global_id : bond.global_ids) {
				buffer << to_string(global_id) << delimiter;
			}
			for (auto& type : bond.bonded_typenames) {
				buffer << type << delimiter;
			}
			buffer << to_string(bond.b0) << delimiter << to_string(bond.kb) << endl;
		}
		buffer << FFPrintHelpers::endBlock();
		file << buffer.str();  // Write the entire string to the file at once
	}


	{
		std::stringstream buffer;  // Create a stringstream to build the data

		buffer << FFPrintHelpers::titleH2("Anglebonds");
		buffer << FFPrintHelpers::titleH3("{Atom-IDs (global & unique) \t Atomtypes \t theta_0 [rad] \t k_theta [J/(mol * rad^2)}");
		buffer << FFPrintHelpers::titleH3("Potential(theta) = 0.5 * k_theta * (theta-theta_0)^2");
		buffer << FFPrintHelpers::parserTitle("anglebonds");

		for (auto& bond : topology.anglebonds) {
			for (auto& global_id : bond.global_ids) {
				buffer << to_string(global_id) << delimiter;
			}
			for (auto& type : bond.bonded_typenames) {
				buffer << type << delimiter;
			}
			buffer << to_string(bond.theta0) << delimiter << to_string(bond.ktheta) << endl;
		}

		buffer << FFPrintHelpers::endBlock();  // Include endBlock in the buffer

		file << buffer.str();  // Write the entire string to the file at once
	}



	{
		std::stringstream buffer;  // Create a stringstream to build the data

		buffer << FFPrintHelpers::titleH2("Dihedrals");
		buffer << FFPrintHelpers::titleH3("{Atom IDs (global & unique) \t Atomtypes \t phi_0 [rad] \t k_phi [J/(mol * rad^2)] \t n}");
		buffer << FFPrintHelpers::titleH3("Potential(phi) = k_phi * (1 + cos(n * phi - phi_0))");
		buffer << FFPrintHelpers::parserTitle("dihedralbonds");

		for (auto& dihedral : topology.dihedralbonds) {
			for (auto& global_id : dihedral.global_ids) {
				buffer << to_string(global_id) << delimiter;
			}
			for (auto& type : dihedral.bonded_typenames) {
				buffer << type << delimiter;
			}
			buffer << to_string(dihedral.phi0) << delimiter << to_string(dihedral.kphi) << delimiter << to_string(dihedral.n) << endl;
		}

		buffer << FFPrintHelpers::endBlock();  // Include endBlock in the buffer

		file << buffer.str();  // Write the entire string to the file at once
	}

	{
		std::stringstream buffer;  // Create a stringstream to build the data

		buffer << FFPrintHelpers::titleH2("ImproperDihedrals");
		buffer << FFPrintHelpers::titleH3("{Atom IDs (global & unique) \t Atomtypes \t psi_0 [rad] \t k_psi [J/(mol * rad^2)]}");
		buffer << FFPrintHelpers::titleH3("Potential(psi) = 0.5 * k_psi * (psi-psi_0)^2");
		buffer << FFPrintHelpers::parserTitle("improperdihedralbonds");

		for (auto improper : topology.improperdeihedralbonds) {
			for (auto& global_id : improper.global_ids) {
				buffer << to_string(global_id) << delimiter;
			}
			for (auto& type : improper.bonded_typenames) {
				buffer << type << delimiter;
			}
			buffer << to_string(improper.psi0) << delimiter << to_string(improper.kpsi) << endl;
		}

		buffer << FFPrintHelpers::endBlock();  // Include endBlock in the buffer

		file << buffer.str();  // Write the entire string to the file at once
	}

	file.close();
}

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


	return files;
}





void LimaForcefieldBuilder::buildForcefield(const std::string& molecule_dir, const std::string& output_dir,
	const ParsedTopologyFile& topol_file, EnvMode envmode)
{
	
	//LimaLogger logger{ LimaLogger::LogMode::compact, envmode, "forcefieldmaker", work_dir };
	const bool verbose(envmode != Headless);
	//const string conf_path = lfs::pathJoin(molecule_dir, conf_name);
	//const string topol_path = lfs::pathJoin(molecule_dir, topol_name);
	//lfs::assertPath(conf_path);
	//lfs::assertPath(topol_path);


	const char ignore_atomtype = IGNORE_HYDROGEN ? 'H' : '.';
	// Check if filtered files already exists, if so return
	BondedTypes forcefield;

	std::vector<std::string> files = getFiles();

	for (auto& file_path : files) {
		const SimpleParsedFile parsedfile = lfs::parsePrmFile(file_path, verbose);
		loadFileIntoForcefield(parsedfile, forcefield);
	}
	
	Topology topology = loadTopology(molecule_dir, topol_file, ignore_atomtype);

	// Filter for the atomtypes used in this simulation and map to them
	const std::vector<NB_Atomtype> atomtypes_filtered = filterAtomtypes(topology, forcefield);
	const std::vector<AtomtypeMapping> atomtype_map = mapGroidsToSimulationspecificAtomtypeids(topology, atomtypes_filtered);

	// Match the topology with the forcefields.
	fillTBondParametersFromForcefield(forcefield, topology);


	printFFNonbonded(lfs::pathJoin(output_dir, "ffnonbonded.lff"), atomtype_map, atomtypes_filtered);
	printFFBonded(lfs::pathJoin(output_dir, "ffbonded.lff"), topology);
	//logger.finishSection("Prepare Forcefield has finished");
}