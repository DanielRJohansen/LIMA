#include <vector>

#include "ForcefieldTypes.h"
#include "Filehandling.h"
#include "ForcefieldMerger.h"


#include "ForcefieldMaker.h"
#include <filesystem>
#include <assert.h>
#include <format>


using std::string;

//#ifdef __linux__
//	string sim_path = "../../Simulation";
//	//string sim_path = "/home/lima/Desktop/LIMA/Simulation";
//#else
//	string sim_path = "C:\\PROJECTS\\QUANTOM\\Simulation";
//#endif

//string mol_path = FileHelpers::pathJoin(sim_path, "Molecule");
//string forcefield_path = FileHelpers::pathJoin(sim_path, "Forcefield");

ForcefieldMaker::ForcefieldMaker( string workdir, string default_ff_dir, string conf_file, string topol_file) :
	molecule_dir(FileHelpers::pathJoin(workdir, "molecule")),
	forcefield_dir(default_ff_dir)
{
	ff_bonded_path = FileHelpers::pathJoin(default_ff_dir, "LIMA_ffbonded.txt");
	ff_nonbonded_path = FileHelpers::pathJoin(default_ff_dir, "LIMA_ffnonbonded.txt");
	assertPath(ff_bonded_path);
	assertPath(ff_nonbonded_path);

	conf_path = FileHelpers::pathJoin(molecule_dir, conf_file);
	topol_path = FileHelpers::pathJoin(molecule_dir, topol_file);
	assertPath(conf_path);
	assertPath(topol_path);
}


vector<NB_Atomtype> ForcefieldMaker::makeFilteredNonbondedFF(Map* map) {
//	vector<vector<string>> simconf_rows = readFile("C:\\PROJECTS\\Quantom\\Molecules\\t4lys_full\\conf.gro");
	vector<vector<string>> simconf_rows = Reader::readFile(conf_path);
	vector<string> simconf = FTHelpers::parseConf(simconf_rows);

	vector<vector<string>> ffnonbonded_rows = Reader::readFile(ff_nonbonded_path, {'/'}, true);
	vector<NB_Atomtype> ffnonbonded = NB_Atomtype::parseNonbonded(ffnonbonded_rows);

	return NB_Atomtype::filterUnusedTypes(ffnonbonded, simconf, map);
}

vector<Atom> makeTopologyAtoms(vector<vector<string>> topology_rows, vector<NB_Atomtype>* ff_nonbonded_active, Map* map) {
	vector<Atom> atoms = Atom::parseTopolAtoms(topology_rows);
	Atom::assignAtomtypeIDs(&atoms, ff_nonbonded_active, map);
	return atoms;
}

vector<Bondtype> makeTopologyBonds(vector<vector<string>>* ffbonded_rows, vector<vector<string>>* topology_rows, vector<Atom>* atoms) {
	vector<Bondtype> forcefield = Bondtype::parseFFBondtypes(*ffbonded_rows);
	vector<Bondtype> topology_bonds = Bondtype::parseTopolBondtypes(*topology_rows);

	Bondtype::assignTypesFromAtomIDs(&topology_bonds, *atoms);
	//Bondtype::assignFFParametersFromBondtypes(&topology_bonds, &ffbondtypes);
	FTHelpers::assignForceVariablesFromForcefield(&topology_bonds, &forcefield);

	return topology_bonds;
}

vector<Angletype> makeTopologyAngles(vector<vector<string>>* ffbonded_rows, vector<vector<string>>* topology_rows, vector<Atom>* atoms) {
	vector<Angletype> forcefield = Angletype::parseFFAngletypes(*ffbonded_rows);
	vector<Angletype> topology_angles = Angletype::parseTopolAngletypes(*topology_rows);

	Angletype::assignTypesFromAtomIDs(&topology_angles, *atoms);
	FTHelpers::assignForceVariablesFromForcefield(&topology_angles, &forcefield);

	return topology_angles;
}

vector<Dihedraltype> makeTopologyDihedrals(vector<vector<string>> ffbonded_rows, vector<vector<string>> topology_rows, vector<Atom> atoms) {
	vector<Dihedraltype> forcefield = Dihedraltype::parseFFDihedraltypes(ffbonded_rows);
	vector<Dihedraltype> topology_dihedrals = Dihedraltype::parseTopolDihedraltypes(topology_rows);

	Dihedraltype::assignTypesFromAtomIDs(&topology_dihedrals, atoms);
	FTHelpers::assignForceVariablesFromForcefield(&topology_dihedrals, &forcefield);

	return topology_dihedrals;
}


void ForcefieldMaker::prepSimulationForcefield() {
	// 95% of runtime is spent matching topologies to forcefield. 
	// To speedup, split the topology_x_types and run multithreaded?

	// TODO: Do a check whether we are building a ff on the same base files, if so, then return

	Map map;
	vector<NB_Atomtype> ff_nonbonded_active = makeFilteredNonbondedFF(&map);		// The map is made here, so this function must come first



	// These two vectors contain information about all types, but are parsed individually. This is very slightly slower, but MUCH more readable!
//	vector<vector<string>> ffbonded_rows = readFile("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\ffbonded.itp");
	//vector<vector<string>> ffbonded_rows = readFile(sim_path + "/Forcefield/ffbonded.itp");
	vector<vector<string>> ffbonded_rows = Reader::readFile(ff_bonded_path, { '/' }, true);
	//	vector<vector<string>> topology_rows = readFile("C:\\PROJECTS\\Quantom\\Molecules\\t4lys_full\\topol.top");
		//vector<vector<string>> topology_rows = readFile(sim_path + "/Molecule/topol.top");
	vector<vector<string>> topology_rows = Reader::readFile(topol_path);


	vector<Atom> atoms = makeTopologyAtoms(topology_rows, &ff_nonbonded_active, &map);



	vector<Bondtype> topology_bonds = makeTopologyBonds(&ffbonded_rows, &topology_rows, &atoms);

	vector<Angletype> topology_angles = makeTopologyAngles(&ffbonded_rows, &topology_rows, &atoms);

	vector<Dihedraltype> topology_dihedrals = makeTopologyDihedrals(ffbonded_rows, topology_rows, atoms);



	Printer::printForcefieldSummary(FileHelpers::pathJoin(molecule_dir, "LIMA_ffnonbonded_filtered.txt"), ff_nonbonded_active, &map);


	Printer::printForcefield(
		FileHelpers::pathJoin(molecule_dir, "LIMA_ffbonded_filtered.txt"),
		atoms,
		topology_bonds,
		topology_angles,
		topology_dihedrals
	);
}

//void mergeForcefieldFiles() {
//	vector<string> files;
//
//	// Some files are commented out because it is NOT clear whether they are using rmin or rmin/2
//
//	//files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all35_ethers.prm");
//	//files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_carb.prm");
//#ifdef __linux__
//	files.push_back(FileHelpers::pathJoin(forcefield_path, "par_all36_lipid.prm"));
//	files.push_back(FileHelpers::pathJoin(forcefield_path, "par_all36_na.prm"));
//	files.push_back(FileHelpers::pathJoin(forcefield_path, "par_all36m_prot.prm"));
//#else
//	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_lipid.prm");
//	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_na.prm");
//	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36m_prot.prm");
//#endif
//
//
//	//files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_cgenff.prm");	// CHARMM wants this to be read last for some reason??
//
//
//	ForcefieldMerger FM;
//	FM.mergeForcefields(files);
//	Printer::printFFNonbonded(FileHelpers::pathJoin(forcefield_path, "LIMA_ffnonbonded.txt"), FM.ff_nonbonded);
//	Printer::printFFBonded(FileHelpers::pathJoin(forcefield_path, "LIMA_ffbonded.txt"), FM.ff_bondtypes, FM.ff_angletypes, FM.ff_dihedraltypes);
//}


