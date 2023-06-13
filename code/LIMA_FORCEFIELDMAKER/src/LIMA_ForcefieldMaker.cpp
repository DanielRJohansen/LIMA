#include <vector>

#include "LIMA_FORCEFIELDMAKER/include/ForcefieldTypes.h"
#include "LIMA_FORCEFIELDMAKER/include/Filehandling.h"
#include "LIMA_FORCEFIELDMAKER/include/ForcefieldMerger.h"
#include "LIMA_FORCEFIELDMAKER/include/ForcefieldMaker.h"

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

ForcefieldMaker::ForcefieldMaker(const string& workdir, EnvMode envmode, const string& default_ff_dir, const string& conf_file, const string& topol_file) :
	molecule_dir(FileHelpers::pathJoin(workdir, "molecule")),
	forcefield_dir(default_ff_dir), 
	logger(LimaLogger::LogMode::compact, envmode, "forcefieldmaker", workdir),
	m_verbose(envmode != Headless)
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
	vector<vector<string>> simconf_rows = Reader::readFile(conf_path, logger);
	vector<string> simconf = FTHelpers::parseConf(simconf_rows);
	logger.print(std::format("{} atom types in conf\n", simconf.size()));


	vector<vector<string>> ffnonbonded_rows = Reader::readFile(ff_nonbonded_path, logger, {'/'}, true);


	auto forcefieldInsertFunction = [](const std::vector<string>& row, vector<NB_Atomtype>& records) {
		if (records.empty()) {
			records.push_back(NB_Atomtype("WATER", 0, water_mass, water_sigma, water_epsilon));		// Solvent type always first!
		}

		records.push_back(NB_Atomtype(
			row[0],
			stof(row[1]),
			stof(row[2]),
			stof(row[3])
		));
	};
	vector<NB_Atomtype> ffnonbonded = FTHelpers::parseFFBondtypes<NB_Atomtype, FTHelpers::STATE::FF_NONBONDED>(ffnonbonded_rows, forcefieldInsertFunction);
	//vector<NB_Atomtype> ffnonbonded = NB_Atomtype::parseNonbonded(ffnonbonded_rows);
	logger.print(std::format("{} atom types read from file\n", ffnonbonded.size()));

	return NB_Atomtype::filterUnusedTypes(ffnonbonded, simconf, *map, logger, false);
}

vector<Atom> makeTopologyAtoms(vector<vector<string>> topology_rows, vector<NB_Atomtype>* ff_nonbonded_active, Map* map, bool verbose) {
	vector<Atom> atoms = Atom::parseTopolAtoms(topology_rows, verbose);
	Atom::assignAtomtypeIDs(atoms, *ff_nonbonded_active, *map);
	return atoms;
}

vector<Singlebondtype> makeTopologyBonds(vector<vector<string>>* ffbonded_rows, vector<vector<string>>* topology_rows, vector<Atom>* atoms, bool verbose) {
	auto forcefieldInsertFunction = [](const std::vector<string>& row, vector<Singlebondtype>& records) {
		std::array<string, 2> bonded_typenames{ row[0], row[1] };
		records.push_back(Singlebondtype(bonded_typenames, stof(row[2]), stof(row[3])));
	};
	vector<Singlebondtype> forcefield = FTHelpers::parseFFBondtypes<Singlebondtype, FTHelpers::STATE::FF_BONDTYPES>(*ffbonded_rows, forcefieldInsertFunction);

	auto topologyInsertFunction = [](const std::vector<string>& row, vector<Singlebondtype>& records) {
		std::array<int, 2> gro_ids{ stoi(row[0]), stoi(row[1]) };
		records.push_back(Singlebondtype(gro_ids));
	};
	vector<Singlebondtype> topology_bonds = FTHelpers::parseTopolBondtypes<Singlebondtype, FTHelpers::STATE::FF_BONDTYPES>(*topology_rows, topologyInsertFunction);
	


	Singlebondtype::assignTypesFromAtomIDs(topology_bonds, *atoms);
	FTHelpers::assignForceVariablesFromForcefield(&topology_bonds, &forcefield);

	return topology_bonds;
}

vector<Anglebondtype> makeTopologyAngles(vector<vector<string>>* ffbonded_rows, vector<vector<string>>* topology_rows, vector<Atom>* atoms, bool verbose) {
	auto forcefieldInsertFunction = [](const std::vector<string>& row, vector<Anglebondtype>& records) {
		std::array<string, 3> bonded_typenames{ row[0], row[1], row[2]};
		records.push_back(Anglebondtype(bonded_typenames, stof(row[3]), stof(row[4])));
	};
	vector<Anglebondtype> forcefield = FTHelpers::parseFFBondtypes<Anglebondtype, FTHelpers::STATE::FF_ANGLETYPES>(*ffbonded_rows, forcefieldInsertFunction);

	auto topologyInsertFunction = [](const std::vector<string>& row, vector<Anglebondtype>& records) {
		std::array<int, 3> gro_ids{ stoi(row[0]), stoi(row[1]), stoi(row[2]) };
		records.push_back(Anglebondtype(gro_ids));
	};
	vector<Anglebondtype> topology_angles = FTHelpers::parseTopolBondtypes<Anglebondtype, FTHelpers::STATE::FF_ANGLETYPES>(*topology_rows, topologyInsertFunction);

	Anglebondtype::assignTypesFromAtomIDs(topology_angles, *atoms);
	FTHelpers::assignForceVariablesFromForcefield(&topology_angles, &forcefield);

	return topology_angles;
}

vector<Dihedralbondtype> makeTopologyDihedrals(vector<vector<string>> ffbonded_rows, vector<vector<string>> topology_rows, vector<Atom> atoms, bool verbose) {
	auto forcefieldInsertFunction = [](const std::vector<string>& row, vector<Dihedralbondtype>& records) {
		std::array<string, 4> bonded_typenames{ row[0], row[1], row[2], row[3]};
		records.push_back(Dihedralbondtype(bonded_typenames, stof(row[4]), stof(row[5]), stoi(row[6])));
		//records.push_back(Dihedralbondtype(row[0], row[1], row[2], row[3], stof(row[4]), stof(row[5]), stoi(row[6])));
	};
	vector<Dihedralbondtype> forcefield = FTHelpers::parseFFBondtypes<Dihedralbondtype, FTHelpers::STATE::FF_DIHEDRALTYPES>(ffbonded_rows, forcefieldInsertFunction);
	//vector<Dihedralbondtype> forcefield = Dihedralbondtype::parseFFDihedraltypes(ffbonded_rows, verbose);

	auto topologyInsertFunction = [](const std::vector<string>& row, vector<Dihedralbondtype>& records) {
		std::array<int, 4> gro_ids{ stoi(row[0]), stoi(row[1]), stoi(row[2]), stoi(row[3])};
		//records.push_back(Anglebondtype(gro_ids));
		//records.push_back(Dihedralbondtype(stoi(row[0]), stoi(row[1]), stoi(row[2]), stoi(row[3])));
		records.push_back(Dihedralbondtype(gro_ids));
	};
	vector<Dihedralbondtype> topology_dihedrals = FTHelpers::parseTopolBondtypes<Dihedralbondtype, FTHelpers::STATE::FF_DIHEDRALTYPES>(topology_rows, topologyInsertFunction);
	//vector<Dihedralbondtype> topology_dihedrals = Dihedralbondtype::parseTopolDihedraltypes(topology_rows, verbose);

	Dihedralbondtype::assignTypesFromAtomIDs(topology_dihedrals, atoms);
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
	vector<vector<string>> ffbonded_rows = Reader::readFile(ff_bonded_path, logger, { '/' }, true);
	vector<vector<string>> topology_rows = Reader::readFile(topol_path, logger);

	vector<Atom> atoms = makeTopologyAtoms(topology_rows, &ff_nonbonded_active, &map, m_verbose);

	vector<Singlebondtype> topology_bonds = makeTopologyBonds(&ffbonded_rows, &topology_rows, &atoms, m_verbose);

	vector<Anglebondtype> topology_angles = makeTopologyAngles(&ffbonded_rows, &topology_rows, &atoms, m_verbose);

	vector<Dihedralbondtype> topology_dihedrals = makeTopologyDihedrals(ffbonded_rows, topology_rows, atoms, m_verbose);



	Printer::printForcefieldSummary(FileHelpers::pathJoin(molecule_dir, "LIMA_ffnonbonded_filtered.txt"), ff_nonbonded_active, &map);


	Printer::printForcefield(
		FileHelpers::pathJoin(molecule_dir, "LIMA_ffbonded_filtered.txt"),
		atoms,
		topology_bonds,
		topology_angles,
		topology_dihedrals
	);
	logger.finishSection("Prepare Forcefield has finished");
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