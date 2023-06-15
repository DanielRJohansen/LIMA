//#include "ForcefieldMerger.h"

void ForcefieldMaker::mergeForcefieldFiles() {
	vector<string> files;

	// Some files are commented out because it is NOT clear whether they are using rmin or rmin/2

#ifdef __linux__
	throw "Add the other files before trying this on linux";
	files.push_back(FileHelpers::pathJoin(forcefield_path, "par_all36_lipid.prm"));
	files.push_back(FileHelpers::pathJoin(forcefield_path, "par_all36_na.prm"));
	files.push_back(FileHelpers::pathJoin(forcefield_path, "par_all36m_prot.prm"));
#else
	//files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all35_ethers.prm");
	//files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_carb.prm");

	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_lipid.prm");
	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_na.prm");
	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36m_prot.prm");
#endif


	//files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_cgenff.prm");	// CHARMM wants this to be read last for some reason??

	const std::string bonded_path = "C:/Users/Daniel/git_repo/LIMA/resources/Forcefields/charm36/LIMA_ffbonded.txt";
	const std::string nonbonded_path = "C:/Users/Daniel/git_repo/LIMA/resources/Forcefields/charm36/LIMA_ffnonbonded.txt";

	ForcefieldMerger FM;
	FM.mergeForcefields(files);
	Printer::printFFNonbonded(nonbonded_path, FM.ff_nonbonded);
	Printer::printFFBonded(bonded_path, FM.ff_bondtypes, FM.ff_angletypes, FM.ff_dihedraltypes, FM.ff_improperdihedraltypes);
}