#pragma once

#include <string>
#include "Utilities.h"

struct Map;
struct NB_Atomtype;

class ForcefieldMaker {
public:
	ForcefieldMaker(
		const std::string& workdir,
		const std::string& default_ff_dir = "C:/Users/Daniel/git_repo/LIMA/resources/Forcefields/charm36/",
		const std::string& conf_file = "conf.gro",
		const std::string& topol_file = "topol.top"
	);

	void prepSimulationForcefield();
	void mergeForcefiledFiles();	// TODO: implement?



private:
	const std::string molecule_dir;
	const std::string forcefield_dir;

	std::string conf_path = "";
	std::string topol_path = "";
	std::string ff_bonded_path = "";
	std::string ff_nonbonded_path = "";
	

	LimaLogger logger;





	vector<NB_Atomtype> makeFilteredNonbondedFF(Map* map);

};