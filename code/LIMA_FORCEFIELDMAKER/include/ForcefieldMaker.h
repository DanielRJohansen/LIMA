#pragma once

#include <string>
#include "LIMA_BASE/include/Utilities.h"
#include "LIMA_BASE/include/LimaTypes.cuh"

struct Map;
struct NB_Atomtype;


class ForcefieldMaker {
public:
	ForcefieldMaker(
		const std::string& workdir,
		EnvMode envmode,
		const std::string& ff_dir = "C:/Users/Daniel/git_repo/LIMA/resources/Forcefields/charm36/",
		const std::string& conf_file = "conf.gro",
		const std::string& topol_file = "topol.top"
	);

	void prepSimulationForcefield();


	/// <summary>
	/// Not part of normal operation. Call directly from main when needed
	/// </summary>
	static void mergeForcefieldFiles();	



private:
	const std::string molecule_dir;
	const std::string forcefield_dir;

	const bool m_verbose;


	std::string conf_path = "";
	std::string topol_path = "";
	std::string ff_bonded_path = "";
	std::string ff_nonbonded_path = "";
	
	LimaLogger logger;














	vector<NB_Atomtype> makeFilteredNonbondedFF(Map* map);

};