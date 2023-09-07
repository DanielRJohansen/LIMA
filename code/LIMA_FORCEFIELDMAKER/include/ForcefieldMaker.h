#pragma once

#include <string>

#include "Utilities.h"
#include "LimaTypes.cuh"

struct NB_Atomtype;

class ForcefieldMaker {

public:
	ForcefieldMaker(
		const std::string& workdir,
		EnvMode envmode,
		const std::string& conf_file = "conf.gro",
		const std::string& topol_file = "topol.top"
	);

	void prepSimulationForcefield(const char ignored_atomtype);



private:

	const std::string molecule_dir;
	//const std::string forcefield_dir;

	const bool m_verbose;

	int current_chain_id = -1;

	std::string conf_path = "";
	std::string topol_path = "";
	
	LimaLogger logger;
};
