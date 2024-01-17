#pragma once

#include <string>

#include "Utilities.h"
#include "LimaTypes.cuh"
#include "MDFiles.h"

struct NB_Atomtype;

namespace LimaForcefieldBuilder {

	/// <summary>
	/// Create ffbonded.lff and ffnonbonded.lff files.
	/// </summary>
	/// <param name="molecule_dir">Dir where conf and topol are, with all .itp include files</param>
	/// <param name="output_dir">Where .lff files will be created</param>
	/// <param name="conf_name">name of main .gro file in molecule_dir</param>
	/// <param name="topol_name">name of main .top/.itp file in molecule_dir</param>
	/// <param name="envmode"></param>
	void buildForcefield(const std::string& molecule_dir, const std::string& output_dir,
		const ParsedTopologyFile& topol_file, EnvMode envmode);
}