/// This file is used to build simulations, that can be saved to .gro and .top files.
/// This is NOT used for loading a simulation into LIMA in any ways
#pragma once

#include "MDFiles.h"




namespace SimulationBuilder {
	using Filepair = std::pair<ParsedGroFile, ParsedTopologyFile>;

	Filepair buildMembrane(Filepair inputfiles, Float3 box_dims);

	Filepair makeBilayerFromMonolayer(const Filepair& inputfiles, Float3 box_dims);
};