#pragma once

#include <vector>

enum class SecondaryStructure {
	Coil,
	Helix,
	Sheet,
};

// Static, simulation-independent description of one protein backbone point.
// The particle id refers to the C-alpha atom in Box particle numbering.
struct BackbonePoint {
	int particleId = -1;
	SecondaryStructure secondaryStructure = SecondaryStructure::Coil;
};

struct BackboneChain {
	std::vector<BackbonePoint> points;
};

using BackboneChains = std::vector<BackboneChain>;
