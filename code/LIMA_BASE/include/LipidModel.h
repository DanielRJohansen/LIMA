#pragma once

#include <vector>

enum class LipidHeadgroup {
	Choline,		// PC, SM
	Ethanolamine,	// PE
	Serine,			// PS
	Anionic,		// PG, PA, PI, ... any phosphate without nitrogen
	Sterol,			// No phosphate, the head is the first oxygen
};

// Static, simulation-independent description of one lipid for the lipid representation.
// Particle ids are in Box particle numbering.
struct LipidModel {
	int headParticleId = -1;
	int firstParticleId = -1;
	int particleCount = 0;
	LipidHeadgroup headgroup = LipidHeadgroup::Anionic;
	// Heavy atom paths, ordered towards the chain end. The first tail starts at the head, the following start at the
	// atom where they branch off an earlier tail
	std::vector<std::vector<int>> tails;
};

using LipidModels = std::vector<LipidModel>;
