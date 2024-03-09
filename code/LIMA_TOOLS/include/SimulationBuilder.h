/// This file is used to build simulations, that can be saved to .gro and .top files.
/// This is NOT used for loading a simulation into LIMA in any ways
#pragma once

#include "MDFiles.h"


struct LipidSelect {
	LipidSelect(std::string lipidname, int percentage) : lipidname(lipidname), percentage(percentage) {
		if (std::find(valid_lipids.begin(), valid_lipids.end(), lipidname) == valid_lipids.end()) {
			throw std::runtime_error("LipidSelect: Lipid not supported (yet): " + lipidname);
		};
	}
	const std::string lipidname;
	const int percentage;
	std::unique_ptr<ParsedGroFile> grofile;
	std::unique_ptr<ParsedTopologyFile> topfile;

	static const std::array<std::string, 6> valid_lipids;	// Defined in .cpp file
};
// Since this is a vector of structs with unique_ptrs, it can never be copied, or resized
using LipidsSelection = std::vector<LipidSelect>;

struct AtomtypeSelect {
	const ParsedTopologyFile::AtomsEntry atomtype;
	const float percentage;
};
using AtomsSelection = std::vector<AtomtypeSelect>;

namespace SimulationBuilder {
	using Filepair = std::pair<ParsedGroFile, ParsedTopologyFile>;

	Filepair buildMembrane(const LipidsSelection& lipidselection, Float3 box_dims);
	Filepair makeBilayerFromMonolayer(const Filepair& inputfiles, Float3 box_dims);

	void DistributeParticlesInBox(ParsedGroFile& grofile, ParsedTopologyFile& topfile, const AtomsSelection& particles,
		float minDistBetweenAnyParticle=0.1f, float particlesPerNm3=32.f);

};