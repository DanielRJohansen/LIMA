/// This file is used to build simulations, that can be saved to .gro and .top files.
/// This is NOT used for loading a simulation into LIMA in any ways
#pragma once

#include "MDFiles.h"
#include <algorithm>

struct LipidSelect {
	LipidSelect(std::string lipidname, int percentage) : lipidname(lipidname), percentage(percentage) {
		if (std::find(valid_lipids.begin(), valid_lipids.end(), lipidname) == valid_lipids.end()) {
			throw std::runtime_error("LipidSelect: Lipid not supported (yet): " + lipidname);
		};

		const fs::path lipid_path = Filehandler::GetLimaDir() / ("resources/Lipids/" + lipidname);
		grofile = std::make_unique<GroFile>(lipid_path / (lipidname + ".gro"));
		topfile = std::make_unique<TopologyFile>(lipid_path / (lipidname + ".itp"));
	}
	const std::string lipidname;
	const int percentage;
	std::shared_ptr<GroFile> grofile;
	std::shared_ptr<TopologyFile> topfile;

	static const std::array<std::string, 6> valid_lipids;	// Defined in .cpp file
};
// Since this is a vector of structs with unique_ptrs, it can never be copied, or resized
using LipidsSelection = std::vector<LipidSelect>;

struct AtomtypeSelect {
	const TopologyFile::AtomsEntry atomtype;
	const float percentage;
};
using AtomsSelection = std::vector<AtomtypeSelect>;

namespace SimulationBuilder {
	using namespace MDFiles;



	FilePair buildMembrane(const LipidsSelection& lipidselection, Float3 box_dims);
	FilePair makeBilayerFromMonolayer(const FilePair& inputfiles, Float3 box_dims);

	void DistributeParticlesInBox(GroFile& grofile, TopologyFile& topfile, const AtomsSelection& particles,
		float minDistBetweenAnyParticle=0.1f, float particlesPerNm3=32.f);


	void SolvateGrofile(GroFile& grofile);

	/// <summary>
	/// Insert submolecules into a possibly empty simulation. The molecules will be spread randomly around in the box
	/// The input grofile does NOT need to be centered around origo, but the grofile must have the atoms
	/// in the correct position relative to eachother, as it uses raw positions, instead of hyperpos relative to index 0
	/// </summary>
	void InsertSubmoleculesInSimulation(GroFile& targetGrofile, TopologyFile& targetTopol,
		const GroFile& submolGro, const std::shared_ptr<TopologyFile>& submolTop, int nMoleculesToInsert);
	
	void InsertSubmoleculesOnSphere(
		GroFile& targetGrofile,
		TopologyFile& targetTopol,
		LipidsSelection,
		int nMoleculesToInsert,
		float sphereRadius,
		const Float3& sphereCenter
	);
};
