#pragma once

#include <filesystem>
#include <fstream>
#include <algorithm>

#include "TestUtils.h"
#include "Programs.h"

namespace TestMembraneBuilder {
	using namespace TestUtils;
	namespace fs = std::filesystem;



	void thing(LipidsSelection& s) {

	}

	static LimaUnittestResult testBuildmembraneSmall(EnvMode envmode, bool do_em)
	{
		const fs::path work_dir = simulations_dir / "BuildMembraneSmall";
		Environment env{ work_dir.string(), envmode, false};

		env.CreateSimulation(7.f);
		LipidsSelection lipidselection;
		for (const auto& name : LipidSelect::valid_lipids) {
			lipidselection.emplace_back(LipidSelect{ name, name == "POPC" ? 50 : 10});	// 10% of each lipid, except 50% POPC
		}
		Programs::CreateMembrane(env, lipidselection, do_em, 3.5f, true);


		// Test that the output files match the reference output files
		const fs::path mol_dir = work_dir / "molecule";
		std::vector<std::array<std::string, 2>> files = { {"monolayer.gro", "monolayer_reference.gro"}, {"monolayer.top", "monolayer_reference.top"} };

		if (!do_em) {
			// These files are altered by the em, and thus the comparison cannot be made
			files.push_back({ "membrane.gro", "membrane_reference.gro" });
			files.push_back({ "membrane.top", "membrane_reference.top" });				
		}		

		TopologyFile newTop{ mol_dir / "monolayer.top" };
		TopologyFile refTop{ mol_dir / "monolayer_reference.top" };

		if (newTop.GetAllAtoms() != refTop.GetAllAtoms()) {
			return LimaUnittestResult{ false , "Atom mismatch", envmode == Full };
		}
		if (newTop.GetAllSinglebonds() != refTop.GetAllSinglebonds()) {
			return LimaUnittestResult{ false , "Singlebond mismatch", envmode == Full };
		}
		if (newTop.GetAllPairs() != refTop.GetAllPairs()) {
			return LimaUnittestResult{ false , "Pair mismatch", envmode == Full };
		}
		if (newTop.GetAllAnglebonds() != refTop.GetAllAnglebonds()) {
			return LimaUnittestResult{ false , "Anglebond mismatch", envmode == Full };
		}
		if (newTop.GetAllDihedralbonds() != refTop.GetAllDihedralbonds()) {
			return LimaUnittestResult{ false , "Dihedralbond mismatch", envmode == Full };
		}
		if (newTop.GetAllImproperDihedralbonds() != refTop.GetAllImproperDihedralbonds()) {
			return LimaUnittestResult{ false , "Improper mismatch", envmode == Full };
		}

		return LimaUnittestResult{ true , "No error", envmode == Full};
	}

	LimaUnittestResult BuildAndRelaxVesicle(EnvMode envmode) {
		GroFile grofile;
		grofile.box_size = Float3{ 5.f };
		TopologyFile topfile;
		MoleculeHullCollection mhCol = Programs::MakeLipidVesicle(grofile, topfile, { {"POPC", 10}, {"cholesterol", 30}, {"DMPC", 60} }, 0.5, grofile.box_size/2.f, 3);

		const bool overwriteData = false;

		// Compare before relaxation
		{
			std::vector<Facet> facets;
			GenericCopyToHost(mhCol.facets, facets, mhCol.nFacets);

			std::vector<Float3> vertices(mhCol.nFacets * 3);
			for (int i = 0; i < mhCol.nFacets; i++) {
				for (int j = 0; j < 3; j++) {
					vertices[i * 3 + j] = facets[i].vertices[j];
				}
			}

			if (!TestUtils::CompareVecWithFile(vertices, TestUtils::simulations_dir / fs::path{ "etc" } / "buildvesicle.bin", 0.01, overwriteData))
				return LimaUnittestResult{ false , "Before relaxation mismatch", envmode == Full };
		}


		Programs::MoveMoleculesUntillNoOverlap(mhCol, grofile.box_size);
		
		// Compare after relaxation
		{
			std::vector<Facet> facets;
			GenericCopyToHost(mhCol.facets, facets, mhCol.nFacets);

			std::vector<Float3> vertices(mhCol.nFacets * 3);
			for (int i = 0; i < mhCol.nFacets; i++) {
				for (int j = 0; j < 3; j++) {
					vertices[i * 3 + j] = facets[i].vertices[j];
				}
			}
			
			if (!TestUtils::CompareVecWithFile(vertices, TestUtils::simulations_dir / fs::path{ "etc" } / "relaxedvesicle.bin", 0.01, overwriteData))
				return LimaUnittestResult{ false , "After relaxation mismatch", envmode == Full };
		}

		if (overwriteData)
			return LimaUnittestResult{false, "Overwriting data", envmode == Full};



		return LimaUnittestResult{ true , "No error", envmode == Full };
	}
}