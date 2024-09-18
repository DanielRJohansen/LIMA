#pragma once

#include <filesystem>
#include <fstream>
#include <algorithm>

#include "TestUtils.h"
#include "Programs.h"

namespace TestMembraneBuilder {
	using namespace TestUtils;
	namespace fs = std::filesystem;

	static LimaUnittestResult testBuildmembraneSmall(EnvMode envmode, bool do_em)
	{
		const fs::path work_dir = simulations_dir / "BuildMembraneSmall";
		const fs::path mol_dir = work_dir / "molecule";

		LipidsSelection lipidselection;
		const std::array<std::string, 6> lipids = { "POPC", "POPE", "DDPC", "DMPC", "cholesterol", "DOPC" };
		for (const auto& lipidname : lipids) {
			lipidselection.emplace_back(LipidSelect{ lipidname, work_dir, lipidname == "POPC" ? 50. : 10.});	// 10% of each lipid, except 50% POPC
		}
		auto [gro, top] = Programs::CreateMembrane(work_dir, lipidselection, Float3{ 7.f }, 3.5f, envmode);
		gro->printToFile(mol_dir / "membrane.gro");
		top->printToFile(mol_dir / "membrane.top");

		TopologyFile newTop{ mol_dir / "membrane.top" };
		TopologyFile refTop{ mol_dir / "membrane_reference.top" };

		ASSERT(newTop.GetAllAtoms() == refTop.GetAllAtoms(), "Topology Atom Mismatch");
		ASSERT(newTop.GetAllSinglebonds() == refTop.GetAllSinglebonds(), "Topology Singlebond Mismatch");
		ASSERT(newTop.GetAllPairs() == refTop.GetAllPairs(), "Topology Pair Mismatch");
		ASSERT(newTop.GetAllAnglebonds() == refTop.GetAllAnglebonds(), "Topology Anglebond Mismatch");
		ASSERT(newTop.GetAllDihedralbonds() == refTop.GetAllDihedralbonds(), "Topology Dihedralbond Mismatch");
		ASSERT(newTop.GetAllImproperDihedralbonds() == refTop.GetAllImproperDihedralbonds(), "Topology Improper Mismatch");

		GroFile newGro{ mol_dir / "membrane.gro" };
		GroFile refGro{ mol_dir / "membrane_reference.gro" };

		ASSERT(newGro.box_size == refGro.box_size, "Box size mismatch");
		ASSERT(newGro.atoms.size() == refGro.atoms.size(), "Atom count mismatch");
		for (int i = 0; i < newGro.atoms.size(); i++) {
			ASSERT(newGro.atoms[0].position == refGro.atoms[0].position, "Atom position mismatch");
		}

		return LimaUnittestResult{ true , "No error", envmode == Full};
	}

	LimaUnittestResult BuildAndRelaxVesicle(EnvMode envmode) {
		GroFile grofile;
		grofile.box_size = Float3{ 5.f };
		TopologyFile topfile;
		const fs::path workDir = TestUtils::simulations_dir / "etc";
		MoleculeHullCollection mhCol = Programs::MakeLipidVesicle(grofile, topfile, { {"POPC", workDir , 10}, {"cholesterol", workDir , 30}, {"DMPC", workDir , 60} }, 0.5, grofile.box_size/2.f, 3);

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