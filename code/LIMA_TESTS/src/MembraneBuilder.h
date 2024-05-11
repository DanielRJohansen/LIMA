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
		const fs::path work_dir = simulations_dir + "/BuildMembraneSmall";
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

		ParsedTopologyFile newTop{ mol_dir / "monolayer.top" };
		ParsedTopologyFile refTop{ mol_dir / "monolayer_reference.top" };

		if (newTop.GetAtoms() != refTop.GetAtoms()) {
			return LimaUnittestResult{ LimaUnittestResult::FAIL , "Atom mismatch", envmode == Full };
		}
		if (newTop.GetSinglebonds() != refTop.GetSinglebonds()) {
			return LimaUnittestResult{ LimaUnittestResult::FAIL , "Singlebond mismatch", envmode == Full };
		}
		if (newTop.GetPairs() != refTop.GetPairs()) {
			return LimaUnittestResult{ LimaUnittestResult::FAIL , "Pair mismatch", envmode == Full };
		}
		if (newTop.GetAnglebonds() != refTop.GetAnglebonds()) {
			return LimaUnittestResult{ LimaUnittestResult::FAIL , "Anglebond mismatch", envmode == Full };
		}
		if (newTop.GetDihedralbonds() != refTop.GetDihedralbonds()) {
			return LimaUnittestResult{ LimaUnittestResult::FAIL , "Dihedralbond mismatch", envmode == Full };
		}
		if (newTop.GetImproperDihedralbonds() != refTop.GetImproperDihedralbonds()) {
			return LimaUnittestResult{ LimaUnittestResult::FAIL , "Improper mismatch", envmode == Full };
		}

		return LimaUnittestResult{ LimaUnittestResult::SUCCESS , "No error", envmode == Full};
	}
}