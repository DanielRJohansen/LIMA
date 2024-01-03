#pragma once

#include "TestUtils.h"
#include "MoleculeGraph.h"


namespace TestMinorPrograms {
	using namespace TestUtils;

	// This is not ready to be a test, for now i just need it to create some lipids files
	LimaUnittestResult testReorderMoleculeParticles() {
		const std::string to_folder = "C:/Users/Daniel/git_repo/LIMA/resources/Lipids/POPC/";
		const std::string from_folder = "C:/PROJECTS/Quantom/Molecules/Lipids/POPC/";
		LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(from_folder +"POPC.gro", from_folder + "POPC.itp", to_folder + "POPC.gro", to_folder + "POPC.itp");

		return LimaUnittestResult(LimaUnittestResult::SUCCESS, "", true);
	}
}

