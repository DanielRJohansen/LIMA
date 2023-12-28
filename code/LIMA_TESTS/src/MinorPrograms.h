#pragma once

#include "TestUtils.h"
#include "MoleculeGraph.h"


namespace TestMinorPrograms {
	using namespace TestUtils;

	// This is not ready to be a test, for now i just need it to create some lipids files
	LimaUnittestResult testReorderMoleculeParticles() {
		LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains("C:/Users/Daniel/git_repo/LIMA/resources/Lipids/POPC/POPC.gro", 
			"C:/Users/Daniel/git_repo/LIMA/resources/Lipids/POPC/POPC.itp");

		return LimaUnittestResult(LimaUnittestResult::SUCCESS, "", true);
	}
}

