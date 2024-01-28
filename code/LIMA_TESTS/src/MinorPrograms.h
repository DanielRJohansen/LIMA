#pragma once

#include "TestUtils.h"
#include "MoleculeGraph.h"


namespace TestMinorPrograms {
	using namespace TestUtils;
	namespace fs = std::filesystem;

	// This is not ready to be a test, for now i just need it to create some lipids files
	LimaUnittestResult testReorderMoleculeParticles(EnvMode envmode) {
		//const std::string to_folder = "C:/Users/Daniel/git_repo/LIMA/resources/Lipids/POPC/";
		//const std::string from_folder = "C:/PROJECTS/Quantom/Molecules/Lipids/POPC/";
		const std::string to_folder = "C:/Users/Daniel/git_repo/LIMA_data/ReorderPOPC/";
		const std::string from_folder = "C:/PROJECTS/Quantom/Molecules/Lipids/POPC/";
		LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(from_folder +"POPC.gro", from_folder + "POPC.itp", to_folder + "POPC.gro", to_folder + "POPC.itp");

		std::string error = compareFilesBitwise(fs::path(to_folder) / "POPC.gro", fs::path(to_folder) / "POPC_reference.gro");
		if (error != "") {
			return LimaUnittestResult{ LimaUnittestResult::FAIL , error, envmode == Full };
		}
		error = compareFilesBitwise(fs::path(to_folder) / "POPC.itp", fs::path(to_folder) / "POPC_reference.itp");
		if (error != "") {
			return LimaUnittestResult{ LimaUnittestResult::FAIL , error, envmode == Full };
		}

	return LimaUnittestResult{ LimaUnittestResult::SUCCESS , "No error", envmode == Full };
	}
}

