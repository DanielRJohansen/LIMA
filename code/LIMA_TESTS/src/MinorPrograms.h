#pragma once

#include "TestUtils.h"
#include "MoleculeUtils.h"
#include "SimulationBuilder.h"


namespace TestMinorPrograms {
	using namespace TestUtils;
	namespace fs = std::filesystem;

	LimaUnittestResult InsertMoleculesAndDoStaticbodyEM(EnvMode envmode) {
		const fs::path workDir = simulations_dir / "manyt4";
		
		const int boxSize = 25;
		const int nInsertions = 50;

		GroFile groSrc(workDir / "conf.gro");
		auto topSrc = std::make_shared<TopologyFile>(workDir / "topol.top");

		Programs::EnergyMinimize(groSrc, *topSrc, true, workDir, envmode, false);

		GroFile grofile{};
		grofile.box_size = Float3{ boxSize, boxSize, boxSize };
		TopologyFile topfile{};
		topfile.SetSystem("Many T4");

		SimulationBuilder::InsertSubmoleculesInSimulation(grofile, topfile, groSrc, topSrc, nInsertions, true);		
		ASSERT(grofile.atoms.size() == groSrc.atoms.size() * nInsertions, "Number of atoms in grofile is not correct");

		Programs::StaticbodyEnergyMinimize(grofile, topfile, envmode==Full);

		BoundingBox requiredBB{ Float3{} - grofile.box_size * 0.2, grofile.box_size * 1.2 };

		for (const auto& atom : grofile.atoms) {
			ASSERT(requiredBB.pointIsInBox(atom.position), "Atom outside of bounding box");
		}

		Programs::EnergyMinimize(grofile, topfile, false, workDir, envmode, false);

		return LimaUnittestResult{ true , "No error", envmode == Full };
	}




	// This is not ready to be a test, for now i just need it to create some lipids files
	//LimaUnittestResult testReorderMoleculeParticles(EnvMode envmode) {
	//	//const std::string to_folder = "C:/Users/Daniel/git_repo/LIMA/resources/Lipids/POPC/";
	//	//const std::string from_folder = "C:/PROJECTS/Quantom/Molecules/Lipids/POPC/";
	//	const std::string to_folder = "C:/Users/Daniel/git_repo/LIMA_data/ReorderPOPC/";
	//	const std::string from_folder = "C:/PROJECTS/Quantom/Molecules/Lipids/POPC/";
	//	LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(from_folder +"POPC.gro", from_folder + "POPC.itp", to_folder + "POPC.gro", to_folder + "POPC.itp");

	//	std::string error = compareFilesBitwise(fs::path(to_folder) / "POPC.gro", fs::path(to_folder) / "POPC_reference.gro");
	//	if (error != "") {
	//		return LimaUnittestResult{ false , error, envmode == Full };
	//	}
	//	error = compareFilesBitwise(fs::path(to_folder) / "POPC.itp", fs::path(to_folder) / "POPC_reference.itp");
	//	if (error != "") {
	//		return LimaUnittestResult{ false , error, envmode == Full };
	//	}

	//return LimaUnittestResult{ true , "No error", envmode == Full };
	//}
}

