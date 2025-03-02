#pragma once

#include <filesystem>

#include "TestUtils.h"
#include "Programs.h"
#include "MoleculeHull.cuh"

namespace TestMembraneBuilder {
	using namespace TestUtils;
	namespace fs = std::filesystem;

	// This test ensures that the membrane is built identical to the reference membrane, NOT considering EM
	static LimaUnittestResult TestBuildmembraneSmall(EnvMode envmode, bool do_em)
	{		
		const fs::path work_dir = simulations_dir / "BuildMembraneSmall";
		const fs::path mol_dir = work_dir / "molecule";
		TestUtils::CleanDirIfNotContains(mol_dir, "reference");

		Lipids::Selection lipidselection;
		const std::array<std::string, 6> lipids = { "POPC", "POPE", "DDPC", "DMPC", "Cholesterol", "DOPC" };
		for (const auto& lipidname : lipids) {
			lipidselection.emplace_back(Lipids::Select{ lipidname, work_dir, lipidname == "POPC" ? 50. : 10.});	// 10% of each lipid, except 50% POPC
		}

		// Build the membrane, and write it to disk
		auto [gro, top] = SimulationBuilder::CreateMembrane(lipidselection, Float3{ 7.f }, 3.5f);
		gro->printToFile(mol_dir / "membrane.gro");
		top->printToFile(mol_dir / "membrane.top");

		// Test the topology is identical to reference
		TopologyFile newTop{ mol_dir / "membrane.top" };
		TopologyFile refTop{ mol_dir / "membrane_reference.top" };

		ASSERT(std::ranges::equal(newTop.GetAllElements<TopologyFile::AtomsEntry>(), refTop.GetAllElements<TopologyFile::AtomsEntry>()), "Topology Atom Mismatch");
		ASSERT(std::ranges::equal(newTop.GetAllElements<TopologyFile::SingleBond>(), refTop.GetAllElements<TopologyFile::SingleBond>()), "Topology Atom Mismatch");
		ASSERT(std::ranges::equal(newTop.GetAllElements<TopologyFile::PairBond>(), refTop.GetAllElements<TopologyFile::PairBond>()), "Topology Atom Mismatch");
		ASSERT(std::ranges::equal(newTop.GetAllElements<TopologyFile::AngleBond>(), refTop.GetAllElements<TopologyFile::AngleBond>()), "Topology Atom Mismatch");
		ASSERT(std::ranges::equal(newTop.GetAllElements<TopologyFile::DihedralBond>(), refTop.GetAllElements<TopologyFile::DihedralBond>()), "Topology Atom Mismatch");
		ASSERT(std::ranges::equal(newTop.GetAllElements<TopologyFile::ImproperDihedralBond>(), refTop.GetAllElements<TopologyFile::ImproperDihedralBond>()), "Topology Atom Mismatch");
		
		// Test the conf is identical to reference
		GroFile newGro{ mol_dir / "membrane.gro" };
		GroFile refGro{ mol_dir / "membrane_reference.gro" };
		ASSERT(newGro.box_size == refGro.box_size, "Box size mismatch");
		ASSERT(newGro.atoms.size() == refGro.atoms.size(), "Atom count mismatch");
		for (int i = 0; i < newGro.atoms.size(); i++) {
			auto b = (newGro.atoms[i].position - refGro.atoms[i].position).len();
			if (newGro.atoms[i].position != refGro.atoms[i].position)
				int a=0;
			ASSERT(newGro.atoms[i].position == refGro.atoms[i].position, "Atom position mismatch");
		}

		// Finally test if we can stabilize the simulation
		const float emtol = 200.f;
		auto sim = Programs::EnergyMinimize(*gro, *top, true, work_dir, envmode, true, emtol);

		return LimaUnittestResult{ sim->maxForceBuffer.back().second < emtol, std::format("Failed to energy minimize membrane {:.2f}/{:.2f}", sim->maxForceBuffer.back().second, emtol), envmode == Full};
	}

	static LimaUnittestResult TestBuildmembraneWithCustomlipidAndCustomForcefield(EnvMode envmode) {
		const fs::path work_dir = simulations_dir / "BuildMembraneCustom";
		const fs::path mol_dir = work_dir / "molecule";

		//TestUtils::CleanDirectory(mol_dir);
		fs::remove_all(mol_dir);
		fs::create_directory(mol_dir);

		Lipids::Selection lipidselection;
		const std::vector<std::pair<std::string, double>> lipids = { {"POPC", 70.}, {"CUST" , 30.} };
		for (const auto& [lipidname, percentage] : lipids) {
			lipidselection.emplace_back(Lipids::Select{ lipidname, work_dir, percentage });	// 10% of each lipid, except 50% POPC
		}

		auto [gro, top] = SimulationBuilder::CreateMembrane(lipidselection, Float3{ 7.f }, 3.5f);
		Programs::EnergyMinimize(*gro, *top, true, work_dir, envmode, true, 300000.f); // high emtol, because we dont care about EM, we just want to see if the simulation can even start

		gro->printToFile(mol_dir / "membrane.gro");
		top->printToFile(mol_dir / "membrane.top");

		TopologyFile newTop{ mol_dir / "membrane.top" };
		GroFile newGro{ mol_dir / "membrane.gro" };

		ASSERT(std::ranges::equal(newTop.GetAllElements<TopologyFile::AtomsEntry>(), top->GetAllElements<TopologyFile::AtomsEntry>()), "Topology Atom Mismatch");


		SimParams params{};
		params.em_variant = true;
		Environment env(work_dir, envmode);
		env.CreateSimulation(newGro, newTop, params);

		return LimaUnittestResult{ true , "No error", envmode == Full };
	}

	LimaUnittestResult TestAllStockholmlipids(EnvMode envmode) {
		const fs::path work_dir = simulations_dir / "BuildMembraneSmall";

		const fs::path path = FileUtils::GetLimaDir() / "resources/Slipids";
		std::vector<std::string> targets;
		for (const auto& entry : fs::directory_iterator(path)) {
			if (entry.path().extension() == ".gro") {
				std::string base_name = entry.path().stem().string();
				fs::path itp_file = path / (base_name + ".itp");
				if (fs::exists(itp_file)) {
					targets.push_back(base_name);
				}
			}
		}

		Lipids::Selection lipidselection;
		for (const auto& lipidname : targets) {
			lipidselection.emplace_back(Lipids::Select{ lipidname, work_dir, 100. / static_cast<double>(targets.size())});	// 10% of each lipid, except 50% POPC
		}

		// The first test is pretty much just to see if this function throws
		auto [grofile, topfile] = SimulationBuilder::CreateMembrane(lipidselection, Float3{ 10.f }, 5.f);
		/*for (const auto& molecule : topfile->GetSystem().molecules) {
			ASSERT(molecule.moleculetype->readFromCache, "This lipid top should have been read from a cached file");
		}*/

		// The third test is to see if this function throws
		const float emtol = 1000.f;
		auto sim = Programs::EnergyMinimize(*grofile, *topfile, false, work_dir, envmode, true, emtol);

		ASSERT(sim->maxForceBuffer.back().second < emtol, "Failed to energy minimize membrane");

		return LimaUnittestResult{ true , "", envmode == Full };
	}

	LimaUnittestResult BuildAndRelaxVesicle(EnvMode envmode) {
		GroFile grofile;
		grofile.box_size = Float3{ 5.f };
		TopologyFile topfile;
		const fs::path workDir = TestUtils::simulations_dir / "etc";
		MoleculeHullCollection mhCol = Programs::MakeLipidVesicle(grofile, topfile, { {"POPC", workDir , 10}, {"Cholesterol", workDir , 30}, {"DMPC", workDir , 60} }, 0.5, grofile.box_size/2.f, 3);

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


		Programs::MoveMoleculesUntillNoOverlap(mhCol, grofile.box_size, envmode==Full);
		
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