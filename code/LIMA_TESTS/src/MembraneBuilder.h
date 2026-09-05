#pragma once

#include <filesystem>
#include <numeric>

#include "Statistics.h"
#include "TestUtils.h"
#include "Programs.h"
#include "MoleculeHull.cuh"
#include "LiveEditCommands.h"

namespace TestMembraneBuilder {
	using namespace TestUtils;
	namespace fs = std::filesystem;

	static std::vector<Float3> ResidueCenters(const GroFile& grofile) {
		std::map<int, std::pair<Float3, int>> accumulators;
		for (const auto& atom : grofile.atoms) {
			auto& [sum, count] = accumulators[atom.residue_number];
			sum += atom.position;
			++count;
		}

		std::vector<Float3> centers;
		centers.reserve(accumulators.size());
		for (const auto& [_, value] : accumulators)
			centers.push_back(value.first / static_cast<float>(value.second));
		return centers;
	}

	static float NearestNeighborSpacingVariation(const std::vector<Float3>& directions) {
		std::vector<float> nearestDistances;
		nearestDistances.reserve(directions.size());
		for (size_t i = 0; i < directions.size(); ++i) {
			float nearestDistance = FLT_MAX;
			for (size_t j = 0; j < directions.size(); ++j) {
				if (i != j)
					nearestDistance = std::min(nearestDistance, (directions[i] - directions[j]).len());
			}
			nearestDistances.push_back(nearestDistance);
		}
		return Statistics::StdDev(nearestDistances) / Statistics::Mean(nearestDistances);
	}

	static float MinimumNearestNeighborSpacing(const std::vector<Float3>& points) {
		float minimumSpacing = FLT_MAX;
		for (size_t i = 0; i < points.size(); ++i) {
			for (size_t j = i + 1; j < points.size(); ++j)
				minimumSpacing = std::min(minimumSpacing, (points[i] - points[j]).len());
		}
		return minimumSpacing;
	}

	static LimaUnittestResult TestSphericalMembraneBuilder(EnvMode envmode) {
		const fs::path workDir = HeavyTestsDir() / "BuildMembraneSphere";
		Lipids::Selection lipids;
		lipids.emplace_back(Lipids::Select{ "DMPC", workDir, 100. });

		const float minimumRadius = SimulationBuilder::MinimumSphereRadius(lipids);
		bool rejectedSmallSphere = false;
		try {
			GroFile tooSmallGrofile;
			tooSmallGrofile.box_size = Float3{ 16.f };
			TopologyFile tooSmallTopfile;
			tooSmallTopfile.SetSystem("Membrane");
			SimulationBuilder::CreateMembrane(tooSmallGrofile, tooSmallTopfile, lipids,
				MembraneGeometry::Sphere{ Float3{ 8.f }, minimumRadius - 0.01f });
		}
		catch (const std::invalid_argument&) {
			rejectedSmallSphere = true;
		}
		ASSERT(rejectedSmallSphere, "A sphere below the lipid-dependent minimum radius was accepted");

		GroFile grofile;
		grofile.box_size = Float3{ 16.f };
		grofile.title = "Spherical membrane";
		TopologyFile topfile;
		topfile.SetSystem("Membrane");
		SimulationBuilder::CreateMembrane(grofile, topfile, lipids,
			MembraneGeometry::Sphere{ Float3{ 8.f }, minimumRadius + 1.f });
		ASSERT(!grofile.atoms.empty(), "Spherical membrane did not contain any atoms");
		ASSERT(!topfile.GetSystem().molecules.empty(), "Spherical membrane topology did not contain any molecules");
		for (const auto& atom : grofile.atoms) {
			ASSERT(std::isfinite(atom.position.x) && std::isfinite(atom.position.y) && std::isfinite(atom.position.z),
				"Spherical membrane contained a non-finite atom position");
		}
		std::vector<float> outerRadii;
		std::vector<float> innerRadii;
		for (const Float3& center : ResidueCenters(grofile)) {
			const float radius = (center - Float3{ 8.f }).len();
			(radius > minimumRadius + 1.f ? outerRadii : innerRadii).push_back(radius);
		}

		
		ASSERT(outerRadii.size() > 2 && innerRadii.size() > 2,
			"Could not identify both spherical membrane leaflets");
		ASSERT(Statistics::StdDev(outerRadii) > 0.02f && Statistics::StdDev(innerRadii) > 0.02f,
			"Spherical membrane radial roughness was too small");
		ASSERT(Statistics::StdDev(outerRadii) < 0.30f && Statistics::StdDev(innerRadii) < 0.30f,
			"Spherical membrane radial roughness was unreasonably large");
		std::vector<Float3> outerDirections;
		std::vector<Float3> innerDirections;
		for (const Float3& center : ResidueCenters(grofile)) {
			const Float3 offset = center - Float3{ 8.f };
			(offset.len() > minimumRadius + 1.f ? outerDirections : innerDirections)
				.push_back(offset.norm());
		}
		ASSERT(NearestNeighborSpacingVariation(outerDirections) < 0.20f
			&& NearestNeighborSpacingVariation(innerDirections) < 0.20f,
			"Equal-radius closed membrane did not preserve uniform spherical spacing");

		GroFile ellipsoidGrofile;
		ellipsoidGrofile.box_size = Float3{ 20.f };
		TopologyFile ellipsoidTopfile;
		ellipsoidTopfile.SetSystem("Membrane");
		const Float3 ellipsoidCenter{ 10.f };
		const Float3 ellipsoidRadii{ 5.f, 6.f, 7.f };
		SimulationBuilder::CreateMembrane(ellipsoidGrofile, ellipsoidTopfile, lipids,
			MembraneGeometry::Ellipsoid{ ellipsoidCenter, ellipsoidRadii });
		ASSERT(!ellipsoidGrofile.atoms.empty(), "Ellipsoid membrane did not contain any atoms");
		int outerEllipsoidLipids = 0;
		int innerEllipsoidLipids = 0;
		std::vector<Float3> outerEllipsoidCenters;
		std::vector<Float3> innerEllipsoidCenters;
		Float3 maximumCenterExtent{};
		for (const Float3& center : ResidueCenters(ellipsoidGrofile)) {
			const Float3 offset = center - ellipsoidCenter;
			const float implicitSurfaceValue =
				offset.x * offset.x / (ellipsoidRadii.x * ellipsoidRadii.x)
				+ offset.y * offset.y / (ellipsoidRadii.y * ellipsoidRadii.y)
				+ offset.z * offset.z / (ellipsoidRadii.z * ellipsoidRadii.z);
			if (implicitSurfaceValue > 1.f) {
				++outerEllipsoidLipids;
				outerEllipsoidCenters.push_back(center);
			}
			else {
				++innerEllipsoidLipids;
				innerEllipsoidCenters.push_back(center);
			}
			maximumCenterExtent = Float3::ElementwiseMax(maximumCenterExtent, offset.abs());
		}
		ASSERT(outerEllipsoidLipids > 2 && innerEllipsoidLipids > 2,
			"Could not identify both ellipsoid membrane leaflets");
		ASSERT(maximumCenterExtent.z > maximumCenterExtent.y
			&& maximumCenterExtent.y > maximumCenterExtent.x,
			"Ellipsoid membrane did not preserve its unequal axis radii");
		ASSERT(MinimumNearestNeighborSpacing(outerEllipsoidCenters) > 0.35f
			&& MinimumNearestNeighborSpacing(innerEllipsoidCenters) > 0.35f,
			"Ellipsoid membrane placement produced overlapping or clustered lipids");

		bool rejectedTightlyCurvedEllipsoid = false;
		try {
			GroFile invalidGrofile;
			invalidGrofile.box_size = Float3{ 20.f };
			TopologyFile invalidTopfile;
			invalidTopfile.SetSystem("Membrane");
			SimulationBuilder::CreateMembrane(invalidGrofile, invalidTopfile, lipids,
				MembraneGeometry::Ellipsoid{
					Float3{ 10.f }, Float3{ minimumRadius, minimumRadius, minimumRadius * 2.f } });
		}
		catch (const std::invalid_argument&) {
			rejectedTightlyCurvedEllipsoid = true;
		}
		ASSERT(rejectedTightlyCurvedEllipsoid,
			"An ellipsoid with an unreasonably small local curvature radius was accepted");

		GroFile planarGrofile;
		planarGrofile.box_size = Float3{ 8.f };
		TopologyFile planarTopfile;
		planarTopfile.SetSystem("Membrane");
		SimulationBuilder::CreateMembrane(planarGrofile, planarTopfile, lipids,
			MembraneGeometry::Plane{ 4.f });
		std::vector<float> topHeights;
		std::vector<float> bottomHeights;
		for (const Float3& center : ResidueCenters(planarGrofile))
			(center.z > 4.f ? topHeights : bottomHeights).push_back(center.z);
		ASSERT(topHeights.size() > 2 && bottomHeights.size() > 2,
			"Could not identify both planar membrane leaflets");
		ASSERT(Statistics::StdDev(topHeights) > 0.03f && Statistics::StdDev(bottomHeights) > 0.03f,
			"Planar membrane height roughness was too small");
		ASSERT(Statistics::StdDev(topHeights) < 0.30f && Statistics::StdDev(bottomHeights) < 0.30f,
			"Planar membrane height roughness was unreasonably large");

		const auto command = LiveEdit::ParseCommand(
			"buildmembrane -lipids DMPC 100 -sphere 8 8 8 4");
		const auto& buildCommand = std::get<LiveEdit::BuildMembrane>(command);
		ASSERT(buildCommand.geometry.has_value(), "Live-edit sphere geometry was not parsed");
		const auto& parsedSphere = std::get<MembraneGeometry::Sphere>(*buildCommand.geometry);
		ASSERT(parsedSphere.center == Float3{ 8.f } && parsedSphere.radius == 4.f,
			"Live-edit sphere geometry had incorrect values");
		const auto ellipsoidCommand = LiveEdit::ParseCommand(
			"buildmembrane -lipids DMPC 100 -ellipsoid 8 9 10 4 5 6");
		const auto& parsedBuildCommand = std::get<LiveEdit::BuildMembrane>(ellipsoidCommand);
		ASSERT(parsedBuildCommand.geometry.has_value(), "Live-edit ellipsoid geometry was not parsed");
		const auto& parsedEllipsoid = std::get<MembraneGeometry::Ellipsoid>(*parsedBuildCommand.geometry);
		ASSERT((parsedEllipsoid.center == Float3{ 8.f, 9.f, 10.f }
			&& parsedEllipsoid.radii == Float3{ 4.f, 5.f, 6.f }),
			"Live-edit ellipsoid geometry had incorrect values");

		return LimaUnittestResult{ true, "", envmode == Full };
	}

	// This test checks topology compatibility and physically bounded coordinate generation, NOT considering EM.
	static LimaUnittestResult TestBuildmembraneSmall(EnvMode envmode, bool do_em)
	{		
		const fs::path workDir = AutomatedTestsDir() / "BuildMembraneSmall";
		const fs::path mol_dir = workDir / "molecule";
		TestUtils::CleanDirIfNotContains(mol_dir, "reference");

		Lipids::Selection lipidselection;
		const std::array<std::string, 6> lipids = { "POPC", "POPE", "DDPC", "DMPC", "cholesterol", "DOPC" };
		for (const auto& lipidname : lipids) {
			lipidselection.emplace_back(Lipids::Select{ lipidname, workDir, lipidname == "POPC" ? 50. : 10.});	// 10% of each lipid, except 50% POPC
		}

		// Build the membrane, and write it to disk
		GroFile gro;
		gro.box_size = Float3{ 7.f };
		gro.title = "Membrane";
		TopologyFile top;
		top.SetSystem("Membrane");
		SimulationBuilder::CreateMembrane(gro, top, lipidselection, 3.5f);
		gro.printToFile(mol_dir / "membrane.gro");
		top.printToFile(mol_dir / "membrane.top");

		// Test the topology is identical to reference
		TopologyFile newTop{ mol_dir / "membrane.top" };
		TopologyFile refTop{ mol_dir / "membrane_reference.top" };

		//std::ostringstream oss;
		//const auto& newAtoms = newTop.GetAllElements<TopologyFile::AtomsEntry>();
		//const auto& refAtoms = refTop.GetAllElements<TopologyFile::AtomsEntry>();
		//if (auto [a, b] = std::ranges::mismatch(newAtoms, refAtoms);
		//	a != newAtoms.end() || b != refAtoms.end()) {
		//	a->composeString(oss); b->composeString(oss);
		//	std::string str = oss.str();
		//	printf(std::format("Mismatch at {}:\n{}\n ", std::distance(newAtoms.begin(), a), str).c_str());
		//}

		LimaUnittestResult topTestResults = TestUtils::CompareTopologyFiles(newTop, refTop, envmode);
		if (!topTestResults.success)
			return topTestResults;
		
		// Test the conf is identical to reference
		GroFile newGro{ mol_dir / "membrane.gro" };
		GroFile refGro{ mol_dir / "membrane_reference.gro" };
		ASSERT(newGro.box_size == refGro.box_size, "Box size mismatch");
		ASSERT(newGro.atoms.size() == refGro.atoms.size(), "Atom count mismatch");
		for (const auto& atom : newGro.atoms) {
			ASSERT(std::isfinite(atom.position.x) && std::isfinite(atom.position.y) && std::isfinite(atom.position.z),
				"Planar membrane contained a non-finite atom position");
		}
		std::vector<float> topLeafletHeights;
		std::vector<float> bottomLeafletHeights;
		for (const Float3& center : ResidueCenters(newGro))
			(center.z > 3.5f ? topLeafletHeights : bottomLeafletHeights).push_back(center.z);
		ASSERT(topLeafletHeights.size() > 2 && bottomLeafletHeights.size() > 2,
			"Could not identify both planar membrane leaflets");
		ASSERT(Statistics::StdDev(topLeafletHeights) > 0.03f && Statistics::StdDev(bottomLeafletHeights) > 0.03f,
			"Planar membrane height roughness was too small");
		ASSERT(Statistics::StdDev(topLeafletHeights) < 0.30f && Statistics::StdDev(bottomLeafletHeights) < 0.30f,
			"Planar membrane height roughness was unreasonably large");

		// Finally test if we can stabilize the simulation
		const float emtol = 200.f;
		auto sim = Programs::EnergyMinimize(gro, top, true, workDir, envmode, true, emtol);
		float finalMaxForce = sim->maxForceBuffer.back().second;

		return LimaUnittestResult{ finalMaxForce < emtol && finalMaxForce != 0, std::format("Failed to energy minimize membrane {:.2f}/{:.2f}", sim->maxForceBuffer.back().second, emtol), envmode == Full};
	}

	static LimaUnittestResult TestBuildmembraneWithCustomlipidAndCustomForcefield(EnvMode envmode) {
		const fs::path workDir = AutomatedTestsDir() / "BuildMembraneCustom";
		const fs::path mol_dir = workDir / "molecule";

		//TestUtils::CleanDirectory(mol_dir);
		fs::remove_all(mol_dir);
		fs::create_directory(mol_dir);

		Lipids::Selection lipidselection;
		const std::vector<std::pair<std::string, double>> lipids = { {"POPC", 70.}, {"CUST" , 30.} };
		for (const auto& [lipidname, percentage] : lipids) {
			lipidselection.emplace_back(Lipids::Select{ lipidname, workDir, percentage });	// 10% of each lipid, except 50% POPC
		}

		GroFile gro;
		gro.box_size = Float3{ 7.f };
		gro.title = "Membrane";
		TopologyFile top;
		top.SetSystem("Membrane");
		SimulationBuilder::CreateMembrane(gro, top, lipidselection, 3.5f);
		Programs::EnergyMinimize(gro, top, true, workDir, envmode, true, 300000.f); // high emtol, because we dont care about EM, we just want to see if the simulation can even start

		gro.printToFile(mol_dir / "membrane.gro");
		top.printToFile(mol_dir / "membrane.top");

		TopologyFile newTop{ mol_dir / "membrane.top" };
		GroFile newGro{ mol_dir / "membrane.gro" };

		//std::ostringstream oss;
		//const auto& newAtoms = newTop.GetAllElements<TopologyFile::AtomsEntry>();
		//const auto& refAtoms = top.GetAllElements<TopologyFile::AtomsEntry>();
		//if (auto [a, b] = std::ranges::mismatch(newAtoms, refAtoms);
		//	a != newAtoms.end() || b != refAtoms.end()) {
		//	a->composeString(oss); b->composeString(oss);
		//	std::string str = oss.str();
		//	printf(std::format("Mismatch at {}:\n{}\n ", std::distance(newAtoms.begin(), a), str).c_str());
		//}
		ASSERT(std::ranges::equal(newTop.GetAllElements<TopologyFile::AtomsEntry>(), top.GetAllElements<TopologyFile::AtomsEntry>()), "Topology Atom Mismatch");


		SimParams params{};
		params.em_variant = true;
		Environment env(workDir, envmode);
		env.CreateSimulation(newGro, newTop, params);

		return LimaUnittestResult{ true , "No error", envmode == Full };
	}

	LimaUnittestResult TestAllStockholmlipids(EnvMode envmode) {
		const fs::path workDir = AutomatedTestsDir() / "BuildMembraneSmall";

		const fs::path path = FileUtils::GetLimaDir() / "resources/Slipids";
		std::vector<std::string> targets;
		for (const auto& entry : fs::directory_iterator(path)) {
			if (entry.path().extension() == ".gro") {
				const std::string base_name = entry.path().stem().string();
				const fs::path itp_file = path / (base_name + ".itp");
				if (fs::exists(itp_file)) {
					targets.emplace_back(base_name);
				}
			}
		}

		Lipids::Selection lipidselection;
		for (const auto& lipidname : targets) {
			lipidselection.emplace_back(Lipids::Select{ lipidname, workDir, 100. / static_cast<double>(targets.size())});	// 10% of each lipid, except 50% POPC
		}

		// The first test is pretty much just to see if this function throws
		GroFile grofile;
		grofile.box_size = Float3{ 10.f };
		grofile.title = "Membrane";
		TopologyFile topfile;
		topfile.SetSystem("Membrane");
		SimulationBuilder::CreateMembrane(grofile, topfile, lipidselection, 5.f);
		/*for (const auto& molecule : topfile.GetSystem().molecules) {
			ASSERT(molecule.moleculetype->readFromCache, "This lipid top should have been read from a cached file");
		}*/

		// The third test is to see if this function throws
		const float emtol = 1000.f;
		auto sim = Programs::EnergyMinimize(grofile, topfile, false, workDir, envmode, true, emtol);

		ASSERT(sim->maxForceBuffer.back().second < emtol, "Failed to energy minimize membrane");

		return LimaUnittestResult{ true , "", envmode == Full };
	}

	LimaUnittestResult BuildAndRelaxVesicle(EnvMode envmode) {
		GroFile grofile;
		grofile.box_size = Float3{ 5.f };
		TopologyFile topfile;
		const fs::path workDir = TestUtils::HeavyTestsDir() / "etc";
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

			if (!TestUtils::CompareVecWithFile(vertices, TestUtils::HeavyTestsDir() / fs::path{ "etc" } / "buildvesicle.bin", 0.01, overwriteData))
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
			
			if (!TestUtils::CompareVecWithFile(vertices, TestUtils::HeavyTestsDir() / fs::path{ "etc" } / "relaxedvesicle.bin", 0.01, overwriteData))
				return LimaUnittestResult{ false , "After relaxation mismatch", envmode == Full };
		}

		if (overwriteData)
			return LimaUnittestResult{false, "Overwriting data", envmode == Full};



		return LimaUnittestResult{ true , "No error", envmode == Full };
	}
}
