#pragma once

#include "SimulationBuilder.h"

#include <string_view>

class MoleculeHullCollection;

namespace Programs {
	enum class WaterModel { Tip3p, Tip4p, Tips3p, Tip5p, Spc, Spce };

	struct GmxConversionResult {
		GroFile grofile;
		TopologyFile topology;
		// Ordered like topology.GetSystem().molecules. Each file contains the
		// position restraints for the corresponding molecule.
		std::vector<GenericItpFile> positionRestraints;
	};

	WaterModel ParseWaterModel(std::string_view name);

	void GetForcefieldParams(const GroFile&, const TopologyFile&, const fs::path& workdir);

	MoleculeHullCollection MakeLipidVesicle(GroFile&, TopologyFile&, Lipids::Selection, float vesicleRadius, 
		Float3 vesicleCenter, std::optional<int> numLipids=std::nullopt);

	void MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize, bool renderProgress);

	void StaticbodyEnergyMinimize(GroFile&, const TopologyFile&, bool render);

	/// Build in-memory CHARMM27 coordinates, topology, and position restraints
	/// from a protein PDB or mmCIF structure.
	GmxConversionResult ToGmx(const fs::path& structureFile, WaterModel waterModel = WaterModel::Tip3p);
}
