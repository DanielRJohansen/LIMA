#pragma once

#include "SimulationBuilder.h"

class Environment;
class MoleculeHullCollection;
class Simulation;

namespace Programs {
	void GetForcefieldParams(const GroFile&, const TopologyFile&, const fs::path& workdir);

	MoleculeHullCollection MakeLipidVesicle(GroFile&, TopologyFile&, Lipids::Selection, float vesicleRadius, 
		Float3 vesicleCenter, std::optional<int> numLipids=std::nullopt);

	void MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize);

	/// <summary></summary>
	/// <param name="writePositionsToGrofile">If false, the grofile will not be modified</param>
	/// <param name="mayOverlapEdges">If the box contents may spill over the edge, set this to true.
	/// Then we will first run a pre-EM with boxEdgePotential enabled</param>
	/// <returns>Can be discarded if not needed. Only makes sense to discard if overwriting grofile</returns>
	std::unique_ptr<Simulation> EnergyMinimize(GroFile&, const TopologyFile&,
		bool writePositionsToGrofile, const fs::path& workDir, EnvMode, bool mayOverlapEdges, float emtol=100.f);
}