#pragma once

#include "SimulationBuilder.h"

class Environment;
class MoleculeHullCollection;

namespace Programs {
	void SetMoleculeCenter(GroFile& grofile, Float3 targetCenter);

	
	void GetForcefieldParams(const GroFile&, const TopologyFile&, const fs::path& workdir);

	MoleculeHullCollection MakeLipidVesicle(GroFile&, TopologyFile&, Lipids::Selection, float vesicleRadius, 
		Float3 vesicleCenter, std::optional<int> numLipids=std::nullopt);

	void MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize);

	// Load file into a box, optionally solvate it and then run untill energy is at a stable level
	void EnergyMinimize(Environment& env, GroFile& grofile, const TopologyFile& topFile, bool solvate, float boxlenNM);

	// Same as above, but starts with extreme EM and iteratively increases dt
	void EnergyMinimizeMax(GroFile& grofile, const TopologyFile& topfile, const fs::path& workDir, EnvMode envmode);
} // namespace Programs