#pragma once

#include "SimulationBuilder.h"

class Environment;
class MoleculeHullCollection;
class Simulation;

namespace Programs {
	void SetMoleculeCenter(GroFile& grofile, Float3 targetCenter);

	
	void GetForcefieldParams(const GroFile&, const TopologyFile&, const fs::path& workdir);

	MoleculeHullCollection MakeLipidVesicle(GroFile&, TopologyFile&, Lipids::Selection, float vesicleRadius, 
		Float3 vesicleCenter, std::optional<int> numLipids=std::nullopt);

	void MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize);

	// Load file into a box, optionally solvate it and then run untill energy is at a stable level
	void EnergyMinimize(Environment& env, GroFile& grofile, const TopologyFile& topFile, bool solvate, float boxlenNM);

	
	/// <summary>
	/// Same as above, but starts with extreme EM and iteratively increases dt 
	/// </summary>
	/// <param name="grofile"></param>
	/// <param name="topfile"></param>
	/// <param name="workDir"></param>
	/// <param name="envmode"></param>
	/// <returns></returns>
	/// 
	

	/// <summary>
	/// Same as above, but starts with extreme EM and iteratively increases dt 
	/// </summary>
	/// <param name="writePositionsToGrofile">If false, the grofile will not be modified</param>
	/// <returns>Can be discarded if not needed. Only makes sense to discard if overwriting grofile</returns>
	std::unique_ptr<Simulation> EnergyMinimizeMax(GroFile&, const TopologyFile&, 
		bool writePositionsToGrofile, const fs::path& workDir, EnvMode);
} // namespace Programs