#pragma once

#include "SimulationBuilder.h"

class Environment;
class MoleculeHullCollection;

namespace Programs {
	

	MDFiles::FilePair CreateMembrane(Environment&, LipidsSelection&, bool doEM, 
		float centerCoordinate, bool writeFiles);

	void SetMoleculeCenter(GroFile& grofile, Float3 targetCenter);

	// Load file into a box, optionally solvate it and then run untill energy is at a stable level
	void EnergyMinimize(Environment& env, GroFile& grofile, const TopologyFile& topFile, bool solvate, float boxlenNM);

	void GetForcefieldParams(const GroFile&, const TopologyFile&, const fs::path& workdir);

	MoleculeHullCollection MakeLipidVesicle(GroFile& grofile, TopologyFile& topfile, LipidsSelection);

	void MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize);
} // namespace Programs