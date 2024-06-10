#pragma once

#include "SimulationBuilder.h"

class Environment;

namespace Programs {
	

	MDFiles::FilePair CreateMembrane(
		Environment& env, LipidsSelection& lipidselection, bool carryout_em, float centerCoordinate, bool writeFiles);

	void SetMoleculeCenter(GroFile& grofile, Float3 targetCenter);

	// Load file into a box, optionally solvate it and then run untill energy is at a stable level
	void EnergyMinimize(Environment& env, GroFile& grofile, const TopologyFile& topFile, bool solvate, float boxlenNM);

} // namespace Programs