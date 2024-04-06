#pragma once

#include "SimulationBuilder.h"

class Environment;

namespace Programs {
	

	MDFiles::FilePair CreateMembrane(
		Environment& env, LipidsSelection& lipidselection, bool carryout_em, float centerCoordinate, bool writeFiles);

	void SetMoleculeCenter(ParsedGroFile& grofile, Float3 targetCenter);

} // namespace Programs