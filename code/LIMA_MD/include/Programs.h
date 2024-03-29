#pragma once

#include "SimulationBuilder.h"

class Environment;

namespace Programs {
	

	void CreateMembrane(Environment& env, LipidsSelection& lipidselection, bool carryout_em, float centerCoordinate);

} // namespace Programs