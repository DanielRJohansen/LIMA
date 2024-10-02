#pragma once

#include "MDFiles.h"

namespace MoleculeUtils {

	Float3 GeometricCenter(const GroFile& grofile);

	void SetMoleculeCenter(GroFile& grofile, Float3 targetCenter);
}