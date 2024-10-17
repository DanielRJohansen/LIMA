#pragma once

#include "MDFiles.h"

namespace MoleculeUtils {

	Float3 GeometricCenter(const GroFile& grofile);

	float Radius(const GroFile& grofile, const Float3& center);

	void MakeMoleculeWholeAfterPBCFragmentation(GroFile& grofile, const TopologyFile& topfile);

	// Center molecule around targetCenter. Defaults to grofile.boxlen/2
	void CenterMolecule(GroFile& grofile, const TopologyFile& topfile, std::optional<Float3> targetCenter=std::nullopt);
}