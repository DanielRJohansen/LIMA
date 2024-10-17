#pragma once

#include "MDFiles.h"

namespace MoleculeUtils {

	Float3 GeometricCenter(const GroFile& grofile);

	float Radius(const GroFile& grofile, const Float3& center);

	void MakeMoleculeWholeAfterPBCFragmentation(GroFile& grofile, const TopologyFile& topfile);

	// Center molecule around targetCenter. Defaults to grofile.boxlen/2
	void CenterMolecule(GroFile& grofile, const TopologyFile& topfile, std::optional<Float3> targetCenter=std::nullopt);


	/// <summary>
	/// Rotates the molecule around the geometric center. Assumes that it is whole.
	/// Rotates around z, y then x axis.
	/// </summary>
	/// <param name="grofile"></param>
	/// <param name="rotation">[rad]</param>
	void RotateMolecule(GroFile& grofile, Float3 rotation);
}