#pragma once

#include "MDFiles.h"

namespace MoleculeUtils {

	Float3 GeometricCenter(const GroFile& grofile);

	float Radius(const GroFile& grofile, const Float3& center);

	void MakeMoleculeWholeAfterPBCFragmentation(GroFile& grofile, const TopologyFile::Moleculetype& topfile);
	void MakeMoleculeWholeAfterPBCFragmentation(GroFile& grofile, const TopologyFile& topfile);

	// Resize the box to contain all atoms with the requested padding on every side,
	// and translate the atoms so the structure is centered in the new box.
	void FitMoleculeInBox(GroFile& grofile, float padding = 1.f);

	// Center molecule around targetCenter. Defaults to grofile.boxlen/2
	void CenterMolecule(GroFile& grofile, const TopologyFile::Moleculetype& topfile, std::optional<Float3> targetCenter=std::nullopt);


	/// <summary>
	/// Rotates the molecule around the geometric center. Assumes that it is whole.
	/// Rotates around z, y then x axis.
	/// </summary>
	/// <param name="grofile"></param>
	/// <param name="rotation">[rad]</param>
	void RotateMolecule(GroFile& grofile, Float3 rotation);
}
