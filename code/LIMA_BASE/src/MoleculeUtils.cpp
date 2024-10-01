#include "MoleculeUtils.h"
#include "BoundaryConditionPublic.h"

void MoleculeUtils::SetMoleculeCenter(GroFile& grofile, Float3 targetCenter) {

	// First make sure the molecule is not split due to PBC ;; TODO We are not currently checking that
	for (auto& particle : grofile.atoms) {
		//const Float3 center{ 0, 0, 0 };
		const Float3 center = grofile.atoms[0].position;
		//BoundaryConditionPublic::applyHyperposNM(center, particle.position, grofile.box_size.x, BoundaryConditionSelect::PBC);
	}

	Double3 sum = { 0,0,0 };
	for (auto& particle : grofile.atoms)
		sum += particle.position;

	const Double3 currentCenter = sum / static_cast<double>(grofile.atoms.size());
	const Float3 diff = targetCenter - Float3{ currentCenter.x, currentCenter.y, currentCenter.z };

	for (auto& particle : grofile.atoms) {
		particle.position += diff;
	}
}