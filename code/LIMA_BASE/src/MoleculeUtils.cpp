#include "MoleculeUtils.h"
#include "BoundaryConditionPublic.h"



Float3 MoleculeUtils::GeometricCenter(const GroFile& grofile) {
	Float3 bbMin{ FLT_MAX }, bbMax{ -FLT_MAX };

	for (const auto& atom : grofile.atoms) {
		bbMin = Float3::ElementwiseMin(bbMin, atom.position);
		bbMax = Float3::ElementwiseMax(bbMax, atom.position);
	}

	return (bbMin + bbMax) / 2;
}

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

	//TODO This is WRONG, we should not use CoM, but geometric center
	const Double3 currentCenter = sum / static_cast<double>(grofile.atoms.size());
	const Float3 diff = targetCenter - Float3{ currentCenter.x, currentCenter.y, currentCenter.z };

	for (auto& particle : grofile.atoms) {
		particle.position += diff;
	}
}