#include "Bodies.cuh"


SingleBond::SingleBond(std::array<uint8_t, nAtoms> ids, const Parameters& params) : params{ params } {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

AngleUreyBradleyBond::AngleUreyBradleyBond(std::array<uint8_t, nAtoms> ids, const Parameters& params) : params{ params} {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

DihedralBond::DihedralBond(std::array<uint8_t, 4> ids, const Parameters& params) : params{ params } {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

ImproperDihedralBond::ImproperDihedralBond(std::array<uint8_t, nAtoms> ids, const Parameters& params) : params{ params} {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
__host__ Float3 CompoundInterimState::sumForce(int particleIndex) const {
	return Float3{}
		//forceEnergyFarneighborShortrange[particleIndex].force 
		+ forceEnergyImmediateneighborShortrange[particleIndex].force 
		+ forceEnergyBonds[particleIndex].force 
		+ forceEnergyBridge[particleIndex].force;
}
__host__ float CompoundInterimState::sumPotentialenergy(int particleIndex) const {
	return 0.f
		//forceEnergyFarneighborShortrange[particleIndex].potE
		+ forceEnergyImmediateneighborShortrange[particleIndex].potE
		+ forceEnergyBonds[particleIndex].potE
		+ forceEnergyBridge[particleIndex].potE;
}