#include "Bodies.cuh"


SingleBond::SingleBond(std::array<uint8_t, nAtoms> ids, const Parameters& params) : params{ params } {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
SingleBond::Parameters SingleBond::Parameters::CreateFromCharmm(float b0, float kB) {
	return Parameters{ b0, kB * KILO };
}

PairBond::PairBond(std::array<uint8_t, nAtoms> ids, const Parameters& params) : params{ params } {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
PairBond::Parameters PairBond::Parameters::CreateFromCharmm(float sigma, float epsilon) {
	return Parameters{ sigma, epsilon * KILO };
}

AngleUreyBradleyBond::AngleUreyBradleyBond(std::array<uint8_t, nAtoms> ids, const Parameters& params) : params{ params} {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
AngleUreyBradleyBond::Parameters AngleUreyBradleyBond::Parameters::CreateFromCharmm(float theta0, float kTheta, float r0, float kUB) {
	return Parameters{ theta0 * DEG_TO_RAD, kTheta * KILO, r0, kUB * KILO };
}

DihedralBond::DihedralBond(std::array<uint8_t, 4> ids, const Parameters& params) : params{ params } {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
DihedralBond::Parameters DihedralBond::Parameters::CreateFromCharmm(float phi0, float kPhi, int n) {
	return Parameters{ phi0 * DEG_TO_RAD, kPhi * KILO / 2.f, static_cast<float>(n)}; //TODO: Move the /2 from here to the force calculation to save an op there
}

ImproperDihedralBond::ImproperDihedralBond(std::array<uint8_t, nAtoms> ids, const Parameters& params) : params{ params} {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
ImproperDihedralBond::Parameters ImproperDihedralBond::Parameters::CreateFromCharmm(float phi0, float kPhi) {
	return Parameters{ phi0 * DEG_TO_RAD, kPhi * KILO }; //TODO: Move the /2 from here to the force calculation to save an op there? Still need??
}

//__host__ Float3 CompoundInterimState::sumForce(int particleIndex) const {
//	return Float3{}
//		//forceEnergyFarneighborShortrange[particleIndex].force 
//		+ forceEnergyImmediateneighborShortrange[particleIndex].force 
//		+ forceEnergyBonds[particleIndex].force 
//		+ forceEnergyBridge[particleIndex].force;
//}
//__host__ float CompoundInterimState::sumPotentialenergy(int particleIndex) const {
//	return 0.f
//		//forceEnergyFarneighborShortrange[particleIndex].potE
//		+ forceEnergyImmediateneighborShortrange[particleIndex].potE
//		+ forceEnergyBonds[particleIndex].potE
//		+ forceEnergyBridge[particleIndex].potE;
//}