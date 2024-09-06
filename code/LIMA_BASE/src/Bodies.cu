#include "Bodies.cuh"


SingleBond::SingleBond(std::array<uint8_t, nAtoms> ids, float b0, float kb) : params{ b0, kb } {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

AngleUreyBradleyBond::AngleUreyBradleyBond(std::array<uint8_t, nAtoms> ids, float t0, float kT, float ub0, float kUB) : params{ t0, kT, ub0, kUB} {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}


DihedralBond::DihedralBond(std::array<uint8_t, 4> ids, float phi_0, float k_phi, float n) 
	//: params{static_cast<half>(phi_0), static_cast<half>(k_phi), static_cast<half>(n)} {
	: params{ phi_0, k_phi, n } {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

ImproperDihedralBond::ImproperDihedralBond(std::array<uint8_t, nAtoms> ids, float psi_0, float k_psi) : params{ psi_0, k_psi } {
	for (int i = 0; i < nAtoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
