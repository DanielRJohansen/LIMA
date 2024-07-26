#include "Bodies.cuh"


SingleBond::SingleBond(std::array<uint8_t, n_atoms> ids, float b0, float kb) : params{ b0, kb } {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

AngleBond::AngleBond(std::array<uint8_t, n_atoms> ids, float theta_0, float k_theta) : params{ theta_0, k_theta } {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}


DihedralBond::DihedralBond(std::array<uint8_t, 4> ids, float phi_0, float k_phi, float n) 
	: params{static_cast<half>(phi_0), static_cast<half>(k_phi), static_cast<half>(n)} {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

ImproperDihedralBond::ImproperDihedralBond(std::array<uint8_t, n_atoms> ids, float psi_0, float k_psi) : params{ psi_0, k_psi } {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
