#include "Bodies.cuh"


SingleBond::SingleBond(std::array<uint8_t, n_atoms> ids, float b0, float kb) : b0(b0), kb(kb) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
SingleBondFactory::SingleBondFactory(std::array<uint32_t, n_atoms> ids, float b0, float kb) : SingleBond{ {0,0}, b0, kb } {
	for (int i = 0; i < n_atoms; i++) {
		global_atom_indexes[i] = ids[i];
	}
}
AngleBond::AngleBond(std::array<uint8_t, n_atoms> ids, float theta_0, float k_theta) : theta_0(theta_0), k_theta(k_theta) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
AngleBondFactory::AngleBondFactory(std::array<uint32_t, n_atoms> ids, float theta_0, float k_theta) : AngleBond{ {0,0,0}, theta_0, k_theta } {
	for (int i = 0; i < n_atoms; i++) {
		global_atom_indexes[i] = ids[i];
	}
}

DihedralBond::DihedralBond(std::array<uint8_t, 4> ids, float phi_0, float k_phi, float n) : phi_0(phi_0), k_phi(k_phi), n(n) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
DihedralBondFactory::DihedralBondFactory(std::array<uint32_t, 4> ids, float phi_0, float k_phi, float n) : DihedralBond{ {0,0,0,0}, phi_0, k_phi, n } {
	for (int i = 0; i < n_atoms; i++) {
		global_atom_indexes[i] = ids[i];
	}
}

ImproperDihedralBond::ImproperDihedralBond(std::array<uint8_t, n_atoms> ids, float psi_0, float k_psi) : psi_0(psi_0), k_psi(k_psi) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}
ImproperDihedralBondFactory::ImproperDihedralBondFactory(std::array<uint32_t, n_atoms> ids, float psi_0, float k_psi) : ImproperDihedralBond{ {0,0,0,0}, psi_0, k_psi } {
	for (int i = 0; i < n_atoms; i++) {
		global_atom_indexes[i] = ids[i];
	}
}


