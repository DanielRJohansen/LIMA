#include "Bodies.cuh"



SingleBond::SingleBond(std::array<int, n_atoms> ids) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

SingleBond::SingleBond(std::array<uint32_t, n_atoms> ids, float b0, float kb) : b0(b0), kb(kb) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}


AngleBond::AngleBond(std::array<int, n_atoms> ids) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

AngleBond::AngleBond(std::array<uint32_t, n_atoms> ids, float theta_0, float k_theta) : theta_0(theta_0), k_theta(k_theta) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}


DihedralBond::DihedralBond(std::array<int, n_atoms> ids) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

DihedralBond::DihedralBond(std::array<uint32_t, 4> ids, float phi_0, float k_phi, float n) : phi_0(phi_0), k_phi(k_phi), n(n) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}


ImproperDihedralBond::ImproperDihedralBond(std::array<uint32_t, n_atoms> ids) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i]; 
	}
}

ImproperDihedralBond::ImproperDihedralBond(std::array<uint32_t, n_atoms> ids, float psi_0, float k_psi) : psi_0(psi_0), k_psi(k_psi) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}






void NeighborList::addCompound(uint16_t new_id) {
	if (n_compound_neighbors >= NEIGHBORLIST_MAX_COMPOUNDS) {
		printf("\nFailed to insert compound neighbor id %d!\n", new_id);
		throw std::runtime_error("Neighborlist overflow");
	}

	neighborcompound_ids[n_compound_neighbors++] = new_id;
}

void NeighborList::removeCompound(uint16_t neighbor_id) {

	for (int i = 0; i < n_compound_neighbors; i++) {
		if (neighborcompound_ids[i] == neighbor_id) {
			neighborcompound_ids[i] = neighborcompound_ids[n_compound_neighbors - 1];	// Overwrite this positions with compound at the back of the line
			n_compound_neighbors--;
			return;
		}
	}

	printf("\nFailed to locate neighbor compound id: %d of %d compound_ids in current nlist\n", neighbor_id, n_compound_neighbors);
	throw std::runtime_error("Nlist failed to remove neighbor");
}

#ifdef ENABLE_SOLVENTS
__host__ void NeighborList::addGridnode(uint16_t gridnode_id) {
	if (n_gridnodes >= max_gridnodes) { 
		throw std::runtime_error("No room for more nearby gridnodes"); }
	gridnode_ids[n_gridnodes++] = gridnode_id;
}


__host__ void NeighborList::removeGridnode(uint16_t gridnode_id) {
	for (int i = 0; i < n_gridnodes; i++) {
		if (gridnode_ids[i] == gridnode_id) {
			gridnode_ids[i] = gridnode_ids[n_gridnodes - 1];
			n_gridnodes--;
			return;
		}
	}
	throw("Failed to remove gridnode from nlist");
}
#endif

__host__ void CompoundGridNode::addNearbyCompound(int16_t compound_id)
{
	if (n_nearby_compounds >= max_nearby_compounds) {
		throw std::runtime_error("Failed to add compound to CompoundGridNode\n");
	}
	nearby_compound_ids[n_nearby_compounds++] = compound_id;
}
