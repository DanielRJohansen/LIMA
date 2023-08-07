#include "Bodies.cuh"

// Create the solventblockgrid on host, and fill out origo for all blocks for all steps
void SolventBlockHelpers::createSolventblockGrid(SolventBlockGrid** solventblockgrid_circularqueue_device) {	// TODO: rename var to what is it now

	*solventblockgrid_circularqueue_device = new SolventBlockGrid[STEPS_PER_SOLVENTBLOCKTRANSFER];

	for (int i = 0; i < STEPS_PER_SOLVENTBLOCKTRANSFER; i++) {
		for (int z = 0; z < BOXGRID_N_NODES; z++) {
			for (int y = 0; y < BOXGRID_N_NODES; y++) {
				for (int x = 0; x < BOXGRID_N_NODES; x++) {
					NodeIndex origo{ x, y, z };	// Doubles as the 3D index of the block!
					//gridqueue_host[i].getBlockPtr(origo)->origo = origo;
					(*solventblockgrid_circularqueue_device)[i].getBlockPtr(origo)->origo = origo;
				}
			}
		}
	}
}





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








bool NeighborList::addId(uint16_t new_id, NEIGHBOR_TYPE nt) {
	switch (nt)
	{
	case NeighborList::COMPOUND:
		if (n_compound_neighbors < NEIGHBORLIST_MAX_COMPOUNDS) {
			neighborcompound_ids[n_compound_neighbors++] = new_id;
			return true;
		}
		printf("\nFailed to insert compound neighbor id %d!\n", new_id);
		exit(1);
		break;
	default:
		break;
	}
	exit(1);
	return false;
}
bool NeighborList::removeId(uint16_t neighbor_id, NEIGHBOR_TYPE nt) {

	switch (nt) {
	case NeighborList::COMPOUND:
		for (int i = 0; i < n_compound_neighbors; i++) {
			if (neighborcompound_ids[i] == neighbor_id) {
				neighborcompound_ids[i] = neighborcompound_ids[n_compound_neighbors - 1];
				n_compound_neighbors--;
				return true;
			}
		}
		printf("\nFailed to remove neighbor compound %d of %d\n", neighbor_id, n_compound_neighbors);
		break;
	default:
		printf("Faulty neighbortype\n");
		break;
	}

	return false;
}


__host__ void NeighborList::addGridnode(uint16_t gridnode_id) {
	if (n_gridnodes >= max_gridnodes) { 
		throw std::exception("No room for more nearby gridnodes"); }
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


__host__ void CompoundGridNode::addNearbyCompound(int16_t compound_id)
{
	if (n_nearby_compounds >= max_nearby_compounds) {
		throw std::exception("Failed to add compound to CompoundGridNode\n");
	}
	nearby_compound_ids[n_nearby_compounds++] = compound_id;
}