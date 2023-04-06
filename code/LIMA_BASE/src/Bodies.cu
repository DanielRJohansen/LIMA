#include "Bodies.cuh"




void CompoundCoords::copyInitialCoordConfiguration(CompoundCoords* coords, CompoundCoords* coords_prev, CompoundCoords* coordarray_circular_queue) {
	// Move pos_t
	cudaMemcpy(coordarray_circular_queue, coords, sizeof(CompoundCoords) * MAX_COMPOUNDS, cudaMemcpyHostToDevice);

	// Move pos_t - 1
	const int index0_of_prev = (STEPS_PER_LOGTRANSFER - 1) * MAX_COMPOUNDS;
	cudaMemcpy(&coordarray_circular_queue[index0_of_prev], coords_prev, sizeof(CompoundCoords) * MAX_COMPOUNDS, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Error during coord's initial configuration copyToDevice\n");
		exit(1);
	}
}


void SolventCoord::copyInitialCoordConfiguration(SolventCoord* coords, SolventCoord* coords_prev, SolventCoord* solventcoordarray_circular_queue) {
	// Move pos_t
	cudaMemcpy(solventcoordarray_circular_queue, coords, sizeof(SolventCoord) * MAX_SOLVENTS, cudaMemcpyHostToDevice);

	// Move pos_t - 1
	const int index0_of_prev = (STEPS_PER_LOGTRANSFER - 1) * MAX_SOLVENTS;
	cudaMemcpy(&solventcoordarray_circular_queue[index0_of_prev], coords_prev, sizeof(SolventCoord) * MAX_SOLVENTS, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Error during solventcoord's initial configuration copyToDevice\n");
		exit(1);
	}
}



bool SolventBlockHelpers::copyInitialConfiguration(const SolventBlockGrid& grid, const SolventBlockGrid& grid_prev, SolventBlockGrid* grid_circular_queue) {
	// Move pos_t
	cudaMemcpy(grid_circular_queue, &grid, sizeof(SolventBlockGrid), cudaMemcpyHostToDevice);

	// Move pos_t - 1
	const int index0_of_prev = (STEPS_PER_SOLVENTBLOCKTRANSFER - 1);
	cudaMemcpy(&grid_circular_queue[index0_of_prev], &grid_prev, sizeof(SolventBlockGrid), cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Error during solventcoord's initial configuration copyToDevice\n");
		exit(1);
	}
}


void SolventBlockHelpers::createSolventblockGrid(SolventBlockGrid** solventblockgrid_circularqueue_device) {
	const uint64_t n_bytes_solventblockgrids = sizeof(SolventBlockGrid) * STEPS_PER_SOLVENTBLOCKTRANSFER;
	cudaMalloc(solventblockgrid_circularqueue_device, n_bytes_solventblockgrids);

	auto gridqueue_host = new SolventBlockGrid[STEPS_PER_SOLVENTBLOCKTRANSFER];
	for (int i = 0; i < STEPS_PER_SOLVENTBLOCKTRANSFER; i++) {
		for (int z = 0; z < BOXGRID_N_NODES; z++) {
			for (int y = 0; y < BOXGRID_N_NODES; y++) {
				for (int x = 0; x < BOXGRID_N_NODES; x++) {
					NodeIndex origo{ x, y, z };	// Doubles as the 3D index of the block!
					gridqueue_host[i].getBlockPtr(origo)->origo = origo;
				}
			}
		}
	}
	cudaMemcpy(*solventblockgrid_circularqueue_device, gridqueue_host, sizeof(SolventBlockGrid) * STEPS_PER_SOLVENTBLOCKTRANSFER, cudaMemcpyHostToDevice);
	delete[] gridqueue_host;
}

__host__  void SolventBlockHelpers::createSolventblockTransfermodules(SolventBlockTransfermodule** transfermodule_array) {
	const size_t array_bytesize = sizeof(SolventBlockTransfermodule) * SolventBlockGrid::blocks_total;

	// Allocate the memory on device
	cudaMalloc(transfermodule_array, array_bytesize);

	// Initialize the memory on host, then transfer to dest on device
	auto array_host = new SolventBlockTransfermodule[SolventBlockGrid::blocks_total];
	for (int i = 0; i < SolventBlockGrid::blocks_total; i++) {
		array_host[i].n_remain = 0;
	}
	cudaMemcpy(*transfermodule_array, array_host, array_bytesize, cudaMemcpyHostToDevice);

	delete[] array_host;
}

void SolventBlockHelpers::setupBlockMetaOnHost(SolventBlockGrid* grid, SolventBlockGrid* grid_prev) {
	for (int z = 0; z < BOXGRID_N_NODES; z++) {
		for (int y = 0; y < BOXGRID_N_NODES; y++) {
			for (int x = 0; x < BOXGRID_N_NODES; x++) {
				NodeIndex origo{ x, y, z };	// Doubles as the 3D index of the block!
				grid->getBlockPtr(origo)->origo = origo;
				grid_prev->getBlockPtr(origo)->origo = origo;
			}
		}
	}
}


SingleBond::SingleBond(std::array<int, n_atoms> ids) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

SingleBond::SingleBond(int id1, int id2, float b0, float kb) : b0(b0), kb(kb) {
	// This is only for loading the forcefield, so the ID's refers to id's given in .conf file!
	atom_indexes[0] = id1;
	atom_indexes[1] = id2;
}

AngleBond::AngleBond(std::array<int, n_atoms> ids) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

AngleBond::AngleBond(int id1, int id2, int id3, float theta_0, float k_theta) : theta_0(theta_0), k_theta(k_theta) {
	atom_indexes[0] = id1;
	atom_indexes[1] = id2;
	atom_indexes[2] = id3;
}

DihedralBond::DihedralBond(std::array<int, n_atoms> ids) {
	for (int i = 0; i < n_atoms; i++) {
		atom_indexes[i] = ids[i];
	}
}

DihedralBond::DihedralBond(int id1, int id2, int id3, int id4, float phi_0, float k_phi, int n) : phi_0(phi_0), k_phi(k_phi), n(n) {
	atom_indexes[0] = id1;
	atom_indexes[1] = id2;
	atom_indexes[2] = id3;
	atom_indexes[3] = id4;
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

	case NeighborList::SOLVENT:
		if (n_solvent_neighbors < NEIGHBORLIST_MAX_SOLVENTS) {
			neighborsolvent_ids[n_solvent_neighbors++] = new_id;
			//if (associated_id==0 || new_id == 0)
				//printf("%d added id %d\n", associated_id, new_id);
			return true;
		}

		printf("\nFailed to insert solvent neighbor id %d of %d\n", new_id, n_solvent_neighbors);
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
	case NeighborList::SOLVENT:
		for (int i = 0; i < n_solvent_neighbors; i++) {
			if (neighborsolvent_ids[i] == neighbor_id) {
				neighborsolvent_ids[i] = neighborsolvent_ids[n_solvent_neighbors - 1];
				n_solvent_neighbors--;
				return true;
			}
		}
		printf("\n%d Failed to remove neighbor solvent ID: %d of %d total IDs. max+1: %d\n", associated_id, neighbor_id, n_solvent_neighbors, neighborcompound_ids[n_solvent_neighbors]);
		for (int i = 0; i < n_solvent_neighbors; i++) {
			//printf("%d\n", neighborsolvent_ids[i]);
		}
		break;
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
	if (n_gridnodes == max_gridnodes) { throw ("No room for more nearby gridnodes"); }
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
	if (n_nearby_compounds == max_elements) {
		throw std::exception("Failed to add compound to CompoundGridNode\n");
	}
	nearby_compound_ids[n_nearby_compounds++] = compound_id;
}

__host__ void CompoundGridNode::addAssociatedCompound(int16_t compound_id) {
	if (n_associated_compounds == 4) {
		throw std::exception("Failed to add compound to CompoundGridNode\n");
	}
	associated_ids[n_associated_compounds++] = compound_id;
}