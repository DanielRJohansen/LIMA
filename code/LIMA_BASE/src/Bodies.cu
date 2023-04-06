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

//Float3 CompoundCoords::getAbsolutePositionLM(const int particle_index) {
//	return (origo * NANO_TO_LIMA).toFloat3() + rel_positions[particle_index].toFloat3();
//}


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



//SolventCoord SolventCoord::createFromPositionNM(const Float3& solvent_pos) {
//	SolventCoord coord;
//	const Float3 origo_f = solvent_pos.piecewiseRound();
//	const Float3 relpos_f = (solvent_pos - origo_f) * NANO_TO_LIMA;
//	coord.origo = Coord{ origo_f };
//	coord.rel_position = Coord{ relpos_f };
//	return coord;
//}



//bool SolventBlockHelpers::insertSolventcoordInGrid(SolventBlockGrid& grid, const SolventCoord& coord, uint32_t solvent_id)
//{
//	if (coord.origo.x >= 7 || coord.origo.y >= 7 || coord.origo.z >= 7) {
//		printf("");
//	}
//	return grid.getBlockPtr(coord.origo)->addSolvent(coord.rel_position, solvent_id);
//}





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


Float3 CompoundCarrier::calcCOM() {
	Float3 com;
	for (int i = 0; i < n_particles; i++) {
		com += state.positions[i] / static_cast<float>(n_particles);
	}
	return com;
}

__host__ void CompoundCarrier::addParticle(int atomtype_id, Float3 pos) {	// TODO: Delete?
	if (n_particles == MAX_COMPOUND_PARTICLES) {
		printf("ERROR: Cannot add particle to compound!\n");
		exit(1);
	}

	atom_types[n_particles] = atomtype_id;
	//prev_positions[n_particles] = pos * LIMA_SCALE;
	state.positions[n_particles] = pos;
	state_tsub1.positions[n_particles] = state.positions[n_particles];

	n_particles++;
}

__host__ void CompoundCarrier::addParticle(int atomtype_id, Float3 pos, int atomtype_color_id, int global_id) {
	if (n_particles == MAX_COMPOUND_PARTICLES) {
		printf("ERROR: Cannot add particle to compound!\n");
		exit(1);
	}

	state.positions[n_particles] = pos;
	state_tsub1.positions[n_particles] = state.positions[n_particles];

	atom_types[n_particles] = atomtype_id;
	atom_color_types[n_particles] = atomtype_color_id;

#ifdef LIMA_DEBUGMODE
	particle_global_ids[n_particles] = global_id;
	//printf("%d global id\n", global_id);
#endif
	n_particles++;
}

__host__ void CompoundCarrier::calcParticleSphere() {
	Float3 com = calcCOM();

	float furthest = 0.f;
	float closest = FLT_MAX;
	int closest_index = 0;

	for (int i = 0; i < n_particles; i++) {
		float dist = (state.positions[i] - com).len();
		closest_index = dist < closest ? i : closest_index;
		closest = std::min(closest, dist);
		furthest = std::max(furthest, dist);
	}

	key_particle_index = closest_index;
	confining_particle_sphere = furthest;
}



CompoundBridgeBundleCompact::CompoundBridgeBundleCompact(CompoundBridgeBundle* bundle, bool verbose) {
	n_bridges = bundle->n_bridges;
	//printf("Transferring %d %d bridges\n", n_bridges, bundle->n_bridges);
	for (int i = 0; i < n_bridges; i++) {
		compound_bridges[i] = CompoundBridgeCompact(bundle->compound_bridges[i], verbose);
		//printf("bridge %d has %d particles\n\n", i, compound_bridges[i].n_particles);
	}
}

CompoundCollection::CompoundCollection() {
	//compounds = new Compound[MAX_COMPOUNDS];
	bonded_particles_lut_manager = new BondedParticlesLUTManager(0);
	//compound_bridge_bundle = new CompoundBridgeBundleCompact;
}

Float3 CompoundCollection::calcCOM() {
	Float3 com(0.f);
	/*for (int i = 0; i < n_compounds; i++) {
		com += (compounds[i].calcCOM() * (1.f / (float)n_compounds));*/
	for (auto& compound : compounds) {
		com += compound.calcCOM() / static_cast<float>(n_compounds);
	}
	return com;
}












bool CompoundBridge::bondBelongsInBridge(GenericBond* bond) const {
	return (compound_id_left == bond->compound_ids[0] && compound_id_right == bond->compound_ids[1]);
}
bool CompoundBridge::particleAlreadyStored(ParticleRef* p_ref) {
	for (int i = 0; i < n_particles; i++) {
		if (particle_refs[i] == *p_ref) {
			return true;
		}
	}
	return false;
}
void CompoundBridge::addParticle(ParticleRef* particle_ref, CompoundCollection* molecule) {
	if (n_particles == MAX_PARTICLES_IN_BRIDGE) {
		printf("Too many particles in bridge\n");
		exit(0);
	}
	particle_ref->bridge_id = 0;
	particle_ref->local_id_bridge = n_particles;
	atom_types[n_particles] = molecule->compounds[particle_ref->compound_id].atom_types[particle_ref->local_id_compound];
	particle_refs[n_particles] = *particle_ref;
	n_particles++;
	//printf("Adding particle with global id: %d\n", particle_ref->global_id);
}

void CompoundBridge::addBondParticles(GenericBond* bond, CompoundCollection* molecule) {
	for (int p = 0; p < bond->n_particles; p++) {
		if (!particleAlreadyStored(&bond->particles[p])) {
			addParticle(&bond->particles[p], molecule);
		}
	}
}

void CompoundBridge::addGenericBond(SingleBond pb) {
	if (n_singlebonds == MAX_SINGLEBONDS_IN_BRIDGE) {
		printf("Cannot add bond to bridge\n");
		exit(0);
	}
	localizeIDs(&pb, 2);
	singlebonds[n_singlebonds++] = pb;
}
void CompoundBridge::addGenericBond(AngleBond ab) {
	if (n_anglebonds == MAX_ANGLEBONDS_IN_BRIDGE) {
		printf("Cannot add angle to bridge\n");
		exit(0);
	}
	localizeIDs(&ab, 3);
	anglebonds[n_anglebonds++] = ab;
}
void CompoundBridge::addGenericBond(DihedralBond db) {
	if (n_dihedrals == MAX_DIHEDRALBONDS_IN_BRIDGE) {
		printf("Cannot add dihedral to bridge\n");
		exit(0);
	}
	localizeIDs(&db, 4);
	dihedrals[n_dihedrals++] = db;
}









bool CompoundBridgeBundle::addBridge(uint16_t left_c_id, uint16_t right_c_id) {
	if (left_c_id > right_c_id) { std::swap(left_c_id, right_c_id); }

	if (n_bridges == COMPOUNDBRIDGES_IN_BUNDLE) {
		printf("FATAL ERROR: MAXIMUM bridges in bundle reached. Missing implementation of multiple bundles");
		exit(1);
	}

	compound_bridges[n_bridges++] = CompoundBridge(left_c_id, right_c_id);
	return true;
}

CompoundBridge* CompoundBridgeBundle::getBelongingBridge(GenericBond* bond) {
	for (int i = 0; i < n_bridges; i++) {
		CompoundBridge* bridge = &compound_bridges[i];
		if (bridge->bondBelongsInBridge(bond)) {
			return bridge;
		}
	}
	printf("FATAL ERROR: Failed to find belonging bridge");
	exit(1);
}






CompoundBridgeCompact::CompoundBridgeCompact(const CompoundBridge& bridge, bool verbose) :
	compound_id_left{ bridge.compound_id_left },
	compound_id_right{ bridge.compound_id_right }
{
	n_particles = bridge.n_particles;

	for (int i = 0; i < n_particles; i++) {
		particle_refs[i] = ParticleRefCompact(bridge.particle_refs[i]);
		atom_types[i] = bridge.atom_types[i];
	}
	n_singlebonds = bridge.n_singlebonds;
	for (int i = 0; i < n_singlebonds; i++) {
		singlebonds[i] = bridge.singlebonds[i];
	}
	n_anglebonds = bridge.n_anglebonds;
	for (int i = 0; i < n_anglebonds; i++) {
		anglebonds[i] = bridge.anglebonds[i];
	}
	n_dihedrals = bridge.n_dihedrals;
	for (int i = 0; i < n_dihedrals; i++) {
		dihedrals[i] = bridge.dihedrals[i];
	}

	if (verbose) {
		printf("Loading bridge with %d particles %d bonds %d angles %d dihedrals\n", n_particles, n_singlebonds, n_anglebonds, n_dihedrals);
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

GenericBond::GenericBond(ParticleRef* particle_refs, int n) {
	for (int i = 0; i < n; i++) {
		particles[n_particles++] = particle_refs[i];


		if (compound_ids[0] == -1) {
			compound_ids[0] = particles[i].compound_id;
		}
		else if (compound_ids[0] != particles[i].compound_id) {
			compound_ids[1] = particles[i].compound_id;
		}
	}
	if (compound_ids[0] > compound_ids[1] && compound_ids[1] != -1) { std::swap(compound_ids[0], compound_ids[1]); }
}

bool GenericBond::spansTwoCompounds() {	// Simply check if any id has been put into index 1
	return (compound_ids[1] != -1);
}
bool GenericBond::allParticlesExist() {	// Returns false if any particle in bond does not exist, for example when ignoring hydrogens
	for (int i = 0; i < n_particles; i++) {
		if (particles[i].local_id_compound == -1)
			return false;
	}
	return true;
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