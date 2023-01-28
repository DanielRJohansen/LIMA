#include "Bodies.cuh"




//__host__ void Compound::init() {
//	center_of_mass = calcCOM();
//	//printf("")
//	//radius = singlebonds[0].reference_dist * n_particles * 0.5f;
//	//center_of_mass.print('C');
//	//printf("Radius %f\n", radius);
//}




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


static void CoordArrayQueueHelpers::copyInitialCoordConfiguration(SolventCoord* coords, SolventCoord* coords_prev, SolventCoord* solventcoordarray_circular_queue) {
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















Float3 Compound_Carrier::calcCOM() {
	Float3 com;
	for (int i = 0; i < n_particles; i++) {
		com += state.positions[i] / static_cast<float>(n_particles);
	}
	return com;
}

__host__ void Compound_Carrier::addParticle(int atomtype_id, Float3 pos) {	// TODO: Delete?
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

__host__ void Compound_Carrier::addParticle(int atomtype_id, Float3 pos, int atomtype_color_id, int global_id) {
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

__host__ void Compound_Carrier::calcParticleSphere() {
	Float3 com = calcCOM();// calcCOM(compound);

	//float furthest = LONG_MIN;
	float furthest = FLT_MIN;
	//float closest = LONG_MAX;
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
		compound_bridges[i] = CompoundBridgeCompact(&bundle->compound_bridges[i], verbose);
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

void CompoundBridge::addGenericBond(PairBond pb) {
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






CompoundBridgeCompact::CompoundBridgeCompact(CompoundBridge* bridge, bool verbose) :
	compound_id_left{ bridge->compound_id_left },
	compound_id_right{ bridge->compound_id_right }
{
	n_particles = bridge->n_particles;

	for (int i = 0; i < n_particles; i++) {
		particle_refs[i] = ParticleRefCompact(bridge->particle_refs[i]);
		atom_types[i] = bridge->atom_types[i];
	}
	n_singlebonds = bridge->n_singlebonds;
	for (int i = 0; i < n_singlebonds; i++) {
		singlebonds[i] = bridge->singlebonds[i];
	}
	n_anglebonds = bridge->n_anglebonds;
	for (int i = 0; i < n_anglebonds; i++) {
		anglebonds[i] = bridge->anglebonds[i];
	}
	n_dihedrals = bridge->n_dihedrals;
	for (int i = 0; i < n_dihedrals; i++) {
		dihedrals[i] = bridge->dihedrals[i];
	}

	if (verbose) {
		printf("Loading bridge with %d particles %d bonds %d angles %d dihedrals\n", n_particles, n_singlebonds, n_anglebonds, n_dihedrals);
	}

}



PairBond::PairBond(int id1, int id2, float b0, float kb) : b0(b0), kb(kb) {
	// This is only for loading the forcefield, so the ID's refers to id's given in .conf file!
	atom_indexes[0] = id1;
	atom_indexes[1] = id2;
}

AngleBond::AngleBond(int id1, int id2, int id3, float theta_0, float k_theta) : theta_0(theta_0), k_theta(k_theta) {
	atom_indexes[0] = id1;
	atom_indexes[1] = id2;
	atom_indexes[2] = id3;
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