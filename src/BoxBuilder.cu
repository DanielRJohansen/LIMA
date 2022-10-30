#include "BoxBuilder.cuh"



void BoxBuilder::buildBox(Simulation* simulation) {
	printf("Building box...\n");
	simulation->box->compounds = new Compound[MAX_COMPOUNDS];
	simulation->box->solvents = new Solvent[MAX_SOLVENTS];


	simulation->box->compound_state_array = new CompoundState[MAX_COMPOUNDS];
	cudaMalloc(&simulation->box->compound_state_array_next, sizeof(CompoundState) * MAX_COMPOUNDS);
	cudaMalloc(&simulation->box->solvents_next, sizeof(Solvent) * MAX_SOLVENTS);


	simulation->box->solvent_neighborlists = new NeighborList[MAX_SOLVENTS];	
	simulation->box->compound_neighborlists = new NeighborList[MAX_COMPOUNDS];


	simulation->box->dt = simulation->dt;
}

void BoxBuilder::addSingleMolecule(Simulation* simulation, Molecule* molecule) {
	Float3 desired_molecule_center = Float3(BOX_LEN_HALF);
	Float3 offset = desired_molecule_center - molecule->calcCOM();

	printf("Molecule offset for centering: ");
	offset.print(' ');
	most_recent_offset_applied = offset;			// Needed so solvents can be offset identically later. Not needed if making solvent positions using LIMA

	for (int c = 0; c < molecule->n_compounds; c++) {
		Compound* compound = &molecule->compounds[c];
		for (int i = 0; i < compound->n_particles; i++) {
			compound->prev_positions[i] += offset;
		}
		integrateCompound(compound, simulation);
	}

	simulation->total_compound_particles = molecule->n_atoms_total;						// TODO: Unknown behavior, if multiple molecules are added!
	simulation->total_particles += molecule->n_atoms_total;

	simulation->box->bridge_bundle = new CompoundBridgeBundleCompact;
	*simulation->box->bridge_bundle = *molecule->compound_bridge_bundle;					// TODO: Breaks if multiple compounds are added, as only one bridgebundle can exist for now!

	simulation->box->bonded_particles_lut_manager = new BondedParticlesLUTManager{ *molecule->bonded_particles_lut_manager };
	//*simulation->box->bonded_particles_lut_manager = &molecule->bonded_particleslut_manager;
	//exit(0);
	

	printf("Molecule added to box\n");
}

void BoxBuilder::addScatteredMolecules(Simulation* simulation, Compound* molecule, int n_copies)
{
	placeMultipleCompoundsRandomly(simulation, molecule, n_copies);
	printf("Scattered %d Compounds in box\n", simulation->box->n_compounds);
}

void BoxBuilder::finishBox(Simulation* simulation) {
	// Load meta information
	simulation->copyBoxVariables();
	printf("Box contains %d compounds, %d bridges and %d solvents\n\n", simulation->n_compounds, simulation->n_bridges, simulation->n_solvents);


	// Need this variable both on host and device
	simulation->total_particles_upperbound = simulation->box->n_compounds * MAX_COMPOUND_PARTICLES + simulation->box->n_solvents;											// BAD AMBIGUOUS AND WRONG CONSTANTS
	printf("Total particles upperbound: %d\n", simulation->total_particles_upperbound);
	simulation->box->total_particles_upperbound = simulation->total_particles_upperbound;											// BAD AMBIGUOUS AND WRONG CONSTANTS


	cudaMemcpy(simulation->box->compound_state_array_next, simulation->box->compound_state_array, sizeof(CompoundState) * MAX_COMPOUNDS, cudaMemcpyHostToDevice);	// Just make sure they have the same n_particles info...
	cudaMemcpy(simulation->box->solvents_next, simulation->box->solvents, sizeof(Solvent) * MAX_SOLVENTS, cudaMemcpyHostToDevice);

	
	// Permanent Outputs for energy & trajectory analysis
	int n_points = simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER;
	printf("n points %d\n", n_points);
	printf("Malloc %.2f KB on device for data buffers\n",(float) ((sizeof(double) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER + sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER) * 1e-3));
	printf("Malloc %.2f MB on host for data buffers\n", (float) ((sizeof(double) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER + sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER) * 1e-6));
	cudaMallocManaged(&simulation->box->potE_buffer, sizeof(float) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER);	// Can only log molecules of size 3 for now...
	simulation->potE_buffer = new float[simulation->total_particles_upperbound * simulation->n_steps];

	cudaMallocManaged(&simulation->box->traj_buffer, sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER);
	simulation->traj_buffer = new Float3[simulation->total_particles_upperbound * simulation->n_steps];

	simulation->temperature_buffer = new float[SIMULATION_STEPS / STEPS_PER_THERMOSTAT + 1];









	// TRAINING DATA and TEMPRARY OUTPUTS
	int n_loggingdata_device = 10 * STEPS_PER_LOGTRANSFER;
	uint64_t n_traindata_device = N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->n_compounds * STEPS_PER_TRAINDATATRANSFER;
	long double total_bytes = static_cast<long double>(sizeof(float) * static_cast<long double>(n_loggingdata_device) + sizeof(Float3) * n_traindata_device);
	printf("Reserving %.2f MB device mem for logging + training data\n", (float) ((total_bytes) * 1e-6));
	cudaMallocManaged(&simulation->box->outdata, sizeof(float) * 10 * STEPS_PER_LOGTRANSFER);	// 10 data streams for 10k steps. 1 step at a time.

	cudaMallocManaged(&simulation->box->data_GAN, sizeof(Float3) * N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->n_compounds * STEPS_PER_TRAINDATATRANSFER);


	uint64_t n_loggingdata_host = 10 * simulation->n_steps;
	uint64_t n_traindata_host = N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->n_compounds * (uint64_t) simulation->n_steps;
	printf("Reserving %.2f GB host mem for logging + training data\n",(float) (sizeof(Float3) * n_traindata_host + sizeof(float) * n_loggingdata_host) * 1e-9);
	simulation->logging_data = new float[n_loggingdata_host];
	simulation->traindata_buffer = new Float3[n_traindata_host];





	cudaDeviceSynchronize();
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Error during log-data mem. allocation\n");
		exit(1);
	}


	simulation->box->moveToDevice();
	printf("Boxbuild complete!\n\n\n");
}





int BoxBuilder::solvateBox(Simulation* simulation)
{
	simulation->box->solvents = new Solvent[MAX_SOLVENTS];

	const int bodies_per_dim = static_cast<int>(ceil(cbrt((double)N_SOLVATE_MOLECULES)));
	const float dist_between_compounds = (BOX_LEN) / static_cast<float>(bodies_per_dim);	// dist_per_index
	const float base = box_base + dist_between_compounds / 2.f;
	printf("Bodies per dim: %d. Dist per dim: %.3f\n", bodies_per_dim, dist_between_compounds);


	for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
				if (simulation->box->n_solvents == N_SOLVATE_MOLECULES)
					break;

				Float3 solvent_center = Float3(base + dist_between_compounds * static_cast<float>(x_index), base + dist_between_compounds * static_cast<float>(y_index), base + dist_between_compounds * static_cast<float>(z_index));
				//double solvent_radius = 0.2;

				if (spaceAvailable(simulation->box, solvent_center)) {
					simulation->box->solvents[simulation->box->n_solvents++] = createSolvent(solvent_center, simulation->dt);
				}
			}
		}
	}
	simulation->total_particles += simulation->box->n_solvents;
	printf("%d solvents added to box\n", simulation->box->n_solvents);
	return simulation->box->n_solvents;
}

int BoxBuilder::solvateBox(Simulation* simulation, std::vector<Float3>* solvent_positions)	// Accepts the position of the center or Oxygen of a solvate molecule. No checks are made wh
{
	for (Float3 sol_pos : *solvent_positions) {
		if (simulation->box->n_solvents == MAX_SOLVENTS) {
			printf("Too many solvents added!\n\n\n\n");
			exit(1);
		}

		sol_pos += most_recent_offset_applied;			// So solvents are re-aligned with an offsat molecule.

		if (spaceAvailable(simulation->box, sol_pos) && simulation->box->n_solvents < SOLVENT_TESTLIMIT) {						// Should i check? Is this what energy-min is for?
			simulation->box->solvents[simulation->box->n_solvents++] = createSolvent(sol_pos, simulation->dt);
		}
	}

	simulation->total_particles += simulation->box->n_solvents;
	printf("%lu of %lu solvents added to box\n", simulation->box->n_solvents, solvent_positions->size());
	return simulation->box->n_solvents;
}

void BoxBuilder::integrateCompound(Compound* compound, Simulation* simulation)
{
	compound->init();
	CompoundState* state = &simulation->box->compound_state_array[simulation->box->n_compounds];
	Float3 compound_united_vel = Float3(random(), random(), random()).norm() * v_rms * 0.f;			// Giving individual comp in molecule different uniform vels is sub-optimal...

	for (int i = 0; i < compound->n_particles; i++) {
		state->positions[i] = compound->prev_positions[i];
		state->n_particles++;
	}


	for (int i = 0; i < compound->n_particles; i++) {
		Float3 atom_pos_sub1 = state->positions[i] - compound_united_vel * simulation->dt;
//		compound->particles[i].pos_tsub1 = atom_pos_sub1;												// Overwrite prev pos here, since we have assigned the former prev pos to the state buffer.
		compound->prev_positions[i] = atom_pos_sub1;												// Overwrite prev pos here, since we have assigned the former prev pos to the state buffer.
		//compound->prev_positions[i].print('p');
	}

	simulation->box->compounds[simulation->box->n_compounds++] = *compound;
}



Solvent BoxBuilder::createSolvent(Float3 com, float dt) {
	Float3 solvent_vel = Float3(random(), random(), random()).norm() * v_rms * 0.f;		// TODO: I dont know, but i think we need to freeze solvents to avoid unrealisticly large forces at step 1
	return Solvent(com, com - solvent_vel * dt);
}






/*
void BoxBuilder::compoundLinker(Simulation* simulation) {
	for (int i = 0; i < simulation->box->n_compounds; i++) {
		for (int j = i+1; j < simulation->box->n_compounds; j++) {
			simulation->box->compound_neighborlists[i].addId(j, NeighborList::NEIGHBOR_TYPE::COMPOUND);
			simulation->box->compound_neighborlists[j].addId(i, NeighborList::NEIGHBOR_TYPE::COMPOUND);
		}
	}
}

void BoxBuilder::solvateLinker(Simulation* simulation)
{
	for (int i = 0; i < simulation->box->n_solvents; i++) {
		Solvent* self = &simulation->box->solvents[i];
		for (int j = i; j < simulation->box->n_solvents; j++) {														// +1 here!
			Solvent* other = &simulation->box->solvents[j];
			if (i != j) {
				if ((self->pos - other->pos).len() < (CUTOFF)) {
					simulation->box->solvent_neighborlists[i].addId(j, NeighborList::NEIGHBOR_TYPE::SOLVENT);
					simulation->box->solvent_neighborlists[j].addId(i, NeighborList::NEIGHBOR_TYPE::SOLVENT);
				}
			}
		}
	}
}

void BoxBuilder::solvateCompoundCrosslinker(Simulation* simulation)
{
	for (int i = 0; i < simulation->box->n_compounds; i++) {
		Compound* compound = &simulation->box->compounds[i];
		for (int j = 0; j < simulation->box->n_solvents; j++) {
			Solvent* solvent= &simulation->box->solvents[j];
			if ((compound->center_of_mass - solvent->pos).len() < CUTOFF) {
				simulation->box->compound_neighborlists[i].addId(j, NeighborList::NEIGHBOR_TYPE::SOLVENT);
				simulation->box->solvent_neighborlists[j].addId(i, NeighborList::NEIGHBOR_TYPE::COMPOUND);
			}
		}
	}
	printf("Compound n0 solvents: %d\n", simulation->box->compound_neighborlists[0].n_solvent_neighbors);
}

*/

void BoxBuilder::placeMultipleCompoundsRandomly(Simulation* simulation, Compound* template_compound, int n_copies)
{
	int copies_placed = 0;
	while (copies_placed < n_copies) {
		Compound* c = randomizeCompound(template_compound);

		if (spaceAvailable(simulation->box, c)) {
			integrateCompound(c, simulation);
			copies_placed++;
		}
		delete c;
	}


	// Temporary check that no to molecules placed are colliding.
	for (int i = 0; i < simulation->box->n_compounds; i++) {
		Compound* c = &simulation->box->compounds[i];
		for (int ii = 0; ii < simulation->box->n_compounds; ii++) {
			Compound* c2 = &simulation->box->compounds[ii];
			if (ii != i) {
				if (!verifyPairwiseParticleMindist(c, c2)) {
					printf("Illegal compound positioning %d %d\n", i, ii);
					exit(0);
				}				
			}
		}
	}
}

Compound* BoxBuilder::randomizeCompound(Compound* original_compound)
{
	Compound* compound = new Compound;
	*compound = *original_compound;

	Float3 xyz_rot = get3Random() * (2.f*PI);
	//rotateCompound(compound, xyz_rot);


	Float3 xyz_target = (get3Random() * 0.6 + Float3(0.2))* BOX_LEN;
	Float3 xyz_mov = xyz_target - original_compound->calcCOM();// calcCompoundCom(original_compound);
	moveCompound(compound, xyz_mov);

	return compound;
}

void BoxBuilder::moveCompound(Compound* compound, Float3 vector)
{
	for (int i = 0; i < compound->n_particles; i++) {
		compound->prev_positions[i] += vector;
		//compound->particles[i].pos_tsub1 += vector;
	}		
}

void BoxBuilder::rotateCompound(Compound* compound, Float3 xyz_rot)
{
	Float3 vec_to_origo = Float3(0, 0, 0) - compound->calcCOM();
	moveCompound(compound, vec_to_origo);

	for (int i = 0; i < compound->n_particles; i++) {
		compound->prev_positions[i].rotateAroundOrigo(xyz_rot);
		//compound->particles[i].pos_tsub1.rotateAroundOrigo(xyz_rot);
	}
		

	moveCompound(compound, vec_to_origo * -1);
}

BoundingBox BoxBuilder::calcCompoundBoundingBox(Compound* compound)
{
	BoundingBox bb(Float3(9999, 9999, 9999), Float3(-9999, -9999, -9999));
	for (int i = 0; i < compound->n_particles; i++) {
		//Float3 pos = compound->particles[i].pos_tsub1;
		Float3 pos = compound->prev_positions[i];
		for (int dim = 0; dim < 3; dim++) {
			*bb.min.placeAt(dim) = std::min(bb.min.at(dim), pos.at(dim));
			*bb.max.placeAt(dim) = std::max(bb.max.at(dim), pos.at(dim));
		}
	}
	return bb;
}

bool BoxBuilder::spaceAvailable(Box* box, Compound* compound)
{
	BoundingBox bb_a = calcCompoundBoundingBox(compound);
	bb_a.addPadding(MIN_NONBONDED_DIST);
	for (size_t c_index = 0; c_index < box->n_compounds; c_index++) {
		BoundingBox bb_b = calcCompoundBoundingBox(&box->compounds[c_index]);

		/*if (box->n_compounds == 13) {
			printf("\n c index %d\n", c_index);
			bb_a.min.print('a');
			bb_a.max.print('a');
			bb_b.min.print('b');
			bb_b.max.print('b');
		}*/

		if (bb_a.intersects(bb_b)) {
			//printf("Verifying %d %d\n", box->n_compounds, c_index);
			if (!verifyPairwiseParticleMindist(compound, &box->compounds[c_index]))
				return false;			
		}
	}
	return true;
}

float minDist(Compound* compound, Float3 particle_pos) {
	float mindist = 999999;
	for (size_t i = 0; i < compound->n_particles; i++) {
		float dist = (compound->prev_positions[i] - particle_pos).len();
		mindist = std::min(mindist, dist);
	}
	return mindist;
}

bool BoxBuilder::spaceAvailable(Box* box, Float3 particle_center)
{
	for (int c_index = 0; c_index < box->n_compounds; c_index++) {
		if (minDist(&box->compounds[c_index], particle_center) < 0.2f)
			return false;



		//BoundingBox bb = calcCompoundBoundingBox(&box->compounds[c_index]);
		//bb.addPadding(MIN_NONBONDED_DIST);
		//if (bb.pointIsInBox(particle_center)) {
		//	return false;
		//}
	}
	for (int si = 0; si < box->n_solvents; si++) {
		float dist = (box->solvents[si].pos - particle_center).len();
		if (dist < 0.25)
			return false;		
	}

	return true;
}



bool BoxBuilder::verifyPairwiseParticleMindist(Compound* a, Compound* b)
{
	for (int ia = 0; ia < a->n_particles; ia++) {
		for (int ib = 0; ib < b->n_particles; ib++) {
			//Float3 pos_a = a->particles[ia].pos_tsub1;
			//Float3 pos_b = b->particles[ib].pos_tsub1;
			Float3 pos_a = a->prev_positions[ia];
			Float3 pos_b = b->prev_positions[ib];

			float dist = (pos_a - pos_b).len();
			if (dist < MIN_NONBONDED_DIST)
				return false;
		}
	}
	return true;
}











/*
int BoxBuilder::solvateBox(Simulation* simulation)
{
	int bodies_per_dim = ceil(cbrt((double)N_SOLVATE_MOLECULES));
	double dist_between_compounds = (BOX_LEN) / (double)bodies_per_dim;	// dist_per_index
	double base = box_base + dist_between_compounds / 2.f;
	printf("Bodies per dim: %d. Dist per dim: %.3f\n", bodies_per_dim, dist_between_compounds);


	for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
				if (simulation->box->n_compounds == N_SOLVATE_MOLECULES)
					break;

				Float3 compound_center = Float3(base + dist_between_compounds * (double)x_index, base + dist_between_compounds * (double)y_index, base + dist_between_compounds * (double)z_index);
				double compound_radius = 0.2;

				if (spaceAvailable(compound_center, compound_radius)) {
					simulation->box->compounds[simulation->box->n_compounds++] = createSolvent(
						compound_center,
						simulation->box->n_compounds,
						&simulation->box->compound_state_array[simulation->box->n_compounds],
						&simulation->box->compound_neighborlist_array[simulation->box->n_compounds],
						simulation->dt
					);
				}
			}
		}
	}
	return simulation->box->n_compounds;
}
*/
