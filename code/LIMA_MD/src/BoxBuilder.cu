#include "BoxBuilder.cuh"
#include "EngineUtils.cuh"

#include "Printer.h"

using namespace LIMA_Print;

void BoxBuilder::buildBox(Simulation* simulation) {
	//printf("Building box...\n");
	printH1("Building box", true, false);
	simulation->box->compounds = new Compound[MAX_COMPOUNDS];
	//simulation->box->solvents = new Solvent[MAX_SOLVENTS];


//	cudaMalloc(&simulation->box->solvents_next, sizeof(Solvent) * MAX_SOLVENTS);

	//simulation->box->compound_coord_array = new CompoundCoords[MAX_COMPOUNDS];								// Need to manipulate this on host
	//simulation->box->compound_coord_array_prev = new CompoundCoords[MAX_COMPOUNDS];							// Need to manipulate this on host
	//cudaMalloc(&simulation->box->compound_coord_array_next, sizeof(CompoundCoords) * MAX_COMPOUNDS);	// Only need this on device


	// This is where the coords will eventually reside (device)
	const uint64_t n_bytes_compoundcoords = sizeof(CompoundCoords) * MAX_COMPOUNDS * STEPS_PER_LOGTRANSFER;
	cudaMalloc(&simulation->box->coordarray_circular_queue, n_bytes_compoundcoords);

	//const uint64_t n_bytes_solventcoords = sizeof(SolventCoord) * MAX_SOLVENTS * STEPS_PER_LOGTRANSFER;
	//cudaMalloc(&simulation->box->solventcoordarray_circular_queue, n_bytes_solventcoords);

	//const uint64_t n_bytes_solventblockgrids = sizeof(SolventBlockGrid) * STEPS_PER_SOLVENTBLOCKTRANSFER;
	//cudaMalloc(&simulation->box->solventblockgrid_circular_queue, n_bytes_solventblockgrids);
	SolventBlockHelpers::createSolventblockGrid(&simulation->box->solventblockgrid_circular_queue);
	SolventBlockHelpers::createSolventblockTransfermodules(&simulation->box->transfermodule_array);


	// This is very the coords reside while build (host)
	coordarray = new CompoundCoords[MAX_COMPOUNDS];
	coordarray_prev = new CompoundCoords[MAX_COMPOUNDS];

	solventcoords = new SolventCoord[MAX_SOLVENTS];
	solventcoords_prev = new SolventCoord[MAX_SOLVENTS];

	solventblocks = new SolventBlockGrid{};
	solventblocks_prev = new SolventBlockGrid{};
	SolventBlockHelpers::setupBlockMetaOnHost(solventblocks, solventblocks_prev);

	CompoundGrid::createCompoundGrid(&simulation->box->compound_grid);
	//EngineUtils::genericErrorCheck("")

	simulation->box->solvent_neighborlists = new NeighborList[MAX_SOLVENTS];	
	simulation->box->compound_neighborlists = new NeighborList[MAX_COMPOUNDS];

	simulation->box->bridge_bundle = new CompoundBridgeBundleCompact{};
	//simulation->box->bonded_particles_lut_manager = new BondedParticlesLUTManager{};

	simulation->box->dt = simulation->dt;

	cudaDeviceSynchronize();
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Error during buildBox()\n");
		exit(1);
	}
}


void BoxBuilder::addCompoundCollection(Simulation* simulation, const CompoundCollection& compound_collection) {
	for (const CompoundFactory& compound : compound_collection.compounds) {
		integrateCompound(compound, simulation);
	}

	simulation->total_compound_particles = compound_collection.total_compound_particles;						// TODO: Unknown behavior, if multiple molecules are added!
	simulation->total_particles += compound_collection.total_compound_particles;


	*simulation->box->bridge_bundle = compound_collection.bridgebundle;					// TODO: Breaks if multiple compounds are added, as only one bridgebundle can exist for now!

	simulation->box->bonded_particles_lut_manager = compound_collection.bp_lut_manager;	// TODO: release a unique ptr here!

	printf("CompoundCollection added to box\n");
}

void BoxBuilder::finishBox(Simulation* simulation, const ForceField_NB& forcefield) {
	// Load meta information
	simulation->copyBoxVariables();
	printf("Box contains %d compounds, %d bridges and %d solvents\n\n", simulation->n_compounds, simulation->n_bridges, simulation->n_solvents);


	simulation->box->forcefield_device_box = new ForceField_NB{ forcefield };// Copy

	//cudaMemcpy(simulation->box->solvents_next, simulation->box->solvents, sizeof(Solvent) * MAX_SOLVENTS, cudaMemcpyHostToDevice);


	// Move the positions to the appropriate places in the circular queue
	CompoundCoords::copyInitialCoordConfiguration(coordarray, coordarray_prev, simulation->box->coordarray_circular_queue);
	//SolventCoord::copyInitialCoordConfiguration(solventcoords, solventcoords_prev, simulation->box->solventcoordarray_circular_queue);
	SolventBlockHelpers::copyInitialConfiguration(*solventblocks, *solventblocks_prev, simulation->box->solventblockgrid_circular_queue);

	
	// Permanent Outputs for energy & trajectory analysis
	int n_points = simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER;
	printf("n points %d\n", n_points);
	printf("Malloc %.2f MB on device for data buffers\n",(float) ((sizeof(double) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER + sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER) * 1e-6));
	printf("Malloc %.2f MB on host for data buffers\n", (float) ((sizeof(double) * simulation->total_particles_upperbound * simulation->n_steps + sizeof(Float3) * simulation->total_particles_upperbound * simulation->n_steps) * 1e-6));
	cudaMallocManaged(&simulation->box->potE_buffer, sizeof(float) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER);	// Can only log molecules of size 3 for now...
	simulation->potE_buffer = new float[simulation->total_particles_upperbound * simulation->n_steps];

	cudaMallocManaged(&simulation->box->traj_buffer, sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER);
	simulation->traj_buffer = new Float3[simulation->total_particles_upperbound * simulation->n_steps];

	simulation->temperature_buffer = new float[simulation->n_steps / STEPS_PER_THERMOSTAT + 1];

#ifdef USEDEBUGF3
	uint64_t bytes_for_debugf3 = sizeof(Float3) * DEBUGDATAF3_NVARS * simulation->total_particles_upperbound * simulation->n_steps;
	cudaMallocManaged(&simulation->box->debugdataf3, bytes_for_debugf3);
#endif






	// TRAINING DATA and TEMPRARY OUTPUTS
	int n_loggingdata_device = 10 * STEPS_PER_LOGTRANSFER;
	uint64_t n_traindata_device = static_cast<uint64_t>(N_DATAGAN_VALUES) * MAX_COMPOUND_PARTICLES * simulation->n_compounds * STEPS_PER_TRAINDATATRANSFER;
	long double total_bytes = static_cast<long double>(sizeof(float) * static_cast<long double>(n_loggingdata_device) + sizeof(Float3) * n_traindata_device);
	printf("Reserving %.4f MB device mem for logging + training data\n", (float) ((total_bytes) * 1e-6));
	cudaMallocManaged(&simulation->box->outdata, sizeof(float) * 10 * STEPS_PER_LOGTRANSFER);	// 10 data streams for 10k steps. 1 step at a time.

	cudaMallocManaged(&simulation->box->data_GAN, sizeof(Float3) * N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->n_compounds * STEPS_PER_TRAINDATATRANSFER);


	uint64_t n_loggingdata_host = 10 * simulation->n_steps;
	uint64_t n_traindata_host = N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->n_compounds * (uint64_t) simulation->n_steps;
	printf("Reserving %.4f GB host mem for logging + training data\n",(float) (sizeof(Float3) * n_traindata_host + sizeof(float) * n_loggingdata_host) * 1e-9);
	simulation->logging_data = new float[n_loggingdata_host];
	simulation->traindata_buffer = new Float3[n_traindata_host];





	cudaDeviceSynchronize();
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Error during log-data mem. allocation\n");
		exit(1);
	}

	printf("Total particles upperbound: %d\n", simulation->total_particles_upperbound);
	printf("Max particles in compound: %d", MAX_COMPOUND_PARTICLES);

	simulation->box->moveToDevice();
	//printf("Boxbuild complete!\n\n\n");
	printH1("Boxbuild complete!", false, true);
}





int BoxBuilder::solvateBox(Simulation* simulation)
{
	//simulation->box->solvents = new Solvent[MAX_SOLVENTS];
	//
	//// First we calculate how to set up the particles in a perfect grid
	//const int bodies_per_dim = static_cast<int>(ceil(cbrt((double)N_SOLVATE_MOLECULES)));
	//const float dist_between_particles = (BOX_LEN) / static_cast<float>(bodies_per_dim);	// dist_per_index
	//const float base = box_base + dist_between_particles / 2.f;
	//printf("Bodies per dim: %d. Dist per dim: %.3f\n", bodies_per_dim, dist_between_particles);


	//for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
	//	for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
	//		for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
	//			if (simulation->box->n_solvents == N_SOLVATE_MOLECULES)
	//				break;

	//			Float3 solvent_center = Float3(base + dist_between_particles * static_cast<float>(x_index), base + dist_between_particles * static_cast<float>(y_index), base + dist_between_particles * static_cast<float>(z_index));
	//			
	//			// Randomly offset the particle within 80% of the perfect grid
	//			solvent_center += (get3Random() - Float3(0.5f)) * 0.8f * dist_between_particles;

	//			if (spaceAvailable(simulation->box, solvent_center)) {
	//				simulation->box->solvents[simulation->box->n_solvents++] = createSolvent(solvent_center, simulation->dt);
	//			}
	//		}
	//	}
	//}
	//simulation->total_particles += simulation->box->n_solvents;
	//printf("%d solvents added to box\n", simulation->box->n_solvents);
	//return simulation->box->n_solvents;
	return 0;
}

int BoxBuilder::solvateBox(Simulation* simulation, const std::vector<Float3>& solvent_positions)	// Accepts the position of the center or Oxygen of a solvate molecule. No checks are made wh
{
	for (Float3 sol_pos : solvent_positions) {
		if (simulation->box->n_solvents == MAX_SOLVENTS) {
			printf("Too many solvents added!\n\n\n\n");
			exit(1);
		}

		sol_pos += most_recent_offset_applied;			// So solvents are re-aligned with an offsat molecule.

		if (spaceAvailable(simulation->box, sol_pos, true) && simulation->box->n_solvents < SOLVENT_TESTLIMIT) {						// Should i check? Is this what energy-min is for?
			LimaPosition position = LIMAPOSITIONSYSTEM::createLimaPosition(sol_pos);
			const SolventCoord solventcoord = LIMAPOSITIONSYSTEM::createSolventcoordFromAbsolutePosition(position);

			
			solventcoords[simulation->box->n_solvents] = solventcoord;	// REMOVE soon?
			solventcoords_prev[simulation->box->n_solvents] = solventcoord;	// TODO: Add a subtraction here for initial velocity.


			solventblocks->addSolventToGrid(solventcoord, simulation->box->n_solvents);
			solventblocks_prev->addSolventToGrid(solventcoord, simulation->box->n_solvents);

			if (solventcoord.origo == NodeIndex{ 1,0,0 } &&
				(solventblocks->getBlockPtr(1)->n_solvents == 26 || solventblocks->getBlockPtr(1)->n_solvents == 37)
				) {
				int h = 0;
			}

			simulation->box->n_solvents++;
		}
	}

	for (int sbi = 0; sbi < SolventBlockGrid::blocks_total; sbi++) {
		auto sb = solventblocks->getBlockPtr(sbi);
		for (int i = 0; i < sb->n_solvents; i++) {
			auto posi = sb->rel_pos[i];


			for (int sbj = sbi; sbj < SolventBlockGrid::blocks_total; sbj++) {
				auto sb2 = solventblocks->getBlockPtr(sbj);
				for (int j = 0; j < sb2->n_solvents; j++) {

					if (sbi == sbj && i == j) { continue; }	// same solvent

					auto posj = sb2->rel_pos[j];

					auto dist = EngineUtils::calcDistance(sb->origo, posi, sb2->origo, posj);
					if (dist < 0.1) {
						int a = 0;
					}
				}
			}
		}
	}

	simulation->total_particles += simulation->box->n_solvents;
	printf("%lu of %lld solvents added to box\n", simulation->box->n_solvents, solvent_positions.size());
	return simulation->box->n_solvents;
}




/*
* These two funcitons are in charge of normalizing ALL coordinates!!
*/
//void BoxBuilder::integrateCompound(CompoundCarrier* compound, Simulation* simulation)
//{
//	compound->init();
//
//	std::vector<LimaPosition> positions;
//	std::vector<LimaPosition> positions_prev;
//	positions.reserve(MAX_COMPOUND_PARTICLES);
//	positions_prev.reserve(MAX_COMPOUND_PARTICLES);
//
//	Float3 compound_united_vel = Float3(random(), random(), random()).norm() * v_rms * 0.f;			// Giving individual comp in molecule different uniform vels is sub-optimal...
//
//
//
//	for (int i = 0; i < compound->n_particles; i++) {		
//		positions.push_back(LIMAPOSITIONSYSTEM::createLimaPosition(compound->state.positions[i]));
//
//		const Float3 pos_prev_nm = (compound->state.positions[i] - compound_united_vel * simulation->dt);
//		positions_prev.push_back(LIMAPOSITIONSYSTEM::createLimaPosition(pos_prev_nm));
//	}
//
//
//	coordarray[simulation->box->n_compounds] = LIMAPOSITIONSYSTEM::positionCompound(positions, 0);
//	coordarray_prev[simulation->box->n_compounds] = LIMAPOSITIONSYSTEM::positionCompound(positions_prev, 0);
//
//	simulation->box->compounds[simulation->box->n_compounds++] = *compound;
//}

void BoxBuilder::integrateCompound(const CompoundFactory& compound, Simulation* simulation)
{
	std::vector<LimaPosition> positions;
	std::vector<LimaPosition> positions_prev;
	positions.reserve(MAX_COMPOUND_PARTICLES);
	positions_prev.reserve(MAX_COMPOUND_PARTICLES);

	Float3 compound_united_vel = Float3(random(), random(), random()).norm() * v_rms * 0.f;			// Giving individual comp in molecule different uniform vels is sub-optimal...

	for (int i = 0; i < compound.n_particles; i++) {
		const Float3& extern_position = compound.positions[i];
		positions.push_back(LIMAPOSITIONSYSTEM::createLimaPosition(extern_position));

		if (compound.gro_ids[i] == 177 ) {
			int a = 0;
		}

		const Float3 pos_prev_nm = (extern_position - compound_united_vel * simulation->dt);
		positions_prev.push_back(LIMAPOSITIONSYSTEM::createLimaPosition(pos_prev_nm));
	}

	coordarray[simulation->box->n_compounds] = LIMAPOSITIONSYSTEM::positionCompound(positions, 0);
	coordarray_prev[simulation->box->n_compounds] = LIMAPOSITIONSYSTEM::positionCompound(positions_prev, 0);

	simulation->box->compounds[simulation->box->n_compounds++] = Compound{ compound };
}

//Solvent BoxBuilder::createSolvent(Float3 com, float dt) {
//	com = com / NORMALIZER * 1e+6;	// concvert to normalized [fm]
//	Float3 solvent_vel = Float3(random(), random(), random()).norm() * v_rms * VEL_RMS_SCALAR / NORMALIZER;		// TODO: I dont know, but i think we need to freeze solvents to avoid unrealisticly large forces at step 1
//	return Solvent(com, com - solvent_vel * dt);
//}
/*
* These two funcitons are in charge of normalizing ALL coordinates!!
*/



//void BoxBuilder::placeMultipleCompoundsRandomly(Simulation* simulation, Compound* template_compound, int n_copies)
//{
//	int copies_placed = 0;
//	while (copies_placed < n_copies) {
//		Compound* c = randomizeCompound(template_compound);
//
//		if (spaceAvailable(simulation->box, c)) {
//			integrateCompound(c, simulation);
//			copies_placed++;
//		}
//		delete c;
//	}
//
//
//	// Temporary check that no to molecules placed are colliding.
//	for (uint32_t i = 0; i < simulation->box->n_compounds; i++) {
//		Compound* c = &simulation->box->compounds[i];
//		for (uint32_t ii = 0; ii < simulation->box->n_compounds; ii++) {
//			Compound* c2 = &simulation->box->compounds[ii];
//			if (ii != i) {
//				if (!verifyPairwiseParticleMindist(c, c2)) {
//					printf("Illegal compound positioning %d %d\n", i, ii);
//					exit(0);
//				}				
//			}
//		}
//	}
//}

//Compound* BoxBuilder::randomizeCompound(Compound* original_compound)
//{
//	Compound* compound = new Compound;
//	*compound = *original_compound;
//
//	Float3 xyz_rot = get3Random() * (2.f*PI);
//	//rotateCompound(compound, xyz_rot);
//
//
//	Float3 xyz_target = (get3Random() * 0.6f + Float3(0.2f))* BOX_LEN;
//	Float3 xyz_mov = xyz_target - original_compound->calcCOM();// calcCompoundCom(original_compound);
//	moveCompound(compound, xyz_mov);
//
//	return compound;
//}

//void BoxBuilder::moveCompound(Compound* compound, Float3 vector)
//{
//	for (int i = 0; i < compound->n_particles; i++) {
//		compound->prev_positions[i] += vector;
//		//compound->particles[i].pos_tsub1 += vector;
//	}		
//}
//
//void BoxBuilder::rotateCompound(Compound* compound, Float3 xyz_rot)
//{
//	Float3 vec_to_origo = Float3(0, 0, 0) - compound->calcCOM();
//	moveCompound(compound, vec_to_origo);
//
//	for (int i = 0; i < compound->n_particles; i++) {
//		compound->prev_positions[i].rotateAroundOrigo(xyz_rot);
//		//compound->particles[i].pos_tsub1.rotateAroundOrigo(xyz_rot);
//	}
//		
//
//	moveCompound(compound, vec_to_origo * -1);
//}
//
//BoundingBox BoxBuilder::calcCompoundBoundingBox(Compound* compound)
//{
//	BoundingBox bb(Float3(9999, 9999, 9999), Float3(-9999, -9999, -9999));
//	for (int i = 0; i < compound->n_particles; i++) {
//		//Float3 pos = compound->particles[i].pos_tsub1;
//		Float3 pos = compound->prev_positions[i];
//		for (int dim = 0; dim < 3; dim++) {
//			*bb.min.placeAt(dim) = std::min(bb.min.at(dim), pos.at(dim));
//			*bb.max.placeAt(dim) = std::max(bb.max.at(dim), pos.at(dim));
//		}
//	}
//	return bb;
//}

//bool BoxBuilder::spaceAvailable(Box* box, Compound* compound)
//{
//	BoundingBox bb_a = calcCompoundBoundingBox(compound);
//	bb_a.addPadding(MIN_NONBONDED_DIST);
//	for (size_t c_index = 0; c_index < box->n_compounds; c_index++) {
//		BoundingBox bb_b = calcCompoundBoundingBox(&box->compounds[c_index]);
//
//		if (bb_a.intersects(bb_b)) {
//			if (!verifyPairwiseParticleMindist(compound, &box->compounds[c_index]))
//				return false;			
//		}
//	}
//	return true;
//}

float minDist(CompoundState& compoundstate, Float3 particle_pos) {
	float mindist = 999999;
	for (size_t i = 0; i < compoundstate.n_particles; i++) {
		//float dist = EngineUtils::calcHyperDist(&compound->prev_positions[i], &particle_pos);
		float dist = EngineUtils::calcHyperDistNM(&compoundstate.positions[i], &particle_pos);		// Hmmm doesn't fit the namespace...
		mindist = std::min(mindist, dist);
	}
	return mindist;
}

bool BoxBuilder::spaceAvailable(Box* box, Float3 particle_center, bool verbose)
{
	particle_center = particle_center;
	for (uint32_t c_index = 0; c_index < box->n_compounds; c_index++) {
		//if (minDist(&box->compounds[c_index], particle_center) < MIN_NONBONDED_DIST)

		// This no longer works, as box doesn't store compound state arrays!
		/*if (minDist(box->compound_state_array[c_index], particle_center) < MIN_NONBONDED_DIST)
			return false;*/
	}

	// THis also no longer works
	/*for (int si = 0; si < box->n_solvents; si++) {
		float dist = EngineUtils::calcHyperDist(&box->solvents[si].pos, &particle_center);
		if (dist < MIN_NONBONDED_DIST) {
			printf("\tWARNING: Skipping particle with dist %f\n", dist);
			return false;
		}
	}*/

	return true;
}

//bool BoxBuilder::verifyPairwiseParticleMindist(Compound* a, Compound* b)
//{
//	for (int ia = 0; ia < a->n_particles; ia++) {
//		for (int ib = 0; ib < b->n_particles; ib++) {
//			//Float3 pos_a = a->particles[ia].pos_tsub1;
//			//Float3 pos_b = b->particles[ib].pos_tsub1;
//			Float3 pos_a = a->prev_positions[ia];
//			Float3 pos_b = b->prev_positions[ib];
//
//			float dist = (pos_a - pos_b).len();
//			if (dist < MIN_NONBONDED_DIST)
//				return false;
//		}
//	}
//	return true;
//}

