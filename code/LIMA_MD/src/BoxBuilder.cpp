#include "BoxBuilder.cuh"
#include "EngineUtils.cuh"
#include "Printer.h"
//#include "BoundaryCondition.cuh"	// TODO: This should be private to engine
#include <random>
#include <format>

using namespace LIMA_Print;

void BoxBuilder::buildBox(Simulation* simulation, float boxsize_nm) {
	m_logger->startSection("Building box");

	simulation->box_host->boxparams.dims = Float3{ boxsize_nm };

	simulation->box_host->compounds = new Compound[MAX_COMPOUNDS];
	simulation->box_host->coordarray_circular_queue = new CompoundCoords[Box::coordarray_circular_queue_n_elements];

	simulation->box_host->solventblockgrid_circularqueue = SolventBlocksCircularQueue::createQueue();
	//simulation->box_host->transfermodule_array = new SolventBlockTransfermodule[SolventBlocksCircularQueue::blocks_per_grid];

	simulation->box_host->bridge_bundle = new CompoundBridgeBundleCompact{};

	simulation->box_host->owns_members = true; // I guess this sorta requires all ptr's to be allocated in the same scope otherwise the destructor will fail with this param / or not get all members

	cudaDeviceSynchronize();
	if (cudaGetLastError() != cudaSuccess) {
		throw std::runtime_error("Error during buildBox()");		
	}
}


void BoxBuilder::addBoxImage(Simulation* simulation, BoxImage& compound_collection) {
	for (const CompoundFactory& compound : compound_collection.compounds) {		
		//CALL_FUNCTION_WITH_BC(insertCompoundInBox, simulation->simparams_host.constparams.bc_select, compound, *simulation);
		//insertCompoundInBox<BoundaryCondition>(compound, *simulation);
		insertCompoundInBox(compound, *simulation);
	}

	simulation->extraparams.total_compound_particles = compound_collection.total_compound_particles;						// TODO: Unknown behavior, if multiple molecules are added!
	simulation->box_host->boxparams.total_particles += compound_collection.total_compound_particles;


	simulation->box_host->bridge_bundle = compound_collection.bridgebundle.release();					// TODO: Breaks if multiple compounds are added, as only one bridgebundle can exist for now!
	simulation->box_host->boxparams.n_bridges = simulation->box_host->bridge_bundle->n_bridges;

	simulation->box_host->bonded_particles_lut_manager = compound_collection.bp_lut_manager.release();

	m_logger->print("BoxImage added to box\n");
}




void BoxBuilder::setupDataBuffers(Simulation& simulation, const uint64_t n_steps) {
	// Permanent Outputs for energy & trajectory analysis
	const size_t n_datapoints = simulation.boxparams_host.total_particles_upperbound * n_steps / LOG_EVERY_N_STEPS;
	const auto datasize_str = std::to_string((float)((2. * sizeof(float) * n_datapoints + sizeof(Float3) * n_datapoints) * 1e-6));
	m_logger->print("Malloc " + datasize_str + " MB on host for data buffers\n");


	simulation.potE_buffer = std::make_unique<ParticleDataBuffer<float>>(simulation.boxparams_host.total_particles_upperbound, simulation.boxparams_host.n_compounds, n_steps);
	simulation.vel_buffer = std::make_unique<ParticleDataBuffer<float>>(simulation.boxparams_host.total_particles_upperbound, simulation.boxparams_host.n_compounds, n_steps);

#ifndef DONTGETDATA
	simulation.traj_buffer = std::make_unique<ParticleDataBuffer<Float3>>(simulation.boxparams_host.total_particles_upperbound, simulation.boxparams_host.n_compounds, n_steps);
#endif

	simulation.temperature_buffer.reserve(n_steps / STEPS_PER_THERMOSTAT + 1);
}

void BoxBuilder::setupTrainingdataBuffers(Simulation& simulation, const uint64_t n_steps) {
#ifdef GENERATETRAINDATA
	uint64_t n_loggingdata_host = 10 * n_steps;
	uint64_t n_traindata_host = n_steps * N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation.boxparams_host.n_compounds;
	auto datasize_str = std::to_string((float)(sizeof(Float3) * n_traindata_host + sizeof(float) * n_loggingdata_host) * 1e-9);
	m_logger->print("Reserving " + datasize_str + "GB host mem for logging and training data\n");

	simulation.loggingdata.resize(n_loggingdata_host);
	simulation.trainingdata.resize(n_traindata_host);
#endif
}
void BoxBuilder::finishBox(Simulation* simulation) {
	const int compoundparticles_upperbound = simulation->box_host->boxparams.n_compounds * MAX_COMPOUND_PARTICLES;

	simulation->box_host->boxparams.total_particles_upperbound = compoundparticles_upperbound + simulation->box_host->boxparams.n_solvents;
	

	// Load meta information
	simulation->copyBoxVariables();
	m_logger->print("Box contains " + std::to_string(simulation->boxparams_host.n_compounds) + " compounds, " 
		+ std::to_string(simulation->boxparams_host.n_bridges) + " bridges and " + std::to_string(simulation->boxparams_host.n_solvents) + " solvents\n");

	// Copy forcefield to sim
	simulation->box_host->forcefield = new ForceField_NB{ simulation->forcefield->getNBForcefield()};	// Copy

	// Allocate buffers. We need to allocate for atleast 1 step, otherwise the bootstrapping mechanism will fail.
	const auto n_steps = std::max(simulation->simparams_host.n_steps, uint64_t{ 1 });
	setupDataBuffers(*simulation, n_steps);
	setupTrainingdataBuffers(*simulation, n_steps);
	LIMA_UTILS::genericErrorCheck("Error during log-data mem. allocation");

	m_logger->print("Total particles upperbound: " +  std::to_string(simulation->boxparams_host.total_particles_upperbound) + "\n");
	m_logger->print("Max particles in compound: " + std::to_string(MAX_COMPOUND_PARTICLES) + "\n");

	m_logger->finishSection("Boxbuild complete");
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
	const float solvent_mass = simulation->forcefield->getNBForcefield().particle_parameters[ATOMTYPE_SOLVENT].mass;
	const float default_solvent_start_temperature = 310;	// [K]

	for (Float3 sol_pos : solvent_positions) {
		if (simulation->box_host->boxparams.n_solvents == MAX_SOLVENTS) {
			throw std::runtime_error("Solvents surpass MAX_SOLVENT");
		}

		sol_pos += most_recent_offset_applied;			// So solvents are re-aligned with an offsat molecule.

		if (spaceAvailable(*simulation->box_host, sol_pos, true) && simulation->box_host->boxparams.n_solvents < MAX_SOLVENTS) {						// Should i check? Is this what energy-min is for?
			LimaPosition position = LIMAPOSITIONSYSTEM::createLimaPosition(sol_pos);

			const SolventCoord solventcoord = LIMAPOSITIONSYSTEM::createSolventcoordFromAbsolutePosition(
				position, simulation->box_host->boxparams.dims.x, simulation->simparams_host.bc_select);


			simulation->box_host->solventblockgrid_circularqueue->addSolventToGrid(solventcoord, simulation->box_host->boxparams.n_solvents, 0);

			//const Float3 direction = get3RandomSigned().norm();
			//const float velocity = EngineUtils::tempToVelocity(default_solvent_start_temperature, solvent_mass);	// [m/s]
			//const Float3 deltapos_lm = (direction * velocity * (simulation->simparams_host.constparams.dt));

			//SolventCoord solventcoord_prev = solventcoord;
			//solventcoord_prev.rel_position -= Coord{ deltapos_lm };

			//simulation->box_host->solventblockgrid_circularqueue->addSolventToGrid(solventcoord_prev, simulation->box_host->boxparams.n_solvents, SolventBlocksCircularQueue::first_step_prev);

			simulation->box_host->boxparams.n_solvents++;
		}
		else {
			// TODO: I should fill this out
			throw std::runtime_error("No room for solvent");
		}
	}

	// Setup forces and vel's for VVS
	simulation->box_host->solvents = new Solvent[simulation->box_host->boxparams.n_solvents];
	for (int i = 0; i < simulation->box_host->boxparams.n_solvents; i++) {
		simulation->box_host->solvents[i].force_prev = Float3{0};

		// Give a random velocity
		const Float3 direction = get3RandomSigned().norm();
		const float velocity = EngineUtils::tempToVelocity(default_solvent_start_temperature, solvent_mass);
		simulation->box_host->solvents[i].vel_prev = direction * velocity;
	}


	simulation->box_host->boxparams.total_particles += simulation->box_host->boxparams.n_solvents;
	auto a = std::to_string(simulation->box_host->boxparams.n_solvents);
	auto b = std::to_string(solvent_positions.size());
	m_logger->print(a + " of " + b + " solvents added to box\n");
	return simulation->box_host->boxparams.n_solvents;
}

// Do a unit-test that ensures velocities from a EM is correctly carried over to the simulation
void BoxBuilder::copyBoxState(Simulation* simulation, std::unique_ptr<Box> boxsrc, const SimSignals& simparams_src, uint32_t boxsrc_current_step)
{
	if (boxsrc_current_step < 1) { throw std::runtime_error("It is not yet possible to create a new box from an old un-run box"); }

	//simulation->box_host = SimUtils::copyToHost(boxsrc);
	simulation->box_host = std::move(boxsrc);

	// Copy current compoundcoord configuration, and put zeroes everywhere else so we can easily spot if something goes wrong
	{
		// Create temporary storage
		std::vector<CompoundCoords> coords_t0(MAX_COMPOUNDS);
		const size_t bytesize = sizeof(CompoundCoords) * MAX_COMPOUNDS;

		// Copy only the current step to temporary storage
		CompoundCoords* src_t0 = CoordArrayQueueHelpers::getCoordarrayRef(simulation->box_host->coordarray_circular_queue, boxsrc_current_step, 0);
		memcpy(coords_t0.data(), src_t0, bytesize);

		// Clear all of the data
		for (size_t i = 0; i < Box::coordarray_circular_queue_n_elements; i++) {
			simulation->box_host->coordarray_circular_queue[i] = CompoundCoords{};
		}

		// Copy the temporary storage back into the queue
		CompoundCoords* dest_t0 = CoordArrayQueueHelpers::getCoordarrayRef(simulation->box_host->coordarray_circular_queue, 0, 0);
		memcpy(dest_t0, coords_t0.data(), bytesize);
	}

	// Do the same for solvents
	{
		// Create temporary storage
		std::vector<SolventBlock> solvents_t0(SolventBlocksCircularQueue::blocks_per_grid);

		// Copy only the current step to temporary storage
		SolventBlock* src_t0 = simulation->box_host->solventblockgrid_circularqueue->getBlockPtr(0, boxsrc_current_step);
		memcpy(solvents_t0.data(), src_t0, SolventBlocksCircularQueue::grid_bytesize);

		// Clear all of the data
		delete simulation->box_host->solventblockgrid_circularqueue;
		simulation->box_host->solventblockgrid_circularqueue = SolventBlocksCircularQueue::createQueue();


		// Copy the temporary storage back into the queue
		SolventBlock* dest_t0 = simulation->box_host->solventblockgrid_circularqueue->getBlockPtr(0, 0);
		memcpy(dest_t0, solvents_t0.data(), SolventBlocksCircularQueue::grid_bytesize);
	}
}

bool BoxBuilder::verifyAllParticlesIsInsideBox(Simulation& sim, float padding, bool verbose) {
	
	for (int cid = 0; cid < sim.boxparams_host.n_compounds; cid++) {
		for (int pid = 0; pid < sim.box_host->compounds[cid].n_particles; pid++) 
		{
			const int index = LIMALOGSYSTEM::getMostRecentDataentryIndex(sim.simsignals_host.step - 1);

			Float3 pos = sim.traj_buffer->getCompoundparticleDatapointAtIndex(cid, pid, index);
			BoundaryConditionPublic::applyBCNM(pos, sim.boxparams_host.dims.x, sim.simparams_host.bc_select);

			for (int i = 0; i < 3; i++) {
				if (pos.at(i) < padding || pos.at(i) > (sim.boxparams_host.dims.x - padding)) {
					m_logger->print(std::format("Found particle not inside the appropriate pdding of the box {}", pos.toString()));
					return false;
				}
			}
		}
	}

	// Handle solvents somehow

	return true;
}











// ---------------------------------------------------------------- Private Functions ---------------------------------------------------------------- //

	//const float M = SOLVENT_MASS;				// kg/mol
	//const double T = 313.;	// Kelvin
	//const double R = 8.3144;					// J/(Kelvin*mol)
	//const float v_rms = static_cast<float>(sqrt(3 * R * T / M));

	//Float3 compound_united_vel = Float3(random(), random(), random()).norm() * v_rms * 0.f;			// Giving individual comp in molecule different uniform vels is sub-optimal...
void BoxBuilder::insertCompoundInBox(const CompoundFactory& compound, Simulation& simulation, Float3 offset)
{
	std::vector<LimaPosition> positions;
	positions.reserve(MAX_COMPOUND_PARTICLES);

	for (int i = 0; i < compound.n_particles; i++) {
		const Float3& extern_position = compound.positions[i];
		positions.push_back(LIMAPOSITIONSYSTEM::createLimaPosition(extern_position));
	}

	CompoundCoords& coords_now = *CoordArrayQueueHelpers::getCoordarrayRef(simulation.box_host->coordarray_circular_queue, 0, simulation.box_host->boxparams.n_compounds);
	coords_now = LIMAPOSITIONSYSTEM::positionCompound(positions, compound.centerparticle_index, simulation.box_host->boxparams.dims.x, simulation.simparams_host.bc_select);
	if (simulation.simparams_host.bc_select == PBC && !coords_now.origo.isInBox(BOXGRID_N_NODES)) {
		throw std::runtime_error(std::format("Invalid compound origo {}", coords_now.origo.toString()));
	}

	simulation.box_host->compounds[simulation.box_host->boxparams.n_compounds++] = Compound{ compound };	// Cast and copy only the base of the factory
}



















//Solvent BoxBuilder::createSolvent(Float3 com, float dt) {
//	com = com / NORMALIZER * 1e+6;	// concvert to normalized [fm]
//	Float3 solvent_vel = Float3(random(), random(), random()).norm() * v_rms * VEL_RMS_SCALAR / NORMALIZER;		// TODO: I dont know, but i think we need to freeze solvents to avoid unrealisticly large forces at step 1
//	return Solvent(com, com - solvent_vel * dt);
//}
/*
* These two funcitons are in charge of normalizing ALL coordinates!!
*/


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

//float minDist(CompoundState& compoundstate, Float3 particle_pos) {
//	float mindist = 999999;
//	for (size_t i = 0; i < compoundstate.n_particles; i++) {
//		//float dist = EngineUtils::calcHyperDist(&compound->prev_positions[i], &particle_pos);
//		float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM<PeriodicBoundaryCondition>(&compoundstate.positions[i], &particle_pos);		// Hmmm doesn't fit the namespace...
//		mindist = std::min(mindist, dist);
//	}
//	return mindist;
//}

// This funciton is currently blank! TODO: fix
bool BoxBuilder::spaceAvailable(const Box& box, Float3 particle_center, bool verbose)
{
	particle_center = particle_center;
	for (uint32_t c_index = 0; c_index < box.boxparams.n_compounds; c_index++) {
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