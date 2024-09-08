#include "BoxBuilder.cuh"
#include "Printer.h"
#include "PhysicsUtils.cuh"
#include "LimaPositionSystem.cuh"

#include <random>
#include <format>

using namespace LIMA_Print;


std::unique_ptr<Box> BoxBuilder::BuildBox(const SimParams& simparams, BoxImage& boxImage) {
	auto box = std::make_unique<Box>(static_cast<int>(boxImage.box_size));

	for (const CompoundFactory& compound : boxImage.compounds) {
		InsertCompoundInBox(compound, *box, simparams);
	}

	box->boxparams.total_compound_particles = boxImage.total_compound_particles;						// TODO: Unknown behavior, if multiple molecules are added!
	box->boxparams.total_particles += boxImage.total_compound_particles;


	box->bridge_bundle = boxImage.bridgebundle.release();					// TODO: Breaks if multiple compounds are added, as only one bridgebundle can exist for now!
	box->boxparams.n_bridges = box->bridge_bundle->n_bridges;

	box->bonded_particles_lut_manager = boxImage.bp_lut_manager.release();

#ifdef ENABLE_SOLVENTS
	SolvateBox(*box, boxImage.forcefield, simparams, boxImage.solvent_positions);
#endif

	m_logger->print("BoxImage added to box\n");

	return box;
}




void BoxBuilder::setupDataBuffers(Simulation& simulation, const uint64_t n_steps) {
	// Permanent Outputs for energy & trajectory analysis
	const size_t n_datapoints = simulation.boxparams_host.total_particles_upperbound * n_steps / simulation.simparams_host.data_logging_interval;
	const auto datasize_str = std::to_string((float)((2. * sizeof(float) * n_datapoints + sizeof(Float3) * n_datapoints) * 1e-6));
	m_logger->print("Malloc " + datasize_str + " MB on host for data buffers\n");


	simulation.potE_buffer = std::make_unique<ParticleDataBuffer<float>>(simulation.boxparams_host.total_particles_upperbound, simulation.boxparams_host.n_compounds, n_steps, simulation.simparams_host.data_logging_interval);
	simulation.vel_buffer = std::make_unique<ParticleDataBuffer<float>>(simulation.boxparams_host.total_particles_upperbound, simulation.boxparams_host.n_compounds, n_steps, simulation.simparams_host.data_logging_interval);
	simulation.forceBuffer= std::make_unique<ParticleDataBuffer<Float3>>(simulation.boxparams_host.total_particles_upperbound, simulation.boxparams_host.n_compounds, n_steps, simulation.simparams_host.data_logging_interval);

#ifndef DONTGETDATA
	simulation.traj_buffer = std::make_unique<ParticleDataBuffer<Float3>>(simulation.boxparams_host.total_particles_upperbound, simulation.boxparams_host.n_compounds, n_steps, simulation.simparams_host.data_logging_interval);
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
	//simulation->copyBoxVariables();
	simulation->boxparams_host = simulation->box_host->boxparams;
	simulation->compounds_host.resize(simulation->boxparams_host.n_compounds);
	for (int i = 0; i < simulation->boxparams_host.n_compounds; i++)
		simulation->compounds_host[i] = simulation->box_host->compounds[i];

	m_logger->print("Box contains " + std::to_string(simulation->boxparams_host.n_compounds) + " compounds, " 
		+ std::to_string(simulation->boxparams_host.n_bridges) + " bridges and " + std::to_string(simulation->boxparams_host.n_solvents) + " solvents\n");

	// Copy forcefield to sim
	//simulation->box_host->forcefield = new ForceField_NB{ simulation->forcefield};	// Copy

	// Allocate buffers. We need to allocate for atleast 1 step, otherwise the bootstrapping mechanism will fail.
	const auto n_steps = std::max(simulation->simparams_host.n_steps, uint64_t{ 1 });
	setupDataBuffers(*simulation, n_steps);
	setupTrainingdataBuffers(*simulation, n_steps);
	LIMA_UTILS::genericErrorCheck("Error during log-data mem. allocation");

	m_logger->print("Total particles upperbound: " +  std::to_string(simulation->boxparams_host.total_particles_upperbound) + "\n");
	m_logger->print("Max particles in compound: " + std::to_string(MAX_COMPOUND_PARTICLES) + "\n");

	m_logger->finishSection("Boxbuild complete");
}


int BoxBuilder::SolvateBox(Box& box, const ForceField_NB& forcefield, const SimParams& simparams, const std::vector<Float3>& solvent_positions)	// Accepts the position of the center or Oxygen of a solvate molecule. No checks are made wh
{
	for (Float3 sol_pos : solvent_positions) {
		if (box.boxparams.n_solvents == MAX_SOLVENTS) {
			throw std::runtime_error("Solvents surpass MAX_SOLVENT");
		}

		sol_pos += most_recent_offset_applied;			// So solvents are re-aligned with an offsat molecule.

		if (spaceAvailable(box, sol_pos, true) && box.boxparams.n_solvents < MAX_SOLVENTS) {						// Should i check? Is this what energy-min is for?

			const SolventCoord solventcoord = LIMAPOSITIONSYSTEM::createSolventcoordFromAbsolutePosition(
				sol_pos, static_cast<float>(box.boxparams.boxSize), simparams.bc_select);

			box.solventblockgrid_circularqueue->addSolventToGrid(solventcoord, box.boxparams.n_solvents, 0, box.boxparams.boxSize);
			box.boxparams.n_solvents++;
		}
		else {
			// TODO: I should fill this out
			throw std::runtime_error("No room for solvent");
		}
	}

	// Setup forces and vel's for VVS
	const float solvent_mass = forcefield.particle_parameters[ATOMTYPE_SOLVENT].mass;
	const float default_solvent_start_temperature = 310;	// [K]
	box.solvents = new Solvent[box.boxparams.n_solvents];
	for (int i = 0; i < box.boxparams.n_solvents; i++) {
		box.solvents[i].force_prev = Float3{0};

		// Give a random velocity
		const Float3 direction = get3RandomSigned().norm();
		const float velocity = PhysicsUtils::tempToVelocity(default_solvent_start_temperature, solvent_mass);
		box.solvents[i].vel_prev = direction * velocity;
	}


	box.boxparams.total_particles += box.boxparams.n_solvents;
	auto a = std::to_string(box.boxparams.n_solvents);
	auto b = std::to_string(solvent_positions.size());
	m_logger->print(a + " of " + b + " solvents added to box\n");
	return box.boxparams.n_solvents;
}

// Do a unit-test that ensures velocities from a EM is correctly carried over to the simulation
void BoxBuilder::copyBoxState(Simulation* simulation, std::unique_ptr<Box> boxsrc, const SimSignals& simparams_src, uint32_t boxsrc_current_step)
{
	if (boxsrc_current_step < 1) { throw std::runtime_error("It is not yet possible to create a new box from an old un-run box"); }

	simulation->box_host = std::move(boxsrc);

	// Copy current compoundcoord configuration, and put zeroes everywhere else so we can easily spot if something goes wrong
	{
		// Create temporary storage
		std::vector<CompoundCoords> coords_t0(MAX_COMPOUNDS);
		const size_t bytesize = sizeof(CompoundCoords) * MAX_COMPOUNDS;

		// Copy only the current step to temporary storage
		CompoundCoords* src_t0 = simulation->box_host->compoundcoordsCircularQueue->getCoordarrayRef(boxsrc_current_step, 0);
		memcpy(coords_t0.data(), src_t0, bytesize);

		// Clear all of the data
		simulation->box_host->compoundcoordsCircularQueue->Flush();

		// Copy the temporary storage back into the queue
		CompoundCoords* dest_t0 = simulation->box_host->compoundcoordsCircularQueue->getCoordarrayRef(0, 0);
		memcpy(dest_t0, coords_t0.data(), bytesize);
	}

	// Do the same for solvents
	{
		// Create temporary storage
		const int blocksInGrid = BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(simulation->box_host->boxparams.boxSize));
		std::vector<SolventBlock> solvents_t0(blocksInGrid);


		//TODO: This is jsut temp:
		const int solventBlocksGridBytesize = sizeof(SolventBlock) * blocksInGrid;

		// Copy only the current step to temporary storage
		SolventBlock* src_t0 = simulation->box_host->solventblockgrid_circularqueue->getBlockPtr(0, boxsrc_current_step);
		memcpy(solvents_t0.data(), src_t0, solventBlocksGridBytesize);

		// Clear all of the data
		delete simulation->box_host->solventblockgrid_circularqueue;
		simulation->box_host->solventblockgrid_circularqueue = SolventBlocksCircularQueue::createQueue(simulation->box_host->boxparams.boxSize);


		// Copy the temporary storage back into the queue
		SolventBlock* dest_t0 = simulation->box_host->solventblockgrid_circularqueue->getBlockPtr(0, 0);
		memcpy(dest_t0, solvents_t0.data(), solventBlocksGridBytesize);
	}
}

bool BoxBuilder::verifyAllParticlesIsInsideBox(Simulation& sim, float padding, bool verbose) {
	
	for (int cid = 0; cid < sim.boxparams_host.n_compounds; cid++) {
		for (int pid = 0; pid < sim.box_host->compounds[cid].n_particles; pid++) 
		{
			const int index = LIMALOGSYSTEM::getMostRecentDataentryIndex(sim.simsignals_host.step - 1, sim.simparams_host.data_logging_interval);

			Float3 pos = sim.traj_buffer->getCompoundparticleDatapointAtIndex(cid, pid, index);
			BoundaryConditionPublic::applyBCNM(pos, (float) sim.boxparams_host.boxSize, sim.simparams_host.bc_select);

			for (int i = 0; i < 3; i++) {
				if (pos[i] < padding || pos[i] > (static_cast<float>(sim.boxparams_host.boxSize) - padding)) {
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

void BoxBuilder::InsertCompoundInBox(const CompoundFactory& compound, Box& box, const SimParams& simparams, Float3 offset)
{
	std::vector<Float3> positions;
	positions.reserve(MAX_COMPOUND_PARTICLES);

	for (int i = 0; i < compound.n_particles; i++) {
		const Float3& extern_position = compound.positions[i];
		positions.push_back(extern_position);
	}

	CompoundCoords& coords_now = *box.compoundcoordsCircularQueue->getCoordarrayRef(0, box.boxparams.n_compounds);
	coords_now = LIMAPOSITIONSYSTEM::positionCompound(positions, compound.centerparticle_index, static_cast<float>(box.boxparams.boxSize), simparams.bc_select);
	if (simparams.bc_select == PBC && !coords_now.origo.isInBox(BoxGrid::NodesPerDim(box.boxparams.boxSize))) {
		throw std::runtime_error(std::format("Invalid compound origo {}", coords_now.origo.toString()));
	}

	box.compounds[box.boxparams.n_compounds++] = Compound{ compound };	// Cast and copy only the base of the factory
}
















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

