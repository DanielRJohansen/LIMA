#include "BoxBuilder.cuh"
#include "Printer.h"
#include "PhysicsUtils.cuh"
#include "LimaPositionSystem.cuh"

#include <random>
#include <format>

using namespace LIMA_Print;




// ---------------------------------------------------------------- Private Functions ---------------------------------------------------------------- //

void InsertCompoundInBox(const CompoundFactory& compound, Box& box, const SimParams& simparams, Float3 offset = Float3{})
{
	if (box.compounds.size() >= MAX_COMPOUNDS) {
		throw std::runtime_error("Compounds surpass MAX_COMPOUNDS");
	}
	std::vector<Float3> positions;
	positions.reserve(MAX_COMPOUND_PARTICLES);
	for (int i = 0; i < compound.n_particles; i++) {
		const Float3& extern_position = compound.positions[i];
		positions.push_back(extern_position);
	}

	/*CompoundCoords& coords_now = *box.compoundcoordsCircularQueue->getCoordarrayRef(0, box.boxparams.n_compounds);
	coords_now = */
	box.compoundCoordsBuffer.emplace_back(LIMAPOSITIONSYSTEM::positionCompound(positions, compound.centerparticle_index, static_cast<float>(box.boxparams.boxSize), simparams.bc_select));
	if (simparams.bc_select == PBC && !box.compoundCoordsBuffer.back().origo.isInBox(BoxGrid::NodesPerDim(box.boxparams.boxSize))) {
		throw std::runtime_error(std::format("Invalid compound origo {}", box.compoundCoordsBuffer.back().origo.toString()));
	}

	CompoundInterimState compoundState{};
	memset(&compoundState, 0, sizeof(CompoundInterimState));
	for (int i = 0; i < compound.n_particles; i++) {		
		compoundState.coords[i] = Coord(box.compoundCoordsBuffer.back().rel_positions[i]);
	}

	box.compounds.emplace_back(Compound{compound});	// Cast and copy only the base of the factory
	box.compoundInterimStates.emplace_back(compoundState);
	box.boxparams.n_compounds++;

	if (!simparams.enable_electrostatics)
		memset(box.compounds.back().atom_charges, 0, sizeof(half) * MAX_COMPOUND_PARTICLES);
}

Float3 get3Random() {	// Returns 3 numbers between 0-1
	return Float3(
		(float) (rand() % RAND_MAX / (double) RAND_MAX),
		(float) (rand() % RAND_MAX / (double) RAND_MAX),
		(float) (rand() % RAND_MAX / (double) RAND_MAX)
	);
}
Float3 get3RandomSigned() {	// returns 3 numbers between -0.5->0.5
	return get3Random() - Float3(0.5f);
}









// ---------------------------------------------------------------- Public Functions ---------------------------------------------------------------- //

std::unique_ptr<Box> BoxBuilder::BuildBox(const SimParams& simparams, BoxImage& boxImage) {
	srand(290128309);

	auto box = std::make_unique<Box>(static_cast<int>(boxImage.grofile.box_size.x));

	box->compounds.reserve(boxImage.compounds.size());
	box->compoundInterimStates.reserve(boxImage.compounds.size());
	box->compoundCoordsBuffer.reserve(boxImage.compounds.size());
	for (const CompoundFactory& compound : boxImage.compounds) {
		InsertCompoundInBox(compound, *box, simparams);
	}	

	box->boxparams.total_compound_particles = boxImage.total_compound_particles;						// TODO: Unknown behavior, if multiple molecules are added!
	box->boxparams.total_particles += boxImage.total_compound_particles;


	box->bridge_bundle = std::move(boxImage.bridgebundle);					// TODO: Breaks if multiple compounds are added, as only one bridgebundle can exist for now!
	box->boxparams.n_bridges = box->bridge_bundle->n_bridges;

	box->bpLutCollection = std::move(boxImage.bpLutCollection);

#ifdef ENABLE_SOLVENTS
	SolvateBox(*box, boxImage.forcefield, simparams, boxImage.solvent_positions);
#endif

	const int compoundparticles_upperbound = box->boxparams.n_compounds * MAX_COMPOUND_PARTICLES;
	box->boxparams.total_particles_upperbound = compoundparticles_upperbound + box->boxparams.n_solvents;

	// Ndof = 3*nParticles - nConstraints - nCOM : https://manual.gromacs.org/current/reference-manual/algorithms/molecular-dynamics.html eq:24
	box->boxparams.degreesOfFreedom = box->boxparams.total_particles * 3 - 0 - 3;

	return box;
}





int BoxBuilder::SolvateBox(Box& box, const ForceField_NB& forcefield, const SimParams& simparams, const std::vector<Float3>& solvent_positions)	// Accepts the position of the center or Oxygen of a solvate molecule. No checks are made wh
{
	for (Float3 sol_pos : solvent_positions) {
		if (box.boxparams.n_solvents == MAX_SOLVENTS) {
			throw std::runtime_error("Solvents surpass MAX_SOLVENT");
		}



		const SolventCoord solventcoord = LIMAPOSITIONSYSTEM::createSolventcoordFromAbsolutePosition(
			sol_pos, static_cast<float>(box.boxparams.boxSize), simparams.bc_select);

		box.solventblockgrid_circularqueue->addSolventToGrid(solventcoord, box.boxparams.n_solvents, 0, box.boxparams.boxSize);
		box.boxparams.n_solvents++;

	}

	// Setup forces and vel's for VVS
	const float solvent_mass = forcefield.particle_parameters[ATOMTYPE_SOLVENT].mass;
	const float default_solvent_start_temperature = 310;	// [K]
	box.solvents.reserve(box.boxparams.n_solvents);
	for (int i = 0; i < box.boxparams.n_solvents; i++) {		
		// Give a random velocity
		const Float3 direction = get3RandomSigned().norm();
		const float velocity = PhysicsUtils::tempToVelocity(default_solvent_start_temperature, solvent_mass);

		box.solvents.emplace_back(Solvent{ direction * velocity, Float3{} });
	}


	box.boxparams.total_particles += box.boxparams.n_solvents;
	return box.boxparams.n_solvents;
}

// Do a unit-test that ensures velocities from a EM is correctly carried over to the simulation
void BoxBuilder::copyBoxState(Simulation& simulation, std::unique_ptr<Box> boxsrc, uint32_t boxsrc_current_step)
{
	if (boxsrc_current_step < 1) { throw std::runtime_error("It is not yet possible to create a new box from an old un-run box"); }

	simulation.box_host = std::move(boxsrc);

	// Copy current compoundcoord configuration, and put zeroes everywhere else so we can easily spot if something goes wrong
	{
		//simulation.box_host->compoundCoordsBuffer = boxsrc->compoundCoordsBuffer;

		//// Create temporary storage
		//std::vector<CompoundCoords> coords_t0(MAX_COMPOUNDS);
		//const size_t bytesize = sizeof(CompoundCoords) * MAX_COMPOUNDS;

		//// Copy only the current step to temporary storage
		//CompoundCoords* src_t0 = simulation.box_host->compoundcoordsCircularQueue->getCoordarrayRef(boxsrc_current_step, 0);
		//memcpy(coords_t0.data(), src_t0, bytesize);

		//// Clear all of the data
		//simulation.box_host->compoundcoordsCircularQueue->Flush();

		//// Copy the temporary storage back into the queue
		//for (int i = 0; i < 3; i++) {
		//	CompoundCoords* dest_t0 = simulation.box_host->compoundcoordsCircularQueue->getCoordarrayRef(i, 0);
		//	memcpy(dest_t0, coords_t0.data(), bytesize);
		//}

		// TODO ERROR: we dont copy CompoundInterimState, so it is not a true state copy
	}

	// Do the same for solvents
	{
		// Create temporary storage
		const int blocksInGrid = BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(simulation.box_host->boxparams.boxSize));
		std::vector<SolventBlock> solvents_t0(blocksInGrid);


		//TODO: This is jsut temp:
		const int solventBlocksGridBytesize = sizeof(SolventBlock) * blocksInGrid;

		// Copy only the current step to temporary storage
		SolventBlock* src_t0 = simulation.box_host->solventblockgrid_circularqueue->getBlockPtr(0, boxsrc_current_step);
		memcpy(solvents_t0.data(), src_t0, solventBlocksGridBytesize);

		// Clear all of the data
		//delete simulation.box_host->solventblockgrid_circularqueue;
		simulation.box_host->solventblockgrid_circularqueue = SolventBlocksCircularQueue::createQueue(simulation.box_host->boxparams.boxSize);


		// Copy the temporary storage back into the queue
		SolventBlock* dest_t0 = simulation.box_host->solventblockgrid_circularqueue->getBlockPtr(0, 0);
		memcpy(dest_t0, solvents_t0.data(), solventBlocksGridBytesize);
	}
}

bool BoxBuilder::verifyAllParticlesIsInsideBox(Simulation& sim, float padding, bool verbose) {
	
	for (int cid = 0; cid < sim.box_host->boxparams.n_compounds; cid++) {
		for (int pid = 0; pid < sim.box_host->compounds[cid].n_particles; pid++) 
		{
			const int index = LIMALOGSYSTEM::getMostRecentDataentryIndex(sim.getStep() - 1, sim.simparams_host.data_logging_interval);

			Float3 pos = sim.traj_buffer->getCompoundparticleDatapointAtIndex(cid, pid, index);
			BoundaryConditionPublic::applyBCNM(pos, (float) sim.box_host->boxparams.boxSize, sim.simparams_host.bc_select);

			for (int i = 0; i < 3; i++) {
				if (pos[i] < padding || pos[i] > (static_cast<float>(sim.box_host->boxparams.boxSize) - padding)) {
					//m_logger->print(std::format("Found particle not inside the appropriate pdding of the box {}", pos.toString()));
					return false;
				}
			}
		}
	}

	// Handle solvents somehow

	return true;
}









