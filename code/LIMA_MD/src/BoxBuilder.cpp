#include "BoxBuilder.cuh"
#include "Printer.h"
#include "PhysicsUtils.cuh"
#include "LimaPositionSystem.cuh"

#include <random>
#include <format>
#include <numeric>

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

int SolvateBox(Box& box, const ForcefieldTinymol& forcefield, const SimParams& simparams, const std::vector<TinyMolFactory>& tinyMols)	// Accepts the position of the center or Oxygen of a solvate molecule. No checks are made wh
{
	for (const auto& tinyMol : tinyMols) {
		if (box.boxparams.nTinymolParticles + tinyMol.nParticles >= MAX_SOLVENTS) {
			throw std::runtime_error("Solvents surpass MAX_SOLVENT");
		}

		auto [nodeIndexOfTinymol, _] = LIMAPOSITIONSYSTEM::absolutePositionPlacement(tinyMol.positions[0], static_cast<float>(box.boxparams.boxSize), simparams.bc_select);
		SolventBlock& solventBlock = SolventBlocksCircularQueue::GetBlockRef(box.solventblockgrid_circularqueue, nodeIndexOfTinymol, 0, box.boxparams.boxSize);
		
		std::vector<Coord> relPos(tinyMol.nParticles);
		std::vector<uint32_t> ids(tinyMol.nParticles);
		std::vector<uint8_t> atomtypeIds(tinyMol.nParticles);
		std::vector<TinyMolParticleState> states(tinyMol.nParticles);
		for (int i = 0; i < tinyMol.nParticles; i++) {
			Float3 hyperPos = tinyMol.positions[i];
			BoundaryConditionPublic::applyHyperposNM(tinyMol.positions[0], hyperPos, static_cast<float>(box.boxparams.boxSize), PBC);
			//auto relposFloat = hyperPos - nodeIndexOfTinymol.toFloat3();
			relPos[i] = LIMAPOSITIONSYSTEM::getRelativeCoord(hyperPos, nodeIndexOfTinymol, 1, box.boxparams.boxSize, PBC);
			//relPos[i] = Coord{ hyperPos - nodeIndexOfTinymol.toFloat3()};
			ids[i] = box.boxparams.nTinymolParticles + i; // TODO: THese should've been made in compoundbuilder
			atomtypeIds[i] = tinyMol.states[i].tinymolTypeIndex;
			states[i] = tinyMol.states[i];
		}

		solventBlock.addSolvent(relPos, ids, atomtypeIds, tinyMol.bondgroup, states);
		box.boxparams.nTinymolParticles += tinyMol.nParticles;
		box.boxparams.nTinymols++;
	}

	std::mt19937 gen(1238971);
	std::uniform_real_distribution<float> distribution(-1.f, 1.f); // TODO: GROMACS COMPARISON: This is why we dont match gromacs in RMSD

	// Setup forces and vel's for VVS
	box.tinyMolParticlesState.resize(0);
	box.tinyMolParticlesState.reserve(box.boxparams.nTinymols);
	for (int i = 0; i < box.boxparams.nTinymols; i++) {
		
		// Give a random velocity. This seems.. odd, but accoring to chatGPT this is what GROMACS does

		const float moleculeMass = std::accumulate(tinyMols[i].states.begin(), tinyMols[i].states.begin() + tinyMols[i].nParticles, 0.f, 
			[&forcefield](float sum, const TinyMolParticleState& state) {return sum + forcefield.types[state.tinymolTypeIndex].mass; }
		);
		const Float3 direction = Float3{ distribution(gen), distribution(gen), distribution(gen) }.norm();
		const float velocity = PhysicsUtils::tempToVelocity(DEFAULT_TINYMOL_START_TEMPERATURE, moleculeMass);

		for (int j = 0; j < tinyMols[i].nParticles; j++) {
			box.tinyMolParticlesState.emplace_back() = tinyMols[i].states[j];
			box.tinyMolParticlesState.back().vel_prev = direction * velocity;
		}

		

		//box.tinyMols.emplace_back(TinyMolParticleState{ direction * velocity, Float3{}, tinyMols[i].state.tinymolTypeIndex });
	}    
	box.boxparams.total_particles += box.boxparams.nTinymolParticles;
	return box.boxparams.nTinymolParticles;
}


// ---------------------------------------------------------------- Public Functions ---------------------------------------------------------------- //

std::unique_ptr<Box> BoxBuilder::BuildBox(const SimParams& simparams, BoxImage& boxImage) {
	auto box = std::make_unique<Box>(boxImage.grofile.box_size);

	box->compounds.reserve(boxImage.compounds.size());
	box->compoundInterimStates.reserve(boxImage.compounds.size());
	box->compoundCoordsBuffer.reserve(boxImage.compounds.size());
	for (const CompoundFactory& compound : boxImage.compounds) {
		InsertCompoundInBox(compound, *box, simparams);
	}	

	box->boxparams.total_compound_particles = boxImage.total_compound_particles;
	box->boxparams.total_particles += boxImage.total_compound_particles;


	box->bpLutCollection = std::move(boxImage.bpLutCollection);

	box->bondgroups = boxImage.bondgroups;// Honestly maybe have these as smart ptrs to avoid copy?

#ifdef ENABLE_SOLVENTS
	SolvateBox(*box, boxImage.tinymolTypes, simparams, boxImage.solvent_positions);
#endif

	const int compoundparticles_upperbound = box->boxparams.n_compounds * MAX_COMPOUND_PARTICLES;
	box->boxparams.total_particles_upperbound = compoundparticles_upperbound + box->boxparams.nTinymolParticles; // Compounds often read/write uncompressed, while tinymols always read/write compressed

	// Ndof = 3*nParticles - nConstraints - nCOM : https://manual.gromacs.org/current/reference-manual/algorithms/molecular-dynamics.html eq:24
	box->boxparams.degreesOfFreedom = box->boxparams.total_particles * 3 - 0 - 3;

	return box;
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
		SolventBlock* src_t0 = SolventBlocksCircularQueue::getBlockPtr(simulation.box_host->solventblockgrid_circularqueue.data(), simulation.box_host->boxparams.boxSize, 0, boxsrc_current_step);
		memcpy(solvents_t0.data(), src_t0, solventBlocksGridBytesize);

		// Clear all of the data
		//delete simulation.box_host->solventblockgrid_circularqueue;
		simulation.box_host->solventblockgrid_circularqueue = SolventBlocksCircularQueue::createQueue(simulation.box_host->boxparams.boxSize);


		// Copy the temporary storage back into the queue
		SolventBlock* dest_t0 = SolventBlocksCircularQueue::getBlockPtr(simulation.box_host->solventblockgrid_circularqueue.data(), simulation.box_host->boxparams.boxSize, 0, 0);
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









