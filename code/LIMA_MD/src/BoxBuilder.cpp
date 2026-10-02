#include "BoxBuilder.cuh"
#include "Printer.h"
#include "PhysicsUtils.cuh"
#include "LimaPositionSystem.cuh"

#include <random>
#include <format>
#include <numeric>

using namespace LIMA_Print;


// ---------------------------------------------------------------- Public Functions ---------------------------------------------------------------- //

std::unique_ptr<Box> BoxBuilder::BuildBox(const SimParams& simparams, BoxImage& boxImage, EnvMode envmode) {
	auto box = std::make_unique<Box>(boxImage.grofile.box_size);

	// All pclusters start with zero previous forces and velocities
	box->pclusterInterimStates.resize(boxImage.persistentClusters.size());

	/*box->boxparams.total_compound_particles = boxImage.total_compound_particles;
	box->boxparams.total_particles += boxImage.total_compound_particles;*/
	box->boxparams.totalParticles = boxImage.totalParticles;

	//box->bpLutCollection = std::move(boxImage.bpLutCollection);

	// This is a bit dirty, consider having it as a shared_ptr instead.
	box->bondgroups = std::move(boxImage.bondgroups);// Honestly maybe have these as smart ptrs to avoid copy?
	boxImage.bondgroups.clear();

	// Ndof = 3*nParticles - nConstraints - nCOM : https://manual.gromacs.org/current/reference-manual/algorithms/molecular-dynamics.html eq:24
	box->boxparams.degreesOfFreedom = box->boxparams.totalParticles * 3 - 0 - 3;

	// I dont like doing this here..
	if (!simparams.enable_electrostatics) {
		for (auto& pc : boxImage.persistentClusters) {
			for (int i = 0; i < PersistentCluster::maxParticles; i++) {
				pc.pqd[i].params.charge = 0;
			}
		}
	}

	box->persistentClusters = boxImage.persistentClusters;
	box->persistentClustersMetadata = boxImage.persistentClustersMetadata;
	box->particlesBondedToParticle = std::move(boxImage.particleBondedToParticle);
	box->pclustersBondedToPcluster = std::move(boxImage.pclusterBondedToPcluster);

	// Only display uses backbones, so if no display we dont need backbone
	if (envmode == EnvMode::Full)
		box->backboneChains = InterpretBackboneChains(boxImage.grofile);
	//box->particleToCompoundOrSolventMapping = boxImage.particleToCompoundOrSolventMapping;

	return box;
}







// Do a unit-test that ensures velocities from a EM is correctly carried over to the simulation
void BoxBuilder::copyBoxState(Simulation& simulation, std::unique_ptr<Box> boxsrc, uint32_t boxsrc_current_step)
{
	if (boxsrc_current_step < 1) { throw std::runtime_error("It is not yet possible to create a new box from an old un-run box"); }

	simulation.box = std::move(boxsrc);

	// Copy current compoundcoord configuration, and put zeroes everywhere else so we can easily spot if something goes wrong
	{
		//simulation.box->compoundCoordsBuffer = boxsrc->compoundCoordsBuffer;

		//// Create temporary storage
		//std::vector<CompoundCoords> coords_t0(MAX_COMPOUNDS);
		//const size_t bytesize = sizeof(CompoundCoords) * MAX_COMPOUNDS;

		//// Copy only the current step to temporary storage
		//CompoundCoords* src_t0 = simulation.box->compoundcoordsCircularQueue->getCoordarrayRef(boxsrc_current_step, 0);
		//memcpy(coords_t0.data(), src_t0, bytesize);

		//// Clear all of the data
		//simulation.box->compoundcoordsCircularQueue->Flush();

		//// Copy the temporary storage back into the queue
		//for (int i = 0; i < 3; i++) {
		//	CompoundCoords* dest_t0 = simulation.box->compoundcoordsCircularQueue->getCoordarrayRef(i, 0);
		//	memcpy(dest_t0, coords_t0.data(), bytesize);
		//}

		// TODO ERROR: we dont copy CompoundInterimState, so it is not a true state copy
	}

}

bool BoxBuilder::verifyAllParticlesIsInsideBox(Simulation& sim, float padding, bool verbose) {
//TODO!	
	//for (int cid = 0; cid < sim.box->boxparams.n_compounds; cid++) {
	//	for (int pid = 0; pid < sim.box->compounds[cid].n_particles; pid++) 
	//	{
	//		const int index = LIMALOGSYSTEM::getMostRecentDataentryIndex(sim.getStep() - 1, sim.simParams.data_logging_interval);

	//		Float3 pos = sim.traj_buffer->getCompoundparticleDatapointAtIndex(cid, pid, index);
	//		BoundaryConditionPublic::applyBCNM(pos, sim.box->boxparams.BoxSizeFloat(), sim.simParams.bc_select);

	//		for (int i = 0; i < 3; i++) {
	//			if (pos[i] < padding || pos[i] > (sim.box->boxparams.BoxSizeFloat()[i] - padding)) {
	//				//m_logger->print(std::format("Found particle not inside the appropriate pdding of the box {}", pos.toString()));
	//				return false;
	//			}
	//		}
	//	}
	//}

	// Handle solvents somehow

	return true;
}









