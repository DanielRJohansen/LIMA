
#include "Engine.cuh"
#include "ForceComputations.cuh"
#include "Utilities.h"
#include "KernelWarnings.cuh"
#include "EngineUtils.cuh"

#include "SimulationDevice.cuh"
#include "BoundaryCondition.cuh"
#include "SupernaturalForces.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/memcpy_async.h>
//#include <cuda/pipeline>

#pragma warning(push)
#pragma warning(disable: E0020)
#pragma warning(push)
#pragma warning(disable: 20054)

#pragma diag_suppress 20054

// ----------------------------------------------------------------------------------- FILE-SPECIFIC FORCEFIELD -------------------------------------------------------------------------------------------//

// TODO?: Pre-calculate a solvent-X paired forcefield, to save ALOT of calc in kernels
__constant__ ForceField_NB forcefield_device;

void Engine::setDeviceConstantMemory() {
	const int forcefield_bytes = sizeof(ForceField_NB);
	cudaMemcpyToSymbol(forcefield_device, &simulation->forcefield->getNBForcefieldRef(), sizeof(ForceField_NB), 0, cudaMemcpyHostToDevice);	// So there should not be a & before the device __constant__
	cudaDeviceSynchronize();
	LIMA_UTILS::genericErrorCheck("Error while moving forcefield to device\n");
}

// __constant__ mem version
__device__ inline float calcSigma(uint8_t atomtype1, uint8_t atomtype2) {
	return (forcefield_device.particle_parameters[atomtype1].sigma + forcefield_device.particle_parameters[atomtype2].sigma) * 0.5f;
}
// __shared__ mem version
__device__ inline float calcSigma(uint8_t atomtype1, uint8_t atomtype2, const ForceField_NB& forcefield) {
	return (forcefield.particle_parameters[atomtype1].sigma + forcefield.particle_parameters[atomtype2].sigma) * 0.5f;
}

__device__ inline float calcEpsilon(uint8_t atomtype1, uint8_t atomtype2) {
	return sqrtf(forcefield_device.particle_parameters[atomtype1].epsilon * forcefield_device.particle_parameters[atomtype2].epsilon);
}
__device__ inline float calcEpsilon(uint8_t atomtype1, uint8_t atomtype2, const ForceField_NB& forcefield) {
	//return 1.4f;
	return __fsqrt_rn(forcefield.particle_parameters[atomtype1].epsilon * forcefield.particle_parameters[atomtype2].epsilon);
}

// ------------------------------------------------------------------------------------------- DEVICE FUNCTIONS -------------------------------------------------------------------------------------------//

// Must be called with ALL threads in a block
__device__ Coord getRandomCoord(int lcg_seed) {
	Coord randcoord;
	for (int i = 0; i < blockDim.x; i++) {
		if (threadIdx.x == i) {
			randcoord = Coord{
				EngineUtils::genPseudoRandomNum(lcg_seed),
				EngineUtils::genPseudoRandomNum(lcg_seed),
				EngineUtils::genPseudoRandomNum(lcg_seed)
			};
		}
		__syncthreads();
	}
	return randcoord;
}

// Two variants of this exists, with and without lut
__device__ Float3 computeCompoundCompoundLJForces(const Float3& self_pos, uint8_t atomtype_self, float& potE_sum,
	const Float3* const neighbor_positions, int neighbor_n_particles, const uint8_t* const atom_types, 
	const BondedParticlesLUT* const bonded_particles_lut, LimaForcecalc::CalcLJOrigin ljorigin, const ForceField_NB& forcefield)
{
	Float3 force(0.f);

	for (int neighborparticle_id = 0; neighborparticle_id < neighbor_n_particles; neighborparticle_id++) {

		// If thread's assoc. particle is bonded to the particle in neighborcompound, continue
		if (bonded_particles_lut->get(threadIdx.x, neighborparticle_id)) { continue; }

		const int neighborparticle_atomtype = atom_types[neighborparticle_id];

		const Float3 diff = (neighbor_positions[neighborparticle_id] - self_pos);
		const float dist_sq_reciprocal = 1.f / diff.lenSquared();
		if (EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) { continue; }

		force += LimaForcecalc::calcLJForceOptim(diff, dist_sq_reciprocal, potE_sum,
			calcSigma(atomtype_self, neighborparticle_atomtype, forcefield), calcEpsilon(atomtype_self, neighborparticle_atomtype, forcefield),
			ljorigin,
			threadIdx.x,neighborparticle_id
			//global_id_self, neighbor_compound->particle_global_ids[neighborparticle_id]
		);
	}
	return force * 24.f;
}

__device__ Float3 computeCompoundCompoundLJForces(const Float3& self_pos, uint8_t atomtype_self, float& potE_sum,
	const Float3* const neighbor_positions, int neighbor_n_particles, const uint8_t* const atom_types, const ForceField_NB& forcefield) 
{
	Float3 force(0.f);

	for (int neighborparticle_id = 0; neighborparticle_id < neighbor_n_particles; neighborparticle_id++) {
		const int neighborparticle_atomtype = atom_types[neighborparticle_id];

		const Float3 diff = (neighbor_positions[neighborparticle_id] - self_pos);
		const float dist_sq_reciprocal = 1.f / diff.lenSquared();
		if (EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) { continue; }

		force += LimaForcecalc::calcLJForceOptim(diff, dist_sq_reciprocal, potE_sum,
			calcSigma(atomtype_self, neighborparticle_atomtype, forcefield), calcEpsilon(atomtype_self, neighborparticle_atomtype, forcefield),
			LimaForcecalc::CalcLJOrigin::ComComInter
			//global_id_self, neighbor_compound->particle_global_ids[neighborparticle_id]
			);
	}
	return force * 24.f;
}


__device__ Float3 computeSolventToSolventLJForces(const Float3& relpos_self, const Float3* const relpos_others, int n_elements, const bool exclude_own_index, float& potE_sum) {	// Specific to solvent kernel
	Float3 force{};

	for (int i = 0; i < n_elements; i++) {

		// If computing within block, dont compute force against thread's solvent
		if (exclude_own_index && threadIdx.x == i) { continue; }

		const Float3 diff = (relpos_others[i] - relpos_self);
		const float dist_sq_reciprocal = 1.f / diff.lenSquared();
		if (EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) { continue; }

		force += LimaForcecalc::calcLJForceOptim(diff, dist_sq_reciprocal, potE_sum,
			forcefield_device.particle_parameters[ATOMTYPE_SOLVENT].sigma,
			forcefield_device.particle_parameters[ATOMTYPE_SOLVENT].epsilon,
			exclude_own_index ? LimaForcecalc::CalcLJOrigin::SolSolIntra : LimaForcecalc::CalcLJOrigin::SolSolInter,
			threadIdx.x, i
		);
	}
	return force * 24.f;
}

__device__ Float3 computeSolventToCompoundLJForces(const Float3& self_pos, const int n_particles, const Float3* const positions, float& potE_sum, const uint8_t atomtype_self,
const ForceField_NB& forcefield) {	// Specific to solvent kernel
	Float3 force{};
	for (int i = 0; i < n_particles; i++) {

		const Float3 diff = (positions[i] - self_pos);
		const float dist_sq_reciprocal = 1.f / diff.lenSquared();
		if (EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) { continue; }

		force += LimaForcecalc::calcLJForceOptim(diff, dist_sq_reciprocal, potE_sum,
			calcSigma(atomtype_self, ATOMTYPE_SOLVENT, forcefield),
			calcEpsilon(atomtype_self, ATOMTYPE_SOLVENT, forcefield),
			LimaForcecalc::CalcLJOrigin::SolCom,
			atomtype_self, ATOMTYPE_SOLVENT
		);
	}
	return force * 24.f;
}

__device__ Float3 computeCompoundToSolventLJForces(const Float3& self_pos, const int n_particles, const Float3* const positions,
	float& potE_sum, const uint8_t* atomtypes_others, const int sol_id) 
{
	Float3 force(0.f);
	for (int i = 0; i < n_particles; i++) {

		const Float3 diff = (positions[i] - self_pos);
		const float dist_sq_reciprocal = 1.f / diff.lenSquared();
		if (EngineUtils::isOutsideCutoff(dist_sq_reciprocal)) { continue; }

		force += LimaForcecalc::calcLJForceOptim(diff, dist_sq_reciprocal, potE_sum,
			calcSigma(ATOMTYPE_SOLVENT, atomtypes_others[i]),
			calcEpsilon(ATOMTYPE_SOLVENT, atomtypes_others[i]),
			LimaForcecalc::CalcLJOrigin::ComSol,
			sol_id, -1
		);
	}
	return force * 24.f;
}





template <typename BoundaryCondition>
__device__ void getCompoundHyperpositionsAsFloat3(const NodeIndex& origo_self, const CompoundCoords* const querycompound, 
	void* output_buffer, Float3& utility_float3, const int n_particles) 
{
	if (threadIdx.x == 0) {
		const NodeIndex querycompound_hyperorigo = LIMAPOSITIONSYSTEM::getHyperNodeIndex<BoundaryCondition>(origo_self, querycompound->origo);
		KernelHelpersWarnings::assertHyperorigoIsValid(querycompound_hyperorigo, origo_self);

		// calc Relative LimaPosition Shift from the origo-shift
		utility_float3 = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, origo_self).toFloat3();
	}
//	__syncthreads();

	auto block = cooperative_groups::this_thread_block();
	cooperative_groups::memcpy_async(block, (Coord*)output_buffer, querycompound->rel_positions, sizeof(Coord) * n_particles);
	cooperative_groups::wait(block);
	__syncthreads();

	// Eventually i could make it so i only copy the active particles in the compound
	if (threadIdx.x < n_particles) {
		const Coord queryparticle_coord = ((Coord*)output_buffer)[threadIdx.x];
		((Float3*)output_buffer)[threadIdx.x] = queryparticle_coord.toFloat3() + utility_float3;
	}
	__syncthreads();
}

__device__ void getCompoundHyperpositionsAsFloat3Async(const CompoundCoords* const querycompound,
	void* output_buffer, const int n_particles, const Float3& relshift)
{
	auto block = cooperative_groups::this_thread_block();
	cooperative_groups::memcpy_async(block, (Coord*)output_buffer, querycompound->rel_positions, sizeof(Coord) * n_particles);
	cooperative_groups::wait(block);

	// Eventually i could make it so i only copy the active particles in the compound
	if (threadIdx.x < n_particles) {
		const Coord queryparticle_coord = ((Coord*)output_buffer)[threadIdx.x];
		((Float3*)output_buffer)[threadIdx.x] = queryparticle_coord.toFloat3() + relshift;
	}
}

/// <summary>
/// </summary>
/// <param name="onehot_remainers">Must be size of MAX_SOLVENTS_IN_BLOCK</param>
/// <param name="utility">Must be large enough to store temporary sums for blelloch sum</param>
/// <returns></returns>
__device__ void doBlellochPrefixSum(uint8_t* onehot_remainers, uint8_t* utility) {
	// Forward scan
	//for (int leap = 1; leap < MAX_SOLVENTS_IN_BLOCK - 1; leap *= 2) {
	//	if (threadIdx.x % leap == 0) {
	//		int index = threadIdx.x + leap;
	//	}
	//}
}

// SLOW - Returns sum of actives before, thus must be -1 for 0-based index :)
__device__ void doSequentialPrefixSum(uint8_t* onehot_remainers, int n_elements) {
	for (int i = 1; i < n_elements; i++) {
		if (threadIdx.x == i) {
			onehot_remainers[i] += onehot_remainers[i - 1];
			KernelHelpersWarnings::verifyOnehotRemaindersIsValid(onehot_remainers, i);
		}
		__syncthreads();
	}
}

__device__ uint8_t computePrefixSum(const bool remain, uint8_t* utility_buffer, int n_elements) {
	utility_buffer[threadIdx.x] = static_cast<uint8_t>(remain);
	__syncthreads();

	doSequentialPrefixSum(utility_buffer, n_elements);	
	//doBlellochPrefixSum

	const uint8_t solventindex_new = utility_buffer[threadIdx.x] - 1; // Underflow here doesn't matter, as the underflowing threads wont remain anyways :)
	return solventindex_new;	
}







/// <summary></summary>
/// <param name="solventblock">In shared memory</param>
/// <param name="transferqueues">In shared memory</param>
/// <param name="relpos_next">Register</param>
/// <param name="transfermodules">Global memory</param>
template <typename BoundaryCondition>
__device__ void transferOut(const NodeIndex& transfer_dir, const SolventBlock& solventblock_current_local, const int new_blockid, 
	const Coord& relpos_next, STransferQueue* transferqueues, SolventBlockTransfermodule* transfermodules, const NodeIndex& blockId3d) {

	// Sequential insertion in shared memory
	for (int i = 0; i < solventblock_current_local.n_solvents; i++) {
		if (threadIdx.x == i && new_blockid != blockIdx.x) {	// Only handle non-remain solvents
			const int queue_index = SolventBlockTransfermodule::getQueueIndex(transfer_dir);

			KernelHelpersWarnings::transferoutVerifyQueueIndex(queue_index, transfer_dir);

			const bool insertionSuccess = transferqueues[queue_index].addElement(relpos_next, solventblock_current_local.rel_pos[threadIdx.x], solventblock_current_local.ids[threadIdx.x]);
			KernelHelpersWarnings::transferoutVerifyInsertion(insertionSuccess);
		}
		__syncthreads();
	}

	// Coaslescing copying to global memory
	for (int queue_index = 0; queue_index < 6; queue_index++) {
		const STransferQueue& queue_local = transferqueues[queue_index];

		if (threadIdx.x < queue_local.n_elements) {

			const NodeIndex transferdir_queue = LIMAPOSITIONSYSTEM::getTransferDirection(queue_local.rel_positions[0]);		// Maybe use a utility-coord a precompute by thread0, or simply hardcode...
			const int blockid_global = EngineUtils::getNewBlockId<BoundaryCondition>(transferdir_queue, blockId3d);
			KernelHelpersWarnings::assertValidBlockId(blockid_global);

			STransferQueue* queue_global = &transfermodules[blockid_global].transfer_queues[queue_index];

			queue_global->fastInsert(
				queue_local.rel_positions[threadIdx.x] - LIMAPOSITIONSYSTEM::nodeIndexToCoord(transferdir_queue),
				queue_local.ids[threadIdx.x]);

			KernelHelpersWarnings::transferOutDebug(queue_global, queue_local, transferdir_queue, queue_index);

			queue_global->n_elements = queue_local.n_elements;

			// Only set n_elements if we get here, meaning atleast 1 new element. Otherwise it will just remain 0
			if (threadIdx.x == 0) { queue_global->n_elements = queue_local.n_elements; }
		}
	}
	__syncthreads();
}


/// <summary>
/// Must be run AFTER transferOut, as it erases all information about the transferring solvents
/// </summary>
/// <param name="solventblock_current">Solventblock in shared memory</param>
/// <param name="solventblock_next">Solventblock belong to cudablock at next step in global memory</param>
/// <param name="relpos_next"></param>
/// <param name="utility_buffer">Buffer of min size MAX_SOLVENTS_IN_BLOCK, maybe more for computing prefix sum</param>
/// <param name="remain_transfermodule">Transfermodule belonging to cudablock</param>
__device__ void compressRemainers(const SolventBlock& solventblock_current_local, SolventBlock* solventblock_next_global,
	const Coord& relpos_next, uint8_t* utility_buffer, SolventBlockTransfermodule* remain_transfermodule, const bool remain) {

	// Compute prefix sum to find new index of solvent belonging to thread
	const uint8_t solventindex_new = computePrefixSum(remain, utility_buffer, solventblock_current_local.n_solvents);


	if (remain) {
		// solventindex_new is only valid for those who remain, the rest *might* have an index of -1
		solventblock_next_global->rel_pos[solventindex_new] = relpos_next;
		solventblock_next_global->ids[solventindex_new] = solventblock_current_local.ids[threadIdx.x];
	}


	const int nsolventsinblock_next = __syncthreads_count(remain);
	if (threadIdx.x == 0) {
		remain_transfermodule->n_remain = nsolventsinblock_next;
		solventblock_next_global->n_solvents = nsolventsinblock_next;	// Doesn't matter, since the transfer kernel handles this. Enabled for debugging now..
	}
}

template <typename BoundaryCondition>
__device__ void transferOutAndCompressRemainders(const SolventBlock& solventblock_current_local, SolventBlock* solventblock_next_global,
	const Coord& relpos_next, uint8_t* utility_buffer, SolventBlockTransfermodule* transfermodules_global, STransferQueue* transferqueues_local) {

	const NodeIndex blockId3d = SolventBlocksCircularQueue::get3dIndex(blockIdx.x);
	const NodeIndex transfer_dir = threadIdx.x < solventblock_current_local.n_solvents 
		? LIMAPOSITIONSYSTEM::getTransferDirection(relpos_next) 
		: NodeIndex{};

	const int new_blockid = EngineUtils::getNewBlockId<BoundaryCondition>(transfer_dir, blockId3d);
	const bool remain = (blockIdx.x == new_blockid) && threadIdx.x < solventblock_current_local.n_solvents;

	
	transferOut<BoundaryCondition>(transfer_dir, solventblock_current_local, new_blockid, relpos_next, transferqueues_local, transfermodules_global, blockId3d);

	SolventBlockTransfermodule* remain_transfermodule = &transfermodules_global[blockIdx.x];
	compressRemainers(solventblock_current_local, solventblock_next_global, relpos_next, utility_buffer, remain_transfermodule, remain);
}






// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//
template <typename BoundaryCondition>
__global__ void compoundBondsAndIntegrationKernel(SimulationDevice* sim) {
	__shared__ CompoundCompact compound;				// Mostly bond information
	__shared__ Float3 compound_positions[THREADS_PER_COMPOUNDBLOCK];
	__shared__ Float3 utility_buffer_f3[THREADS_PER_COMPOUNDBLOCK];
	__shared__ float utility_buffer_f[THREADS_PER_COMPOUNDBLOCK];
	__shared__ NodeIndex compound_origo;

	// Buffer to be cast to different datatypes. This is dangerous!
	__shared__ char utility_buffer[cbkernel_utilitybuffer_size];

	Box* box = sim->box;
	SimParams& simparams = *sim->params;
	SimSignals* signals = sim->signals;
	const Compound* const compound_global = &box->compounds[blockIdx.x];

	if (threadIdx.x == 0) {
		compound.loadMeta(&box->compounds[blockIdx.x]);
	}
	__syncthreads();
	compound.loadData(&box->compounds[blockIdx.x]);

	{
		static_assert(cbkernel_utilitybuffer_size >= sizeof(CompoundCoords), "Utilitybuffer not large enough for CompoundCoords");
		CompoundCoords* compound_coords = (CompoundCoords*)utility_buffer;

		const CompoundCoords& compoundcoords_global = *CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, sim->signals->step, blockIdx.x);
		compound_coords->loadData(compoundcoords_global);
		__syncthreads();

		if (threadIdx.x == 0)
			compound_origo = compound_coords->origo;


		compound_positions[threadIdx.x] = compound_coords->rel_positions[threadIdx.x].toFloat3();
		__syncthreads();
	}

	float potE_sum{};
	Float3 force{};

	// ------------------------------------------------------------ Intracompound Operations ------------------------------------------------------------ //
	{
		{
			__syncthreads();
			static_assert(cbkernel_utilitybuffer_size >= sizeof(SingleBond) * MAX_SINGLEBONDS_IN_COMPOUND, "Utilitybuffer not large enough for single bonds");
			SingleBond* singlebonds = (SingleBond*)utility_buffer;

			auto block = cooperative_groups::this_thread_block();
			cooperative_groups::memcpy_async(block, singlebonds, box->compounds[blockIdx.x].singlebonds, sizeof(SingleBond) * compound.n_singlebonds);
			cooperative_groups::wait(block);

			force += LimaForcecalc::computeSinglebondForces(singlebonds, compound.n_singlebonds, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum, 0);
			__syncthreads();
		}
		{
			__syncthreads();
			static_assert(cbkernel_utilitybuffer_size >= sizeof(AngleBond) * MAX_ANGLEBONDS_IN_COMPOUND, "Utilitybuffer not large enough for angle bonds");
			AngleBond* anglebonds = (AngleBond*)utility_buffer;

			auto block = cooperative_groups::this_thread_block();
			cooperative_groups::memcpy_async(block, anglebonds, box->compounds[blockIdx.x].anglebonds, sizeof(AngleBond) * compound.n_anglebonds);
			cooperative_groups::wait(block);

			force += LimaForcecalc::computeAnglebondForces(anglebonds, compound.n_anglebonds, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);
			__syncthreads();
		}
		{
			__syncthreads();
			static_assert(cbkernel_utilitybuffer_size >= sizeof(DihedralBond) * MAX_DIHEDRALBONDS_IN_COMPOUND, "Utilitybuffer not large enough for dihedrals");
			DihedralBond* dihedrals = (DihedralBond*)utility_buffer;

			auto block = cooperative_groups::this_thread_block();
			cooperative_groups::memcpy_async(block, dihedrals, box->compounds[blockIdx.x].dihedrals, sizeof(DihedralBond) * compound.n_dihedrals);
			cooperative_groups::wait(block);

			force += LimaForcecalc::computeDihedralForces(dihedrals, compound.n_dihedrals, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);
			__syncthreads();
		}
		{
			__syncthreads();
			static_assert(cbkernel_utilitybuffer_size >= sizeof(ImproperDihedralBond) * MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND, "Utilitybuffer not large enough for improper dihedrals");
			ImproperDihedralBond* impropers = (ImproperDihedralBond*)utility_buffer;

			auto block = cooperative_groups::this_thread_block();
			cooperative_groups::memcpy_async(block, impropers, box->compounds[blockIdx.x].impropers, sizeof(ImproperDihedralBond) * compound.n_improperdihedrals);
			cooperative_groups::wait(block);

			force += LimaForcecalc::computeImproperdihedralForces(impropers, compound.n_improperdihedrals, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);

			__syncthreads();
		}
	}

	// Fetch interims from other kernels
	if (threadIdx.x < compound.n_particles) {
		force += box->compounds[blockIdx.x].forces_interim[threadIdx.x];
		potE_sum += box->compounds[blockIdx.x].potE_interim[threadIdx.x];
	}
	


	// ------------------------------------------------------------ Supernatural Forces --------------------------------------------------------------- //	
	if (simparams.snf_select == HorizontalSqueeze) {
		const float mass = forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass;
		SupernaturalForces::applyHorizontalSqueeze(utility_buffer_f3, utility_buffer_f, utility_buffer, compound_positions, compound.n_particles, compound_origo, force, mass);
	}



	// -------------------------------------------------------------- Integration & PBC --------------------------------------------------------------- //	
	{
		__syncthreads();

		static_assert(clj_utilitybuffer_bytes >= sizeof(CompoundCoords), "Utilitybuffer not large enough for CompoundCoords");
		CompoundCoords* compound_coords = (CompoundCoords*)utility_buffer;
		if (threadIdx.x == 0) {
			compound_coords->origo = compound_origo;
		}
		//compound_coords->rel_positions[threadIdx.x] = Coord{ compound_positions[threadIdx.x] };	// Alternately load from global again? Will be more precise

		const CompoundCoords& compoundcoords_global = *CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, signals->step, blockIdx.x);
		compound_coords->rel_positions[threadIdx.x] = compoundcoords_global.rel_positions[threadIdx.x];
		__syncthreads();

		float speed = 0.f;
		
		if (threadIdx.x < compound.n_particles) {
			//force.print('F');
			const float mass = forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass;

			const Float3 force_prev = box->compounds[blockIdx.x].forces_prev[threadIdx.x];	// OPTIM: make ref?
			const Float3 vel_prev = box->compounds[blockIdx.x].vels_prev[threadIdx.x];
			const Float3 vel_now = EngineUtils::integrateVelocityVVS(vel_prev, force_prev, force, simparams.dt, mass);
			const Coord pos_now = EngineUtils::integratePositionVVS(compound_coords->rel_positions[threadIdx.x], vel_now, force, mass, simparams.dt);

			const Float3 vel_scaled = vel_now * signals->thermostat_scalar;

			box->compounds[blockIdx.x].forces_prev[threadIdx.x] = force;
			box->compounds[blockIdx.x].vels_prev[threadIdx.x] = vel_scaled;

			speed = vel_scaled.len();

			// Save pos locally, but only push to box as this kernel ends
			compound_coords->rel_positions[threadIdx.x] = pos_now;
		}

		__syncthreads();

		// PBC //
		{
			__shared__ Coord shift_lm;	// Use utility coord for this?
			if (threadIdx.x == 0) {
				shift_lm = LIMAPOSITIONSYSTEM::shiftOrigo(*compound_coords, compound_global->centerparticle_index);
			}
			__syncthreads();

			LIMAPOSITIONSYSTEM_HACK::shiftRelPos(*compound_coords, shift_lm);
			__syncthreads();
		}
		LIMAPOSITIONSYSTEM_HACK::applyBC<BoundaryCondition>(*compound_coords);
		__syncthreads();

		Float3 force_LJ_sol{};	// temp
		EngineUtils::LogCompoundData(compound, box, *compound_coords, &potE_sum, force, force_LJ_sol, simparams, *signals, sim->databuffers, speed);


		// Push positions for next step
		auto* coordarray_next_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, signals->step + 1, blockIdx.x);
		coordarray_next_ptr->loadData(*compound_coords);
	}
}
template __global__ void compoundBondsAndIntegrationKernel<PeriodicBoundaryCondition>(SimulationDevice* sim);
template __global__ void compoundBondsAndIntegrationKernel<NoBoundaryCondition>(SimulationDevice* sim);



#define compound_index blockIdx.x
template <typename BoundaryCondition>
__global__ void compoundLJKernel(SimulationDevice* sim) {
	__shared__ CompoundCompact compound;				// Mostly bond information
	__shared__ Float3 compound_positions[THREADS_PER_COMPOUNDBLOCK];
	__shared__ NeighborList neighborlist;		
	__shared__ Float3 utility_buffer_f3[THREADS_PER_COMPOUNDBLOCK];
	__shared__ Float3 utility_float3;
	//__shared__ int utility_int;
	// Buffer to be cast to different datatypes. This is dangerous!
	__shared__ char utility_buffer[clj_utilitybuffer_bytes];

	__shared__ ForceField_NB forcefield_shared;

	Box* box = sim->box;
	SimParams& simparams = *sim->params;
	SimSignals* signals = sim->signals;

	// Load positions
	if (threadIdx.x == 0) {
		compound.loadMeta(&box->compounds[blockIdx.x]);
		neighborlist.loadMeta(&sim->compound_neighborlists[blockIdx.x]);
	}
	__syncthreads();
	compound.loadData(&box->compounds[blockIdx.x]);
	NodeIndex compound_origo{};
	{
		static_assert(clj_utilitybuffer_bytes >= sizeof(CompoundCoords), "Utilitybuffer not large enough for CompoundCoords");
		CompoundCoords* compound_coords = (CompoundCoords*)utility_buffer;

		const CompoundCoords& compoundcoords_global = *CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, signals->step, blockIdx.x);
		compound_coords->loadData(compoundcoords_global);
		//neighborlist.loadData(&box->compound_neighborlists[blockIdx.x]);
		neighborlist.loadData(&sim->compound_neighborlists[blockIdx.x]);
		__syncthreads();

		compound_origo = compound_coords->origo;
		compound_positions[threadIdx.x] = compound_coords->rel_positions[threadIdx.x].toFloat3();
		__syncthreads();
	}

	// Load Forcefield
	if (threadIdx.x < MAX_ATOM_TYPES)
		forcefield_shared.particle_parameters[threadIdx.x] = forcefield_device.particle_parameters[threadIdx.x];

	
	float potE_sum{};
	Float3 force{};

	// Important to wipe these, or the bondkernel will add again to next step	- or is it still now?
	box->compounds[blockIdx.x].potE_interim[threadIdx.x] = float{};
	box->compounds[blockIdx.x].forces_interim[threadIdx.x] = Float3{};

	// ------------------------------------------------------------ Intracompound Operations ------------------------------------------------------------ //
	{
		__syncthreads();
		static_assert(clj_utilitybuffer_bytes >= sizeof(BondedParticlesLUT), "Utilitybuffer not large enough for BondedParticlesLUT");

		BondedParticlesLUT* bplut_global = box->bonded_particles_lut_manager->get(compound_index, compound_index);
		BondedParticlesLUT* bonded_particles_lut = (BondedParticlesLUT*)utility_buffer;

		bonded_particles_lut->load(*bplut_global);	// A lut always exists within a compound
		__syncthreads();

		if (threadIdx.x < compound.n_particles) {
			// Having this inside vs outside the context makes impact the resulting VC, but it REALLY SHOULD NOT
			force += computeCompoundCompoundLJForces(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum, compound_positions, compound.n_particles,
				compound.atom_types, bonded_particles_lut, LimaForcecalc::CalcLJOrigin::ComComIntra, forcefield_shared);
		}

		__syncthreads();
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //

	// --------------------------------------------------------------- Intercompound forces --------------------------------------------------------------- //
	{
		const int batchsize = 32;
		__shared__ Float3 relshifts[batchsize];
		__shared__ int neighbor_n_particles[batchsize];
		__shared__ CompoundCoords* coords_ptrs[batchsize];

		const int utilitybuffer_reserved_size = sizeof(uint8_t) * MAX_COMPOUND_PARTICLES;
		uint8_t* atomtypes = (uint8_t*)utility_buffer;

		static_assert(clj_utilitybuffer_bytes >= sizeof(BondedParticlesLUT) + utilitybuffer_reserved_size, "Utilitybuffer not large enough for BondedParticlesLUT");
		BondedParticlesLUT* bonded_particles_lut = (BondedParticlesLUT*)(&utility_buffer[utilitybuffer_reserved_size]);

		// This part is scary, but it also takes up by far the majority of compute time. We use the utilitybuffer twice simultaneously, so be careful when making changes
		__syncthreads();
		int batch_index = batchsize;
		for (int i = 0; i < neighborlist.n_compound_neighbors; i++)
		{
			__syncthreads();
			// First check if we need to load a new batch of relshifts & n_particles for the coming 32 compounds
			if (batch_index == batchsize) {
				if (threadIdx.x < batchsize && threadIdx.x + i < neighborlist.n_compound_neighbors) {
					const int query_compound_id = neighborlist.neighborcompound_ids[i + threadIdx.x];
					const CompoundCoords* const querycompound = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, signals->step, query_compound_id);
					const NodeIndex querycompound_hyperorigo = LIMAPOSITIONSYSTEM::getHyperNodeIndex<BoundaryCondition>(compound_origo, querycompound->origo);
					KernelHelpersWarnings::assertHyperorigoIsValid(querycompound_hyperorigo, compound_origo);

					// calc Relative LimaPosition Shift from the origo-shift
					relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, compound_origo).toFloat3();

					neighbor_n_particles[threadIdx.x] = box->compounds[query_compound_id].n_particles;

					coords_ptrs[threadIdx.x] = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, signals->step, query_compound_id);
				}
				batch_index = 0;
				__syncthreads();
			}

			const uint16_t neighborcompound_id = neighborlist.neighborcompound_ids[i];
			const int n_particles_neighbor = neighbor_n_particles[batch_index];

			// Load the current neighbors atomtypes and positions
			if (threadIdx.x < n_particles_neighbor) {
				atomtypes[threadIdx.x] = box->compounds[neighborcompound_id].atom_types[threadIdx.x];
			}
			getCompoundHyperpositionsAsFloat3Async(coords_ptrs[batch_index], utility_buffer_f3, n_particles_neighbor, relshifts[batch_index]);
			__syncthreads();


			// The bonded compounds always comes first in the list
			if (i < compound.n_bonded_compounds)
			{
				BondedParticlesLUT* compoundpair_lut_global = box->bonded_particles_lut_manager->get(compound_index, neighborcompound_id);
				bonded_particles_lut->load(*compoundpair_lut_global);
				__syncthreads();

				if (threadIdx.x < compound.n_particles) {
					force += computeCompoundCompoundLJForces(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum,
						utility_buffer_f3, n_particles_neighbor, atomtypes, bonded_particles_lut, LimaForcecalc::CalcLJOrigin::ComComInter, forcefield_shared);
				}
			}
			else {
				if (threadIdx.x < compound.n_particles) {
					force += computeCompoundCompoundLJForces(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum,
						utility_buffer_f3, n_particles_neighbor, atomtypes, forcefield_shared);
				}
			}
			batch_index++;
		}
	}
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //

	// --------------------------------------------------------------- Solvation forces --------------------------------------------------------------- //
#ifdef ENABLE_SOLVENTS
	for (int i = 0; i < neighborlist.n_gridnodes; i++) {
		const int solventblock_id = neighborlist.gridnode_ids[i];
		const NodeIndex solventblock_hyperorigo = LIMAPOSITIONSYSTEM::getHyperNodeIndex<BoundaryCondition>(compound_origo, SolventBlocksCircularQueue::get3dIndex(solventblock_id));

		const Float3 relpos_shift = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(solventblock_hyperorigo, compound_origo).toFloat3();	// TODO: Only t0 needs to do this

		const SolventBlock* solventblock = box->solventblockgrid_circularqueue->getBlockPtr(solventblock_id, signals->step);
		const int nsolvents_neighbor = solventblock->n_solvents;

		

		// There are many more solvents in a block, than threads in this kernel, so we need to loop in strides of blockdim.x
		__syncthreads();	// Dont load buffer before all are finished with the previous iteration
		for (uint32_t offset = 0; offset < nsolvents_neighbor; offset += blockDim.x) {
			const uint32_t solvent_index = offset + threadIdx.x;
			const int n_elements_this_stride = CPPD::min(nsolvents_neighbor - offset, blockDim.x);

			// Load the positions and add rel shift
			if (solvent_index < nsolvents_neighbor) {
				utility_buffer_f3[threadIdx.x] = solventblock->rel_pos[solvent_index].toFloat3() + relpos_shift;
			}
			__syncthreads();

			if (threadIdx.x < compound.n_particles) {
				force += computeSolventToCompoundLJForces(compound_positions[threadIdx.x], n_elements_this_stride, utility_buffer_f3, potE_sum, compound.atom_types[threadIdx.x], forcefield_shared);
			}
			__syncthreads();
		}
	}
#endif
	// ------------------------------------------------------------------------------------------------------------------------------------------------ //




	// This is the first kernel, so we overwrite
	if (threadIdx.x < compound.n_particles) {
		sim->box->compounds[blockIdx.x].potE_interim[threadIdx.x] = potE_sum;
		sim->box->compounds[blockIdx.x].forces_interim[threadIdx.x] = force;
	}
}
template  __global__ void compoundLJKernel<PeriodicBoundaryCondition>(SimulationDevice* sim);
template __global__ void compoundLJKernel<NoBoundaryCondition>(SimulationDevice* sim);
#undef compound_index



#define solvent_active (threadIdx.x < solventblock.n_solvents)
#define solvent_mass (forcefield_device.particle_parameters[ATOMTYPE_SOLVENT].mass)
static_assert(SolventBlock::MAX_SOLVENTS_IN_BLOCK >= MAX_COMPOUND_PARTICLES, "solventForceKernel was about to reserve an insufficient amount of memory");
template <typename BoundaryCondition>
__global__ void solventForceKernel(SimulationDevice* sim) {
	__shared__ Float3 utility_buffer[SolventBlock::MAX_SOLVENTS_IN_BLOCK];
	__shared__ uint8_t utility_buffer_small[SolventBlock::MAX_SOLVENTS_IN_BLOCK];
	__shared__ SolventBlock solventblock;
	__shared__ SolventTransferqueue<SolventBlockTransfermodule::max_queue_size> transferqueues[6];		// TODO: Use template to make identical kernel, so the kernel with transfer is slower and larger, and the rest remain fast!!!!
	__shared__ int utility_int;
	__shared__ Coord utility_coord;
	__shared__ Float3 utility_float3;

	// Doubles as block_index_3d!
	const NodeIndex block_origo = SolventBlocksCircularQueue::get3dIndex(blockIdx.x);

	Box* box = sim->box;
	SimParams& simparams = *sim->params;
	SimSignals* signals = sim->signals;
	SolventBlock* solventblock_ptr = box->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, signals->step);

	// Init queue, otherwise it will contain wierd values // TODO: only do this on transferstep?
	if (threadIdx.x < 6) { 
		transferqueues[threadIdx.x] = SolventTransferqueue<SolventBlockTransfermodule::max_queue_size>{};
	}




	// temp
	utility_buffer[threadIdx.x] = Float3{0};
	utility_buffer_small[threadIdx.x] = 0;


	if (threadIdx.x == 0) {
		solventblock.loadMeta(*solventblock_ptr);
	}
	__syncthreads();
	solventblock.loadData(*solventblock_ptr);
	__syncthreads();	


	Float3 force{};
	float potE_sum{};
	const Float3 relpos_self = solventblock.rel_pos[threadIdx.x].toFloat3();

	// --------------------------------------------------------------- Molecule Interactions --------------------------------------------------------------- //	
	{
		// Thread 0 finds n nearby compounds
		const CompoundGridNode* compoundgridnode = sim->compound_grid->getBlockPtr(blockIdx.x);
		if (threadIdx.x == 0) { utility_int = compoundgridnode->n_nearby_compounds; }
		__syncthreads();



		for (int i = 0; i < utility_int; i++) {
			const int16_t neighborcompound_index = compoundgridnode->nearby_compound_ids[i];
			const Compound* neighborcompound = &box->compounds[neighborcompound_index];
			const int n_compound_particles = neighborcompound->n_particles;

			// All threads help loading the molecule
			// First load particles of neighboring compound
			const CompoundCoords* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, signals->step, neighborcompound_index);
			getCompoundHyperpositionsAsFloat3<BoundaryCondition>(solventblock.origo, coordarray_ptr, utility_buffer, utility_float3, n_compound_particles);


			// Then load atomtypes of neighboring compound
			if (threadIdx.x < n_compound_particles) {
				utility_buffer_small[threadIdx.x] = neighborcompound->atom_types[threadIdx.x];
			}
			__syncthreads();

			//  We can optimize here by loading and calculate the paired sigma and eps, jsut remember to loop threads, if there are many aomttypes.
			if (solvent_active) {
				force += computeCompoundToSolventLJForces(relpos_self, n_compound_particles, utility_buffer, potE_sum, utility_buffer_small, solventblock.ids[threadIdx.x]);
			}
			__syncthreads();
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //



	// --------------------------------------------------------------- Intrablock Solvent Interactions ----------------------------------------------------- //
	{
		__syncthreads(); // Sync since use of utility
		if (solvent_active) {
			utility_buffer[threadIdx.x] = relpos_self;
		}
		__syncthreads();
		if (solvent_active) {
			force += computeSolventToSolventLJForces(relpos_self, utility_buffer, solventblock.n_solvents, true, potE_sum);
		}
		__syncthreads(); // Sync since use of utility
	}	
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //

	// --------------------------------------------------------------- Interblock Solvent Interactions ----------------------------------------------------- //
	const int query_range = 2;
	for (int x = -query_range; x <= query_range; x++) {
		for (int y = -query_range; y <= query_range; y++) {
			for (int z = -query_range; z <= query_range; z++) {
				const NodeIndex dir{ x,y,z };
				if (dir.sum() > 3) { continue; }
				if (dir.isZero()) { continue; }

				const int blockindex_neighbor = EngineUtils::getNewBlockId<BoundaryCondition>(dir, block_origo);
				KernelHelpersWarnings::assertValidBlockId(blockindex_neighbor);

				const SolventBlock* solventblock_neighbor = box->solventblockgrid_circularqueue->getBlockPtr(blockindex_neighbor, signals->step);
				const int nsolvents_neighbor = solventblock_neighbor->n_solvents;
				const Float3 origoshift_offset = LIMAPOSITIONSYSTEM::nodeIndexToCoord(dir).toFloat3();

				// All threads help loading the solvent, and shifting it's relative position reletive to this solventblock
				__syncthreads();
				if (threadIdx.x < nsolvents_neighbor) {
					utility_buffer[threadIdx.x] = solventblock_neighbor->rel_pos[threadIdx.x].toFloat3() + origoshift_offset;
				}
				__syncthreads();

				if (solvent_active) {
					force += computeSolventToSolventLJForces(relpos_self, utility_buffer, nsolvents_neighbor, false, potE_sum);
				}
				__syncthreads();
			}
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //


	Coord relpos_next{};
	if (solvent_active) {
		const float mass = forcefield_device.particle_parameters[ATOMTYPE_SOLVENT].mass;
		Solvent& solventdata_ref = box->solvents[solventblock.ids[threadIdx.x]];	// Solvent private data, for VVS

		const Float3 vel_now = EngineUtils::integrateVelocityVVS(solventdata_ref.vel_prev, solventdata_ref.force_prev, force, simparams.dt, mass);
		const Coord pos_now = EngineUtils::integratePositionVVS(solventblock.rel_pos[threadIdx.x], vel_now, force, mass, simparams.dt);

		solventdata_ref.vel_prev = vel_now * signals->thermostat_scalar;
		solventdata_ref.force_prev = force;

		// Save pos locally, but only push to box as this kernel ends
		relpos_next = pos_now;

		EngineUtils::LogSolventData(box, potE_sum, solventblock, solvent_active, force, vel_now, signals->step, sim->databuffers);
	}



	// Push new SolventCoord to global mem
	SolventBlock* solventblock_next_ptr = box->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, signals->step + 1);

	if (SolventBlocksCircularQueue::isTransferStep(signals->step)) {
		//transferOutAndCompressRemainders<BoundaryCondition>(solventblock, solventblock_next_ptr, relpos_next, utility_buffer_small, box->transfermodule_array, transferqueues);
		transferOutAndCompressRemainders<BoundaryCondition>(solventblock, solventblock_next_ptr, relpos_next, utility_buffer_small, sim->transfermodule_array, transferqueues);
	}
	else {
		solventblock_next_ptr->rel_pos[threadIdx.x] = relpos_next;
		solventblock_next_ptr->ids[threadIdx.x] = solventblock.ids[threadIdx.x];
		if (threadIdx.x == 0) {
			solventblock_next_ptr->n_solvents = solventblock.n_solvents;
		}
	}
}
template __global__ void solventForceKernel<PeriodicBoundaryCondition>(SimulationDevice* sim);
template __global__ void solventForceKernel<NoBoundaryCondition>(SimulationDevice* sim);

#undef solvent_index
#undef solvent_mass
#undef solvent_active
#undef solventblock_ptr




// This is run before step.inc(), but will always publish results to the first array in grid!
template <typename BoundaryCondition>
__global__ void solventTransferKernel(SimulationDevice* sim) {
	Box* box = sim->box;
	SimParams& simparams = *sim->params;

	SolventBlockTransfermodule* transfermodule = &sim->transfermodule_array[blockIdx.x];
	
	SolventBlock* solventblock_current = box->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, sim->signals->step);
	SolventBlock* solventblock_next = box->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, sim->signals->step + 1);

	SolventTransferWarnings::assertSolventsEqualNRemain(*solventblock_next, *transfermodule);

	// Handling incoming transferring solvents
	int n_solvents_next = transfermodule->n_remain;
	for (int queue_index = 0; queue_index < SolventBlockTransfermodule::n_queues; queue_index++) {
		auto* queue = &transfermodule->transfer_queues[queue_index];
		if (threadIdx.x < queue->n_elements) {
			const int incoming_index = n_solvents_next + threadIdx.x;

			solventblock_next->rel_pos[incoming_index] = queue->rel_positions[threadIdx.x];
			solventblock_next->ids[incoming_index] = queue->ids[threadIdx.x];
		}
		n_solvents_next += queue->n_elements;

		// Signal that all elements of the queues have been moved
		__syncthreads();
		if (threadIdx.x == 0) {
			queue->n_elements = 0;
		}
	}

	SolventTransferWarnings::assertMaxPlacedSolventsIsWithinLimits(n_solvents_next, sim->signals->critical_error_encountered);

	// Finally update the solventblock_next with how many solvents it now contains
	if (threadIdx.x == 0) {
		solventblock_next->n_solvents = n_solvents_next;
	}
}
template __global__ void solventTransferKernel<PeriodicBoundaryCondition>(SimulationDevice* sim);
template __global__ void solventTransferKernel<NoBoundaryCondition>(SimulationDevice* sim);



#define particle_id_bridge threadIdx.x
template <typename BoundaryCondition>
__global__ void compoundBridgeKernel(SimulationDevice* sim) {
	__shared__ CompoundBridge bridge;
	__shared__ Float3 positions[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Float3 utility_buffer[MAX_PARTICLES_IN_BRIDGE];
	__shared__ float utility_buffer_f[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Coord utility_coord[MAX_COMPOUNDS_IN_BRIDGE];

	SimParams& simparams = *sim->params;
	Box* box = sim->box;

	if (threadIdx.x == 0) {
		bridge.loadMeta(&box->bridge_bundle->compound_bridges[blockIdx.x]);

		
	}
	__syncthreads();

	// TODO: we dont need to do this for the first compound, as it will always be 0,0,0
	if (threadIdx.x < bridge.n_compounds) {
		// Calculate necessary shift in relative positions for right, so right share the origo with left.
		utility_coord[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelativeShiftBetweenCoordarrays<BoundaryCondition>(box->coordarray_circular_queue, sim->signals->step, bridge.compound_ids[0], bridge.compound_ids[threadIdx.x]);
	}


	bridge.loadData(&box->bridge_bundle->compound_bridges[blockIdx.x]);
	__syncthreads();

	if (particle_id_bridge < bridge.n_particles) {
		ParticleReference& p_ref = bridge.particle_refs[particle_id_bridge];

		BridgeWarnings::verifyPRefValid(p_ref, bridge);

		const CompoundCoords* coordarray = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, sim->signals->step, p_ref.compound_id);

		Coord relpos = coordarray->rel_positions[p_ref.local_id_compound];
		relpos += utility_coord[p_ref.compoundid_local_to_bridge];
		positions[threadIdx.x] = relpos.toFloat3();
	}
	__syncthreads();

	float potE_sum = 0;
	Float3 force{};

	// ------------------------------------------------------------ Intercompund Operations ------------------------------------------------------------ //
	{											// So for the very first step, these ´should all be 0, but they are not??										TODO: Look into this at some point!!!! 
		force += LimaForcecalc::computeSinglebondForces(bridge.singlebonds, bridge.n_singlebonds, positions, utility_buffer, utility_buffer_f, &potE_sum, 1);
		force += LimaForcecalc::computeAnglebondForces(bridge.anglebonds, bridge.n_anglebonds, positions, utility_buffer, utility_buffer_f, &potE_sum);
		force += LimaForcecalc::computeDihedralForces(bridge.dihedrals, bridge.n_dihedrals, positions, utility_buffer, utility_buffer_f, &potE_sum);
		force += LimaForcecalc::computeImproperdihedralForces(bridge.impropers, bridge.n_improperdihedrals, positions, utility_buffer, utility_buffer_f, &potE_sum);
	}
	__syncthreads();
	// --------------------------------------------------------------------------------------------------------------------------------------------------- //

	// This is 2nd kernel so we add to the interims
	if (particle_id_bridge < bridge.n_particles) {
		ParticleReference* p_ref = &bridge.particle_refs[particle_id_bridge];
		box->compounds[p_ref->compound_id].forces_interim[p_ref->local_id_compound] += force;
		sim->box->compounds[p_ref->compound_id].potE_interim[p_ref->local_id_compound] += potE_sum;
	}
}
template __global__ void compoundBridgeKernel<PeriodicBoundaryCondition>(SimulationDevice* sim);
template __global__ void compoundBridgeKernel<NoBoundaryCondition>(SimulationDevice* sim);
#pragma warning (pop)
#pragma warning (pop)
