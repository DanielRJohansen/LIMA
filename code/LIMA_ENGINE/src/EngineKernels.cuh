
#include "Engine.cuh"
#include "ForceComputations.cuh"
#include "Utilities.h"
#include "KernelWarnings.cuh"
#include "EngineUtils.cuh"

#include "SimulationDevice.cuh"
#include "BoundaryCondition.cuh"
#include "SupernaturalForces.cuh"
#include "SolventBlockTransfers.cuh"
#include "DeviceAlgorithms.cuh"

//#include <cuda/pipeline>
#include "KernelConstants.cuh"

#pragma warning(push)
#pragma warning(disable:E0020)
#pragma warning(push)
#pragma warning(disable: 20054)

#pragma diag_suppress 20054



#include "LennardJonesInteractions.cuh"


template <typename BondType, int max_bondtype_in_compound>
__device__ BondType* LoadBonds(char* utility_buffer, BondType* source, int nBondsToLoad) {
	static_assert(cbkernel_utilitybuffer_size >= sizeof(BondType) * max_bondtype_in_compound, "Utilitybuffer not large enough for bondtype");
	BondType* bonds = (BondType*)utility_buffer;

	auto block = cooperative_groups::this_thread_block();
	cooperative_groups::memcpy_async(block, bonds, source, sizeof(BondType) * nBondsToLoad);
	cooperative_groups::wait(block);

	return bonds;
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

		const CompoundCoords& compoundcoords_global = *box->compoundcoordsCircularQueue->getCoordarrayRef(sim->signals->step, blockIdx.x);
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
		SingleBond* singlebonds = LoadBonds<SingleBond, MAX_SINGLEBONDS_IN_COMPOUND>(utility_buffer, box->compounds[blockIdx.x].singlebonds, compound.n_singlebonds);
		if (simparams.em_variant)
			force += LimaForcecalc::computeSinglebondForces<true>(singlebonds, compound.n_singlebonds, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum, 0);
		else
			force += LimaForcecalc::computeSinglebondForces<false>(singlebonds, compound.n_singlebonds, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum, 0);
			
		AngleBond* anglebonds = LoadBonds<AngleBond, MAX_ANGLEBONDS_IN_COMPOUND>(utility_buffer, box->compounds[blockIdx.x].anglebonds, compound.n_anglebonds);
		force += LimaForcecalc::computeAnglebondForces(anglebonds, compound.n_anglebonds, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);
				
		DihedralBond* dihedrals = LoadBonds<DihedralBond, MAX_DIHEDRALBONDS_IN_COMPOUND>(utility_buffer, box->compounds[blockIdx.x].dihedrals, compound.n_dihedrals);
		force += LimaForcecalc::computeDihedralForces(dihedrals, compound.n_dihedrals, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);
			
		ImproperDihedralBond* impropers = LoadBonds<ImproperDihedralBond, MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND>(utility_buffer, box->compounds[blockIdx.x].impropers, compound.n_improperdihedrals);
		force += LimaForcecalc::computeImproperdihedralForces(impropers, compound.n_improperdihedrals, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);	
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
	if (simparams.snf_select == HorizontalChargeField && threadIdx.x < compound.n_particles) {
		force += box->uniformElectricField.GetForce(box->compounds[blockIdx.x].atom_charges[threadIdx.x]);
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

		const CompoundCoords& compoundcoords_global = *box->compoundcoordsCircularQueue->getCoordarrayRef(signals->step, blockIdx.x);
		compound_coords->rel_positions[threadIdx.x] = compoundcoords_global.rel_positions[threadIdx.x];
		__syncthreads();

		float speed = 0.f;
		
		if (threadIdx.x < compound.n_particles) {
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
		auto* coordarray_next_ptr = box->compoundcoordsCircularQueue->getCoordarrayRef(signals->step + 1, blockIdx.x);
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

	// Buffer to be cast to different datatypes. UNSAFE!
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

		const CompoundCoords& compoundcoords_global = *box->compoundcoordsCircularQueue->getCoordarrayRef(signals->step, blockIdx.x);
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
			force += LJ::computeCompoundCompoundLJForces(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum, compound_positions, compound.n_particles,
				compound.atom_types, bonded_particles_lut, LJ::CalcLJOrigin::ComComIntra, forcefield_shared);
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
			//if (i > 10) break;
			__syncthreads();
			// First check if we need to load a new batch of relshifts & n_particles for the coming 32 compounds
			if (batch_index == batchsize) {
				if (threadIdx.x < batchsize && threadIdx.x + i < neighborlist.n_compound_neighbors) {
					const int query_compound_id = neighborlist.neighborcompound_ids[i + threadIdx.x];
					const CompoundCoords* const querycompound = box->compoundcoordsCircularQueue->getCoordarrayRef(signals->step, query_compound_id);
					const NodeIndex querycompound_hyperorigo = LIMAPOSITIONSYSTEM::getHyperNodeIndex<BoundaryCondition>(compound_origo, querycompound->origo);
					KernelHelpersWarnings::assertHyperorigoIsValid(querycompound_hyperorigo, compound_origo);

					// calc Relative LimaPosition Shift from the origo-shift
					relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, compound_origo).toFloat3();

					neighbor_n_particles[threadIdx.x] = box->compounds[query_compound_id].n_particles;
				
					coords_ptrs[threadIdx.x] = box->compoundcoordsCircularQueue->getCoordarrayRef(signals->step, query_compound_id);
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
			EngineUtils::getCompoundHyperpositionsAsFloat3Async(coords_ptrs[batch_index], utility_buffer_f3, n_particles_neighbor, relshifts[batch_index]);
			__syncthreads();


			// The bonded compounds always comes first in the list
			if (i < compound.n_bonded_compounds)
			{
				BondedParticlesLUT* compoundpair_lut_global = box->bonded_particles_lut_manager->get(compound_index, neighborcompound_id);
				bonded_particles_lut->load(*compoundpair_lut_global);
				__syncthreads();

				if (threadIdx.x < compound.n_particles) {
					force += LJ::computeCompoundCompoundLJForces(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum,
						utility_buffer_f3, n_particles_neighbor, atomtypes, bonded_particles_lut, LJ::CalcLJOrigin::ComComInter, forcefield_shared);
				}
			}
			else {
				if (threadIdx.x < compound.n_particles) {
					force += LJ::computeCompoundCompoundLJForces(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum,
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
		const NodeIndex solventblock_hyperorigo = LIMAPOSITIONSYSTEM::getHyperNodeIndex<BoundaryCondition>(compound_origo, BoxGrid::Get3dIndex(solventblock_id, boxSize_device.boxSizeNM_i));

		const Float3 relpos_shift = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(solventblock_hyperorigo, compound_origo).toFloat3();	// TODO: Only t0 needs to do this

		const SolventBlock* solventblock = box->solventblockgrid_circularqueue->getBlockPtr(solventblock_id, signals->step);
		const int nsolvents_neighbor = solventblock->n_solvents;



		// There are many more solvents in a block, than threads in this kernel, so we need to loop in strides of blockdim.x
		__syncthreads();	// Dont load buffer before all are finished with the previous iteration
		for (uint32_t offset = 0; offset < nsolvents_neighbor; offset += blockDim.x) {
			const uint32_t solvent_index = offset + threadIdx.x;
			const int n_elements_this_stride = LAL::min(nsolvents_neighbor - offset, blockDim.x);

			// Load the positions and add rel shift
			if (solvent_index < nsolvents_neighbor) {
				utility_buffer_f3[threadIdx.x] = solventblock->rel_pos[solvent_index].toFloat3() + relpos_shift;
			}
			__syncthreads();

			if (threadIdx.x < compound.n_particles) {
				force += LJ::computeSolventToCompoundLJForces(compound_positions[threadIdx.x], n_elements_this_stride, utility_buffer_f3, potE_sum, compound.atom_types[threadIdx.x], forcefield_shared);
			}
			__syncthreads();
		}
	}
	
#endif
	// ------------------------------------------------------------------------------------------------------------------------------------------------ //


	// Electrostatic
#ifdef ENABLE_ELECTROSTATICS
	if ( simparams.enable_electrostatics && threadIdx.x < compound.n_particles && SCA::DoRecalc(signals->step)) {
		Float3 abspos = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compound_origo, Coord{ compound_positions[threadIdx.x] });
		PeriodicBoundaryCondition::applyBCNM(abspos);	// TODO: Use generic BC
		sim->charge_octtree->pushChargeToLeaf(abspos, box->compounds[blockIdx.x].atom_charges[threadIdx.x]);
	}
#endif



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
	const NodeIndex block_origo = BoxGrid::Get3dIndex(blockIdx.x, boxSize_device.boxSizeNM_i);

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
			const CompoundCoords* coordarray_ptr = box->compoundcoordsCircularQueue->getCoordarrayRef(signals->step, neighborcompound_index);
			EngineUtils::getCompoundHyperpositionsAsFloat3<BoundaryCondition>(solventblock.origo, coordarray_ptr, utility_buffer, utility_float3, n_compound_particles);


			// Then load atomtypes of neighboring compound
			if (threadIdx.x < n_compound_particles) {
				utility_buffer_small[threadIdx.x] = neighborcompound->atom_types[threadIdx.x];
			}
			__syncthreads();

			//  We can optimize here by loading and calculate the paired sigma and eps, jsut remember to loop threads, if there are many aomttypes.
			if (solvent_active) {
				force += LJ::computeCompoundToSolventLJForces(relpos_self, n_compound_particles, utility_buffer, potE_sum, utility_buffer_small, solventblock.ids[threadIdx.x]);
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
			force += LJ::computeSolventToSolventLJForces(relpos_self, utility_buffer, solventblock.n_solvents, true, potE_sum);
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
					force += LJ::computeSolventToSolventLJForces(relpos_self, utility_buffer, nsolvents_neighbor, false, potE_sum);
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

		EngineUtils::LogSolventData(box, potE_sum, solventblock, solvent_active, force, vel_now, signals->step, sim->databuffers, simparams.data_logging_interval);
	}



	// Push new SolventCoord to global mem
	SolventBlock* solventblock_next_ptr = box->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, signals->step + 1);

	if (SolventBlocksCircularQueue::isTransferStep(signals->step)) {
		SolventBlockTransfers::transferOutAndCompressRemainders<BoundaryCondition>(solventblock, solventblock_next_ptr, relpos_next, utility_buffer_small, sim->transfermodule_array, transferqueues);
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
		utility_coord[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelativeShiftBetweenCoordarrays<BoundaryCondition>(box->compoundcoordsCircularQueue, sim->signals->step, bridge.compound_ids[0], bridge.compound_ids[threadIdx.x]);
	}


	bridge.loadData(&box->bridge_bundle->compound_bridges[blockIdx.x]);
	__syncthreads();

	if (particle_id_bridge < bridge.n_particles) {
		ParticleReference& p_ref = bridge.particle_refs[particle_id_bridge];

		BridgeWarnings::verifyPRefValid(p_ref, bridge);

		const CompoundCoords* coordarray = box->compoundcoordsCircularQueue->getCoordarrayRef(sim->signals->step, p_ref.compound_id);

		Coord relpos = coordarray->rel_positions[p_ref.local_id_compound];
		relpos += utility_coord[p_ref.compoundid_local_to_bridge];
		positions[threadIdx.x] = relpos.toFloat3();
	}
	__syncthreads();

	float potE_sum = 0;
	Float3 force{};

	// ------------------------------------------------------------ Intercompund Operations ------------------------------------------------------------ //
	{											// So for the very first step, these �should all be 0, but they are not??										TODO: Look into this at some point!!!! 
		if (simparams.em_variant)
			force += LimaForcecalc::computeSinglebondForces<true>(bridge.singlebonds, bridge.n_singlebonds, positions, utility_buffer, utility_buffer_f, &potE_sum, 1);
		else 
			force += LimaForcecalc::computeSinglebondForces<false>(bridge.singlebonds, bridge.n_singlebonds, positions, utility_buffer, utility_buffer_f, &potE_sum, 1);

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