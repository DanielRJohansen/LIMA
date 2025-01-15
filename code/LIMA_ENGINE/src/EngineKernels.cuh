//#pragma once - this file must NOT be included multiple times

#include "Engine.cuh"
#include "ForceComputations.cuh"
#include "KernelWarnings.cuh"
#include "EngineUtils.cuh"

#include "SimulationDevice.cuh"
#include "BoundaryCondition.cuh"
#include "SolventBlockTransfers.cuh"
#include "DeviceAlgorithms.cuh"

//#include <cuda/pipeline>
#include "KernelConstants.cuh"

#include "LennardJonesInteractions.cuh"


#include <cfloat>

#pragma warning(push)
#pragma warning(disable:E0020)
#pragma warning(push)
#pragma warning(disable: 20054)

#pragma diag_suppress 20054



















// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//
template <typename BoundaryCondition, bool energyMinimize, bool computePotE> // We dont compute potE if we dont log data this step
__global__ void compoundFarneighborShortrangeInteractionsKernel(const int64_t step, const BoxState boxState, const BoxConfig boxConfig, const NeighborList* const compoundNeighborlists, 
	bool enableES, ForceEnergy* const forceEnergy, const ForceField_NB::ParticleParameters* const ljParams) {
	__shared__ Float3 compound_positions[MAX_COMPOUND_PARTICLES]; // [nm] // TODO: maybe only keep these in register mem    
    __shared__ uint8_t atomTypes[MAX_COMPOUND_PARTICLES];

    const int nParticles = boxConfig.compounds[blockIdx.x].n_particles;
    const int nNonbondedCompoundNeighbors = compoundNeighborlists[blockIdx.x].nNonbondedNeighbors;

	__shared__ Float3 neighborPositions[MAX_COMPOUND_PARTICLES];

	__shared__ float neighborParticlescharges[MAX_COMPOUND_PARTICLES];

	__shared__ ForceField_NB forcefield_shared;
	__shared__ uint8_t neighborAtomstypes[MAX_COMPOUND_PARTICLES];

    const NodeIndex compoundOrigo = boxState.compoundOrigos[blockIdx.x];

	float potE_sum{};
	Float3 force{};
	//NodeIndex compound_origo{}; // TODO: Figure out why TF everything breaks if i make this __shared__???????
	const float particleCharge = enableES	// TODO: this is temporary
	? static_cast<float>(boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x])
	: 0.f;

    atomTypes[threadIdx.x] = boxConfig.compounds[blockIdx.x].atom_types[threadIdx.x];
	{
		auto block = cooperative_groups::this_thread_block();
		cooperative_groups::memcpy_async(block, compound_positions, &boxState.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES], sizeof(Float3) * MAX_COMPOUND_PARTICLES);
		cooperative_groups::memcpy_async(block, &forcefield_shared, &DeviceConstants::forcefield, sizeof(ForceField_NB));
		cooperative_groups::wait(block);
	}



    const int batchsize = 32;
	static_assert(batchsize <= MAX_COMPOUND_PARTICLES, "Not enough threads to load a full batch");
	__shared__ Float3 relshifts[batchsize];	// [nm]
	__shared__ int neighborIds[batchsize]; // either compoundID or solventblockID
	__shared__ int neighborNParticles[batchsize]; // either particlesInCompound or particlesInSolventblock
	
	ForceField_NB::ParticleParameters myParams = ljParams[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
	__shared__ ForceField_NB::ParticleParameters neighborLjParams[MAX_COMPOUND_PARTICLES];
	__shared__ int util;
	
	// --------------------------------------------------------------- Intercompound forces --------------------------------------------------------------- //
	{
		// This part is scary, but it also takes up by far the majority of compute time. We use the utilitybuffer twice simultaneously, so be careful when making changes
		int indexInBatch = batchsize;
		for (int i = 0; i < nNonbondedCompoundNeighbors; i++) {
			__syncthreads();
			util = 0;
			static_assert(MAX_COMPOUND_PARTICLES >= batchsize);
			// First check if we need to load a new batch of relshifts & n_particles for the coming 32 compounds
			if (indexInBatch == batchsize) {
				if (threadIdx.x < batchsize && threadIdx.x + i < nNonbondedCompoundNeighbors) {
					neighborIds[threadIdx.x] = compoundNeighborlists[blockIdx.x].nonbondedNeighborcompoundIds[i + threadIdx.x];

					const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(compoundOrigo, boxState.compoundOrigos[neighborIds[threadIdx.x]]);
					KernelHelpersWarnings::assertHyperorigoIsValid(querycompound_hyperorigo, compoundOrigo);

					// calc Relative LimaPosition Shift from the origo-shift
					relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, compoundOrigo).ToRelpos();
					neighborNParticles[threadIdx.x] = boxConfig.compounds[neighborIds[threadIdx.x]].n_particles;
				}

				indexInBatch = 0;
				__syncthreads();
			}

			const int currentNeighborNParticles = neighborNParticles[indexInBatch];
			if (threadIdx.x < currentNeighborNParticles) {
				// Load first element in batch and sync
				const int currentNeighborId = neighborIds[indexInBatch];

				neighborAtomstypes[threadIdx.x] = boxConfig.compoundsAtomtypes[currentNeighborId * MAX_COMPOUND_PARTICLES + threadIdx.x];
				neighborParticlescharges[threadIdx.x] = boxConfig.compoundsAtomCharges[currentNeighborId * MAX_COMPOUND_PARTICLES + threadIdx.x];
				neighborPositions[threadIdx.x] = boxState.compoundsRelposNm[currentNeighborId * MAX_COMPOUND_PARTICLES + threadIdx.x] + relshifts[indexInBatch];
				neighborLjParams[threadIdx.x] = ljParams[currentNeighborId * MAX_COMPOUND_PARTICLES + threadIdx.x];
			}
			__syncthreads();

            if (threadIdx.x < nParticles) {
                force += LJ::computeCompoundCompoundLJForces<computePotE, energyMinimize>(compound_positions[threadIdx.x], atomTypes[threadIdx.x], potE_sum,
					neighborPositions, neighborNParticles[indexInBatch], neighborAtomstypes, forcefield_shared, particleCharge, neighborParticlescharges, myParams, neighborLjParams, util);
			}

			if (threadIdx.x == 0) {
                //printf("Hitrate %f\n", static_cast<float>(util) / (MAX_COMPOUND_PARTICLES * MAX_COMPOUND_PARTICLES));
			}

			indexInBatch++;
		}
	}
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //

    if (threadIdx.x < nParticles) {
		forceEnergy[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = ForceEnergy{ force, potE_sum };
	}
}


#define compound_index blockIdx.x
template <typename BoundaryCondition, bool energyMinimize, bool computePotE> // We dont compute potE if we dont log data this step
__global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel(SimulationDevice* sim, const int64_t step, ForceEnergy* const forceEnergy) {
	__shared__ CompoundCompact compound;				// Mostly bond information
	__shared__ Float3 compound_positions[MAX_COMPOUND_PARTICLES]; // [nm]

	__shared__ int nGridnodes;

	__shared__ Float3 utility_buffer_f3[MAX_COMPOUND_PARTICLES * 2];

	__shared__ BondedParticlesLUT bpLUT;
	__shared__ float particleChargesBuffers[MAX_COMPOUND_PARTICLES * 2];

	__shared__ ForceField_NB forcefield_shared;

	__shared__ uint8_t neighborAtomstypes[MAX_COMPOUND_PARTICLES];

	BoxState* const boxState = sim->boxState;
	const BoxConfig& boxConfig = sim->boxConfig;
	const SimParams& simparams = sim->params;

	float potE_sum{};
	Float3 force{};
	NodeIndex compound_origo{}; // TODO: Figure out why TF everything breaks if i make this __shared__???????
	const float particleCharge = simparams.enable_electrostatics	// TODO: this is temporary
		? static_cast<float>(boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x])
		: 0.f;


	{
		auto block = cooperative_groups::this_thread_block();

		compound_origo = boxState->compoundOrigos[blockIdx.x];
		cooperative_groups::memcpy_async(block, compound_positions, &boxState->compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES], sizeof(Coord) * MAX_COMPOUND_PARTICLES);

		if (threadIdx.x == 0) {
			compound.loadMeta(&boxConfig.compounds[blockIdx.x]);
			nGridnodes = sim->compound_neighborlists[blockIdx.x].n_gridnodes;
		}

		cooperative_groups::memcpy_async(block, &forcefield_shared, &DeviceConstants::forcefield, sizeof(ForceField_NB));

		cooperative_groups::wait(block);
	}
	compound.loadData(&boxConfig.compounds[blockIdx.x]);



	// ------------------------------------------------------------ Intracompound Operations ------------------------------------------------------------ //
	{
		const BondedParticlesLUT* const bplut_global = BondedParticlesLUTHelpers::get(sim->boxConfig.bpLUTs, compound_index, compound_index);
		bpLUT.load(*bplut_global);	// A lut always exists within a compound

		float* particleChargesCompound = &particleChargesBuffers[0];
		particleChargesCompound[threadIdx.x] = particleCharge;
		__syncthreads();

		if (threadIdx.x < compound.n_particles) {
			// Having this inside vs outside the context makes impact the resulting VC, but it REALLY SHOULD NOT
			force += LJ::computeCompoundCompoundLJForces<computePotE, energyMinimize>(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum, compound_positions, compound.n_particles,
				compound.atom_types, &bpLUT, LJ::CalcLJOrigin::ComComIntra, forcefield_shared,
				particleCharge, particleChargesCompound);
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //



	// For neighborcompound and neighbor solventblocks we utilize a batching system
	
	const int batchsize = 32;
	static_assert(batchsize <= MAX_COMPOUND_PARTICLES, "Not enough threads to load a full batch");
	__shared__ Float3 relshifts[batchsize];	// [nm]
	__shared__ int neighborIds[batchsize]; // either compoundID or solventblockID // should be uint16_t? Does it make a diff?
	__shared__ int neighborNParticles[batchsize]; // either particlesInCompound or particlesInSolventblock

	// --------------------------------------------------------------- Intercompound forces --------------------------------------------------------------- //
	{
		__shared__ const BondedParticlesLUT* compoundPairLutPtrs[Compound::max_bonded_compounds];

		Float3* neighborPositions = utility_buffer_f3;
		float* neighborParticlescharges = particleChargesBuffers;

		if (threadIdx.x < compound.n_bonded_compounds) {
			const uint16_t neighborId = boxConfig.compounds[compound_index].bonded_compound_ids[threadIdx.x];

			const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(compound_origo, boxState->compoundOrigos[neighborId]);
			relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, compound_origo).ToRelpos();
			compoundPairLutPtrs[threadIdx.x] = BondedParticlesLUTHelpers::get(sim->boxConfig.bpLUTs, compound_index, neighborId);
		}
		__syncthreads();

		for (int i = 0; i < compound.n_bonded_compounds; i++) {
			const uint16_t neighborId = boxConfig.compounds[compound_index].bonded_compound_ids[i];
			const int neighborNParticles = boxConfig.compounds[neighborId].n_particles;
			
			neighborPositions[threadIdx.x] = boxState->compoundsRelposNm[neighborId * MAX_COMPOUND_PARTICLES + threadIdx.x] + relshifts[i];
			neighborAtomstypes[threadIdx.x] = boxConfig.compoundsAtomtypes[neighborId * MAX_COMPOUND_PARTICLES + threadIdx.x];
			neighborParticlescharges[threadIdx.x] = boxConfig.compoundsAtomCharges[neighborId * MAX_COMPOUND_PARTICLES + threadIdx.x];
			bpLUT.load(*compoundPairLutPtrs[i]);
			__syncthreads();

			if (threadIdx.x < compound.n_particles) {
				force += LJ::computeCompoundCompoundLJForces<computePotE, energyMinimize>(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum,
					neighborPositions, neighborNParticles, neighborAtomstypes, &bpLUT, LJ::CalcLJOrigin::ComComInter, forcefield_shared,
					particleCharge, neighborParticlescharges);
			}
			__syncthreads();
		}

	}
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //



	// --------------------------------------------------------------- Solvation forces --------------------------------------------------------------- //
#ifdef ENABLE_SOLVENTS
	__shared__ SolventBlock* solventblockPtrs[batchsize];
	__shared__ ForcefieldTinymol forcefieldTinymol_shared;
	if (threadIdx.x < ForcefieldTinymol::MAX_TYPES)
		forcefieldTinymol_shared.types[threadIdx.x] = DeviceConstants::tinymolForcefield.types[threadIdx.x];

	{
		int indexInBatch = batchsize;
		for (int i = 0; i < nGridnodes; i++) {
			__syncthreads();
			if (indexInBatch == batchsize) {
				if (threadIdx.x < batchsize && threadIdx.x + i < nGridnodes) {
					neighborIds[threadIdx.x] = sim->compound_neighborlists[blockIdx.x].gridnode_ids[i + threadIdx.x];

					solventblockPtrs[threadIdx.x] = SolventBlocksCircularQueue::getBlockPtr(boxState->solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, neighborIds[threadIdx.x], step);
					const NodeIndex solventblock_hyperorigo = BoundaryCondition::applyHyperpos_Return(compound_origo, BoxGrid::Get3dIndex(neighborIds[threadIdx.x], DeviceConstants::boxSize.boxSizeNM_i));
					relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(solventblock_hyperorigo, compound_origo).ToRelpos();
					neighborNParticles[threadIdx.x] = solventblockPtrs[threadIdx.x]->n_solvents;
				}
				indexInBatch = 0;
				__syncthreads();
			}

			// There are many more solvents in a block, than threads in this kernel, so we need to loop in strides of blockdim.x
			__syncthreads();	// Dont load buffer before all are finished with the previous iteration
			for (uint32_t offset = 0; offset < neighborNParticles[indexInBatch]; offset += blockDim.x) {
				const uint32_t queryIndex = offset + threadIdx.x;
				const int n_elements_this_stride = LAL::min(neighborNParticles[indexInBatch] - offset, blockDim.x);

				// Load the positions and add rel shift
				if (queryIndex < neighborNParticles[indexInBatch]) {
					utility_buffer_f3[threadIdx.x] = solventblockPtrs[indexInBatch]->rel_pos[queryIndex].ToRelpos() + relshifts[indexInBatch];
					neighborAtomstypes[threadIdx.x] = solventblockPtrs[indexInBatch]->atomtypeIds[queryIndex];
				}
				__syncthreads();

				if (threadIdx.x < compound.n_particles) {
					force += LJ::computeSolventToCompoundLJForces<computePotE, energyMinimize>(compound_positions[threadIdx.x], n_elements_this_stride, 
						utility_buffer_f3, potE_sum, compound.atom_types[threadIdx.x], forcefield_shared, forcefieldTinymol_shared, neighborAtomstypes);
				}
				__syncthreads();
			}
			indexInBatch++;
		}
	}
#endif
	// ------------------------------------------------------------------------------------------------------------------------------------------------ //

	if (threadIdx.x < compound.n_particles) {
		//sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyImmediateneighborShortrange[threadIdx.x] = ForceEnergy{ force, potE_sum };
		forceEnergy[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = ForceEnergy{ force, potE_sum };
	}
}
#undef compound_index

template <typename BoundaryCondition, bool energyMinimize>
__global__ void CompoundSnfKernel(SimulationDevice* sim, int64_t step, const UniformElectricField uniformElectricField, ForceEnergy* const forceEnergy) {
	__shared__ int nParticles;

	if (threadIdx.x == 0) {
		nParticles = sim->boxConfig.compounds[blockIdx.x].n_particles;
	}
	__syncthreads();

	float potE_sum{};
	Float3 force{};

	// ------------------------------------------------------------ Supernatural Forces --------------------------------------------------------------- //	
	if (sim->params.snf_select == HorizontalChargeField && nParticles) {
		force += uniformElectricField.GetForce(sim->boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x]);
		// No potE, as kinE in this field approaches infinity, potE approaches -infinity.
	}

	// ------------------------------------------------------------ Push Data --------------------------------------------------------------- //	
	if (threadIdx.x < nParticles) {
		forceEnergy[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = ForceEnergy{ force, potE_sum };
	}
}

template<typename BoundaryCondition, bool emvariant>
__global__ void CompoundIntegrationKernel(SimulationDevice* sim, int64_t step, const CompoundForceEnergyInterims forceEnergies, const ForceEnergy* const bondgroupForceenergies, const ForceEnergy* const forceEnergiesPME) {
	
	__shared__ CompoundCoords compound_coords;
	__shared__ uint8_t atom_types[MAX_COMPOUND_PARTICLES];
	const int nParticles = sim->boxConfig.compounds[blockIdx.x].n_particles;
	if (threadIdx.x == 0) {
		compound_coords.origo = sim->boxState->compoundOrigos[blockIdx.x];
	}	
	compound_coords.rel_positions[threadIdx.x] = sim->boxState->compoundsInterimState[blockIdx.x].coords[threadIdx.x];
	atom_types[threadIdx.x] = sim->boxConfig.compounds[blockIdx.x].atom_types[threadIdx.x];

	// Fetch interims from other kernels
	/*Float3 force = sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyFarneighborShortrange[threadIdx.x].force;
	float potE_sum = sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyFarneighborShortrange[threadIdx.x].potE;*/

	ForceEnergy forceEnergy = forceEnergies.Sum(blockIdx.x, threadIdx.x);	

	{
		const Compound::BondgroupRefManager* const bgReferences = &sim->boxConfig.compounds[blockIdx.x].bondgroupReferences[threadIdx.x];
		for (int i = 0; i < bgReferences->nBondgroupApperances; i++) {
			const BondgroupRef bondgroupRef = bgReferences->bondgroupApperances[i];
			forceEnergy = forceEnergy + bondgroupForceenergies[bondgroupRef.bondgroupId * BondGroup::maxParticles + bondgroupRef.localIndexInBondgroup];
		}
	}
	

	
	__syncthreads();

	// ------------------------------------------------------------ LongRange Electrostatics --------------------------------------------------------------- //	
	if constexpr (ENABLE_ES_LR) {
		if (threadIdx.x < nParticles)
			forceEnergy = forceEnergy + forceEnergiesPME[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
	}


	// ------------------------------------------------------------ Integration --------------------------------------------------------------- //	
#ifdef FORCE_NAN_CHECK
	if (isnan(forceEnergy.force.len()))
		printf("NAN\n");
#endif
	float speed = 0.f;
	if (threadIdx.x < nParticles) {
		const float mass = sim->boxConfig.compounds[blockIdx.x].atomMasses[threadIdx.x];

		// Energy minimize
		if constexpr (emvariant) {
			const Float3 safeForce = EngineUtils::ForceActivationFunction(forceEnergy.force);

			AdamState* const adamState = &sim->adamState[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
			const Coord pos_now = EngineUtils::IntegratePositionADAM(compound_coords.rel_positions[threadIdx.x], safeForce, adamState, step);

			compound_coords.rel_positions[threadIdx.x] = pos_now;// Save pos locally, but only push to box as this kernel ends
		}
		else {
			const Float3 force_prev = sim->boxState->compoundsInterimState[blockIdx.x].forces_prev[threadIdx.x];	// OPTIM: make ref?
			const Float3 vel_prev = sim->boxState->compoundsInterimState[blockIdx.x].vels_prev[threadIdx.x];
			const Float3 vel_now = EngineUtils::integrateVelocityVVS(vel_prev, force_prev, forceEnergy.force, sim->params.dt, mass);
			const Coord pos_now = EngineUtils::integratePositionVVS(compound_coords.rel_positions[threadIdx.x], vel_now, forceEnergy.force, mass, sim->params.dt);
			compound_coords.rel_positions[threadIdx.x] = pos_now;// Save pos locally, but only push to box as this kernel ends

			Float3 velScaled;
			velScaled = vel_now * DeviceConstants::thermostatScalar;

			sim->boxState->compoundsInterimState[blockIdx.x].forces_prev[threadIdx.x] = forceEnergy.force;
			sim->boxState->compoundsInterimState[blockIdx.x].vels_prev[threadIdx.x] = velScaled;

			speed = velScaled.len();
		}
	}
	__syncthreads();

	// ------------------------------------------------------------ Boundary Condition --------------------------------------------------------------- //	
	{
		__shared__ Coord shift_lm;	// Use utility coord for this?
		if (threadIdx.x == 0) {
			shift_lm = LIMAPOSITIONSYSTEM::shiftOrigo(compound_coords, sim->boxConfig.compounds[blockIdx.x].centerparticle_index);
		}
		__syncthreads();

		LIMAPOSITIONSYSTEM_HACK::shiftRelPos(compound_coords, shift_lm);
		__syncthreads();
	}
	LIMAPOSITIONSYSTEM_HACK::applyBC<BoundaryCondition>(compound_coords);
	__syncthreads();

	Float3 force_LJ_sol{};	// temp
	EngineUtils::LogCompoundData(sim->boxConfig.compounds[blockIdx.x], sim->boxparams.total_particles_upperbound, compound_coords, &forceEnergy.potE, forceEnergy.force,
		force_LJ_sol, sim->params, *sim->signals, sim->potE_buffer, sim->traj_buffer, sim->vel_buffer, sim->forceBuffer, speed, step);


	// Push positions for next step
	if (threadIdx.x == 0)
		sim->boxState->compoundOrigos[blockIdx.x] = compound_coords.origo;
	sim->boxState->compoundsInterimState[blockIdx.x].coords[threadIdx.x] = compound_coords.rel_positions[threadIdx.x];
	sim->boxState->compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = compound_coords.rel_positions[threadIdx.x].ToRelpos();
}


static_assert(SolventBlock::MAX_SOLVENTS_IN_BLOCK >= MAX_COMPOUND_PARTICLES, "solventForceKernel was about to reserve an insufficient amount of memory");
template <typename BoundaryCondition, bool energyMinimize>
__global__ void TinymolCompoundinteractionsKernel(BoxState boxState, const BoxConfig boxConfig, const CompoundGridNode* const compoundGrid, int64_t step) {
	__shared__ Float3 utility_buffer[SolventBlock::MAX_SOLVENTS_IN_BLOCK];
	__shared__ uint8_t utility_buffer_small[SolventBlock::MAX_SOLVENTS_IN_BLOCK];
	__shared__ int neighborblockNumElements;
	__shared__ Float3 utility_float3;
	__shared__ ForcefieldTinymol forcefieldTinymol_shared;
	__shared__ int nElementsInBlock;

	// Doubles as block_index_3d!
	const NodeIndex block_origo = BoxGrid::Get3dIndex(blockIdx.x, DeviceConstants::boxSize.boxSizeNM_i);

	const SolventBlock* const solventblock_ptr = SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step);


	if (threadIdx.x < ForcefieldTinymol::MAX_TYPES) {
		forcefieldTinymol_shared.types[threadIdx.x] = DeviceConstants::tinymolForcefield.types[threadIdx.x]; // TODO Check that im using this and not the __constant??
	}

	if (threadIdx.x == 0) {
		nElementsInBlock = solventblock_ptr->n_solvents;
	}
	__syncthreads();


	const bool threadActive = threadIdx.x < nElementsInBlock;

	Float3 force{};
	float potE_sum{};
	const Float3 relpos_self = solventblock_ptr->rel_pos[threadIdx.x].ToRelpos();
	const uint8_t tinymolTypeId = solventblock_ptr->atomtypeIds[threadIdx.x];
	const uint32_t idSelf = solventblock_ptr->ids[threadIdx.x];
	// --------------------------------------------------------------- Molecule Interactions --------------------------------------------------------------- //	
	{
		// Thread 0 finds n nearby compounds
		const CompoundGridNode* compoundgridnode = BoxGrid::GetNodePtr(compoundGrid, blockIdx.x);
		if (threadIdx.x == 0) { neighborblockNumElements = compoundgridnode->n_nearby_compounds; }
		__syncthreads();



		for (int i = 0; i < neighborblockNumElements; i++) {
			const uint16_t neighborcompound_index = compoundgridnode->compoundidsWithinLjCutoff[i];
			const Compound* neighborcompound = &boxConfig.compounds[neighborcompound_index];
			const int n_compound_particles = neighborcompound->n_particles;

			// All threads help loading the molecule
			// First load particles of neighboring compound
			EngineUtils::getCompoundHyperpositionsAsFloat3<BoundaryCondition>(block_origo, boxState.compoundOrigos[neighborcompound_index], &boxState.compoundsRelposNm[neighborcompound_index * MAX_COMPOUND_PARTICLES],
				utility_buffer, utility_float3, n_compound_particles);


			// Then load atomtypes of neighboring compound
			if (threadIdx.x < n_compound_particles) {
				utility_buffer_small[threadIdx.x] = neighborcompound->atom_types[threadIdx.x];
			}
			__syncthreads();

			//  We can optimize here by loading and calculate the paired sigma and eps, jsut remember to loop threads, if there are many aomttypes.
			if (threadActive) {// TODO: use computePote template param here
				force += LJ::computeCompoundToSolventLJForces<true, energyMinimize>(relpos_self, n_compound_particles, utility_buffer, potE_sum,
					utility_buffer_small, idSelf, DeviceConstants::tinymolForcefield, tinymolTypeId);
			}
			__syncthreads();
		}
	}

	SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step)->forceEnergiesCompoundinteractions[threadIdx.x] = ForceEnergy{ force, potE_sum };
}

static_assert(SolventBlock::MAX_SOLVENTS_IN_BLOCK >= MAX_COMPOUND_PARTICLES, "solventForceKernel was about to reserve an insufficient amount of memory");
template <typename BoundaryCondition, bool energyMinimize>
__global__ void solventForceKernel(BoxState boxState, const BoxConfig boxConfig, const CompoundGridNode* const compoundGrid, int64_t step) {
	__shared__ Float3 utility_buffer[SolventBlock::MAX_SOLVENTS_IN_BLOCK];
	__shared__ uint8_t utility_buffer_small[SolventBlock::MAX_SOLVENTS_IN_BLOCK];
	__shared__ ForcefieldTinymol forcefieldTinymol_shared;
	__shared__ int nElementsInBlock;

	// Doubles as block_index_3d!
	const NodeIndex block_origo = BoxGrid::Get3dIndex(blockIdx.x, DeviceConstants::boxSize.boxSizeNM_i);

	const SolventBlock* const solventblock_ptr = SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step);


	if (threadIdx.x < ForcefieldTinymol::MAX_TYPES) {
		forcefieldTinymol_shared.types[threadIdx.x] = DeviceConstants::tinymolForcefield.types[threadIdx.x]; // TODO Check that im using this and not the __constant??
	}

	if (threadIdx.x == 0) {
		nElementsInBlock = solventblock_ptr->n_solvents;
	}
	__syncthreads();
	
	const bool threadActive = threadIdx.x < nElementsInBlock;

	Float3 force{};
	float potE_sum{};
	const Float3 relpos_self = solventblock_ptr->rel_pos[threadIdx.x].ToRelpos();
	const uint8_t tinymolTypeId = solventblock_ptr->atomtypeIds[threadIdx.x];


	// --------------------------------------------------------------- Intrablock TinyMolState Interactions ----------------------------------------------------- //
	{		
		if (threadActive) {
			utility_buffer[threadIdx.x] = relpos_self;
			utility_buffer_small[threadIdx.x] = tinymolTypeId;
		}
		__syncthreads();
		if (threadActive) {
			force += LJ::computeSolventToSolventLJForces<true, energyMinimize>(relpos_self, tinymolTypeId, utility_buffer, nElementsInBlock, true, potE_sum, forcefieldTinymol_shared, utility_buffer_small);
		}
	}	
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //

	// --------------------------------------------------------------- Interblock TinyMolState Interactions ----------------------------------------------------- //
	__shared__ BoxGrid::TinymolBlockAdjacency::BlockRef nearbyBlock[BoxGrid::TinymolBlockAdjacency::nNearbyBlocks];
	if (threadIdx.x < BoxGrid::TinymolBlockAdjacency::nNearbyBlocks) {
		nearbyBlock[threadIdx.x] = BoxGrid::TinymolBlockAdjacency::GetPtrToNearbyBlockids(blockIdx.x, boxConfig.tinymolNearbyBlockIds)[threadIdx.x];
	}
	__syncthreads();


	for (int i = 0; i < BoxGrid::TinymolBlockAdjacency::nNearbyBlocks; i++) {
		const int blockindex_neighbor = nearbyBlock[i].blockId;

		const SolventBlock* solventblock_neighbor = SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockindex_neighbor, step);
		const int nsolvents_neighbor = solventblock_neighbor->n_solvents;

		// All threads help loading the solvent, and shifting it's relative position reletive to this solventblock
		__syncthreads();
		if (threadIdx.x < nsolvents_neighbor) {
			utility_buffer[threadIdx.x] = solventblock_neighbor->rel_pos[threadIdx.x].ToRelpos() + nearbyBlock[i].relShift;
			utility_buffer_small[threadIdx.x] = solventblock_neighbor->atomtypeIds[threadIdx.x];
		}
		__syncthreads();

		if (threadActive) {
			force += LJ::computeSolventToSolventLJForces<true, energyMinimize>(relpos_self, tinymolTypeId, utility_buffer, nsolvents_neighbor, false, potE_sum, forcefieldTinymol_shared, utility_buffer_small);
		}
	}

	// Finally push force and potE for next kernel
	SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step)->forceEnergiesTinymolinteractions[threadIdx.x] = ForceEnergy{ force, potE_sum };
}


template <typename BoundaryCondition, bool energyMinimize, bool transferOutThisStep>
__global__ void TinymolIntegrationLoggingAndTransferout(SimulationDevice* sim, int64_t step) {
	__shared__ SolventBlock solventblock;
	__shared__ uint8_t utility_buffer_small[SolventBlock::MAX_SOLVENTS_IN_BLOCK];

	// Doubles as block_index_3d!
	const NodeIndex block_origo = BoxGrid::Get3dIndex(blockIdx.x, DeviceConstants::boxSize.boxSizeNM_i);

	BoxState* boxState = sim->boxState;
	const SimParams& simparams = sim->params;
	SolventBlock* solventblock_ptr = SolventBlocksCircularQueue::getBlockPtr(boxState->solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step);

	const Float3 force = solventblock_ptr->forceEnergiesCompoundinteractions[threadIdx.x].force + solventblock_ptr->forceEnergiesTinymolinteractions[threadIdx.x].force;
	const float potE = solventblock_ptr->forceEnergiesCompoundinteractions[threadIdx.x].potE + solventblock_ptr->forceEnergiesTinymolinteractions[threadIdx.x].potE;

	if (threadIdx.x == 0) {
		solventblock.loadMeta(*solventblock_ptr);
	}
	__syncthreads();
	const bool solventActive = threadIdx.x < solventblock.n_solvents;
	solventblock.loadData(*solventblock_ptr);
	__syncthreads();



	Coord relpos_next{};
	if (solventActive) {
		TinyMolState& tinyMols_ref = boxState->tinyMols[solventblock.ids[threadIdx.x]];	// TinyMolState private data, for VVS
		const float mass = DeviceConstants::tinymolForcefield.types[tinyMols_ref.tinymolTypeIndex].mass;

		if constexpr (energyMinimize) {
			const Float3 safeForce = EngineUtils::ForceActivationFunction(force);
			AdamState* const adamStatePtr = &sim->adamState[sim->boxparams.n_compounds * MAX_COMPOUND_PARTICLES + solventblock.ids[threadIdx.x]];
			const Coord pos_now = EngineUtils::IntegratePositionADAM(solventblock.rel_pos[threadIdx.x], safeForce, adamStatePtr, step);

			relpos_next = pos_now;
			EngineUtils::LogSolventData(sim->boxparams, potE, block_origo, solventblock.ids[threadIdx.x], solventblock.rel_pos[threadIdx.x], solventActive,
				force, Float3{}, step, sim->potE_buffer, sim->traj_buffer, sim->vel_buffer, simparams.data_logging_interval);
		}
		else {
			if (force.isNan() || force.lenSquared() >= FLT_MAX)
				force.print('S');

			Float3 vel_now = EngineUtils::integrateVelocityVVS(tinyMols_ref.vel_prev, tinyMols_ref.force_prev, force, simparams.dt, mass);
			const Coord pos_now = EngineUtils::integratePositionVVS(solventblock.rel_pos[threadIdx.x], vel_now, force, mass, simparams.dt);

			vel_now = vel_now * DeviceConstants::thermostatScalar;

			tinyMols_ref.vel_prev = vel_now;
			tinyMols_ref.force_prev = force;

			// Save pos locally, but only push to box as this kernel ends
			relpos_next = pos_now;
			EngineUtils::LogSolventData(sim->boxparams, potE, block_origo, solventblock.ids[threadIdx.x], solventblock.rel_pos[threadIdx.x], solventActive,
				force, vel_now, step, sim->potE_buffer, sim->traj_buffer, sim->vel_buffer, simparams.data_logging_interval);
		}
	}



	// Push new SolventCoord to global mem
	SolventBlock* const solventblock_next_ptr = SolventBlocksCircularQueue::getBlockPtr(boxState->solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step + 1);

	if constexpr (transferOutThisStep) {
		__shared__ SolventTransferqueue<SolventBlockTransfermodule::max_queue_size> transferqueues[6];		// TODO: Use template to make identical kernel, so the kernel with transfer is slower and larger, and the rest remain fast!!!!
		// Init queue, otherwise it will contain wierd values 
		if (threadIdx.x < 6) {
			transferqueues[threadIdx.x] = SolventTransferqueue<SolventBlockTransfermodule::max_queue_size>{};
		}
		utility_buffer_small[threadIdx.x] = 0;
		__syncthreads();
		SolventBlockTransfers::transferOutAndCompressRemainders<BoundaryCondition>(solventblock, solventblock_next_ptr, relpos_next, utility_buffer_small, sim->transfermodule_array, transferqueues);
	}
	else {
		solventblock_next_ptr->rel_pos[threadIdx.x] = relpos_next;
		solventblock_next_ptr->ids[threadIdx.x] = solventblock.ids[threadIdx.x];
		solventblock_next_ptr->atomtypeIds[threadIdx.x] = solventblock.atomtypeIds[threadIdx.x];
		if (threadIdx.x == 0) {
			solventblock_next_ptr->n_solvents = solventblock.n_solvents;
		}
	}
}


// This is run before step.inc(), but will always publish results to the first array in grid!
template <typename BoundaryCondition>
__global__ void solventTransferKernel(SimulationDevice* sim, int64_t step) {
	BoxState* boxState = sim->boxState;

	SolventBlockTransfermodule* transfermodule = &sim->transfermodule_array[blockIdx.x];
	
	SolventBlock* solventblock_current = SolventBlocksCircularQueue::getBlockPtr(boxState->solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step);
	SolventBlock* solventblock_next = SolventBlocksCircularQueue::getBlockPtr(boxState->solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step + 1);

	SolventTransferWarnings::assertSolventsEqualNRemain(*solventblock_next, *transfermodule);

	// Handling incoming transferring solvents
	int n_solvents_next = transfermodule->n_remain;
	for (int queue_index = 0; queue_index < SolventBlockTransfermodule::n_queues; queue_index++) {
		auto* queue = &transfermodule->transfer_queues[queue_index];
		if (threadIdx.x < queue->n_elements) {
			const int incoming_index = n_solvents_next + threadIdx.x;

			solventblock_next->rel_pos[incoming_index] = queue->rel_positions[threadIdx.x];
			solventblock_next->ids[incoming_index] = queue->ids[threadIdx.x];
			solventblock_next->atomtypeIds[incoming_index] = queue->atomtypeIds[threadIdx.x];
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


static const int THREADS_PER_BONDSGROUPSKERNEL = BondGroup::maxParticles;
template <typename BoundaryCondition, bool emVariant>
__global__ void BondgroupsKernel(const BondGroup* const bondGroups, const BoxState boxState, ForceEnergy* const forceEnergiesOut) {
	__shared__ Float3 positions[BondGroup::maxParticles];


	__shared__ Float3 forcesInterrim[BondGroup::maxParticles];
	__shared__ float potEInterrim[BondGroup::maxParticles];

	static const int batchSize = THREADS_PER_BONDSGROUPSKERNEL;
	static const int largestBondBytesize = std::max(sizeof(AngleUreyBradleyBond), sizeof(DihedralBond));
	__shared__ char _bondsBuffer[largestBondBytesize * batchSize];	
	__shared__ NodeIndex origo;

	const BondGroup* const bondGroup = &bondGroups[blockIdx.x];
	const BondGroup::ParticleRef pRef = bondGroup->particles[threadIdx.x];

	if (threadIdx.x == 0) {
		origo = boxState.compoundOrigos[pRef.compoundId];
	}
	__syncthreads();

	if (threadIdx.x < bondGroup->nParticles) {
		// Calculate necessary shift in relative positions for right, so right share the origo with left.
		const NodeIndex myNodeindex = BoundaryCondition::applyHyperpos_Return(origo, boxState.compoundOrigos[pRef.compoundId]);
		//KernelHelpersWarnings::assertHyperorigoIsValid(querycompound_hyperorigo, compoundOrigo);

		// calc Relative LimaPosition Shift from the origo-shift
		const Float3 relShift = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(myNodeindex, origo).ToRelpos();

		positions[threadIdx.x] = boxState.compoundsRelposNm[pRef.compoundId * MAX_COMPOUND_PARTICLES + pRef.localIdInCompound] + relShift;
		
	}
//	__syncthreads();


	Float3 force{};
	float potE{};

	
	{
		SingleBond* bondsBuffer = reinterpret_cast<SingleBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nSinglebonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nSinglebonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->singlebonds[bondIndex];
			}
			__syncthreads();

			force += LimaForcecalc::computeSinglebondForces<emVariant>(bondsBuffer, std::min(batchSize, bondGroup->nSinglebonds - batchStart), positions, forcesInterrim, potEInterrim, &potE, 0);
		}
	}

	{
		AngleUreyBradleyBond* bondsBuffer = reinterpret_cast<AngleUreyBradleyBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nAnglebonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nAnglebonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->anglebonds[bondIndex];
			}
			__syncthreads();

			force += LimaForcecalc::computeAnglebondForces(bondsBuffer, std::min(batchSize, bondGroup->nAnglebonds - batchStart), positions, forcesInterrim, potEInterrim, &potE);
		}
	}

	{
		DihedralBond* bondsBuffer = reinterpret_cast<DihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nDihedralbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nDihedralbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->dihedralbonds[bondIndex];
			}
			__syncthreads();

			force += LimaForcecalc::computeDihedralForces(bondsBuffer, std::min(batchSize, bondGroup->nDihedralbonds - batchStart), positions, forcesInterrim, potEInterrim, &potE);
		}
	}

	{
		ImproperDihedralBond* bondsBuffer = reinterpret_cast<ImproperDihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nImproperdihedralbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nImproperdihedralbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->improperdihedralbonds[bondIndex];
			}
			__syncthreads();

			force += LimaForcecalc::computeImproperdihedralForces(bondsBuffer, std::min(batchSize, bondGroup->nImproperdihedralbonds - batchStart), positions, forcesInterrim, potEInterrim, &potE);
		}
	}


	{
		// TODO: i have no clue if pairbonds should also compute SR electrostatics?
		PairBond* bondsBuffer = reinterpret_cast<PairBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nPairbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nPairbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->pairbonds[bondIndex];
			}
			__syncthreads();

			force += LimaForcecalc::computePairbondForces(bondsBuffer, std::min(batchSize, bondGroup->nPairbonds - batchStart), positions, forcesInterrim, potEInterrim, &potE);
		}
	}

	forceEnergiesOut[blockIdx.x * BondGroup::maxParticles + threadIdx.x] = ForceEnergy{ force, potE };
}

#pragma warning (pop)
#pragma warning (pop)
