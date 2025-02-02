//#pragma once - this file must NOT be included multiple times

#include "Engine.cuh"
#include "ForceComputations.cuh"
#include "KernelWarnings.cuh"
#include "EngineUtils.cuh"

#include "SimulationDevice.cuh"
#include "BoundaryCondition.cuh"
#include "SolventBlockTransfers.cuh"
#include "DeviceAlgorithms.cuh"
#include "Neighborlists.cuh"

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
__global__ void compoundFarneighborShortrangeInteractionsKernel(bool enableES, ForceEnergy* const forceEnergy, const CompoundQuickData* const compoundQuickDataBuffer,
                    const uint16_t* const compoundsNNeighborNonbondedCompounds, const NeighborList::IdAndRelshift* const compoundsNeighborNonbondedCompounds, const uint8_t* const nParticlesInCompoundsBuffer)
{
    const int batchsize = 32;

    __shared__ CompoundQuickData compoundQuickData; // 768 bytes
    __shared__ NeighborList::IdAndRelshift neighborCompounds[batchsize]; // 512 bytes

    __shared__ int nNonbondedCompoundNeighbors;    
    __shared__ int nParticles;

    {
        auto block = cooperative_groups::this_thread_block();
        cooperative_groups::memcpy_async(block, &compoundQuickData, &compoundQuickDataBuffer[blockIdx.x], sizeof(CompoundQuickData));
        cooperative_groups::wait(block);
    }

    if (threadIdx.x == 0) {
        nNonbondedCompoundNeighbors = compoundsNNeighborNonbondedCompounds[blockIdx.x];
        nParticles = nParticlesInCompoundsBuffer[blockIdx.x];

        //if (blockIdx.x == 726) {
        //    printf("\n\n\n new\n\n\n");
        //    for (int i = 0; i < nNonbondedCompoundNeighbors; i++) {
        //        printf("%d\n", compoundsNeighborNonbondedCompounds[blockIdx.x * NeighborList::compoundsMaxNearbyCompounds + i].id);
        //    }
        //}
    }

    const Float3 myPos = compoundQuickData.relPos[threadIdx.x];
    const float myCharge = enableES ? compoundQuickData.charges[threadIdx.x] : 0.f;
    const ForceField_NB::ParticleParameters myParams = compoundQuickData.ljParams[threadIdx.x];

    static_assert(batchsize <= MAX_COMPOUND_PARTICLES, "Not enough threads to load a full batch");

    float potE_sum{};
    Float3 force{};
    // --------------------------------------------------------------- Intercompound forces --------------------------------------------------------------- //
    {
        auto block = cooperative_groups::this_thread_block();

        for (int batchIndex = 0; batchIndex < nNonbondedCompoundNeighbors/batchsize + 1; batchIndex++) {
            __syncthreads();

            cooperative_groups::memcpy_async(block, neighborCompounds, &compoundsNeighborNonbondedCompounds[blockIdx.x * NeighborList::compoundsMaxNearbyCompounds + batchIndex * batchsize], sizeof(NeighborList::IdAndRelshift) * batchsize);
            cooperative_groups::wait(block);

            const int nElementsInBatch = min(nNonbondedCompoundNeighbors-batchIndex*batchsize, batchsize);
            for (int indexInBatch = 0; indexInBatch < nElementsInBatch; indexInBatch++) {
                cooperative_groups::memcpy_async(block, &compoundQuickData, &compoundQuickDataBuffer[neighborCompounds[indexInBatch].id], sizeof(CompoundQuickData));
                cooperative_groups::wait(block);

                if (threadIdx.x < nParticles) {
                    force += LJ::computeCompoundCompoundLJForces<computePotE, energyMinimize>(myPos - neighborCompounds[indexInBatch].relShift, potE_sum,
                        compoundQuickData.relPos, neighborCompounds[indexInBatch].nParticles, myCharge, compoundQuickData.charges, myParams, compoundQuickData.ljParams);
                }
                __syncthreads();
            }
        }
    }
    // ------------------------------------------------------------------------------------------------------------------------------------------------------ //

    static_assert(sizeof(CompoundQuickData) >= sizeof(ForceEnergy) * MAX_COMPOUND_PARTICLES);
    ForceEnergy* forceEnergyOut = (ForceEnergy*)&compoundQuickData;
    __syncthreads();
    forceEnergyOut[threadIdx.x] = ForceEnergy{ force, potE_sum };
    __syncthreads();

    auto block = cooperative_groups::this_thread_block();
    cooperative_groups::memcpy_async(block, &forceEnergy[blockIdx.x * MAX_COMPOUND_PARTICLES], forceEnergyOut, sizeof(ForceEnergy) * MAX_COMPOUND_PARTICLES);
}


#define compound_index blockIdx.x
template <typename BoundaryCondition, bool energyMinimize, bool computePotE> // We dont compute potE if we dont log data this step
__global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel(SimulationDevice* sim, const int64_t step, ForceEnergy* const forceEnergy, NeighborList::Buffers nlistBuffers) {
	__shared__ CompoundCompact compound;				// Mostly bond information
	__shared__ Float3 compound_positions[MAX_COMPOUND_PARTICLES]; // [nm]

	__shared__ int nGridnodes;

	__shared__ Float3 utility_buffer_f3[MAX_COMPOUND_PARTICLES * 2];

	__shared__ BondedParticlesLUT bpLUT;
	__shared__ float particleChargesBuffers[MAX_COMPOUND_PARTICLES * 2];

	__shared__ ForceField_NB forcefield_shared;

	__shared__ uint8_t neighborAtomstypes[MAX_COMPOUND_PARTICLES];

	const BoxState& boxState = sim->boxState;
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

		compound_origo = boxState.compoundOrigos[blockIdx.x];
		cooperative_groups::memcpy_async(block, compound_positions, &boxState.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES], sizeof(Coord) * MAX_COMPOUND_PARTICLES);

		if (threadIdx.x == 0) {
			compound.loadMeta(&boxConfig.compounds[blockIdx.x]);
			nGridnodes = nlistBuffers.compoundsNNearbyGridnodes[blockIdx.x];//  sim->compound_neighborlists[blockIdx.x].n_gridnodes;
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

			const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(compound_origo, boxState.compoundOrigos[neighborId]);
			relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::GetRelShiftFromOrigoShift_Float3(querycompound_hyperorigo, compound_origo);
			compoundPairLutPtrs[threadIdx.x] = BondedParticlesLUTHelpers::get(sim->boxConfig.bpLUTs, compound_index, neighborId);
		}
		__syncthreads();

		for (int i = 0; i < compound.n_bonded_compounds; i++) {
			const uint16_t neighborId = boxConfig.compounds[compound_index].bonded_compound_ids[i];
			const int neighborNParticles = boxConfig.compounds[neighborId].n_particles;
			
			neighborPositions[threadIdx.x] = boxState.compoundsRelposNm[neighborId * MAX_COMPOUND_PARTICLES + threadIdx.x] + relshifts[i];
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
					//printf("Load index %d\n", blockIdx.x * NeighborList::compoundsMaxNearbyGridnodes + i + threadIdx.x);
					neighborIds[threadIdx.x] = nlistBuffers.compoundsNearbyGridnodes[blockIdx.x * NeighborList::compoundsMaxNearbyGridnodes + i + threadIdx.x];
						//rbyBlocks  sim->compound_neighborlists[blockIdx.x].gridnode_ids[i + threadIdx.x];

					solventblockPtrs[threadIdx.x] = SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, neighborIds[threadIdx.x], step);
					const NodeIndex solventblock_hyperorigo = BoundaryCondition::applyHyperpos_Return(compound_origo, BoxGrid::Get3dIndex(neighborIds[threadIdx.x], DeviceConstants::boxSize.boxSizeNM_i));
					relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::GetRelShiftFromOrigoShift_Float3(solventblock_hyperorigo, compound_origo);
					neighborNParticles[threadIdx.x] = solventblockPtrs[threadIdx.x]->nParticles;
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
		forceEnergy[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = ForceEnergy{ force, potE_sum };
	}
}
#undef compound_index

template <typename BoundaryCondition, bool energyMinimize>
__global__ void CompoundSnfKernel(SimulationDevice* sim, const UniformElectricField uniformElectricField, ForceEnergy* const forceEnergy) {
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

// Instead of updating origo every step (which doesnt really give us precision, as the relpos can stretch very far from origo no problem), we only update 
// just before a new NLIST update. This means, nlists can precompute all compound-compound interactions relshift once
template<typename BoundaryCondition, bool emvariant>
__global__ void CompoundIntegrationKernel(SimulationDevice* sim, int64_t step, const ForceEnergyInterims forceEnergies, CompoundQuickData* const compoundQuickData, bool updateOrigo) {

	__shared__ CompoundCoords compound_coords;
	__shared__ uint8_t atom_types[MAX_COMPOUND_PARTICLES];
	const int nParticles = sim->boxConfig.compounds[blockIdx.x].n_particles;
	if (threadIdx.x == 0) {
		compound_coords.origo = sim->boxState.compoundOrigos[blockIdx.x];
	}
	compound_coords.rel_positions[threadIdx.x] = sim->boxState.compoundsInterimState[blockIdx.x].coords[threadIdx.x];
	atom_types[threadIdx.x] = sim->boxConfig.compounds[blockIdx.x].atom_types[threadIdx.x];

	// Fetch interims from other kernels
	ForceEnergy forceEnergy = forceEnergies.SumCompound(blockIdx.x, threadIdx.x);
	{
		const Compound::BondgroupRefManager* const bgReferences = &sim->boxConfig.compounds[blockIdx.x].bondgroupReferences[threadIdx.x];
		for (int i = 0; i < bgReferences->nBondgroupApperances; i++) {
			const BondgroupRef bondgroupRef = bgReferences->bondgroupApperances[i];
			forceEnergy = forceEnergy + forceEnergies.forceEnergiesBondgroups[bondgroupRef.bondgroupId * BondGroup::maxParticles + bondgroupRef.localIndexInBondgroup];
		}
	}



	__syncthreads();

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
			const Float3 force_prev = sim->boxState.compoundsInterimState[blockIdx.x].forces_prev[threadIdx.x];	// OPTIM: make ref?
			const Float3 vel_prev = sim->boxState.compoundsInterimState[blockIdx.x].vels_prev[threadIdx.x];
			const Float3 vel_now = EngineUtils::integrateVelocityVVS(vel_prev, force_prev, forceEnergy.force, sim->params.dt, mass);
			const Coord pos_now = EngineUtils::integratePositionVVS(compound_coords.rel_positions[threadIdx.x], vel_now, forceEnergy.force, mass, sim->params.dt);
			compound_coords.rel_positions[threadIdx.x] = pos_now;// Save pos locally, but only push to box as this kernel ends

			Float3 velScaled;
			velScaled = vel_now * DeviceConstants::thermostatScalar;

			sim->boxState.compoundsInterimState[blockIdx.x].forces_prev[threadIdx.x] = forceEnergy.force;
			sim->boxState.compoundsInterimState[blockIdx.x].vels_prev[threadIdx.x] = velScaled;

			speed = velScaled.len();
		}
	}
	__syncthreads();

	// ------------------------------------------------------------ Boundary Condition --------------------------------------------------------------- //	
	if (updateOrigo) {
		__shared__ Coord shift_lm;	// Use utility coord for this?
		if (threadIdx.x == 0) {
			shift_lm = LIMAPOSITIONSYSTEM::shiftOrigo(compound_coords, sim->boxConfig.compounds[blockIdx.x].centerparticle_index);
		}
		__syncthreads();

		LIMAPOSITIONSYSTEM_HACK::shiftRelPos(compound_coords, shift_lm);
		__syncthreads();

		LIMAPOSITIONSYSTEM_HACK::applyBC<BoundaryCondition>(compound_coords);
		__syncthreads();
	}

	Float3 force_LJ_sol{};	// temp
	EngineUtils::LogCompoundData(sim->boxConfig.compounds[blockIdx.x], sim->boxparams.total_particles_upperbound, compound_coords, &forceEnergy.potE, forceEnergy.force,
		force_LJ_sol, sim->params, *sim->signals, sim->potE_buffer, sim->traj_buffer, sim->vel_buffer, sim->forceBuffer, speed, step);


	// Push positions for next step
	if (threadIdx.x == 0)
		sim->boxState.compoundOrigos[blockIdx.x] = compound_coords.origo;
	sim->boxState.compoundsInterimState[blockIdx.x].coords[threadIdx.x] = compound_coords.rel_positions[threadIdx.x];
	sim->boxState.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = compound_coords.rel_positions[threadIdx.x].ToRelpos();
	compoundQuickData[blockIdx.x].relPos[threadIdx.x] = compound_coords.rel_positions[threadIdx.x].ToRelpos();
}


static_assert(SolventBlock::MAX_SOLVENTS_IN_BLOCK >= MAX_COMPOUND_PARTICLES, "solventForceKernel was about to reserve an insufficient amount of memory");
template <typename BoundaryCondition, bool energyMinimize>
__global__ void TinymolCompoundinteractionsKernel(BoxState boxState, const BoxConfig boxConfig, const NeighborList::Buffers nlistBuffers, int64_t step, ForceEnergy* const forceEnergies) {
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
		nElementsInBlock = solventblock_ptr->nParticles;
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
		const uint16_t* const neighborcompound_indices = &nlistBuffers.gridnodesNearbyCompounds[blockIdx.x * NeighborList::gridnodesMaxNearbyCompounds];
		if (threadIdx.x == 0) { neighborblockNumElements = nlistBuffers.gridnodesNNearbyCompounds[blockIdx.x]; }
		__syncthreads();



		for (int i = 0; i < neighborblockNumElements; i++) {
			const uint16_t neighborcompound_index = neighborcompound_indices[i];
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

	forceEnergies[blockIdx.x * SolventBlock::MAX_SOLVENTS_IN_BLOCK + threadIdx.x] = ForceEnergy{ force, potE_sum };
}

static_assert(SolventBlock::MAX_SOLVENTS_IN_BLOCK >= MAX_COMPOUND_PARTICLES, "solventForceKernel was about to reserve an insufficient amount of memory");
template <typename BoundaryCondition, bool energyMinimize>
__global__ void solventForceKernel(BoxState boxState, const BoxConfig boxConfig, int64_t step, ForceEnergy* const forceEnergies) {
    //__shared__ Float3 utility_buffer[SolventBlock::MAX_SOLVENTS_IN_BLOCK];

    static const int utilityBufferBytesize = sizeof(ForceEnergy) * SolventBlock::MAX_SOLVENTS_IN_BLOCK;
    static_assert(utilityBufferBytesize >= sizeof(Float3) * SolventBlock::MAX_SOLVENTS_IN_BLOCK + sizeof(ForcefieldTinymol));
    static_assert(sizeof(Coord) == sizeof(Float3));
    __shared__ uint8_t utilityBuffer[utilityBufferBytesize];
    Coord* positionsBuffer_coord = (Coord*)&utilityBuffer[0];
    Float3* positionsBuffer_relpos = (Float3*)&utilityBuffer[0];
    ForcefieldTinymol* forcefieldTinymolShared = (ForcefieldTinymol*)&utilityBuffer[sizeof(Float3) * SolventBlock::MAX_SOLVENTS_IN_BLOCK];
    ForceEnergy* forceEnergyOut = (ForceEnergy*)utilityBuffer; // Overlaps all of the buffers above

	__shared__ uint8_t utility_buffer_small[SolventBlock::MAX_SOLVENTS_IN_BLOCK];
	__shared__ int nElementsInBlock;

	const SolventBlock* const solventblock_ptr = SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step);


	if (threadIdx.x < ForcefieldTinymol::MAX_TYPES) {
        forcefieldTinymolShared->types[threadIdx.x] = DeviceConstants::tinymolForcefield.types[threadIdx.x]; // TODO Check that im using this and not the __constant??
	}

	if (threadIdx.x == 0) {
		nElementsInBlock = solventblock_ptr->nParticles;
	}
	__syncthreads();
	
	const bool threadActive = threadIdx.x < nElementsInBlock;

	Float3 force{};
	float potE_sum{};
	const Float3 relpos_self = solventblock_ptr->rel_pos[threadIdx.x].ToRelpos();
	const uint8_t tinymolTypeId = solventblock_ptr->atomtypeIds[threadIdx.x];


	// --------------------------------------------------------------- Intrablock TinyMolParticleState Interactions ----------------------------------------------------- //
	{		
		if (threadActive) {
            positionsBuffer_relpos[threadIdx.x] = relpos_self;
			utility_buffer_small[threadIdx.x] = tinymolTypeId;
		}
		__syncthreads();
		if (threadActive) {
            force += LJ::computeSolventToSolventLJForces<true, energyMinimize, true>
				(relpos_self, tinymolTypeId, positionsBuffer_relpos, nElementsInBlock, potE_sum, *forcefieldTinymolShared, utility_buffer_small, solventblock_ptr->particlesBondgroupIds);
		}
	}	
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //

	// --------------------------------------------------------------- Interblock TinyMolParticleState Interactions ----------------------------------------------------- //
	__shared__ BoxGrid::TinymolBlockAdjacency::BlockRef nearbyBlock[BoxGrid::TinymolBlockAdjacency::nNearbyBlocks];
    static_assert(SolventBlock::MAX_SOLVENTS_IN_BLOCK >= BoxGrid::TinymolBlockAdjacency::nNearbyBlocks);
	if (threadIdx.x < BoxGrid::TinymolBlockAdjacency::nNearbyBlocks) {
		nearbyBlock[threadIdx.x] = BoxGrid::TinymolBlockAdjacency::GetPtrToNearbyBlockids(blockIdx.x, boxConfig.tinymolNearbyBlockIds)[threadIdx.x];
	}
	__syncthreads();


	for (int i = 0; i < BoxGrid::TinymolBlockAdjacency::nNearbyBlocks; i++) {
		const int blockindex_neighbor = nearbyBlock[i].blockId;

		const SolventBlock* solventblock_neighbor = SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockindex_neighbor, step);
		const int nsolvents_neighbor = solventblock_neighbor->nParticles;

		// All threads help loading the solvent, and shifting it's relative position reletive to this solventblock
        __syncthreads();

        //auto block = cooperative_groups::this_thread_block();
        //cooperative_groups::memcpy_async(block, positionsBuffer_coord, solventblock_neighbor->rel_pos, sizeof(Coord) * SolventBlock::MAX_SOLVENTS_IN_BLOCK);
        //cooperative_groups::memcpy_async(block, utility_buffer_small, solventblock_neighbor->atomtypeIds, sizeof(uint8_t) * SolventBlock::MAX_SOLVENTS_IN_BLOCK);
        //cooperative_groups::wait(block);
        //positionsBuffer_relpos[threadIdx.x] = positionsBuffer_coord[threadIdx.x].ToRelpos();

        if (threadIdx.x < nsolvents_neighbor) {
            positionsBuffer_relpos[threadIdx.x] = solventblock_neighbor->rel_pos[threadIdx.x].ToRelpos();
            utility_buffer_small[threadIdx.x] = solventblock_neighbor->atomtypeIds[threadIdx.x];
        }
		__syncthreads();


		if (threadActive) {
            force += LJ::computeSolventToSolventLJForces<true, energyMinimize, false>
				(relpos_self - nearbyBlock[i].relShift, tinymolTypeId, positionsBuffer_relpos, nsolvents_neighbor, potE_sum, *forcefieldTinymolShared, utility_buffer_small, nullptr);
		}
	}

    // Finally push force and potE for next kernel
    __syncthreads();
    forceEnergyOut[threadIdx.x] = ForceEnergy{ force, potE_sum };
    __syncthreads();

    auto block = cooperative_groups::this_thread_block();
    cooperative_groups::memcpy_async(block, &forceEnergies[blockIdx.x * SolventBlock::MAX_SOLVENTS_IN_BLOCK], forceEnergyOut, sizeof(ForceEnergy) * nElementsInBlock);
}




static_assert(BondgroupTinymol::maxSinglebonds <= BondgroupTinymol::maxParticles, "Not enough threads to load all singlebonds");
static_assert(BondgroupTinymol::maxAnglebonds <= BondgroupTinymol::maxParticles, "Not enough threads to load all anglebonds");
// Spawn 1 block per solventblock, blockDim(SOlventblock::MAX_BONDGROUPS, BondgroupTinymol::maxParticles)
template <bool emvariant>
__global__ void TinymolBondgroupsKernel(const SimulationDevice* const sim, const int16_t step, ForceEnergy* const forceEnergies) {
	//__shared__ ForceEnergy
	static_assert(sizeof(Coord) == sizeof(Float3), "Coord and Float3 must be the same size");
	__shared__ Float3 positions[SolventBlock::MAX_SOLVENTS_IN_BLOCK];

	const int solventblockId = blockIdx.x;

	SolventBlock* const solventblockGlobalPtr = SolventBlocksCircularQueue::getBlockPtr(sim->boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, solventblockId, step);
	const int nParticlesInBlock = solventblockGlobalPtr->nParticles;
	const int nBondgroupsInBlock = solventblockGlobalPtr->nBondgroups;

	if (nParticlesInBlock > SolventBlock::MAX_SOLVENTS_IN_BLOCK || nParticlesInBlock < 0 || nBondgroupsInBlock > SolventBlock::maxBondgroups || nBondgroupsInBlock < 0) {
		printf("Illegal count\n");		
	}

	// Load positions
	{
		Coord* positionsAsCoord = (Coord*)positions;
		auto block = cooperative_groups::this_thread_block();
		cooperative_groups::memcpy_async(block, positionsAsCoord, solventblockGlobalPtr->rel_pos, sizeof(Coord) * nParticlesInBlock);
		cooperative_groups::wait(block);

		if (threadIdx.y * blockDim.x + threadIdx.x < SolventBlock::MAX_SOLVENTS_IN_BLOCK) {
			positions[threadIdx.y * blockDim.x + threadIdx.x] = positionsAsCoord[threadIdx.y * blockDim.x + threadIdx.x].ToRelpos();
		}
	}

	__shared__ ForceEnergy forceEnergyInterrimsShared[SolventBlock::MAX_SOLVENTS_IN_BLOCK];	
	static_assert(SolventBlock::maxBondgroups * BondgroupTinymol::maxParticles >= SolventBlock::MAX_SOLVENTS_IN_BLOCK, "Not enough threads");
	if (threadIdx.y * blockDim.x + threadIdx.x < BondgroupTinymol::maxParticles) {
		forceEnergyInterrimsShared[threadIdx.y * blockDim.x + threadIdx.x] = ForceEnergy{};
	}

	// Singlebonds
	{
		BondgroupTinymol* bondgroupPtr = nullptr;
		int bondgroupsFirstAtomIndexInSolventblock = -1;
		float potE{};
		Float3 forces[SingleBond::nAtoms];
		SingleBond singlebond;

		if (threadIdx.x < nBondgroupsInBlock) {
			//bondgroupPtr = &bondgroupsBuffer[solventblockId * SolventBlock::maxBondgroups + threadIdx.x];
			bondgroupPtr = &solventblockGlobalPtr->bondgroups[threadIdx.x];
			bondgroupsFirstAtomIndexInSolventblock = solventblockGlobalPtr->bondgroupsFirstAtomindexInSolventblock[threadIdx.x];
			//printf("N singlebonds %d Indices %d %d\n", , singlebond.atom_indexes[0], singlebond.atom_indexes[1]);
			//printf("Abs indices %d %d\n", singlebond.atom_indexes[0] + bondgroupsFirstAtomIndexInSolventblock, singlebond.atom_indexes[1] + bondgroupsFirstAtomIndexInSolventblock);

			if (threadIdx.y < bondgroupPtr->nSinglebonds) {
				singlebond = bondgroupPtr->singlebonds[threadIdx.y];				
				/*printf("blockidx %d N singlebonds %d Indices %d %d\n", 
					blockIdx.x,
					bondgroupPtr->nSinglebonds, 
					singlebond.atom_indexes[0]+ bondgroupsFirstAtomIndexInSolventblock, 
					singlebond.atom_indexes[1]+ bondgroupsFirstAtomIndexInSolventblock
				);*/


				if constexpr (!LIMA_PUSH) {
					for (int i = 0; i < 2; i++)
						if (singlebond.atom_indexes[i] + bondgroupsFirstAtomIndexInSolventblock < 0 || singlebond.atom_indexes[i] + bondgroupsFirstAtomIndexInSolventblock >= nParticlesInBlock)
							printf("Illegal index %d\n", singlebond.atom_indexes[i]);
				}

				// Sets force and pot, not adding
				LimaForcecalc::calcSinglebondForces<emvariant>(
					positions[singlebond.atom_indexes[0] + bondgroupsFirstAtomIndexInSolventblock],
					positions[singlebond.atom_indexes[1] + bondgroupsFirstAtomIndexInSolventblock],
					singlebond.params, forces, potE, false);
			}
		}

		// We can safely assume there's no overlap between the bondgroups used particles, thus 1 thread per bondgroup can write to shared mem at a time
		for (int i = 0; i < BondgroupTinymol::maxParticles; i++) {
			if (threadIdx.x < nBondgroupsInBlock && threadIdx.y == i && threadIdx.y < bondgroupPtr->nSinglebonds) {
				forceEnergyInterrimsShared[singlebond.atom_indexes[0] + bondgroupsFirstAtomIndexInSolventblock] += ForceEnergy{ forces[0], potE * 0.5f };
				forceEnergyInterrimsShared[singlebond.atom_indexes[1] + bondgroupsFirstAtomIndexInSolventblock] += ForceEnergy{ forces[1], potE * 0.5f };
			}
			__syncthreads();
		}
	}

	// Anglebonds
	{
		//BondgroupTinymol* bondgroupPtr = nullptr;
		//int bondgroupsFirstAtomIndexInSolventblock = -1;
		//float potE{};
		//Float3 forces[AngleUreyBradleyBond::nAtoms];
		//AngleUreyBradleyBond anglebond;

		//if (threadIdx.x < nBondgroupsInBlock) {
		//	bondgroupPtr = &solventblockGlobalPtr->bondgroups[threadIdx.x];
		//	bondgroupsFirstAtomIndexInSolventblock = solventblockGlobalPtr->bondgroupsFirstAtomindexInSolventblock[threadIdx.x];

		//	if (threadIdx.y < bondgroupPtr->nAnglebonds) {
		//		anglebond = bondgroupPtr->anglebonds[threadIdx.y];

		//		if constexpr (!LIMA_PUSH) {
		//			for (int i = 0; i < 3; i++)
		//				if (anglebond.atom_indexes[i] < 0 || anglebond.atom_indexes[i] >= nParticlesInBlock)
		//					printf("Illegal index %d\n", anglebond.atom_indexes[i]);
		//		}

		//		// Sets force and pot, not adding
		//		LimaForcecalc::calcAnglebondForces(
		//			positions[anglebond.atom_indexes[0] + bondgroupsFirstAtomIndexInSolventblock],
		//			positions[anglebond.atom_indexes[1] + bondgroupsFirstAtomIndexInSolventblock],
		//			positions[anglebond.atom_indexes[2] + bondgroupsFirstAtomIndexInSolventblock],
		//			anglebond, forces, potE);
		//	}
//		}

		//// We can safely assume there's no overlap between the bondgroups used particles, thus 1 thread per bondgroup can write to shared mem at a time
		//for (int i = 0; i < BondgroupTinymol::maxParticles; i++) {
		//	if (threadIdx.x < nBondgroupsInBlock && threadIdx.y == i && threadIdx.y < bondgroupPtr->nAnglebonds) {
		//		forceEnergyInterrimsShared[anglebond.atom_indexes[0] + bondgroupsFirstAtomIndexInSolventblock] += ForceEnergy{ forces[0], potE / 3.f }; // OPTIM mul with 0.333?
		//		forceEnergyInterrimsShared[anglebond.atom_indexes[1] + bondgroupsFirstAtomIndexInSolventblock] += ForceEnergy{ forces[1], potE / 3.f };
		//		forceEnergyInterrimsShared[anglebond.atom_indexes[2] + bondgroupsFirstAtomIndexInSolventblock] += ForceEnergy{ forces[2], potE / 3.f };
		//	}
		//	__syncthreads();
		//}
	}

	//if (threadIdx.y == 0 && threadIdx.x < nParticlesInBlock) {
	//	printf("\nid %d force %f %f %f\n", threadIdx.x, forceEnergyInterrimsShared[threadIdx.x].force.x, forceEnergyInterrimsShared[threadIdx.x].force.y, forceEnergyInterrimsShared[threadIdx.x].force.z);
	//}


	// Write forceenergy to global mem
	{
		auto block = cooperative_groups::this_thread_block();
		cooperative_groups::memcpy_async(block, &forceEnergies[solventblockId * SolventBlock::MAX_SOLVENTS_IN_BLOCK], forceEnergyInterrimsShared, sizeof(ForceEnergy) * nParticlesInBlock);
	}	
}

template <typename BoundaryCondition, bool energyMinimize>
__global__ void TinymolIntegrateAndLogKernel(SimulationDevice* sim, int64_t step, const ForceEnergyInterims forceEnergies) {
	__shared__ SolventBlock solventblock;
	__shared__ uint8_t utility_buffer_small[SolventBlock::MAX_SOLVENTS_IN_BLOCK];

	__shared__ Coord relPositionsNext[SolventBlock::MAX_SOLVENTS_IN_BLOCK];


	// Doubles as block_index_3d!
	const NodeIndex block_origo = BoxGrid::Get3dIndex(blockIdx.x, DeviceConstants::boxSize.boxSizeNM_i);

	const BoxState& boxState = sim->boxState;
	const SimParams& simparams = sim->params;
	SolventBlock* solventblock_ptr = SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step);

	const ForceEnergy myForceEnergy = forceEnergies.forceEnergiesCompoundinteractions[blockIdx.x * SolventBlock::MAX_SOLVENTS_IN_BLOCK + threadIdx.x]
		+ forceEnergies.forceEnergiesTinymolinteractions[blockIdx.x * SolventBlock::MAX_SOLVENTS_IN_BLOCK + threadIdx.x]
		+ forceEnergies.forceEnergiesTinymolBondgroups[blockIdx.x * SolventBlock::MAX_SOLVENTS_IN_BLOCK + threadIdx.x];

	//const ForceEnergy myForceEnergy = forceEnergies.forceEnergiesTinymolBondgroups[blockIdx.x * SolventBlock::MAX_SOLVENTS_IN_BLOCK + threadIdx.x];

	if (threadIdx.x == 0) {
		solventblock.loadMeta(*solventblock_ptr);
	}
	__syncthreads();
	const bool solventActive = threadIdx.x < solventblock.nParticles;
	solventblock.loadData(*solventblock_ptr);
	__syncthreads();



	//Coord relpos_next{};
	if (solventActive) {
		TinyMolParticleState& tinyMols_ref = boxState.tinyMolParticlesState[solventblock.ids[threadIdx.x]];	// TinyMolParticleState private data, for VVS
		const float mass = DeviceConstants::tinymolForcefield.types[tinyMols_ref.tinymolTypeIndex].mass;

		if constexpr (energyMinimize) {
			const Float3 safeForce = EngineUtils::ForceActivationFunction(myForceEnergy.force);
			AdamState* const adamStatePtr = &sim->adamState[sim->boxparams.n_compounds * MAX_COMPOUND_PARTICLES + solventblock.ids[threadIdx.x]];
			const Coord pos_now = EngineUtils::IntegratePositionADAM(solventblock.rel_pos[threadIdx.x], safeForce, adamStatePtr, step);

			relPositionsNext[threadIdx.x] = pos_now;
			EngineUtils::LogSolventData(sim->boxparams, myForceEnergy.potE, block_origo, solventblock.ids[threadIdx.x], solventblock.rel_pos[threadIdx.x], solventActive,
				myForceEnergy.force, Float3{}, step, sim->potE_buffer, sim->traj_buffer, sim->vel_buffer, simparams.data_logging_interval);
		}
		else {
            if constexpr (!LIMA_PUSH) {
                if (myForceEnergy.force.isNan() || myForceEnergy.force.lenSquared() >= FLT_MAX)
					myForceEnergy.force.print('S');
            }

			Float3 vel_now = EngineUtils::integrateVelocityVVS(tinyMols_ref.vel_prev, tinyMols_ref.force_prev, myForceEnergy.force, simparams.dt, mass);
			const Coord pos_now = EngineUtils::integratePositionVVS(solventblock.rel_pos[threadIdx.x], vel_now, myForceEnergy.force, mass, simparams.dt);

			vel_now = vel_now * DeviceConstants::thermostatScalar;

			tinyMols_ref.vel_prev = vel_now;
			tinyMols_ref.force_prev = myForceEnergy.force;

			// Save pos locally, but only push to box as this kernel ends
			relPositionsNext[threadIdx.x] = pos_now;
			EngineUtils::LogSolventData(sim->boxparams, myForceEnergy.potE, block_origo, solventblock.ids[threadIdx.x], solventblock.rel_pos[threadIdx.x], solventActive,
				myForceEnergy.force, vel_now, step, sim->potE_buffer, sim->traj_buffer, sim->vel_buffer, simparams.data_logging_interval);
		}
	}

	// TODO: LONG: stop using a circular queue for solvents, just have the data 1 place now that we sync before integration anyways

	// Push new SolventCoord to global mem
	SolventBlock* const solventblock_next_ptr = SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step + 1);
	solventblock_next_ptr->rel_pos[threadIdx.x] = relPositionsNext[threadIdx.x];
	solventblock_next_ptr->ids[threadIdx.x] = solventblock.ids[threadIdx.x];
	solventblock_next_ptr->atomtypeIds[threadIdx.x] = solventblock.atomtypeIds[threadIdx.x];
	if (threadIdx.x == 0) {
		solventblock_next_ptr->nParticles = solventblock.nParticles;
		solventblock_next_ptr->nBondgroups = solventblock.nBondgroups;
	}
	if (threadIdx.x < solventblock.nBondgroups) {
		solventblock_next_ptr->bondgroups[threadIdx.x] = solventblock.bondgroups[threadIdx.x];
		solventblock_next_ptr->bondgroupsFirstAtomindexInSolventblock[threadIdx.x] = solventblock.bondgroupsFirstAtomindexInSolventblock[threadIdx.x];
	}
	
}

template <typename BoundaryCondition>
__global__ void SolventPretransferKernel(SimulationDevice* sim, int64_t _step, const TinymolTransferModule tinymolTransferModule) {
	const NodeIndex directions[6]{
		{1, 0, 0},
		{-1, 0, 0},
		{0, 1, 0},
		{0, -1, 0},
		{0, 0, 1},
		{0, 0, -1}
	};

	__shared__ int directionIndexOfBondgroups[SolventBlock::maxBondgroups]; // -1 for stay
	if (threadIdx.x < SolventBlock::maxBondgroups) {
		directionIndexOfBondgroups[threadIdx.x] = -1;
	}

	const int solventblockId = blockIdx.x;
	const int stepToLoadFrom = _step + 1;
	SolventBlock* const solventblockGlobalPtr = SolventBlocksCircularQueue::getBlockPtr(sim->boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, solventblockId, stepToLoadFrom);
	const int nBondgroupsInBlock = solventblockGlobalPtr->nBondgroups;

	__shared__ int nParticlesInBondgroups[SolventBlock::maxBondgroups];
	if (threadIdx.x < nBondgroupsInBlock) {
		nParticlesInBondgroups[threadIdx.x] = solventblockGlobalPtr->bondgroups[threadIdx.x].nParticles;
	}

	__shared__ int bondgroupsIndexOfFirstParticleInSolventblock[SolventBlock::maxBondgroups];
	if (threadIdx.x < nBondgroupsInBlock) {
		bondgroupsIndexOfFirstParticleInSolventblock[threadIdx.x] = solventblockGlobalPtr->bondgroupsFirstAtomindexInSolventblock[threadIdx.x];
	}

	// Each thread is responsible for where it's bondgroup go
	if (threadIdx.x < nBondgroupsInBlock) {
		const int indexInBlockOfFirstParticle = bondgroupsIndexOfFirstParticleInSolventblock[threadIdx.x];
		const NodeIndex direction = LIMAPOSITIONSYSTEM::getTransferDirection(solventblockGlobalPtr->rel_pos[indexInBlockOfFirstParticle]);
		for (int directionIndex = 0; directionIndex < 6; directionIndex++) {
			if (direction == directions[directionIndex]) {
				directionIndexOfBondgroups[threadIdx.x] = directionIndex;
				break;
			}
		}
	}
	__syncthreads;

	// First 6 threads are responsible for marking a direction
	__shared__ int nBondgroupsThisDirection[6];
	__shared__ int nParticlesThisDirection[6];
	__shared__ int bondgroupIdsThisDirection[TinymolTransferModule::maxOutgoingBondgroups * 6];
	__shared__ int bondgroupsFirstAtomindexInThisDirection[TinymolTransferModule::maxOutgoingBondgroups * 6];
	if (threadIdx.x < 6) {
		int myBondgroupCount = 0;
		int myParticleCount = 0;
		for (int bondgroupIndex = 0; bondgroupIndex < nBondgroupsInBlock; bondgroupIndex++) {
			if (directionIndexOfBondgroups[bondgroupIndex] == threadIdx.x) {
				bondgroupIdsThisDirection[threadIdx.x * TinymolTransferModule::maxOutgoingBondgroups + myBondgroupCount] = bondgroupIndex;
				bondgroupsFirstAtomindexInThisDirection[threadIdx.x * TinymolTransferModule::maxOutgoingBondgroups + myBondgroupCount] = myParticleCount;
				myBondgroupCount++;
				myParticleCount += nParticlesInBondgroups[bondgroupIndex];
			}
		}

		if (myBondgroupCount > TinymolTransferModule::maxOutgoingBondgroups)
			printf("Too many bondgroups in one direction");

		nBondgroupsThisDirection[threadIdx.x] = myBondgroupCount;
		nParticlesThisDirection[threadIdx.x] = myParticleCount;
	}
	__syncthreads();

	// Now all threads loop over the direction, and if they have a particle, they push it directy to the incoming queue in global memory
	for (int directionIndex = 0; directionIndex < 6; directionIndex++) {
		const NodeIndex direction = directions[directionIndex];
		const NodeIndex blockOrigo = BoxGrid::Get3dIndex(blockIdx.x, DeviceConstants::boxSize.boxSizeNM_i);
		const NodeIndex targetBlock = BoundaryCondition::applyBC(blockOrigo + direction, DeviceConstants::boxSize.blocksPerDim);
		const int targetBlockId = BoxGrid::Get1dIndex(targetBlock, DeviceConstants::boxSize.boxSizeNM_i);
		const Coord relposShift = Coord{ -direction.toFloat3()};

		if (targetBlockId >= DeviceConstants::boxSize.blocksPerDim * DeviceConstants::boxSize.blocksPerDim * DeviceConstants::boxSize.blocksPerDim)
			printf("Target block was out of bounds");

		// Write results directly to global mem
		if (threadIdx.x == 0) {
			tinymolTransferModule.nIncomingParticles[targetBlockId * 6 + directionIndex] = nParticlesThisDirection[directionIndex];
			tinymolTransferModule.nIncomingBondgroups[targetBlockId * 6 + directionIndex] = nBondgroupsThisDirection[directionIndex];
		}
		if (threadIdx.x < nBondgroupsThisDirection[directionIndex]) {
			const int bondgroupIndex = bondgroupIdsThisDirection[directionIndex * TinymolTransferModule::maxOutgoingBondgroups + threadIdx.x];
			const int bondgroupTargetIndex = targetBlockId * 6 * TinymolTransferModule::maxOutgoingBondgroups
				+ directionIndex * TinymolTransferModule::maxOutgoingBondgroups
				+ threadIdx.x;
			tinymolTransferModule.incomingBondgroups[bondgroupTargetIndex] = solventblockGlobalPtr->bondgroups[bondgroupIndex];
			tinymolTransferModule.incomingBondgroupsParticlesOffset[bondgroupTargetIndex] = bondgroupsFirstAtomindexInThisDirection[directionIndex * TinymolTransferModule::maxOutgoingBondgroups + threadIdx.x];
			
			const int destOffsetOfFirstParticle = bondgroupsFirstAtomindexInThisDirection[directionIndex * TinymolTransferModule::maxOutgoingBondgroups + threadIdx.x];
			const int sourceOffsetOfFirstParticle = solventblockGlobalPtr->bondgroupsFirstAtomindexInSolventblock[bondgroupIndex];
			for (int particleIndexInBondgroup = 0; particleIndexInBondgroup < nParticlesInBondgroups[bondgroupIndex]; particleIndexInBondgroup++) {
				const int particleDestIndex = targetBlockId * 6 * TinymolTransferModule::maxOutgoingParticles
					+ directionIndex * TinymolTransferModule::maxOutgoingParticles
					+ destOffsetOfFirstParticle + particleIndexInBondgroup;
				const int particleSourceIndex = sourceOffsetOfFirstParticle + particleIndexInBondgroup;

				tinymolTransferModule.incomingPositions[particleDestIndex] = solventblockGlobalPtr->rel_pos[particleSourceIndex] + relposShift;
				tinymolTransferModule.incomingIds[particleDestIndex] = solventblockGlobalPtr->ids[particleSourceIndex];
				tinymolTransferModule.incomingAtomtypeIds[particleDestIndex] = solventblockGlobalPtr->atomtypeIds[particleSourceIndex];
				tinymolTransferModule.incomingBondgroupIds[particleDestIndex] = threadIdx.x;
			}
		}
	}
	__syncthreads();


	// Finally compress remainders. Reuse the direction-of-particle buffer as prefixsum buffer
	__shared__ int bondgroupRemains[SolventBlock::maxBondgroups];
	__shared__ int nBondgroupsRemaining;
	__shared__ int nParticlesRemaining;
	if (threadIdx.x < SolventBlock::maxBondgroups) {
		bondgroupRemains[threadIdx.x] = directionIndexOfBondgroups[threadIdx.x] == -1;
	}
	
	__syncthreads();
	if (threadIdx.x == 0) {
		nBondgroupsRemaining = 0;
		nParticlesRemaining = 0;

		for (int bondgroupIndex = 0; bondgroupIndex < nBondgroupsInBlock; bondgroupIndex++) {
			if (bondgroupRemains[bondgroupIndex]) {


				solventblockGlobalPtr->bondgroups[nBondgroupsRemaining] = solventblockGlobalPtr->bondgroups[bondgroupIndex];
				solventblockGlobalPtr->bondgroupsFirstAtomindexInSolventblock[nBondgroupsRemaining] = nParticlesRemaining;
				nBondgroupsRemaining++;

				for (int particleIndex = 0; particleIndex < nParticlesInBondgroups[bondgroupIndex]; particleIndex++) {
					const int srcIndex = bondgroupsIndexOfFirstParticleInSolventblock[bondgroupIndex] + particleIndex;
					solventblockGlobalPtr->rel_pos[nParticlesRemaining] = solventblockGlobalPtr->rel_pos[srcIndex];
					solventblockGlobalPtr->ids[nParticlesRemaining] = solventblockGlobalPtr->ids[srcIndex];
					solventblockGlobalPtr->atomtypeIds[nParticlesRemaining] = solventblockGlobalPtr->atomtypeIds[srcIndex];
					solventblockGlobalPtr->particlesBondgroupIds[nParticlesRemaining] = nBondgroupsRemaining - 1;
					nParticlesRemaining++;
				}
			}
		}

		solventblockGlobalPtr->nBondgroups = nBondgroupsRemaining;
		solventblockGlobalPtr->nParticles = nParticlesRemaining;
	}
}

__global__ void SolventTransferKernel(SimulationDevice* sim, int64_t _step, const TinymolTransferModule tinymolTransferModule)  {	
	__shared__ int nParticlesInBlock;
	__shared__ int nBondgroupsInBlock;

	const int solventblockId = blockIdx.x;
	const int stepToLoadFrom = _step + 1;
	SolventBlock* const solventblockGlobalPtr = SolventBlocksCircularQueue::getBlockPtr(sim->boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, solventblockId, stepToLoadFrom);

	if (threadIdx.x == 0) {
		nParticlesInBlock = solventblockGlobalPtr->nParticles;
		nBondgroupsInBlock = solventblockGlobalPtr->nBondgroups;
	}
	__syncthreads;

	// All threads loop over the directions. If it has an incoming particle, it copies it directy to global mem
	for (int directionIndex = 0; directionIndex < 6; directionIndex++) {
		const int nIncomingParticles = tinymolTransferModule.nIncomingParticles[blockIdx.x * 6 + directionIndex];
		const int nIncomingBondgroups = tinymolTransferModule.nIncomingBondgroups[blockIdx.x * 6 + directionIndex];
		
		if (threadIdx.x < nIncomingParticles) {
			const int srcIndex = solventblockId * 6 * TinymolTransferModule::maxOutgoingParticles
				+ directionIndex * TinymolTransferModule::maxOutgoingParticles
				+ threadIdx.x;
			solventblockGlobalPtr->rel_pos[nParticlesInBlock + threadIdx.x] = tinymolTransferModule.incomingPositions[srcIndex];
			solventblockGlobalPtr->ids[nParticlesInBlock + threadIdx.x] = tinymolTransferModule.incomingIds[srcIndex];
			solventblockGlobalPtr->atomtypeIds[nParticlesInBlock + threadIdx.x] = tinymolTransferModule.incomingAtomtypeIds[srcIndex];
			solventblockGlobalPtr->particlesBondgroupIds[nParticlesInBlock + threadIdx.x] = tinymolTransferModule.incomingBondgroupIds[srcIndex] + nBondgroupsInBlock;
		}
		if (threadIdx.x < nIncomingBondgroups) {
			const int srcIndex = solventblockId * 6 * TinymolTransferModule::maxOutgoingBondgroups
				+ directionIndex * TinymolTransferModule::maxOutgoingBondgroups
				+ threadIdx.x;
			solventblockGlobalPtr->bondgroups[nBondgroupsInBlock + threadIdx.x] = tinymolTransferModule.incomingBondgroups[srcIndex];
			solventblockGlobalPtr->bondgroupsFirstAtomindexInSolventblock[nBondgroupsInBlock + threadIdx.x] = nParticlesInBlock + tinymolTransferModule.incomingBondgroupsParticlesOffset[srcIndex];
		}

		__syncthreads();
		if (threadIdx.x == 0) {
			nParticlesInBlock += nIncomingParticles;
			nBondgroupsInBlock += nIncomingBondgroups;
		}
		__syncthreads();
	}

	// Finally we write to global mem how many particles there are in the block now
	if (threadIdx.x == 0) {
		solventblockGlobalPtr->nParticles = nParticlesInBlock;
		solventblockGlobalPtr->nBondgroups = nBondgroupsInBlock;
	}
}

//
//
//// This is run before step.inc(), but will always publish results to the first array in grid!
//__global__ void SolventTransferKernel(SimulationDevice* sim, int64_t step) {
//	const BoxState& boxState = sim->boxState;
//
//	SolventBlockTransfermodule* transfermodule = &sim->transfermodule_array[blockIdx.x];
//	
//	SolventBlock* solventblock_current = SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step);
//	SolventBlock* solventblock_next = SolventBlocksCircularQueue::getBlockPtr(boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, blockIdx.x, step + 1);
//
//	SolventTransferWarnings::assertSolventsEqualNRemain(*solventblock_next, *transfermodule);
//
//	// Handling incoming transferring solvents
//	int n_solvents_next = transfermodule->n_remain;
//	for (int queue_index = 0; queue_index < SolventBlockTransfermodule::n_queues; queue_index++) {
//		auto* queue = &transfermodule->transfer_queues[queue_index];
//		if (threadIdx.x < queue->n_elements) {
//			const int incoming_index = n_solvents_next + threadIdx.x;
//
//			solventblock_next->rel_pos[incoming_index] = queue->rel_positions[threadIdx.x];
//			solventblock_next->ids[incoming_index] = queue->ids[threadIdx.x];
//			solventblock_next->atomtypeIds[incoming_index] = queue->atomtypeIds[threadIdx.x];
//		}
//		n_solvents_next += queue->n_elements;
//
//		// Signal that all elements of the queues have been moved
//		__syncthreads();
//		if (threadIdx.x == 0) {
//			queue->n_elements = 0;
//		}
//	}
//
//	SolventTransferWarnings::assertMaxPlacedSolventsIsWithinLimits(n_solvents_next, sim->signals->critical_error_encountered);
//
//	// Finally update the solventblock_next with how many solvents it now contains
//	if (threadIdx.x == 0) {
//		solventblock_next->n_solvents = n_solvents_next;
//	}
//}


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
		const Float3 relShift = LIMAPOSITIONSYSTEM_HACK::GetRelShiftFromOrigoShift_Float3(myNodeindex, origo);

		positions[threadIdx.x] = boxState.compoundsRelposNm[pRef.compoundId * MAX_COMPOUND_PARTICLES + pRef.localIdInCompound] + relShift;
		
	}


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
