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


#include "Electrostatics.cuh"


#pragma warning(push)
#pragma warning(disable:E0020)
#pragma warning(push)
#pragma warning(disable: 20054)

#pragma diag_suppress 20054



#include "LennardJonesInteractions.cuh"

// TODO: move to engine utils
template <typename BondType, int max_bondtype_in_compound>
__device__ BondType* LoadBonds(char* utility_buffer, const BondType* const source, int nBondsToLoad) {
	static_assert(cbkernel_utilitybuffer_size >= sizeof(BondType) * max_bondtype_in_compound, "Utilitybuffer not large enough for bondtype");
	BondType* bonds = (BondType*)utility_buffer;

	auto block = cooperative_groups::this_thread_block();
	cooperative_groups::memcpy_async(block, bonds, source, sizeof(BondType) * nBondsToLoad);
	cooperative_groups::wait(block);

	return bonds;
}















// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//
#define compound_index blockIdx.x
template <typename BoundaryCondition, bool energyMinimize, bool computePotE> // We dont compute potE if we dont log data this step
__global__ void compoundFarneighborShortrangeInteractionsKernel(SimulationDevice* sim, const int64_t step) {
	__shared__ CompoundCompact compound;				// Mostly bond information
	__shared__ Float3 compound_positions[MAX_COMPOUND_PARTICLES]; // [lm]

	//Neighborlist
	__shared__ int nNonbondedCompoundNeighbors;
	__shared__ int nGridnodes;

	__shared__ Float3 utility_buffer_f3[MAX_COMPOUND_PARTICLES*2];
	__shared__ Float3 utility_float3;

	__shared__ half particleChargesBuffers[MAX_COMPOUND_PARTICLES * 2];

	__shared__ ForceField_NB forcefield_shared;

	BoxState* const boxState = sim->boxState;
	const BoxConfig& boxConfig = sim->boxConfig;
	const SimParams& simparams = sim->params;
	SimSignals* signals = sim->signals;

	float potE_sum{};
	Float3 force{};
	NodeIndex compound_origo{}; // TODO: Figure out why TF everything breaks if i make this __shared__???????
	const float particleCharge = simparams.enable_electrostatics	// TODO: this is temporary
		? static_cast<float>(boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x])
		: 0.f;

	{
		auto block = cooperative_groups::this_thread_block();

		const CompoundCoords* const compoundcoords_global = CompoundcoordsCircularQueueUtils::getCoordarrayRef(boxState->compoundcoordsCircularQueue, step, blockIdx.x);
		compound_origo = compoundcoords_global->origo;
		cooperative_groups::memcpy_async(block, (Coord*)compound_positions, compoundcoords_global->rel_positions, sizeof(Coord) * MAX_COMPOUND_PARTICLES);
			
		if (threadIdx.x == 0) {
			compound.loadMeta(&boxConfig.compounds[blockIdx.x]);
			nNonbondedCompoundNeighbors = sim->compound_neighborlists[blockIdx.x].nNonbondedNeighbors;
			nGridnodes = sim->compound_neighborlists[blockIdx.x].n_gridnodes;
		}

		cooperative_groups::memcpy_async(block, &forcefield_shared, &forcefield_device, sizeof(ForceField_NB));

		cooperative_groups::wait(block);
		compound_positions[threadIdx.x] = ((Coord*)compound_positions)[threadIdx.x].toFloat3();
	}
	compound.loadData(&boxConfig.compounds[blockIdx.x]);




	// For neighborcompound and neighbor solventblocks we utilize a batching system
	const int batchsize = 64;
	__shared__ Float3 relshifts[batchsize];	// [lm]
	__shared__ int neighborIds[batchsize]; // either compoundID or solventblockID
	__shared__ int neighborNParticles[batchsize]; // either particlesInCompound or particlesInSolventblock
	__shared__ void* neighborPtrs[batchsize]; // Either CompoundCoords* or SolventBlock*
	// --------------------------------------------------------------- Intercompound forces --------------------------------------------------------------- //
	{
		__shared__ const BondedParticlesLUT* compoundPairLUTs[batchsize];		
		__shared__ uint8_t atomtypesBuffer[MAX_COMPOUND_PARTICLES*2];
		auto block = cooperative_groups::this_thread_block();		

		// This part is scary, but it also takes up by far the majority of compute time. We use the utilitybuffer twice simultaneously, so be careful when making changes
		int indexInBatch = batchsize;
		for (int i = 0; i < nNonbondedCompoundNeighbors; i++) {
			__syncthreads();

			static_assert(sizeof utility_buffer_f3 >= sizeof(Float3) * MAX_COMPOUND_PARTICLES * 2, "Utilitybuffer not large enough for neighbor positions");
			Float3* neighborPositionsCurrent = ((Float3*)utility_buffer_f3) + (MAX_COMPOUND_PARTICLES * (i & 1));
			Float3* neighborPositionsNext = ((Float3*)utility_buffer_f3) + (MAX_COMPOUND_PARTICLES * !(i & 1));
			uint8_t* neighborAtomstypesCurrent = atomtypesBuffer + (MAX_COMPOUND_PARTICLES * (i & 1));
			uint8_t* neighborAtomstypesNext = atomtypesBuffer + (MAX_COMPOUND_PARTICLES * !(i & 1));
			half* neighborParticleschargesCurrent = particleChargesBuffers + (MAX_COMPOUND_PARTICLES * (i & 1));
			half* neighborParticleschargesNext = particleChargesBuffers + (MAX_COMPOUND_PARTICLES * !(i & 1));

			static_assert(MAX_COMPOUND_PARTICLES >= batchsize);
			// First check if we need to load a new batch of relshifts & n_particles for the coming 32 compounds
			if (indexInBatch == batchsize) {
				if (threadIdx.x < batchsize && threadIdx.x + i < nNonbondedCompoundNeighbors) {
					neighborIds[threadIdx.x] = sim->compound_neighborlists[blockIdx.x].nonbondedNeighborcompoundIds[i + threadIdx.x];

					neighborPtrs[threadIdx.x] = (void*)CompoundcoordsCircularQueueUtils::getCoordarrayRef(boxState->compoundcoordsCircularQueue, step, neighborIds[threadIdx.x]);
					const CompoundCoords* const querycompound = CompoundcoordsCircularQueueUtils::getCoordarrayRef(boxState->compoundcoordsCircularQueue, step, neighborIds[threadIdx.x]);
					const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(compound_origo, ((CompoundCoords*)neighborPtrs[threadIdx.x])->origo);
					KernelHelpersWarnings::assertHyperorigoIsValid(querycompound_hyperorigo, compound_origo);

					// calc Relative LimaPosition Shift from the origo-shift
					relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, compound_origo).toFloat3();
					neighborNParticles[threadIdx.x] = boxConfig.compounds[neighborIds[threadIdx.x]].n_particles;
					
					compoundPairLUTs[threadIdx.x] = BondedParticlesLUTHelpers::get(sim->boxConfig.bpLUTs, compound_index, neighborIds[threadIdx.x]);
					
				}
				indexInBatch = 0;
				__syncthreads();

				{
					// Load first element in batch and sync				
					const int currentNeighborId = neighborIds[indexInBatch];
					const int currentNeighborNParticles = neighborNParticles[indexInBatch];
					
					cooperative_groups::memcpy_async(block, neighborAtomstypesCurrent, boxConfig.compoundsAtomtypes + currentNeighborId * MAX_COMPOUND_PARTICLES, sizeof(uint8_t) * MAX_COMPOUND_PARTICLES);
					cooperative_groups::memcpy_async(block, (Coord*)neighborPositionsCurrent, ((CompoundCoords*)neighborPtrs[indexInBatch])->rel_positions, sizeof(Coord) * currentNeighborNParticles);
					cooperative_groups::memcpy_async(block, neighborParticleschargesCurrent, boxConfig.compoundsAtomCharges + currentNeighborId * MAX_COMPOUND_PARTICLES, sizeof(half) * currentNeighborNParticles);
					cooperative_groups::wait(block); // Joins all threads, waits for all copies to complete		

					const Coord queryparticle_coord = ((Coord*)neighborPositionsCurrent)[threadIdx.x];
					neighborPositionsCurrent[threadIdx.x] = queryparticle_coord.toFloat3() + relshifts[indexInBatch];
					__syncthreads();
				}
			}

			if (i + 1 < nNonbondedCompoundNeighbors && indexInBatch+1 < batchsize) {
				static_assert(sizeof(Coord) == sizeof(Float3));
				const int nextNeighborId = neighborIds[indexInBatch + 1];
				const int nextNeighborNParticles = neighborNParticles[indexInBatch + 1];

				cooperative_groups::memcpy_async(block, neighborAtomstypesNext, boxConfig.compoundsAtomtypes + nextNeighborId * MAX_COMPOUND_PARTICLES, sizeof(uint8_t)* MAX_COMPOUND_PARTICLES);
				cooperative_groups::memcpy_async(block, (Coord*)neighborPositionsNext, ((CompoundCoords*)neighborPtrs[indexInBatch+1])->rel_positions, sizeof(Coord) * nextNeighborNParticles);
				cooperative_groups::memcpy_async(block, neighborParticleschargesNext, boxConfig.compoundsAtomCharges + nextNeighborId * MAX_COMPOUND_PARTICLES, sizeof(half) * nextNeighborNParticles);
			}
			
			if (threadIdx.x < compound.n_particles) {
				force += LJ::computeCompoundCompoundLJForces<computePotE>(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum,
					neighborPositionsCurrent, neighborNParticles[indexInBatch], neighborAtomstypesCurrent, forcefield_shared, particleCharge, neighborParticleschargesCurrent);
			}

			cooperative_groups::wait(block); // Joins all threads, waits for all copies to complete		
			// Process the positions
			if (indexInBatch + 1 < batchsize) {
				const Coord queryparticle_coord = ((Coord*)neighborPositionsNext)[threadIdx.x];
				neighborPositionsNext[threadIdx.x] = queryparticle_coord.toFloat3() + relshifts[indexInBatch + 1];
				__syncthreads();
			}	
			indexInBatch++;
		}
	}
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //



	// This is the first kernel, so we overwrite
	if (threadIdx.x < compound.n_particles || true) { // TEMP
		if constexpr (computePotE) {
			sim->boxState->compoundsInterimState[blockIdx.x].potE_interim[threadIdx.x] = potE_sum;
		}
		sim->boxState->compoundsInterimState[blockIdx.x].forces_interim[threadIdx.x] = force;
	}
}
template  __global__ void compoundFarneighborShortrangeInteractionsKernel<PeriodicBoundaryCondition, true, true>(SimulationDevice* sim, int64_t step);
template  __global__ void compoundFarneighborShortrangeInteractionsKernel<PeriodicBoundaryCondition, true, false>(SimulationDevice* sim, int64_t step);
template  __global__ void compoundFarneighborShortrangeInteractionsKernel<PeriodicBoundaryCondition, false, true>(SimulationDevice* sim, int64_t step);
template  __global__ void compoundFarneighborShortrangeInteractionsKernel<PeriodicBoundaryCondition, false, false>(SimulationDevice* sim, int64_t step);
template __global__ void compoundFarneighborShortrangeInteractionsKernel<NoBoundaryCondition, true, true>(SimulationDevice* sim, int64_t step);
template __global__ void compoundFarneighborShortrangeInteractionsKernel<NoBoundaryCondition, true, false>(SimulationDevice* sim, int64_t step);
template __global__ void compoundFarneighborShortrangeInteractionsKernel<NoBoundaryCondition, false, true>(SimulationDevice* sim, int64_t step);
template __global__ void compoundFarneighborShortrangeInteractionsKernel<NoBoundaryCondition, false, false>(SimulationDevice* sim, int64_t step);
#undef compound_index


#define compound_index blockIdx.x
template <typename BoundaryCondition, bool energyMinimize, bool computePotE> // We dont compute potE if we dont log data this step
__global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel(SimulationDevice* sim, const int64_t step) {
	__shared__ CompoundCompact compound;				// Mostly bond information
	__shared__ Float3 compound_positions[MAX_COMPOUND_PARTICLES]; // [lm]
	//Neighborlist
	//__shared__ int nCompoundNeighbors;
	__shared__ int nGridnodes;

	__shared__ Float3 utility_buffer_f3[MAX_COMPOUND_PARTICLES * 2];
	__shared__ Float3 utility_float3;

	__shared__ BondedParticlesLUT bpLUT;
	__shared__ half particleChargesBuffers[MAX_COMPOUND_PARTICLES * 2];

	__shared__ ForceField_NB forcefield_shared;

	BoxState* const boxState = sim->boxState;
	const BoxConfig& boxConfig = sim->boxConfig;
	const SimParams& simparams = sim->params;
	SimSignals* signals = sim->signals;

	float potE_sum{};
	Float3 force{};
	NodeIndex compound_origo{}; // TODO: Figure out why TF everything breaks if i make this __shared__???????
	const float particleCharge = simparams.enable_electrostatics	// TODO: this is temporary
		? static_cast<float>(boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x])
		: 0.f;


	{
		auto block = cooperative_groups::this_thread_block();

		const CompoundCoords* const compoundcoords_global = CompoundcoordsCircularQueueUtils::getCoordarrayRef(boxState->compoundcoordsCircularQueue, step, blockIdx.x);
		compound_origo = compoundcoords_global->origo;
		cooperative_groups::memcpy_async(block, (Coord*)compound_positions, compoundcoords_global->rel_positions, sizeof(Coord) * MAX_COMPOUND_PARTICLES);

		if (threadIdx.x == 0) {
			compound.loadMeta(&boxConfig.compounds[blockIdx.x]);
			nGridnodes = sim->compound_neighborlists[blockIdx.x].n_gridnodes;
		}

		cooperative_groups::memcpy_async(block, &forcefield_shared, &forcefield_device, sizeof(ForceField_NB));

		cooperative_groups::wait(block);
		compound_positions[threadIdx.x] = ((Coord*)compound_positions)[threadIdx.x].toFloat3();
	}
	compound.loadData(&boxConfig.compounds[blockIdx.x]);



	// ------------------------------------------------------------ Intracompound Operations ------------------------------------------------------------ //
	{
		static_assert(clj_utilitybuffer_bytes >= sizeof(BondedParticlesLUT), "Utilitybuffer not large enough for BondedParticlesLUT");

		const BondedParticlesLUT* const bplut_global = BondedParticlesLUTHelpers::get(sim->boxConfig.bpLUTs, compound_index, compound_index);
		bpLUT.load(*bplut_global);	// A lut always exists within a compound

		static_assert(clj_utilitybuffer_bytes >= sizeof(BondedParticlesLUT) + sizeof(half) * MAX_COMPOUND_PARTICLES,
			"Utilitybuffer not large enough for neighbor charges");
		half* particleChargesCompound = &particleChargesBuffers[0];
		particleChargesCompound[threadIdx.x] = particleCharge;
		__syncthreads();

		if (threadIdx.x < compound.n_particles) {
			// Having this inside vs outside the context makes impact the resulting VC, but it REALLY SHOULD NOT
			force += LJ::computeCompoundCompoundLJForces<computePotE>(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum, compound_positions, compound.n_particles,
				compound.atom_types, &bpLUT, LJ::CalcLJOrigin::ComComIntra, forcefield_shared,
				particleCharge, particleChargesCompound);
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //



	// For neighborcompound and neighbor solventblocks we utilize a batching system
	const int batchsize = 64;
	__shared__ Float3 relshifts[batchsize];	// [lm]
	__shared__ int neighborIds[batchsize]; // either compoundID or solventblockID // should be uint16_t? Does it make a diff?
	__shared__ int neighborNParticles[batchsize]; // either particlesInCompound or particlesInSolventblock
	__shared__ void* neighborPtrs[batchsize]; // Either CompoundCoords* or SolventBlock*

	// --------------------------------------------------------------- Intercompound forces --------------------------------------------------------------- //
	{
		__shared__ uint8_t neighborAtomstypes[MAX_COMPOUND_PARTICLES];
		__shared__ const BondedParticlesLUT* compoundPairLutPtrs[Compound::max_bonded_compounds];

		Float3* neighborPositions = utility_buffer_f3;
		half* neighborParticlescharges = particleChargesBuffers;

		if (threadIdx.x < compound.n_bonded_compounds) {
			const uint16_t neighborId = boxConfig.compounds[compound_index].bonded_compound_ids[threadIdx.x];

			const CompoundCoords* neighborCoordsptr = CompoundcoordsCircularQueueUtils::getCoordarrayRef(boxState->compoundcoordsCircularQueue, step, neighborId);
			neighborPtrs[threadIdx.x] = (void*)neighborCoordsptr;
			const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(compound_origo, neighborCoordsptr->origo);
			relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, compound_origo).toFloat3();
			compoundPairLutPtrs[threadIdx.x] = BondedParticlesLUTHelpers::get(sim->boxConfig.bpLUTs, compound_index, neighborId);
		}
		__syncthreads();

		for (int i = 0; i < compound.n_bonded_compounds; i++) {
			const uint16_t neighborId = boxConfig.compounds[compound_index].bonded_compound_ids[i];
			const int neighborNParticles = boxConfig.compounds[neighborId].n_particles;
			
			neighborPositions[threadIdx.x] = ((CompoundCoords*)neighborPtrs[i])->rel_positions[threadIdx.x].toFloat3() + relshifts[i];
			neighborAtomstypes[threadIdx.x] = boxConfig.compoundsAtomtypes[neighborId * MAX_COMPOUND_PARTICLES + threadIdx.x];
			neighborParticlescharges[threadIdx.x] = boxConfig.compoundsAtomCharges[neighborId * MAX_COMPOUND_PARTICLES + threadIdx.x];
			bpLUT.load(*compoundPairLutPtrs[i]);
			__syncthreads();

			if (threadIdx.x < compound.n_particles) {
				force += LJ::computeCompoundCompoundLJForces<computePotE>(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum,
					neighborPositions, neighborNParticles, neighborAtomstypes, &bpLUT, LJ::CalcLJOrigin::ComComInter, forcefield_shared,
					particleCharge, neighborParticlescharges);
			}
			__syncthreads();
		}

	}
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //



	// --------------------------------------------------------------- Solvation forces --------------------------------------------------------------- //
#ifdef ENABLE_SOLVENTS
	{
		int indexInBatch = batchsize;
		for (int i = 0; i < nGridnodes; i++) {
			__syncthreads();
			if (indexInBatch == batchsize) {
				if (threadIdx.x < batchsize && threadIdx.x + i < nGridnodes) {
					neighborIds[threadIdx.x] = sim->compound_neighborlists[blockIdx.x].gridnode_ids[i + threadIdx.x];

					neighborPtrs[threadIdx.x] = (void*)boxState->solventblockgrid_circularqueue->getBlockPtr(neighborIds[threadIdx.x], step);
					const NodeIndex solventblock_hyperorigo = BoundaryCondition::applyHyperpos_Return(compound_origo, BoxGrid::Get3dIndex(neighborIds[threadIdx.x], boxSize_device.boxSizeNM_i));
					relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(solventblock_hyperorigo, compound_origo).toFloat3();
					neighborNParticles[threadIdx.x] = boxState->solventblockgrid_circularqueue->getBlockPtr(neighborIds[threadIdx.x], step)->n_solvents;
				}
				indexInBatch = 0;
				__syncthreads();
			}

			// There are many more solvents in a block, than threads in this kernel, so we need to loop in strides of blockdim.x
			__syncthreads();	// Dont load buffer before all are finished with the previous iteration
			for (uint32_t offset = 0; offset < neighborNParticles[indexInBatch]; offset += blockDim.x) {
				const uint32_t solvent_index = offset + threadIdx.x;
				const int n_elements_this_stride = LAL::min(neighborNParticles[indexInBatch] - offset, blockDim.x);

				// Load the positions and add rel shift
				if (solvent_index < neighborNParticles[indexInBatch]) {
					utility_buffer_f3[threadIdx.x] = ((SolventBlock*)neighborPtrs[indexInBatch])->rel_pos[solvent_index].toFloat3() + relshifts[indexInBatch];
				}
				__syncthreads();

				if (threadIdx.x < compound.n_particles) {
					force += LJ::computeSolventToCompoundLJForces(compound_positions[threadIdx.x], n_elements_this_stride, utility_buffer_f3, potE_sum, compound.atom_types[threadIdx.x], forcefield_shared);
				}
				__syncthreads();
			}
			indexInBatch++;
		}
	}
#endif
	// ------------------------------------------------------------------------------------------------------------------------------------------------ //


	// -------------------------------------------------------------- Distribute charges --------------------------------------------------------------- //	
	if constexpr (ENABLE_ES_LR) {
		if (simparams.enable_electrostatics) {
			__syncthreads();
			char* utilityBuffer = (char*)(&bpLUT);
			Electrostatics::DistributeChargesToChargegrid(compound_origo, compound_positions[threadIdx.x], 
				sim->boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x], sim->chargeGrid, compound.n_particles, utilityBuffer);
		}
	}


	// This is the first kernel, so we overwrite
	if (threadIdx.x < compound.n_particles) {
		if constexpr (computePotE) {
			sim->boxState->compoundsInterimState[blockIdx.x].potE_interim[threadIdx.x] += potE_sum;
		}
		sim->boxState->compoundsInterimState[blockIdx.x].forces_interim[threadIdx.x] += force;
	}
}
template  __global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel<PeriodicBoundaryCondition, true, true>(SimulationDevice* sim, int64_t step);
template  __global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel<PeriodicBoundaryCondition, true, false>(SimulationDevice* sim, int64_t step);
template  __global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel<PeriodicBoundaryCondition, false, true>(SimulationDevice* sim, int64_t step);
template  __global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel<PeriodicBoundaryCondition, false, false>(SimulationDevice* sim, int64_t step);
template __global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel<NoBoundaryCondition, true, true>(SimulationDevice* sim, int64_t step);
template __global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel<NoBoundaryCondition, true, false>(SimulationDevice* sim, int64_t step);
template __global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel<NoBoundaryCondition, false, true>(SimulationDevice* sim, int64_t step);
template __global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel<NoBoundaryCondition, false, false>(SimulationDevice* sim, int64_t step);
#undef compound_index

template <typename BoundaryCondition, bool energyMinimize>
__global__ void compoundBondsAndIntegrationKernel(SimulationDevice* sim, int64_t step) {
	__shared__ CompoundCompact compound;				// Mostly bond information
	__shared__ Float3 compound_positions[THREADS_PER_COMPOUNDBLOCK];
	__shared__ Float3 utility_buffer_f3[THREADS_PER_COMPOUNDBLOCK];
	__shared__ float utility_buffer_f[THREADS_PER_COMPOUNDBLOCK];
	__shared__ NodeIndex compound_origo;

	// Buffer to be cast to different datatypes. This is dangerous!
	__shared__ char utility_buffer[cbkernel_utilitybuffer_size];

	BoxState* boxState = sim->boxState;
	const SimParams& simparams = sim->params;
	SimSignals* signals = sim->signals;
	const Compound* const compound_global = &sim->boxConfig.compounds[blockIdx.x];

	if (threadIdx.x == 0) {
		compound.loadMeta(&sim->boxConfig.compounds[blockIdx.x]);
	}
	__syncthreads();
	compound.loadData(&sim->boxConfig.compounds[blockIdx.x]);

	{
		static_assert(cbkernel_utilitybuffer_size >= sizeof(CompoundCoords), "Utilitybuffer not large enough for CompoundCoords");
		CompoundCoords* compound_coords = (CompoundCoords*)utility_buffer;

		const CompoundCoords& compoundcoords_global = *CompoundcoordsCircularQueueUtils::getCoordarrayRef(boxState->compoundcoordsCircularQueue, step, blockIdx.x);
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
		SingleBond* singlebonds = LoadBonds<SingleBond, MAX_SINGLEBONDS_IN_COMPOUND>(utility_buffer, sim->boxConfig.compounds[blockIdx.x].singlebonds, compound.n_singlebonds);
		force += LimaForcecalc::computeSinglebondForces<energyMinimize>(singlebonds, compound.n_singlebonds, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum, 0);

		AngleUreyBradleyBond* anglebonds = LoadBonds<AngleUreyBradleyBond, MAX_ANGLEBONDS_IN_COMPOUND>(utility_buffer, sim->boxConfig.compounds[blockIdx.x].anglebonds, compound.n_anglebonds);
		force += LimaForcecalc::computeAnglebondForces(anglebonds, compound.n_anglebonds, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);

		DihedralBond* dihedrals = LoadBonds<DihedralBond, MAX_DIHEDRALBONDS_IN_COMPOUND>(utility_buffer, sim->boxConfig.compounds[blockIdx.x].dihedrals, compound.n_dihedrals);
		force += LimaForcecalc::computeDihedralForces(dihedrals, compound.n_dihedrals, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);

		ImproperDihedralBond* impropers = LoadBonds<ImproperDihedralBond, MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND>(utility_buffer, sim->boxConfig.compounds[blockIdx.x].impropers, compound.n_improperdihedrals);
		force += LimaForcecalc::computeImproperdihedralForces(impropers, compound.n_improperdihedrals, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);
	}


	// Fetch interims from other kernels
	if (threadIdx.x < compound.n_particles) {
		force += boxState->compoundsInterimState[blockIdx.x].forces_interim[threadIdx.x];
		potE_sum += boxState->compoundsInterimState[blockIdx.x].potE_interim[threadIdx.x];
	}



	// ------------------------------------------------------------ LongRange Electrostatics --------------------------------------------------------------- //	
	if constexpr (ENABLE_ES_LR) {
		if (simparams.enable_electrostatics && threadIdx.x < compound.n_particles) {
			NodeIndex nodeindex = compound_origo + LIMAPOSITIONSYSTEM::PositionToNodeIndex(compound_positions[threadIdx.x]);
			BoundaryCondition::applyBC(nodeindex);
			const float myCharge = compound_global->atom_charges[threadIdx.x];
			//printf("F %f ES %f\n", force.len(), BoxGrid::GetNodePtr(sim->chargeGridOutputForceAndPot, nodeindex)->force.len());
			force += BoxGrid::GetNodePtr(sim->chargeGridOutputForceAndPot, nodeindex)->forcePart * myCharge;
			potE_sum += BoxGrid::GetNodePtr(sim->chargeGridOutputForceAndPot, nodeindex)->potentialPart * myCharge;
		}
	}



	// ------------------------------------------------------------ Supernatural Forces --------------------------------------------------------------- //	
	__syncthreads();

	//if (simparams.snf_select == HorizontalSqueeze) {
	//	const float mass = forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass;
	//	SupernaturalForces::applyHorizontalSqueeze(utility_buffer_f3, utility_buffer_f, utility_buffer, compound_positions, compound.n_particles, compound_origo, force, mass);
	//
	//}



	if (simparams.snf_select == HorizontalChargeField && threadIdx.x < compound.n_particles) {
		force += sim->boxConfig.uniformElectricField.GetForce(sim->boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x]);
	}
	__syncthreads();


	// -------------------------------------------------------------- Integration & PBC --------------------------------------------------------------- //	
	{
		__syncthreads();

		static_assert(clj_utilitybuffer_bytes >= sizeof(CompoundCoords), "Utilitybuffer not large enough for CompoundCoords");
		CompoundCoords* compound_coords = (CompoundCoords*)utility_buffer;
		if (threadIdx.x == 0) {
			compound_coords->origo = compound_origo;
		}
		//compound_coords->rel_positions[threadIdx.x] = Coord{ compound_positions[threadIdx.x] };	// Alternately load from global again? Will be more precise

		const CompoundCoords& compoundcoords_global = *CompoundcoordsCircularQueueUtils::getCoordarrayRef(boxState->compoundcoordsCircularQueue, step, blockIdx.x);
		compound_coords->rel_positions[threadIdx.x] = compoundcoords_global.rel_positions[threadIdx.x];
		__syncthreads();

		float speed = 0.f;

		float forceLen = force.len();

		__syncthreads();

		if (threadIdx.x < compound.n_particles) {
			const float mass = forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass;

			// Energy minimize
			if constexpr (energyMinimize) {
				const float progress = static_cast<float>(step) / static_cast<float>(simparams.n_steps);
				const Float3 safeForce = EngineUtils::ForceActivationFunction(force);
				const Float3 deltaPosPrev = boxState->compoundsInterimState[blockIdx.x].vels_prev[threadIdx.x];
				const Coord pos_now = EngineUtils::IntegratePositionEM(compound_coords->rel_positions[threadIdx.x], safeForce, mass, simparams.dt, progress, deltaPosPrev);

				// So these dont represent the same as they do in normal MD, but we can use the same buffers
				const Float3 deltaPos = (pos_now - compound_coords->rel_positions[threadIdx.x]).toFloat3();
				//speed = deltaPos.len();
				boxState->compoundsInterimState[blockIdx.x].vels_prev[threadIdx.x] = deltaPos;
				compound_coords->rel_positions[threadIdx.x] = pos_now;// Save pos locally, but only push to box as this kernel ends

			}
			else {
				const Float3 force_prev = boxState->compoundsInterimState[blockIdx.x].forces_prev[threadIdx.x];	// OPTIM: make ref?
				const Float3 vel_prev = boxState->compoundsInterimState[blockIdx.x].vels_prev[threadIdx.x];
				const Float3 vel_now = EngineUtils::integrateVelocityVVS(vel_prev, force_prev, force, simparams.dt, mass);
				const Coord pos_now = EngineUtils::integratePositionVVS(compound_coords->rel_positions[threadIdx.x], vel_now, force, mass, simparams.dt);
				compound_coords->rel_positions[threadIdx.x] = pos_now;// Save pos locally, but only push to box as this kernel ends

				Float3 velScaled;
				velScaled = vel_now * thermostatScalar_device;

				boxState->compoundsInterimState[blockIdx.x].forces_prev[threadIdx.x] = force;
				boxState->compoundsInterimState[blockIdx.x].vels_prev[threadIdx.x] = velScaled;

				speed = velScaled.len();
			}
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
		EngineUtils::LogCompoundData(compound, sim->boxConfig.boxparams.total_particles_upperbound, *compound_coords, &potE_sum, force, force_LJ_sol, simparams, *signals, sim->potE_buffer, sim->traj_buffer, sim->vel_buffer, sim->forceBuffer, speed, step);


		// Push positions for next step
		auto* coordarray_next_ptr = CompoundcoordsCircularQueueUtils::getCoordarrayRef(boxState->compoundcoordsCircularQueue, step + 1, blockIdx.x);
		coordarray_next_ptr->loadData(*compound_coords);
	}
}
template __global__ void compoundBondsAndIntegrationKernel<PeriodicBoundaryCondition, true>(SimulationDevice* sim, int64_t step);
template __global__ void compoundBondsAndIntegrationKernel<PeriodicBoundaryCondition, false>(SimulationDevice* sim, int64_t step);
template __global__ void compoundBondsAndIntegrationKernel<NoBoundaryCondition, true>(SimulationDevice* sim, int64_t step);
template __global__ void compoundBondsAndIntegrationKernel<NoBoundaryCondition, false>(SimulationDevice* sim, int64_t step);





#define solvent_active (threadIdx.x < solventblock.n_solvents)
#define solvent_mass (forcefield_device.particle_parameters[ATOMTYPE_SOLVENT].mass)
static_assert(SolventBlock::MAX_SOLVENTS_IN_BLOCK >= MAX_COMPOUND_PARTICLES, "solventForceKernel was about to reserve an insufficient amount of memory");
template <typename BoundaryCondition, bool energyMinimize>
__global__ void solventForceKernel(SimulationDevice* sim, int64_t step) {
	__shared__ Float3 utility_buffer[SolventBlock::MAX_SOLVENTS_IN_BLOCK];
	__shared__ uint8_t utility_buffer_small[SolventBlock::MAX_SOLVENTS_IN_BLOCK];
	__shared__ SolventBlock solventblock;
	__shared__ SolventTransferqueue<SolventBlockTransfermodule::max_queue_size> transferqueues[6];		// TODO: Use template to make identical kernel, so the kernel with transfer is slower and larger, and the rest remain fast!!!!
	__shared__ int utility_int;
	__shared__ Coord utility_coord;
	__shared__ Float3 utility_float3;

	// Doubles as block_index_3d!
	const NodeIndex block_origo = BoxGrid::Get3dIndex(blockIdx.x, boxSize_device.boxSizeNM_i);

	BoxState* boxState = sim->boxState;
	const SimParams& simparams = sim->params;
	SimSignals* signals = sim->signals;
	SolventBlock* solventblock_ptr = boxState->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, step);

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
		const CompoundGridNode* compoundgridnode = BoxGrid::GetNodePtr(sim->compound_grid, blockIdx.x);
		if (threadIdx.x == 0) { utility_int = compoundgridnode->n_nearby_compounds; }
		__syncthreads();



		for (int i = 0; i < utility_int; i++) {
			const uint16_t neighborcompound_index = compoundgridnode->compoundidsWithinLjCutoff[i];
			const Compound* neighborcompound = &sim->boxConfig.compounds[neighborcompound_index];
			const int n_compound_particles = neighborcompound->n_particles;

			// All threads help loading the molecule
			// First load particles of neighboring compound
			const CompoundCoords* coordarray_ptr = CompoundcoordsCircularQueueUtils::getCoordarrayRef(boxState->compoundcoordsCircularQueue, step, neighborcompound_index);
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

				const SolventBlock* solventblock_neighbor = boxState->solventblockgrid_circularqueue->getBlockPtr(blockindex_neighbor, step);
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
		Solvent& solventdata_ref = boxState->solvents[solventblock.ids[threadIdx.x]];	// Solvent private data, for VVS

		if constexpr (energyMinimize) {
			const float progress = static_cast<float>(step) / static_cast<float>(simparams.n_steps);
			const Float3 safeForce = EngineUtils::ForceActivationFunction(force);
			const Coord pos_now = EngineUtils::IntegratePositionEM(solventblock.rel_pos[threadIdx.x], force, mass, simparams.dt, progress, solventdata_ref.vel_prev);

			const Float3 deltaPos = (pos_now - solventblock.rel_pos[threadIdx.x]).toFloat3();
			solventdata_ref.vel_prev = deltaPos;

			relpos_next = pos_now;
			EngineUtils::LogSolventData(sim->boxConfig.boxparams, potE_sum, solventblock, solvent_active, force, Float3{}, step, sim->potE_buffer, sim->traj_buffer, sim->vel_buffer, simparams.data_logging_interval);
		}
		else {
			Float3 vel_now = EngineUtils::integrateVelocityVVS(solventdata_ref.vel_prev, solventdata_ref.force_prev, force, simparams.dt, mass);
			const Coord pos_now = EngineUtils::integratePositionVVS(solventblock.rel_pos[threadIdx.x], vel_now, force, mass, simparams.dt);
			
			vel_now = vel_now * thermostatScalar_device;

			solventdata_ref.vel_prev = vel_now;
			solventdata_ref.force_prev = force;

			// Save pos locally, but only push to box as this kernel ends
			relpos_next = pos_now;
			EngineUtils::LogSolventData(sim->boxConfig.boxparams, potE_sum, solventblock, solvent_active, force, vel_now, step, sim->potE_buffer, sim->traj_buffer, sim->vel_buffer, simparams.data_logging_interval);
		}
	}



	// Push new SolventCoord to global mem
	SolventBlock* solventblock_next_ptr = boxState->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, step + 1);

	if (SolventBlocksCircularQueue::isTransferStep(step)) {
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
template __global__ void solventForceKernel<PeriodicBoundaryCondition, true>(SimulationDevice*, int64_t);
template __global__ void solventForceKernel<PeriodicBoundaryCondition, false>(SimulationDevice*, int64_t);
template __global__ void solventForceKernel<NoBoundaryCondition, true>(SimulationDevice*, int64_t);
template __global__ void solventForceKernel<NoBoundaryCondition, false>(SimulationDevice*, int64_t);

#undef solvent_index
#undef solvent_mass
#undef solvent_active
#undef solventblock_ptr




// This is run before step.inc(), but will always publish results to the first array in grid!
template <typename BoundaryCondition>
__global__ void solventTransferKernel(SimulationDevice* sim, int64_t step) {
	BoxState* boxState = sim->boxState;

	SolventBlockTransfermodule* transfermodule = &sim->transfermodule_array[blockIdx.x];
	
	SolventBlock* solventblock_current = boxState->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, step);
	SolventBlock* solventblock_next = boxState->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, step + 1);

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
template __global__ void solventTransferKernel<PeriodicBoundaryCondition>(SimulationDevice* sim, int64_t step);
template __global__ void solventTransferKernel<NoBoundaryCondition>(SimulationDevice* sim, int64_t step);



#define particle_id_bridge threadIdx.x
template <typename BoundaryCondition>
__global__ void compoundBridgeKernel(SimulationDevice* sim, int64_t step) {
	__shared__ CompoundBridge bridge;
	__shared__ Float3 positions[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Float3 utility_buffer[MAX_PARTICLES_IN_BRIDGE];
	__shared__ float utility_buffer_f[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Coord utility_coord[MAX_COMPOUNDS_IN_BRIDGE];

	const SimParams& simparams = sim->params;
	BoxState* boxState = sim->boxState;

	if (threadIdx.x == 0) {
		bridge.loadMeta(&sim->boxConfig.bridge_bundle->compound_bridges[blockIdx.x]);

		
	}
	__syncthreads();

	// TODO: we dont need to do this for the first compound, as it will always be 0,0,0
	if (threadIdx.x < bridge.n_compounds) {
		// Calculate necessary shift in relative positions for right, so right share the origo with left.
		utility_coord[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelativeShiftBetweenCoordarrays<BoundaryCondition>(boxState->compoundcoordsCircularQueue, step, bridge.compound_ids[0], bridge.compound_ids[threadIdx.x]);
	}


	bridge.loadData(&sim->boxConfig.bridge_bundle->compound_bridges[blockIdx.x]);
	__syncthreads();

	if (particle_id_bridge < bridge.n_particles) {
		ParticleReference& p_ref = bridge.particle_refs[particle_id_bridge];

		BridgeWarnings::verifyPRefValid(p_ref, bridge);

		const CompoundCoords* coordarray = CompoundcoordsCircularQueueUtils::getCoordarrayRef(boxState->compoundcoordsCircularQueue, step, p_ref.compound_id);

		Coord relpos = coordarray->rel_positions[p_ref.local_id_compound];
		relpos += utility_coord[p_ref.compoundid_local_to_bridge];
		positions[threadIdx.x] = relpos.toFloat3();
	}
	__syncthreads();

	float potE_sum = 0;
	Float3 force{};

	// ------------------------------------------------------------ Intercompund Operations ------------------------------------------------------------ //
	{											// So for the very first step, these ´should all be 0, but they are not??										TODO: Look into this at some point!!!! 
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
		boxState->compoundsInterimState[p_ref->compound_id].forces_interim[p_ref->local_id_compound] += force;
		boxState->compoundsInterimState[p_ref->compound_id].potE_interim[p_ref->local_id_compound] += potE_sum;
	}
}
template __global__ void compoundBridgeKernel<PeriodicBoundaryCondition>(SimulationDevice* sim, int64_t step);
template __global__ void compoundBridgeKernel<NoBoundaryCondition>(SimulationDevice* sim, int64_t step);
#pragma warning (pop)
#pragma warning (pop)
