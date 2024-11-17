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
#include "LennardJonesInteractions.cuh"


#include <cfloat>

#pragma warning(push)
#pragma warning(disable:E0020)
#pragma warning(push)
#pragma warning(disable: 20054)

#pragma diag_suppress 20054



















// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//
__global__ void DistributeCompoundchargesToGridKernel(SimulationDevice* sim) {
	NodeIndex compound_origo = sim->boxState->compoundOrigos[blockIdx.x];
	const Float3 relPos = sim->boxState->compoundsRelposLm[blockIdx.x * MAX_COMPOUND_PARTICLES];
	const int nParticles = sim->boxConfig.compounds[blockIdx.x].n_particles;
	char utilityBuffer[sizeof(int) * (27 * 2 + MAX_COMPOUND_PARTICLES)];


	__syncthreads();
	Electrostatics::DistributeChargesToChargegrid(compound_origo, relPos,
		sim->boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x], sim->chargeGrid, nParticles, utilityBuffer);	
}
// #define compound_index blockIdx.x
// template <typename BoundaryCondition, bool energyMinimize, bool computePotE> // We dont compute potE if we dont log data this step
// __global__ void compoundFarneighborShortrangeInteractionsKernel(const int64_t step, const BoxState boxState, const BoxConfig boxConfig, const NeighborList* const compoundNeighborlists, bool enableES) {
// 	__shared__ CompoundCompact compound;				// Mostly bond information
// 	__shared__ Float3 compound_positions[MAX_COMPOUND_PARTICLES]; // [lm]
//
// 	//Neighborlist
// 	__shared__ int nNonbondedCompoundNeighbors;
//
// 	__shared__ Float3 utility_buffer_f3[MAX_COMPOUND_PARTICLES*2];
//
// 	__shared__ half particleChargesBuffers[MAX_COMPOUND_PARTICLES * 2];
//
// 	__shared__ ForceField_NB forcefield_shared;
//
// 	NodeIndex compoundOrigo;
//
// 	float potE_sum{};
// 	Float3 force{};
// 	//NodeIndex compound_origo{}; // TODO: Figure out why TF everything breaks if i make this __shared__???????
// 	const float particleCharge = enableES	// TODO: this is temporary
// 		? static_cast<float>(boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x])
// 		: 0.f;
//
// 	{
// 		auto block = cooperative_groups::this_thread_block();
//
// 		compoundOrigo = boxState.compoundOrigos[blockIdx.x];
// 		cooperative_groups::memcpy_async(block, compound_positions, &boxState.compoundsRelposLm[blockIdx.x * MAX_COMPOUND_PARTICLES], sizeof(Float3) * MAX_COMPOUND_PARTICLES);
// 		if (threadIdx.x == 0) {
// 			compound.loadMeta(&boxConfig.compounds[blockIdx.x]);
// 			nNonbondedCompoundNeighbors = compoundNeighborlists[blockIdx.x].nNonbondedNeighbors;
// 		}
//
// 		cooperative_groups::memcpy_async(block, &forcefield_shared, &forcefield_device, sizeof(ForceField_NB));
//
// 		cooperative_groups::wait(block);
//
//
// 		/*if (compoundcoords_global->origo != compoundOrigo)
// 			printf("wtf");*/
// 	}
// 	compound.loadData(&boxConfig.compounds[blockIdx.x]);
//
//
//
// 	const int batchsize = 32;
// 	__shared__ Float3 relshifts[batchsize];	// [lm]
// 	__shared__ int neighborIds[batchsize]; // either compoundID or solventblockID
// 	__shared__ int neighborNParticles[batchsize]; // either particlesInCompound or particlesInSolventblock
// 	// --------------------------------------------------------------- Intercompound forces --------------------------------------------------------------- //
// 	{
// 		__shared__ uint8_t atomtypesBuffer[MAX_COMPOUND_PARTICLES*2];
// 		auto block = cooperative_groups::this_thread_block();
// 		// This part is scary, but it also takes up by far the majority of compute time. We use the utilitybuffer twice simultaneously, so be careful when making changes
// 		int indexInBatch = batchsize;
// 		for (int i = 0; i < nNonbondedCompoundNeighbors; i++) {
// 			__syncthreads();
// 			static_assert(sizeof utility_buffer_f3 >= sizeof(Float3) * MAX_COMPOUND_PARTICLES * 2, "Utilitybuffer not large enough for neighbor positions");
// 			Float3* neighborPositionsCurrent = ((Float3*)utility_buffer_f3) + (MAX_COMPOUND_PARTICLES * (i & 1));
// 			Float3* neighborPositionsNext = ((Float3*)utility_buffer_f3) + (MAX_COMPOUND_PARTICLES * !(i & 1));
// 			uint8_t* neighborAtomstypesCurrent = atomtypesBuffer + (MAX_COMPOUND_PARTICLES * (i & 1));
// 			uint8_t* neighborAtomstypesNext = atomtypesBuffer + (MAX_COMPOUND_PARTICLES * !(i & 1));
// 			half* neighborParticleschargesCurrent = particleChargesBuffers + (MAX_COMPOUND_PARTICLES * (i & 1));
// 			half* neighborParticleschargesNext = particleChargesBuffers + (MAX_COMPOUND_PARTICLES * !(i & 1));
//
// 			static_assert(MAX_COMPOUND_PARTICLES >= batchsize);
// 			// First check if we need to load a new batch of relshifts & n_particles for the coming 32 compounds
// 			if (indexInBatch == batchsize) {
// 				if (threadIdx.x < batchsize && threadIdx.x + i < nNonbondedCompoundNeighbors) {
// 					neighborIds[threadIdx.x] = compoundNeighborlists[blockIdx.x].nonbondedNeighborcompoundIds[i + threadIdx.x];
//
// 					const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(compoundOrigo, boxState.compoundOrigos[neighborIds[threadIdx.x]]);
// 					KernelHelpersWarnings::assertHyperorigoIsValid(querycompound_hyperorigo, compoundOrigo);
//
// 					// calc Relative LimaPosition Shift from the origo-shift
// 					relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, compoundOrigo).toFloat3();
// 					neighborNParticles[threadIdx.x] = boxConfig.compounds[neighborIds[threadIdx.x]].n_particles;
// 				}
// 				indexInBatch = 0;
// 				__syncthreads();
//
// 				{
// 					// Load first element in batch and sync
// 					const int currentNeighborId = neighborIds[indexInBatch];
// 					const int currentNeighborNParticles = neighborNParticles[indexInBatch];
//
// 					if (threadIdx.x < currentNeighborNParticles) {
// 						neighborAtomstypesCurrent[threadIdx.x] = boxConfig.compoundsAtomtypes[currentNeighborId * MAX_COMPOUND_PARTICLES + threadIdx.x];
// 						neighborParticleschargesCurrent[threadIdx.x] = boxConfig.compoundsAtomCharges[currentNeighborId * MAX_COMPOUND_PARTICLES + threadIdx.x];
// 						neighborPositionsCurrent[threadIdx.x] = boxState.compoundsRelposLm[currentNeighborId * MAX_COMPOUND_PARTICLES + threadIdx.x] + relshifts[indexInBatch];
// 					}
// 					__syncthreads();
// 				}
// 			}
//
// 			if (i + 1 < nNonbondedCompoundNeighbors && indexInBatch+1 < batchsize) {
// 				static_assert(sizeof(Coord) == sizeof(Float3));
// 				const int nextNeighborId = neighborIds[indexInBatch + 1];
// 				//const int nextNeighborNParticles = neighborNParticles[indexInBatch + 1];
//
// 				//neighborAtomstypesNext[threadIdx.x] = boxConfig.compoundsAtomtypes[neighborIds[indexInBatch+1] * MAX_COMPOUND_PARTICLES + threadIdx.x];
// 				//neighborParticleschargesNext[threadIdx.x] = boxConfig.compoundsAtomCharges[neighborIds[indexInBatch+1] * MAX_COMPOUND_PARTICLES + threadIdx.x];
//
// 				neighborAtomstypesNext[threadIdx.x] = boxConfig.compoundsAtomtypes[nextNeighborId * MAX_COMPOUND_PARTICLES + threadIdx.x];
// 				neighborParticleschargesNext[threadIdx.x] = boxConfig.compoundsAtomCharges[nextNeighborId * MAX_COMPOUND_PARTICLES + threadIdx.x];
// 				neighborPositionsNext[threadIdx.x] = boxState.compoundsRelposLm[nextNeighborId * MAX_COMPOUND_PARTICLES + threadIdx.x] + relshifts[indexInBatch + 1];
//
// 				//cooperative_groups::memcpy_async(block, neighborAtomstypesNext, boxConfig.compoundsAtomtypes + nextNeighborId * MAX_COMPOUND_PARTICLES, sizeof(uint8_t)* MAX_COMPOUND_PARTICLES);
// 				//cooperative_groups::memcpy_async(block, neighborParticleschargesNext, boxConfig.compoundsAtomCharges + nextNeighborId * MAX_COMPOUND_PARTICLES, sizeof(half) * MAX_COMPOUND_PARTICLES);
// 				//cooperative_groups::memcpy_async(block, neighborPositionsNext, &boxState.compoundsRelposLm[nextNeighborId*MAX_COMPOUND_PARTICLES], sizeof(Coord)* MAX_COMPOUND_PARTICLES);
// 			}
//
// 			if (threadIdx.x < compound.n_particles) {
// 				force += LJ::computeCompoundCompoundLJForces<computePotE, energyMinimize>(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum,
// 					neighborPositionsCurrent, neighborNParticles[indexInBatch], neighborAtomstypesCurrent, forcefield_shared, particleCharge, neighborParticleschargesCurrent);
// 			}
//
// 			// Process the positions
// 			if (indexInBatch + 1 < batchsize) {
// 				//cooperative_groups::wait(block); // Joins all threads, waits for all copies to complete
// 				//neighborPositionsNext[threadIdx.x] += relshifts[indexInBatch + 1];
// 			}
// 			indexInBatch++;
// 		}
// 	}
// 	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //
//
// 	// This is the first kernel, so we overwrite
// 	if (threadIdx.x < compound.n_particles) {
//
// 		//if (isinf(force.lenSquared()))
// 		//	printf("Far nei %d\n", threadIdx.x < compound.n_particles);
// 		if constexpr (computePotE) {
// 			boxState.compoundsInterimState[blockIdx.x].potE_interim[threadIdx.x] = potE_sum;
// 		}
// 		boxState.compoundsInterimState[blockIdx.x].forces_interim[threadIdx.x] = force;
// 	}
// }
// #undef compound_index

template <typename BoundaryCondition, bool energyMinimize, bool computePotE> // We dont compute potE if we dont log data this step
__global__ void compoundFarneighborShortrangeInteractionsKernel(const int64_t step, const BoxState boxState, const BoxConfig boxConfig, const NeighborList* const compoundNeighborlists, bool enableES) {
	__shared__ CompoundCompact compound;				// Mostly bond information
	__shared__ Float3 compound_positions[MAX_COMPOUND_PARTICLES]; // [lm] // TODO: maybe only keep these in register mem

	//Neighborlist
	__shared__ int nNonbondedCompoundNeighbors;

	__shared__ Float3 neighborPositions[MAX_COMPOUND_PARTICLES];

	__shared__ float neighborParticlescharges[MAX_COMPOUND_PARTICLES];

	__shared__ ForceField_NB forcefield_shared;
	__shared__ uint8_t neighborAtomstypes[MAX_COMPOUND_PARTICLES];

	NodeIndex compoundOrigo;

	float potE_sum{};
	Float3 force{};
	//NodeIndex compound_origo{}; // TODO: Figure out why TF everything breaks if i make this __shared__???????
	const float particleCharge = enableES	// TODO: this is temporary
	? static_cast<float>(boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x])
	: 0.f;

	{
		auto block = cooperative_groups::this_thread_block();

		compoundOrigo = boxState.compoundOrigos[blockIdx.x];
		cooperative_groups::memcpy_async(block, compound_positions, &boxState.compoundsRelposLm[blockIdx.x * MAX_COMPOUND_PARTICLES], sizeof(Float3) * MAX_COMPOUND_PARTICLES);
		if (threadIdx.x == 0) {
			compound.loadMeta(&boxConfig.compounds[blockIdx.x]);
			nNonbondedCompoundNeighbors = compoundNeighborlists[blockIdx.x].nNonbondedNeighbors;
		}

		cooperative_groups::memcpy_async(block, &forcefield_shared, &forcefield_device, sizeof(ForceField_NB));

		cooperative_groups::wait(block);
	}
	compound.loadData(&boxConfig.compounds[blockIdx.x]);



    const int batchsize = 32;
	__shared__ Float3 relshifts[batchsize];	// [lm]
	__shared__ int neighborIds[batchsize]; // either compoundID or solventblockID
	__shared__ int neighborNParticles[batchsize]; // either particlesInCompound or particlesInSolventblock
	// --------------------------------------------------------------- Intercompound forces --------------------------------------------------------------- //
	{
		//auto block = cooperative_groups::this_thread_block();
		// This part is scary, but it also takes up by far the majority of compute time. We use the utilitybuffer twice simultaneously, so be careful when making changes
		int indexInBatch = batchsize;
		for (int i = 0; i < nNonbondedCompoundNeighbors; i++) {
			__syncthreads();

			static_assert(MAX_COMPOUND_PARTICLES >= batchsize);
			// First check if we need to load a new batch of relshifts & n_particles for the coming 32 compounds
			if (indexInBatch == batchsize) {
				if (threadIdx.x < batchsize && threadIdx.x + i < nNonbondedCompoundNeighbors) {
					neighborIds[threadIdx.x] = compoundNeighborlists[blockIdx.x].nonbondedNeighborcompoundIds[i + threadIdx.x];

					const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(compoundOrigo, boxState.compoundOrigos[neighborIds[threadIdx.x]]);
					KernelHelpersWarnings::assertHyperorigoIsValid(querycompound_hyperorigo, compoundOrigo);

					// calc Relative LimaPosition Shift from the origo-shift
					relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, compoundOrigo).toFloat3();
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
				neighborPositions[threadIdx.x] = boxState.compoundsRelposLm[currentNeighborId * MAX_COMPOUND_PARTICLES + threadIdx.x] + relshifts[indexInBatch];
			}
			__syncthreads();

			if (threadIdx.x < compound.n_particles) {
				force += LJ::computeCompoundCompoundLJForces<computePotE, energyMinimize>(compound_positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum,
					neighborPositions, neighborNParticles[indexInBatch], neighborAtomstypes, forcefield_shared, particleCharge, neighborParticlescharges);
			}

			indexInBatch++;
		}
	}
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //

	// This is the first kernel, so we overwrite
	if (threadIdx.x < compound.n_particles) {

		boxState.compoundsInterimState[blockIdx.x].forceEnergyFarneighborShortrange[threadIdx.x] = ForceEnergy{ force, potE_sum };
		//if (isinf(force.lenSquared()))
		//	printf("Far nei %d\n", threadIdx.x < compound.n_particles);
		/*if constexpr (computePotE) {
			boxState.compoundsInterimState[blockIdx.x].potE_interim[threadIdx.x] = potE_sum;
		}
		boxState.compoundsInterimState[blockIdx.x].forces_interim[threadIdx.x] = force;*/
	}
}


#define compound_index blockIdx.x
template <typename BoundaryCondition, bool energyMinimize, bool computePotE> // We dont compute potE if we dont log data this step
__global__ void compoundImmediateneighborAndSelfShortrangeInteractionsKernel(SimulationDevice* sim, const int64_t step) {
	__shared__ CompoundCompact compound;				// Mostly bond information
	__shared__ Float3 compound_positions[MAX_COMPOUND_PARTICLES]; // [lm]

	__shared__ int nGridnodes;

	__shared__ Float3 utility_buffer_f3[MAX_COMPOUND_PARTICLES * 2];
	__shared__ Float3 utility_float3;

	__shared__ BondedParticlesLUT bpLUT;
	__shared__ float particleChargesBuffers[MAX_COMPOUND_PARTICLES * 2];

	__shared__ ForceField_NB forcefield_shared;

	__shared__ uint8_t neighborAtomstypes[MAX_COMPOUND_PARTICLES];

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

		compound_origo = boxState->compoundOrigos[blockIdx.x];
		cooperative_groups::memcpy_async(block, compound_positions, &boxState->compoundsRelposLm[blockIdx.x * MAX_COMPOUND_PARTICLES], sizeof(Coord) * MAX_COMPOUND_PARTICLES);

		if (threadIdx.x == 0) {
			compound.loadMeta(&boxConfig.compounds[blockIdx.x]);
			nGridnodes = sim->compound_neighborlists[blockIdx.x].n_gridnodes;
		}

		cooperative_groups::memcpy_async(block, &forcefield_shared, &forcefield_device, sizeof(ForceField_NB));

		cooperative_groups::wait(block);
	}
	compound.loadData(&boxConfig.compounds[blockIdx.x]);



	// ------------------------------------------------------------ Intracompound Operations ------------------------------------------------------------ //
	{
		static_assert(clj_utilitybuffer_bytes >= sizeof(BondedParticlesLUT), "Utilitybuffer not large enough for BondedParticlesLUT");

		const BondedParticlesLUT* const bplut_global = BondedParticlesLUTHelpers::get(sim->boxConfig.bpLUTs, compound_index, compound_index);
		bpLUT.load(*bplut_global);	// A lut always exists within a compound

		static_assert(clj_utilitybuffer_bytes >= sizeof(BondedParticlesLUT) + sizeof(float) * MAX_COMPOUND_PARTICLES,
			"Utilitybuffer not large enough for neighbor charges");
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
	const int batchsize = 64;
	__shared__ Float3 relshifts[batchsize];	// [lm]
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
			relshifts[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, compound_origo).toFloat3();
			compoundPairLutPtrs[threadIdx.x] = BondedParticlesLUTHelpers::get(sim->boxConfig.bpLUTs, compound_index, neighborId);
		}
		__syncthreads();

		for (int i = 0; i < compound.n_bonded_compounds; i++) {
			const uint16_t neighborId = boxConfig.compounds[compound_index].bonded_compound_ids[i];
			const int neighborNParticles = boxConfig.compounds[neighborId].n_particles;
			
			neighborPositions[threadIdx.x] = boxState->compoundsRelposLm[neighborId * MAX_COMPOUND_PARTICLES + threadIdx.x] + relshifts[i];
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
		forcefieldTinymol_shared.types[threadIdx.x] = tinymolForcefield_device.types[threadIdx.x];

	{
		int indexInBatch = batchsize;
		for (int i = 0; i < nGridnodes; i++) {
			__syncthreads();
			if (indexInBatch == batchsize) {
				if (threadIdx.x < batchsize && threadIdx.x + i < nGridnodes) {
					neighborIds[threadIdx.x] = sim->compound_neighborlists[blockIdx.x].gridnode_ids[i + threadIdx.x];

					solventblockPtrs[threadIdx.x] = boxState->solventblockgrid_circularqueue->getBlockPtr(neighborIds[threadIdx.x], step);
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
				const uint32_t queryIndex = offset + threadIdx.x;
				const int n_elements_this_stride = LAL::min(neighborNParticles[indexInBatch] - offset, blockDim.x);

				// Load the positions and add rel shift
				if (queryIndex < neighborNParticles[indexInBatch]) {
					utility_buffer_f3[threadIdx.x] = solventblockPtrs[indexInBatch]->rel_pos[queryIndex].toFloat3() + relshifts[indexInBatch];
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
		sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyImmediateneighborShortrange[threadIdx.x] = ForceEnergy{ force, potE_sum };
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
__global__ void compoundBondsKernel(SimulationDevice* sim, int64_t step, const UniformElectricField uniformElectricField) {
	__shared__ CompoundCompact compound;				// Mostly bond information
	__shared__ Float3 compound_positions[THREADS_PER_COMPOUNDBLOCK];
	__shared__ Float3 utility_buffer_f3[THREADS_PER_COMPOUNDBLOCK];
	__shared__ float utility_buffer_f[THREADS_PER_COMPOUNDBLOCK];
	//__shared__ NodeIndex compound_origo;

	// Buffer to be cast to different datatypes. This is dangerous!
	__shared__ char utility_buffer[cbkernel_utilitybuffer_size];

	BoxState* boxState = sim->boxState;

	if (threadIdx.x == 0) {
		compound.loadMeta(&sim->boxConfig.compounds[blockIdx.x]);
	}
	__syncthreads();
	compound.loadData(&sim->boxConfig.compounds[blockIdx.x]);

	{
		static_assert(cbkernel_utilitybuffer_size >= sizeof(CompoundCoords), "Utilitybuffer not large enough for CompoundCoords");
		__syncthreads();
		compound_positions[threadIdx.x] = boxState->compoundsRelposLm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
		__syncthreads();
	}

	float potE_sum{};
	Float3 force{};

	// ------------------------------------------------------------ Intracompound Operations ------------------------------------------------------------ //
	{
		SingleBond* singlebonds = EngineUtils::LoadBonds<SingleBond, MAX_SINGLEBONDS_IN_COMPOUND>(utility_buffer, sim->boxConfig.compounds[blockIdx.x].singlebonds, compound.n_singlebonds);
		force += LimaForcecalc::computeSinglebondForces<energyMinimize>(singlebonds, compound.n_singlebonds, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum, 0);

		AngleUreyBradleyBond* anglebonds = EngineUtils::LoadBonds<AngleUreyBradleyBond, MAX_ANGLEBONDS_IN_COMPOUND>(utility_buffer, sim->boxConfig.compounds[blockIdx.x].anglebonds, compound.n_anglebonds);
		force += LimaForcecalc::computeAnglebondForces(anglebonds, compound.n_anglebonds, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);

		DihedralBond* dihedrals = EngineUtils::LoadBonds<DihedralBond, MAX_DIHEDRALBONDS_IN_COMPOUND>(utility_buffer, sim->boxConfig.compounds[blockIdx.x].dihedrals, compound.n_dihedrals);
		force += LimaForcecalc::computeDihedralForces(dihedrals, compound.n_dihedrals, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);

		ImproperDihedralBond* impropers = EngineUtils::LoadBonds<ImproperDihedralBond, MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND>(utility_buffer, sim->boxConfig.compounds[blockIdx.x].impropers, compound.n_improperdihedrals);
		force += LimaForcecalc::computeImproperdihedralForces(impropers, compound.n_improperdihedrals, compound_positions, utility_buffer_f3, utility_buffer_f, &potE_sum);
	}

	// ------------------------------------------------------------ Supernatural Forces --------------------------------------------------------------- //	
	if (sim->params.snf_select == HorizontalChargeField && threadIdx.x < compound.n_particles) {		
		force += uniformElectricField.GetForce(sim->boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x]);
	}

	// ------------------------------------------------------------ Push Data --------------------------------------------------------------- //	
	if (threadIdx.x < compound.n_particles) {
		boxState->compoundsInterimState[blockIdx.x].forceEnergyBonds[threadIdx.x] = ForceEnergy{ force, potE_sum };
	}
}

template<typename BoundaryCondition, bool emvariant>
__global__ void CompoundIntegrationKernel(SimulationDevice* sim, int64_t step) {
	
	__shared__ CompoundCoords compound_coords;
	__shared__ uint8_t atom_types[MAX_COMPOUND_PARTICLES];
	const int nParticles = sim->boxConfig.compounds[blockIdx.x].n_particles;
	if (threadIdx.x == 0) {
		compound_coords.origo = sim->boxState->compoundOrigos[blockIdx.x];
	}	
	compound_coords.rel_positions[threadIdx.x] = sim->boxState->compoundsInterimState[blockIdx.x].coords[threadIdx.x];
	atom_types[threadIdx.x] = sim->boxConfig.compounds[blockIdx.x].atom_types[threadIdx.x];

	// Fetch interims from other kernels
	Float3 force = sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyFarneighborShortrange[threadIdx.x].force;
	float potE_sum = sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyFarneighborShortrange[threadIdx.x].potE;

	force += sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyImmediateneighborShortrange[threadIdx.x].force;
	potE_sum += sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyImmediateneighborShortrange[threadIdx.x].potE;

	force += sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyBonds[threadIdx.x].force;
	potE_sum += sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyBonds[threadIdx.x].potE;

	force += sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyBridge[threadIdx.x].force;
	potE_sum += sim->boxState->compoundsInterimState[blockIdx.x].forceEnergyBridge[threadIdx.x].potE;
	


	
	__syncthreads();

	// ------------------------------------------------------------ LongRange Electrostatics --------------------------------------------------------------- //	
	if constexpr (ENABLE_ES_LR) {
		if (sim->params.enable_electrostatics && threadIdx.x < nParticles) {
			NodeIndex nodeindex = compound_coords.origo + LIMAPOSITIONSYSTEM::PositionToNodeIndex(compound_coords.rel_positions[threadIdx.x].toFloat3());
			BoundaryCondition::applyBC(nodeindex);
			const float myCharge = sim->boxConfig.compounds[blockIdx.x].atom_charges[threadIdx.x];
			//printf("F %f ES %f\n", force.len(), BoxGrid::GetNodePtr(sim->chargeGridOutputForceAndPot, nodeindex)->force.len());			
			if (BoxGrid::GetNodePtr(sim->chargeGridOutputForceAndPot, nodeindex) == nullptr) {
				//printf("nullptr 0");
				auto a = compound_coords.origo + LIMAPOSITIONSYSTEM::PositionToNodeIndex(compound_coords.rel_positions[threadIdx.x].toFloat3());
				printf("abs %d %d %d hyper %d %d %d  BPD %d\n", a.x, a.y, a.z, nodeindex.x, nodeindex.y, nodeindex.z, boxSize_device.blocksPerDim);
			}


			force += BoxGrid::GetNodePtr(sim->chargeGridOutputForceAndPot, nodeindex)->forcePart * myCharge;
			potE_sum += BoxGrid::GetNodePtr(sim->chargeGridOutputForceAndPot, nodeindex)->potentialPart * myCharge;
		}
	}


	// ------------------------------------------------------------ Integration --------------------------------------------------------------- //	
#ifdef FORCE_NAN_CHECK
	if (isnan(force.len()))
		printf("NAN\n");
#endif
	float speed = 0.f;
	if (threadIdx.x < nParticles) {
		const float mass = sim->boxConfig.compounds[blockIdx.x].atomMasses[threadIdx.x];

		// Energy minimize
		if constexpr (emvariant) {
			const float progress = static_cast<float>(step) / static_cast<float>(sim->params.n_steps);
			const Float3 safeForce = EngineUtils::ForceActivationFunction(force);

			AdamState* const adamState = &sim->adamState[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
			const Coord pos_now = EngineUtils::IntegratePositionADAM(compound_coords.rel_positions[threadIdx.x], safeForce, adamState, step);

			compound_coords.rel_positions[threadIdx.x] = pos_now;// Save pos locally, but only push to box as this kernel ends
		}
		else {
			const Float3 force_prev = sim->boxState->compoundsInterimState[blockIdx.x].forces_prev[threadIdx.x];	// OPTIM: make ref?
			const Float3 vel_prev = sim->boxState->compoundsInterimState[blockIdx.x].vels_prev[threadIdx.x];
			const Float3 vel_now = EngineUtils::integrateVelocityVVS(vel_prev, force_prev, force, sim->params.dt, mass);
			const Coord pos_now = EngineUtils::integratePositionVVS(compound_coords.rel_positions[threadIdx.x], vel_now, force, mass, sim->params.dt);
			compound_coords.rel_positions[threadIdx.x] = pos_now;// Save pos locally, but only push to box as this kernel ends

			Float3 velScaled;
			velScaled = vel_now * thermostatScalar_device;

			sim->boxState->compoundsInterimState[blockIdx.x].forces_prev[threadIdx.x] = force;
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
	EngineUtils::LogCompoundData(sim->boxConfig.compounds[blockIdx.x], sim->boxparams.total_particles_upperbound, compound_coords, &potE_sum, force,
		force_LJ_sol, sim->params, *sim->signals, sim->potE_buffer, sim->traj_buffer, sim->vel_buffer, sim->forceBuffer, speed, step);


	// Push positions for next step
	if (threadIdx.x == 0)
		sim->boxState->compoundOrigos[blockIdx.x] = compound_coords.origo;
	sim->boxState->compoundsInterimState[blockIdx.x].coords[threadIdx.x] = compound_coords.rel_positions[threadIdx.x];
	sim->boxState->compoundsRelposLm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = compound_coords.rel_positions[threadIdx.x].toFloat3();
}
template __global__ void CompoundIntegrationKernel<PeriodicBoundaryCondition, true>(SimulationDevice* sim, int64_t step);
template __global__ void CompoundIntegrationKernel<PeriodicBoundaryCondition, false>(SimulationDevice* sim, int64_t step);
template __global__ void CompoundIntegrationKernel<NoBoundaryCondition, true>(SimulationDevice* sim, int64_t step);
template __global__ void CompoundIntegrationKernel<NoBoundaryCondition, false>(SimulationDevice* sim, int64_t step);



static_assert(SolventBlock::MAX_SOLVENTS_IN_BLOCK >= MAX_COMPOUND_PARTICLES, "solventForceKernel was about to reserve an insufficient amount of memory");
template <typename BoundaryCondition, bool energyMinimize>
__global__ void solventForceKernel(BoxState boxState, const BoxConfig boxConfig, const CompoundGridNode* const compoundGrid, int64_t step) {
	__shared__ Float3 utility_buffer[SolventBlock::MAX_SOLVENTS_IN_BLOCK];
	__shared__ uint8_t utility_buffer_small[SolventBlock::MAX_SOLVENTS_IN_BLOCK];
	__shared__ int neighborblockNumElements;
	__shared__ Float3 utility_float3;
	__shared__ ForcefieldTinymol forcefieldTinymol_shared;
	__shared__ int nElementsInBlock;

	// Doubles as block_index_3d!
	const NodeIndex block_origo = BoxGrid::Get3dIndex(blockIdx.x, boxSize_device.boxSizeNM_i);

	const SolventBlock* const solventblock_ptr = boxState.solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, step);


	if (threadIdx.x < ForcefieldTinymol::MAX_TYPES) {
		forcefieldTinymol_shared.types[threadIdx.x] = tinymolForcefield_device.types[threadIdx.x]; // TODO Check that im using this and not the __constant??
	}

	if (threadIdx.x == 0) {
		nElementsInBlock = solventblock_ptr->n_solvents;
	}
	__syncthreads();
	

	const bool threadActive = threadIdx.x < nElementsInBlock;

	Float3 force{};
	float potE_sum{};
	const Float3 relpos_self = solventblock_ptr->rel_pos[threadIdx.x].toFloat3();
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
			EngineUtils::getCompoundHyperpositionsAsFloat3<BoundaryCondition>(block_origo, boxState.compoundOrigos[neighborcompound_index], &boxState.compoundsRelposLm[neighborcompound_index * MAX_COMPOUND_PARTICLES],
				utility_buffer, utility_float3, n_compound_particles);


			// Then load atomtypes of neighboring compound
			if (threadIdx.x < n_compound_particles) {
				utility_buffer_small[threadIdx.x] = neighborcompound->atom_types[threadIdx.x];
			}
			__syncthreads();

			//  We can optimize here by loading and calculate the paired sigma and eps, jsut remember to loop threads, if there are many aomttypes.
			if (threadActive) {// TODO: use computePote template param here
				force += LJ::computeCompoundToSolventLJForces<true, energyMinimize>(relpos_self, n_compound_particles, utility_buffer, potE_sum, 
					utility_buffer_small, idSelf, tinymolForcefield_device, tinymolTypeId);
			}
			__syncthreads();
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //



	// --------------------------------------------------------------- Intrablock TinyMolState Interactions ----------------------------------------------------- //
	{		
		__syncthreads(); // Sync since use of utility
		if (threadActive) {
			utility_buffer[threadIdx.x] = relpos_self;
			utility_buffer_small[threadIdx.x] = tinymolTypeId;
		}
		__syncthreads();
		if (threadActive) {
			force += LJ::computeSolventToSolventLJForces<true, energyMinimize>(relpos_self, tinymolTypeId, utility_buffer, nElementsInBlock, true, potE_sum, forcefieldTinymol_shared, utility_buffer_small);
		}
		__syncthreads(); // Sync since use of utility
	}	
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //

	// --------------------------------------------------------------- Interblock TinyMolState Interactions ----------------------------------------------------- //
	const int query_range = 2;
	for (int x = -query_range; x <= query_range; x++) {
		for (int y = -query_range; y <= query_range; y++) {
			for (int z = -query_range; z <= query_range; z++) {
				const NodeIndex dir{ x,y,z };
				if (dir.sum() > 3) { continue; }
				if (dir.isZero()) { continue; }

				const int blockindex_neighbor = EngineUtils::getNewBlockId<BoundaryCondition>(dir, block_origo);
				KernelHelpersWarnings::assertValidBlockId(blockindex_neighbor);

				const SolventBlock* solventblock_neighbor = boxState.solventblockgrid_circularqueue->getBlockPtr(blockindex_neighbor, step);
				const int nsolvents_neighbor = solventblock_neighbor->n_solvents;
				const Float3 origoshift_offset = LIMAPOSITIONSYSTEM::nodeIndexToCoord(dir).toFloat3();

				// All threads help loading the solvent, and shifting it's relative position reletive to this solventblock
				__syncthreads();
				if (threadIdx.x < nsolvents_neighbor) {
					utility_buffer[threadIdx.x] = solventblock_neighbor->rel_pos[threadIdx.x].toFloat3() + origoshift_offset;
					utility_buffer_small[threadIdx.x] = solventblock_neighbor->atomtypeIds[threadIdx.x];
				}
				__syncthreads();

				if (threadActive) {
					force += LJ::computeSolventToSolventLJForces<true, energyMinimize>(relpos_self, tinymolTypeId, utility_buffer, nsolvents_neighbor, false, potE_sum, forcefieldTinymol_shared, utility_buffer_small);
				}
				__syncthreads();
			}
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //

	// Finally push force and potE for next kernel
	boxState.solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, step)->forceEnergies[threadIdx.x] = ForceEnergy{ force, potE_sum };
}
template __global__ void solventForceKernel<PeriodicBoundaryCondition, true>(const BoxState, const BoxConfig, const CompoundGridNode* const, int64_t);
template __global__ void solventForceKernel<PeriodicBoundaryCondition, false>(const BoxState, const BoxConfig, const CompoundGridNode* const, int64_t);
template __global__ void solventForceKernel<NoBoundaryCondition, true>(const BoxState, const BoxConfig, const CompoundGridNode* const, int64_t);
template __global__ void solventForceKernel<NoBoundaryCondition, false>(const BoxState, const BoxConfig, const CompoundGridNode* const, int64_t);


template <typename BoundaryCondition, bool energyMinimize, bool transferOutThisStep>
__global__ void TinymolIntegrationLoggingAndTransferout(SimulationDevice* sim, int64_t step) {
	__shared__ SolventBlock solventblock;
	__shared__ uint8_t utility_buffer_small[SolventBlock::MAX_SOLVENTS_IN_BLOCK];

	// Doubles as block_index_3d!
	const NodeIndex block_origo = BoxGrid::Get3dIndex(blockIdx.x, boxSize_device.boxSizeNM_i);

	BoxState* boxState = sim->boxState;
	const SimParams& simparams = sim->params;
	SimSignals* signals = sim->signals;
	SolventBlock* solventblock_ptr = boxState->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, step);

	//const ForceEnergy forceEnergy = solventblock_ptr->forceEnergies[threadIdx.x];
	const Float3 force = solventblock_ptr->forceEnergies[threadIdx.x].force;
	const float potE = solventblock_ptr->forceEnergies[threadIdx.x].potE;
	
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
		const float mass = tinymolForcefield_device.types[tinyMols_ref.tinymolTypeIndex].mass;

		if constexpr (energyMinimize) {
			const float progress = static_cast<float>(step) / static_cast<float>(simparams.n_steps);
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

			vel_now = vel_now * thermostatScalar_device;

			tinyMols_ref.vel_prev = vel_now;
			tinyMols_ref.force_prev = force;

			// Save pos locally, but only push to box as this kernel ends
			relpos_next = pos_now;
			EngineUtils::LogSolventData(sim->boxparams, potE, block_origo, solventblock.ids[threadIdx.x], solventblock.rel_pos[threadIdx.x], solventActive,
				force, vel_now, step, sim->potE_buffer, sim->traj_buffer, sim->vel_buffer, simparams.data_logging_interval);
		}
	}



	// Push new SolventCoord to global mem
	SolventBlock* const solventblock_next_ptr = boxState->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, step + 1);

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
template __global__ void TinymolIntegrationLoggingAndTransferout<PeriodicBoundaryCondition, true, true>(SimulationDevice*, int64_t);
template __global__ void TinymolIntegrationLoggingAndTransferout<PeriodicBoundaryCondition, true, false>(SimulationDevice*, int64_t);
template __global__ void TinymolIntegrationLoggingAndTransferout<PeriodicBoundaryCondition, false, true>(SimulationDevice*, int64_t);
template __global__ void TinymolIntegrationLoggingAndTransferout<PeriodicBoundaryCondition, false, false>(SimulationDevice*, int64_t);
template __global__ void TinymolIntegrationLoggingAndTransferout<NoBoundaryCondition, true, true>(SimulationDevice*, int64_t);
template __global__ void TinymolIntegrationLoggingAndTransferout<NoBoundaryCondition, true, false>(SimulationDevice*, int64_t);
template __global__ void TinymolIntegrationLoggingAndTransferout<NoBoundaryCondition, false, true>(SimulationDevice*, int64_t);
template __global__ void TinymolIntegrationLoggingAndTransferout<NoBoundaryCondition, false, false>(SimulationDevice*, int64_t);


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
		bridge.loadMeta(&sim->boxConfig.compoundBridges[blockIdx.x]);

		
	}
	__syncthreads();

	// TODO: we dont need to do this for the first compound, as it will always be 0,0,0
	if (threadIdx.x < bridge.n_compounds) {
		// Calculate necessary shift in relative positions for right, so right share the origo with left.
		utility_coord[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelativeShiftBetweenCoordarrays<BoundaryCondition>(boxState->compoundOrigos, step, bridge.compound_ids[0], bridge.compound_ids[threadIdx.x]);
	}


	bridge.loadData(&sim->boxConfig.compoundBridges[blockIdx.x]);
	__syncthreads();

	if (particle_id_bridge < bridge.n_particles) {
		ParticleReference& p_ref = bridge.particle_refs[particle_id_bridge];

		BridgeWarnings::verifyPRefValid(p_ref, bridge);

		positions[threadIdx.x] = boxState->compoundsRelposLm[p_ref.compound_id * MAX_COMPOUND_PARTICLES + p_ref.local_id_compound] + utility_coord[p_ref.compoundid_local_to_bridge].toFloat3();
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

#ifdef FORCE_NAN_CHECK
		if (isinf(force.lenSquared()))
			printf("Bridge inf\n");
#endif
		boxState->compoundsInterimState[p_ref->compound_id].forceEnergyBridge[p_ref->local_id_compound] = ForceEnergy{ force, potE_sum };
	}
}
template __global__ void compoundBridgeKernel<PeriodicBoundaryCondition>(SimulationDevice* sim, int64_t step);
template __global__ void compoundBridgeKernel<NoBoundaryCondition>(SimulationDevice* sim, int64_t step);
#pragma warning (pop)
#pragma warning (pop)
