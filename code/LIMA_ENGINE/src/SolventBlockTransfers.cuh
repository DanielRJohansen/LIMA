#pragma once

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "EngineUtils.cuh"
#include "DeviceAlgorithms.cuh"
#include "KernelConstants.cuh"

//namespace SolventBlockTransfers {
//	/// <summary></summary>
//	/// <param name="solventblock">In shared memory</param>
//	/// <param name="transferqueues">In shared memory</param>
//	/// <param name="relpos_next">Register</param>
//	/// <param name="transfermodules">Global memory</param>
//	template <typename BoundaryCondition>
//	__device__ void transferOut(const SolventBlock& solventblock_current_local, const Coord* const relPositionsNext, STransferQueue* transferqueues, 
//		SolventBlockTransfermodule* transfermodules, const NodeIndex& blockId3d, NodeIndex* transferDirectionSharedbuffer) {
//
//		// Parallel insertion in shared memory
//		if (threadIdx.x < 6) {
//			const NodeIndex myTaskDir = SolventblockTransferUtils::indexToDirection[threadIdx.x];
//			for (int i = 0; i < solventblock_current_local.n_solvents; i++) {
//				if (transferDirectionSharedbuffer[i] == myTaskDir) {
//					bool success = transferqueues[threadIdx.x].addElement(relPositionsNext[i], solventblock_current_local.rel_pos[i], solventblock_current_local.ids[i]);
//				}
//			}
//		}
//		__syncthreads();
//
//		// Coaslescing copying to global memory
//		for (int queue_index = 0; queue_index < 6; queue_index++) {
//			const STransferQueue& queue_local = transferqueues[queue_index];
//
//			if (threadIdx.x < queue_local.n_elements) {
//
//				const NodeIndex transferdir_queue = SolventblockTransferUtils::indexToDirection[queue_index];
//				const int blockid_global = EngineUtils::getNewBlockId<BoundaryCondition>(transferdir_queue, blockId3d);
//				KernelHelpersWarnings::assertValidBlockId(blockid_global);
//
//				STransferQueue* queue_global = &transfermodules[blockid_global].transfer_queues[queue_index];
//
//				queue_global->rel_positions[threadIdx.x] = queue_local.rel_positions[threadIdx.x] - Coord(transferdir_queue);
//				queue_global->ids[threadIdx.x] = queue_local.ids[threadIdx.x];
//				queue_global->atomtypeIds[threadIdx.x] = queue_local.atomtypeIds[threadIdx.x];
//
//				KernelHelpersWarnings::transferOutDebug(queue_global, queue_local, transferdir_queue, queue_index);
//
//				// Only set n_elements if we get here, meaning atleast 1 new element. Otherwise it will just remain 0
//				if (threadIdx.x == 0) { queue_global->n_elements = queue_local.n_elements; }
//			}
//		}
//		__syncthreads();
//	}
//
//
//	/// <summary>
//	/// Must be run AFTER transferOut, as it erases all information about the transferring solvents
//	/// </summary>
//	/// <param name="solventblock_current">Solventblock in shared memory</param>
//	/// <param name="solventblock_next">Solventblock belong to cudablock at next step in global memory</param>
//	/// <param name="relpos_next"></param>
//	/// <param name="utility_buffer">Buffer of min size MAX_SOLVENTS_IN_BLOCK, maybe more for computing prefix sum</param>
//	/// <param name="remain_transfermodule">Transfermodule belonging to cudablock</param>
//	__device__ void compressRemainers(const SolventBlock& solventblock_current_local, SolventBlock* solventblock_next_global,
//		const Coord* const relPositionsNext, uint8_t* utility_buffer, SolventBlockTransfermodule* remain_transfermodule, const bool remain) {
//
//		// Compute prefix sum to find new index of solvent belonging to thread
//		const uint8_t solventindex_new = LAL::computePrefixSum(remain, utility_buffer, solventblock_current_local.n_solvents);
//
//
//		if (remain) {
//			SolventBlock::Transfer(solventblock_current_local, solventblock_next_global, solventindex_new, relPositionsNext[threadIdx.x]);
//		}
//
//
//		const int nsolventsinblock_next = __syncthreads_count(remain);
//		if (threadIdx.x == 0) {
//			remain_transfermodule->n_remain = nsolventsinblock_next;
//			solventblock_next_global->n_solvents = nsolventsinblock_next;	// Doesn't matter, since the transfer kernel handles this. Enabled for debugging now..
//		}
//	}
//
//	template <typename BoundaryCondition>
//	__device__ void transferOutAndCompressRemainders(const SolventBlock& solventblock_current_local, SolventBlock* solventblock_next_global,
//		const Coord* const relPositionsNext, uint8_t* utility_buffer, SolventBlockTransfermodule* transfermodules_global, STransferQueue* transferqueues_local, NodeIndex* transferDirectionSharedbuffer) 
//	{
//		transferDirectionSharedbuffer[threadIdx.x] = threadIdx.x < solventblock_current_local.n_solvents
//			? LIMAPOSITIONSYSTEM::getTransferDirection(relPositionsNext[threadIdx.x])
//			: NodeIndex{};
//
//		const NodeIndex blockId3d = BoxGrid::Get3dIndex(blockIdx.x, DeviceConstants::boxSize.boxSizeNM_i);
//		const int new_blockid = EngineUtils::getNewBlockId<BoundaryCondition>(transferDirectionSharedbuffer[threadIdx.x], blockId3d);
//		const bool remain = (blockIdx.x == new_blockid) && threadIdx.x < solventblock_current_local.n_solvents;
//
//
//		__syncthreads();
//
//		transferOut<BoundaryCondition>(solventblock_current_local, relPositionsNext, transferqueues_local, transfermodules_global, blockId3d, transferDirectionSharedbuffer);
//
//		SolventBlockTransfermodule* remain_transfermodule = &transfermodules_global[blockIdx.x];
//		compressRemainers(solventblock_current_local, solventblock_next_global, relPositionsNext, utility_buffer, remain_transfermodule, remain);
//	}
//}


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
	__syncthreads();

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

				if constexpr (INDEXING_CHECKS)
					if (myBondgroupCount >= TinymolTransferModule::maxOutgoingBondgroups)
						printf("Too many bondgroups in one direction");

				bondgroupIdsThisDirection[threadIdx.x * TinymolTransferModule::maxOutgoingBondgroups + myBondgroupCount] = bondgroupIndex;
				bondgroupsFirstAtomindexInThisDirection[threadIdx.x * TinymolTransferModule::maxOutgoingBondgroups + myBondgroupCount] = myParticleCount;
				myBondgroupCount++;
				myParticleCount += nParticlesInBondgroups[bondgroupIndex];
			}
		}

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
		if constexpr (INDEXING_CHECKS)
			if (targetBlockId < 0 || targetBlockId > BoxGrid::BlocksTotal(DeviceConstants::boxSize.blocksPerDim))
				printf("Target block was out of bounds");
		
		const Coord relposShift = Coord{ -direction.toFloat3() };

		if (targetBlockId >= DeviceConstants::boxSize.blocksPerDim * DeviceConstants::boxSize.blocksPerDim * DeviceConstants::boxSize.blocksPerDim)
			printf("Target block was out of bounds");

		// Write results directly to global mem
		if (threadIdx.x == 0) {
			tinymolTransferModule.nIncomingParticles[targetBlockId * 6 + directionIndex] = nParticlesThisDirection[directionIndex];
			tinymolTransferModule.nIncomingBondgroups[targetBlockId * 6 + directionIndex] = nBondgroupsThisDirection[directionIndex];
		}
		// Each thread takes 1 bondgroup
		if (threadIdx.x < nBondgroupsThisDirection[directionIndex]) {
			const int bondgroupIndex = bondgroupIdsThisDirection[directionIndex * TinymolTransferModule::maxOutgoingBondgroups + threadIdx.x];
			const int bondgroupTargetIndex = targetBlockId * 6 * TinymolTransferModule::maxOutgoingBondgroups
				+ directionIndex * TinymolTransferModule::maxOutgoingBondgroups
				+ threadIdx.x;
			tinymolTransferModule.incomingBondgroups[bondgroupTargetIndex] = solventblockGlobalPtr->bondgroups[bondgroupIndex];
			tinymolTransferModule.incomingBondgroupsParticlesOffset[bondgroupTargetIndex] = bondgroupsFirstAtomindexInThisDirection[directionIndex * TinymolTransferModule::maxOutgoingBondgroups + threadIdx.x];

			const int destOffsetOfFirstParticle = bondgroupsFirstAtomindexInThisDirection[directionIndex * TinymolTransferModule::maxOutgoingBondgroups + threadIdx.x];
			const int sourceOffsetOfFirstParticle = solventblockGlobalPtr->bondgroupsFirstAtomindexInSolventblock[bondgroupIndex];
			for (int particleIndexInBondgroup = 0; particleIndexInBondgroup < nParticlesInBondgroups[bondgroupIndex]; particleIndexInBondgroup++) 
			{
				const int particleDestIndex = 
					targetBlockId * 6 * TinymolTransferModule::maxOutgoingParticles
					+ directionIndex * TinymolTransferModule::maxOutgoingParticles
					+ destOffsetOfFirstParticle + particleIndexInBondgroup;
				const int particleSourceIndex = sourceOffsetOfFirstParticle + particleIndexInBondgroup;

				if constexpr (INDEXING_CHECKS) {
					if (particleSourceIndex >= SolventBlock::MAX_SOLVENTS_IN_BLOCK || particleSourceIndex < 0)
						printf("Source index out of bounds\n");
					if (particleDestIndex > (targetBlockId+1) * 6 * TinymolTransferModule::maxOutgoingParticles || particleDestIndex < 0)
						printf("Dest index out of bounds {%d}\n", particleDestIndex);
				}

				tinymolTransferModule.incomingPositions[particleDestIndex] = solventblockGlobalPtr->rel_pos[particleSourceIndex] + relposShift;
				tinymolTransferModule.incomingIds[particleDestIndex] = solventblockGlobalPtr->ids[particleSourceIndex];
				tinymolTransferModule.incomingAtomtypeIds[particleDestIndex] = solventblockGlobalPtr->atomtypeIds[particleSourceIndex];
				tinymolTransferModule.incomingBondgroupIds[particleDestIndex] = threadIdx.x;
				tinymolTransferModule.incomingStates[particleDestIndex] = solventblockGlobalPtr->states[particleSourceIndex];
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
					solventblockGlobalPtr->states[nParticlesRemaining] = solventblockGlobalPtr->states[srcIndex];
					nParticlesRemaining++;
				}
			}
		}

		if constexpr (INDEXING_CHECKS) {
			if (nParticlesRemaining > solventblockGlobalPtr->nParticles)
				printf("More particles remain that we started with\n");
		}

		solventblockGlobalPtr->nBondgroups = nBondgroupsRemaining;
		solventblockGlobalPtr->nParticles = nParticlesRemaining;
	}
}

__global__ void SolventTransferKernel(SimulationDevice* sim, int64_t _step, const TinymolTransferModule tinymolTransferModule) {
	__shared__ int nParticlesInBlock;
	__shared__ int nBondgroupsInBlock;

	const int solventblockId = blockIdx.x;
	const int stepToLoadFrom = _step + 1;
	SolventBlock* const solventblockGlobalPtr = SolventBlocksCircularQueue::getBlockPtr(sim->boxState.solventblockgrid_circularqueue, DeviceConstants::boxSize.boxSizeNM_i, solventblockId, stepToLoadFrom);

	if (threadIdx.x == 0) {
		nParticlesInBlock = solventblockGlobalPtr->nParticles;
		nBondgroupsInBlock = solventblockGlobalPtr->nBondgroups;

		if (nParticlesInBlock > SolventBlock::MAX_SOLVENTS_IN_BLOCK || nBondgroupsInBlock > SolventBlock::maxBondgroups)
			printf("Error %d %d\n", nParticlesInBlock, nBondgroupsInBlock);

	}
	__syncthreads();

	// All threads loop over the directions. If it has an incoming particle, it copies it directy to global mem
	for (int directionIndex = 0; directionIndex < 6; directionIndex++) {
		const int nIncomingParticles = tinymolTransferModule.nIncomingParticles[blockIdx.x * 6 + directionIndex];
		const int nIncomingBondgroups = tinymolTransferModule.nIncomingBondgroups[blockIdx.x * 6 + directionIndex];

		if constexpr (INDEXING_CHECKS) {
			if (nParticlesInBlock + nIncomingParticles >= SolventBlock::MAX_SOLVENTS_IN_BLOCK)
				printf("Inserting too many particles %d %d did %d pid %d\n", nParticlesInBlock, nIncomingParticles, directionIndex, threadIdx.x);
			if (nBondgroupsInBlock + nIncomingBondgroups >= SolventBlock::maxBondgroups)
				printf("Inserting too many bondgroups %d %d did %d pid %d\n", nBondgroupsInBlock, nIncomingBondgroups, directionIndex, threadIdx.x);
		}


		if (threadIdx.x < nIncomingParticles) {
			const int srcIndex = solventblockId * 6 * TinymolTransferModule::maxOutgoingParticles
				+ directionIndex * TinymolTransferModule::maxOutgoingParticles
				+ threadIdx.x;
			solventblockGlobalPtr->rel_pos[nParticlesInBlock + threadIdx.x] = tinymolTransferModule.incomingPositions[srcIndex];
			solventblockGlobalPtr->ids[nParticlesInBlock + threadIdx.x] = tinymolTransferModule.incomingIds[srcIndex];
			solventblockGlobalPtr->atomtypeIds[nParticlesInBlock + threadIdx.x] = tinymolTransferModule.incomingAtomtypeIds[srcIndex];
			solventblockGlobalPtr->particlesBondgroupIds[nParticlesInBlock + threadIdx.x] = tinymolTransferModule.incomingBondgroupIds[srcIndex] + nBondgroupsInBlock;
			solventblockGlobalPtr->states[nParticlesInBlock + threadIdx.x] = tinymolTransferModule.incomingStates[srcIndex];
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
	__syncthreads();

	// Finally we write to global mem how many particles there are in the block now
	if (threadIdx.x == 0) {
		solventblockGlobalPtr->nParticles = nParticlesInBlock;
		solventblockGlobalPtr->nBondgroups = nBondgroupsInBlock;
	}
}