#pragma once

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "EngineUtils.cuh"
#include "DeviceAlgorithms.cuh"
#include "KernelConstants.cuh"

namespace SolventBlockTransfers {
	/// <summary></summary>
	/// <param name="solventblock">In shared memory</param>
	/// <param name="transferqueues">In shared memory</param>
	/// <param name="relpos_next">Register</param>
	/// <param name="transfermodules">Global memory</param>
	template <typename BoundaryCondition>
	__device__ void transferOut(const SolventBlock& solventblock_current_local, const Coord* const relPositionsNext, STransferQueue* transferqueues, 
		SolventBlockTransfermodule* transfermodules, const NodeIndex& blockId3d, NodeIndex* transferDirectionSharedbuffer) {

		// Parallel insertion in shared memory
		if (threadIdx.x < 6) {
			const NodeIndex myTaskDir = SolventblockTransferUtils::indexToDirection[threadIdx.x];
			for (int i = 0; i < solventblock_current_local.n_solvents; i++) {
				if (transferDirectionSharedbuffer[i] == myTaskDir) {
					bool success = transferqueues[threadIdx.x].addElement(relPositionsNext[i], solventblock_current_local.rel_pos[i], solventblock_current_local.ids[i]);
				}
			}
		}
		__syncthreads();

		// Coaslescing copying to global memory
		for (int queue_index = 0; queue_index < 6; queue_index++) {
			const STransferQueue& queue_local = transferqueues[queue_index];

			if (threadIdx.x < queue_local.n_elements) {

				const NodeIndex transferdir_queue = SolventblockTransferUtils::indexToDirection[queue_index];
				const int blockid_global = EngineUtils::getNewBlockId<BoundaryCondition>(transferdir_queue, blockId3d);
				KernelHelpersWarnings::assertValidBlockId(blockid_global);

				STransferQueue* queue_global = &transfermodules[blockid_global].transfer_queues[queue_index];

				queue_global->rel_positions[threadIdx.x] = queue_local.rel_positions[threadIdx.x] - Coord(transferdir_queue);
				queue_global->ids[threadIdx.x] = queue_local.ids[threadIdx.x];
				queue_global->atomtypeIds[threadIdx.x] = queue_local.atomtypeIds[threadIdx.x];

				KernelHelpersWarnings::transferOutDebug(queue_global, queue_local, transferdir_queue, queue_index);

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
		const Coord* const relPositionsNext, uint8_t* utility_buffer, SolventBlockTransfermodule* remain_transfermodule, const bool remain) {

		// Compute prefix sum to find new index of solvent belonging to thread
		const uint8_t solventindex_new = LAL::computePrefixSum(remain, utility_buffer, solventblock_current_local.n_solvents);


		if (remain) {
			SolventBlock::Transfer(solventblock_current_local, solventblock_next_global, solventindex_new, relPositionsNext[threadIdx.x]);
		}


		const int nsolventsinblock_next = __syncthreads_count(remain);
		if (threadIdx.x == 0) {
			remain_transfermodule->n_remain = nsolventsinblock_next;
			solventblock_next_global->n_solvents = nsolventsinblock_next;	// Doesn't matter, since the transfer kernel handles this. Enabled for debugging now..
		}
	}

	template <typename BoundaryCondition>
	__device__ void transferOutAndCompressRemainders(const SolventBlock& solventblock_current_local, SolventBlock* solventblock_next_global,
		const Coord* const relPositionsNext, uint8_t* utility_buffer, SolventBlockTransfermodule* transfermodules_global, STransferQueue* transferqueues_local, NodeIndex* transferDirectionSharedbuffer) 
	{
		transferDirectionSharedbuffer[threadIdx.x] = threadIdx.x < solventblock_current_local.n_solvents
			? LIMAPOSITIONSYSTEM::getTransferDirection(relPositionsNext[threadIdx.x])
			: NodeIndex{};

		const NodeIndex blockId3d = BoxGrid::Get3dIndex(blockIdx.x, DeviceConstants::boxSize.boxSizeNM_i);
		const int new_blockid = EngineUtils::getNewBlockId<BoundaryCondition>(transferDirectionSharedbuffer[threadIdx.x], blockId3d);
		const bool remain = (blockIdx.x == new_blockid) && threadIdx.x < solventblock_current_local.n_solvents;


		__syncthreads();

		transferOut<BoundaryCondition>(solventblock_current_local, relPositionsNext, transferqueues_local, transfermodules_global, blockId3d, transferDirectionSharedbuffer);

		SolventBlockTransfermodule* remain_transfermodule = &transfermodules_global[blockIdx.x];
		compressRemainers(solventblock_current_local, solventblock_next_global, relPositionsNext, utility_buffer, remain_transfermodule, remain);
	}
}