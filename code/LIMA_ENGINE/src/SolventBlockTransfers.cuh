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
					queue_local.ids[threadIdx.x], queue_local.atomtypeIds[threadIdx.x]);

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
		const uint8_t solventindex_new = LAL::computePrefixSum(remain, utility_buffer, solventblock_current_local.n_solvents);


		if (remain) {
			// solventindex_new is only valid for those who remain, the rest *might* have an index of -1
			/*solventblock_next_global->rel_pos[solventindex_new] = relpos_next;
			solventblock_next_global->ids[solventindex_new] = solventblock_current_local.ids[threadIdx.x];*/
			SolventBlock::Transfer(solventblock_current_local, solventblock_next_global, solventindex_new, relpos_next);
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

		const NodeIndex blockId3d = BoxGrid::Get3dIndex(blockIdx.x, boxSize_device.boxSizeNM_i);
		const NodeIndex transfer_dir = threadIdx.x < solventblock_current_local.n_solvents
			? LIMAPOSITIONSYSTEM::getTransferDirection(relpos_next)
			: NodeIndex{};

		const int new_blockid = EngineUtils::getNewBlockId<BoundaryCondition>(transfer_dir, blockId3d);
		const bool remain = (blockIdx.x == new_blockid) && threadIdx.x < solventblock_current_local.n_solvents;


		transferOut<BoundaryCondition>(transfer_dir, solventblock_current_local, new_blockid, relpos_next, transferqueues_local, transfermodules_global, blockId3d);

		SolventBlockTransfermodule* remain_transfermodule = &transfermodules_global[blockIdx.x];
		compressRemainers(solventblock_current_local, solventblock_next_global, relpos_next, utility_buffer, remain_transfermodule, remain);
	}
}