#pragma once

#include "Bodies.cuh"
#include "LimaTypes.cuh"

#include "Constants.h"

// Highest concentration in smtv test is only 48 tinymols, but could be larger with more ions and less solvent
static const int MAX_PARTICLES_IN_BOXGRIDNODE = 64;




namespace BoxGrid {
	static const int blocksizeNM = 1;
	constexpr int NodesPerDim(int boxlenNM) { return boxlenNM; }
	constexpr Int3 NodesPerDim(Int3 boxlenNM) {
		return Int3{ NodesPerDim(boxlenNM.x), NodesPerDim(boxlenNM.y), NodesPerDim(boxlenNM.z) };
	}
	constexpr int BlocksTotal(Int3 blocksPerDim) { return NodesPerDim(blocksPerDim.x) * NodesPerDim(blocksPerDim.y) * NodesPerDim(blocksPerDim.z); }

	constexpr int Get1dIndex(const NodeIndex& index3d, Int3 boxSizeNM) {
		if constexpr (INDEXING_CHECKS) {
			if (index3d.Min() < 0 || index3d.x >= boxSizeNM.x || index3d.y >= boxSizeNM.y || index3d.z >= boxSizeNM.z) {
				printf("Illegal nodexindex fetch, blockindex: %d %d %d\n", index3d.x, index3d.y, index3d.z);
			}
		}
		return index3d.x + index3d.y * NodesPerDim(boxSizeNM.x) + index3d.z * NodesPerDim(boxSizeNM.x) * NodesPerDim(boxSizeNM.y);
	}

	constexpr NodeIndex Get3dIndexWithNNodes(int index1d, Int3 npd) {
		int z = index1d / (npd.x * npd.y);
		index1d -= z * npd.x * npd.y;
		int y = index1d / npd.x;
		index1d -= y * npd.x;
		int x = index1d;
		return NodeIndex{ x, y, z };
	}
	constexpr NodeIndex Get3dIndex(int index1d, const Int3& boxlenNM) {
		return Get3dIndexWithNNodes(index1d, NodesPerDim(boxlenNM));
	}




	/// <summary>
	/// Create a BoxGrid of NodeTyoe. This function assumes the type NodeType does not need to be initialized
	/// </summary>
	template <typename NodeType>
	__host__ static NodeType* MallocOnDevice(int boxSizeNM) {
		const int blocksTotal = BoxGrid::NodesPerDim(boxSizeNM) * BoxGrid::NodesPerDim(boxSizeNM) * BoxGrid::NodesPerDim(boxSizeNM);

		NodeType* grid_dev = nullptr;
		cudaMalloc(&grid_dev, sizeof(NodeType) * blocksTotal);
		cudaMemset(grid_dev, 0, sizeof(NodeType) * blocksTotal);
		return grid_dev;
	}

	// This function assumes the user has used PBC
	template <typename NodeType>
	__device__ __host__ inline NodeType* GetNodePtr(NodeType* grid, const int index1d) {
		return &grid[index1d];
	}

	namespace TinymolBlockAdjacency {
        static const int nNearbyBlocks = 32;

		
		struct BlockRef {
			/*int blockId = -1;
			Float3 relShift{};*/
			uint16_t blockId = 0xFFFFFFFF;
			Float3Compressed relShift{};
		};


		// Constant
		struct NearbyBlocksSequences {
			// Optim Pack this info into a single uint32_t?
			struct Sequence{
				int blockIndexStart = -1;
				int nBlocks = -1;
			};
			static const int maxSequences = 22;
			
			Sequence sequences[maxSequences];
			int nSequences = 0;
		};

		struct NearbyBlocksSequencesParticles {
			struct Sequence {
				int indexOfFirstParticleInSequence = -1;
				int nParticlesInSequence = -1;
			};
			Sequence sequences[NearbyBlocksSequences::maxSequences];
			int nSequences = 0;
		};


		// Returns a cudapointer to the data
        BlockRef* PrecomputeNeabyBlockIds(Int3 boxlenNM, float ljCutoffNm);
		NearbyBlocksSequences* PrecomputeNearbyBlockSequences(Int3 boxlenNM);

		__device__ inline const BlockRef* GetPtrToNearbyBlockids(int blockId, const BlockRef* const nearbyBlockIdsData) {
			return &nearbyBlockIdsData[blockId * nNearbyBlocks];
		}
	}
};


namespace SolventBlocksCircularQueue {
	//static const int STEPS_PER_SOLVENTBLOCKTRANSFER = 5;	// If we go below 2, we might see issue in solventtransfers
	//static const int SOLVENTBLOCK_TRANSFERSTEP = STEPS_PER_SOLVENTBLOCKTRANSFER - 1;
	//static const int queue_len = STEPS_PER_SOLVENTBLOCKTRANSFER;


	//constexpr int nElementsTotal(Int3 boxlenNM) {
	//	return BoxGrid::BlocksTotal(boxlenNM) * queue_len;
	//}

	//static std::vector<SolventBlock> createQueue(Int3 boxlenNM) {
	//	return std::vector<SolventBlock>(nElementsTotal(boxlenNM));
	//}

	//constexpr bool isTransferStep(int64_t step) {
	//	return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == SOLVENTBLOCK_TRANSFERSTEP;
	//}
	//constexpr bool isFirstStepAfterTransfer(int64_t step) {
	//	return (step % STEPS_PER_SOLVENTBLOCKTRANSFER) == 0;
	//}

	//// This function assumes the user has used PBC
	//__device__ __host__ static SolventBlock* getBlockPtr(SolventBlock* queue, Int3 boxlenNM, const size_t index1d, const size_t step) {
	//	const size_t step_offset = (step % queue_len) * BoxGrid::BlocksTotal(boxlenNM);
	//	return &queue[index1d + step_offset];
	//}

	//__host__ static SolventBlock& GetBlockRef(std::vector<SolventBlock>& queue, NodeIndex index3d, const int64_t step, Int3 boxSizeNm) {
	//	if (index3d.x >= boxSizeNm.x || index3d.y >= boxSizeNm.y || index3d.z >= boxSizeNm.z
	//		|| index3d.x < 0 || index3d.y < 0 || index3d.z < 0) {
	//		throw std::runtime_error("Bad 3d index for blockptr\n");
	//	}
	//	return *getBlockPtr(queue.data(), boxSizeNm, BoxGrid::Get1dIndex(index3d, boxSizeNm), step);
	//}
};

