#pragma once

#include "EngineBodies.cuh"
#include "CudaBuffer.h"
#include "CubWrappers.h"

#include <climits>
#include <vector>


class InteractionToken {
	uint32_t data{};
public:
	constexpr InteractionToken() {}
	constexpr InteractionToken(int queryId, bool useNointeractionMatrix) {
		data = ((uint32_t)queryId << 1) | (uint32_t)useNointeractionMatrix;
	}
	constexpr int GetQueryId() const { return (int)(data >> 1); }
	constexpr bool UseNointeractionMatrix() const { return data & 1; }
};

// A supercluster that a supercluster owns an interaction with (queryScId >= its own id), see FindNeighborsKernel
struct ScNeighbor {
	int queryScId;
	uint16_t quarterMask;	// Bit ownQuarter*4+queryQuarter set if that 4x4 block has a pair within the list radius
	uint16_t bonded;		// The superclusters contain bonded particles, or are the same supercluster
};


// Pass to GPU. The interaction lists are only used, and allocated, in EM
struct TaskBuilderControlContents {
	static const int maxTasksPerSc = 256;

	ParticlesBondedToParticle* particlesBondedToParticle = nullptr;
	PclustersBondedToPcluster* pclustersBondedToPcluster = nullptr;

	int* nInteractionsOwned = nullptr;
	InteractionToken* interactionsOwned = nullptr;
	int* nInteractionsNonowned = nullptr;
	int* scIdsQueryNonowned = nullptr;
	int* nResults = nullptr;

	int* nResultsPrefixsum = nullptr;
	int* nQueryBuffersPrefixsum = nullptr;
};

// Keep on CPU
class TaskBuilderControl {

public:
	const int nSuperclustersUpperbound;
	TaskBuilderControlContents contents;

	// Neighbor search, used by both MD and EM
	CudaBuffer<float4> scSpheres;			// Bounding sphere of each supercluster
	CudaBuffer<int> scCells;				// The grid cell each supercluster belongs to
	CudaBuffer<uint16_t> scValidMasks;		// Bit i set if particle i is not padding
	CudaBuffer<float4> cellMin;				// AABB of the particles of each cell's superclusters
	CudaBuffer<float4> cellMax;
	CudaBuffer<ScNeighbor> neighbors;		// maxTasksPerSc per supercluster
	CudaBuffer<int> nNeighbors;
	CudaBuffer<int> overflow;

	// MD
	CudaBuffer<int> entryCounts;			// nSuperclusters+1, for the scan
	CudaBuffer<int> entryStarts;

	TaskBuilderControl(const TaskBuilderControl&) = delete;
	TaskBuilderControl& operator=(const TaskBuilderControl&) = delete;
	TaskBuilderControl(int nSuperclustersUpperbound, const std::vector<ParticlesBondedToParticle>& particlesBondedToParticle,
		const std::vector<PclustersBondedToPcluster>& pclustersBondedToPcluster)
		: nSuperclustersUpperbound(nSuperclustersUpperbound)
	{
		contents.particlesBondedToParticle = GenericCopyToDevice(particlesBondedToParticle);
		contents.pclustersBondedToPcluster = GenericCopyToDevice(pclustersBondedToPcluster);
	}

	// The interaction lists take maxTasksPerSc entries per supercluster, so they are only allocated once EM needs them
	void AllocateEmBuffers() {
		if (contents.interactionsOwned)
			return;
		cudaMalloc(&contents.nInteractionsOwned, sizeof(int) * (nSuperclustersUpperbound + 1));
		cudaMalloc(&contents.interactionsOwned, sizeof(InteractionToken) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		cudaMalloc(&contents.nInteractionsNonowned, sizeof(int) * nSuperclustersUpperbound);
		cudaMalloc(&contents.scIdsQueryNonowned, sizeof(int) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		cudaMalloc(&contents.nResults, sizeof(int) * (nSuperclustersUpperbound + 1));

		cudaMalloc(&contents.nResultsPrefixsum, sizeof(int) * (nSuperclustersUpperbound + 1));
		cudaMalloc(&contents.nQueryBuffersPrefixsum, sizeof(int) * (nSuperclustersUpperbound + 1));
	}

	// EM only
	void Reset(cudaStream_t stream) {
		cudaMemsetAsync(contents.interactionsOwned, 0xFF,
			sizeof(InteractionToken) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound, stream);
		CubWrappers::FillN(contents.scIdsQueryNonowned,
			TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound, INT_MAX, stream);
		cudaMemsetAsync(contents.nInteractionsNonowned, 0, sizeof(int) * nSuperclustersUpperbound, stream);
	}

	~TaskBuilderControl() {
		cudaFree(contents.particlesBondedToParticle);
		cudaFree(contents.pclustersBondedToPcluster);
		cudaFree(contents.nInteractionsOwned);
		cudaFree(contents.interactionsOwned);
		cudaFree(contents.nInteractionsNonowned);
		cudaFree(contents.scIdsQueryNonowned);
		cudaFree(contents.nResults);

		cudaFree(contents.nResultsPrefixsum);
		cudaFree(contents.nQueryBuffersPrefixsum);
	}
};
