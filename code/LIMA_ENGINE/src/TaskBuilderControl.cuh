#pragma once

#include "EngineBodies.cuh"
#include "CubWrappers.h"

#include <array>
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


// Pass to GPU
struct TaskBuilderControlContents {
	static const int maxTasksPerSc = 256;

	std::array<float4, 4>* superclusterPositionSpheres;

	ParticlesBondedToParticle* particlesBondedToParticle;
	PclustersBondedToPcluster* pclustersBondedToPcluster;

	int* nInteractionsOwned;
	InteractionToken* interactionsOwned;
	int* nInteractionsNonowned;
	int* scIdsQueryNonowned;
	int* nResults;

	int* nResultsPrefixsum;
	int* nQueryBuffersPrefixsum;
};

// Keep on CPU
class TaskBuilderControl {

public:
	const int nSuperclustersUpperbound;
	TaskBuilderControlContents contents;

	TaskBuilderControl(const TaskBuilderControl&) = delete;
	TaskBuilderControl& operator=(const TaskBuilderControl&) = delete;
	TaskBuilderControl(int nSuperclustersUpperbound, const std::vector<ParticlesBondedToParticle>& particlesBondedToParticle,
		const std::vector<PclustersBondedToPcluster>& pclustersBondedToPcluster, cudaStream_t stream)
		: nSuperclustersUpperbound(nSuperclustersUpperbound)
	{
		contents.particlesBondedToParticle = GenericCopyToDevice(particlesBondedToParticle);
		contents.pclustersBondedToPcluster = GenericCopyToDevice(pclustersBondedToPcluster);

		cudaMalloc(&contents.superclusterPositionSpheres, sizeof(std::array<float4, 4>) * nSuperclustersUpperbound);
		cudaMalloc(&contents.nInteractionsOwned, sizeof(int) * (nSuperclustersUpperbound + 1));
		cudaMalloc(&contents.interactionsOwned, sizeof(InteractionToken) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		cudaMalloc(&contents.nInteractionsNonowned, sizeof(int) * nSuperclustersUpperbound);
		cudaMalloc(&contents.scIdsQueryNonowned, sizeof(int) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		cudaMalloc(&contents.nResults, sizeof(int) * (nSuperclustersUpperbound + 1));

		cudaMalloc(&contents.nResultsPrefixsum, sizeof(int) * (nSuperclustersUpperbound + 1));
		cudaMalloc(&contents.nQueryBuffersPrefixsum, sizeof(int) * (nSuperclustersUpperbound + 1));

		Reset(stream);
	}

	void Reset(cudaStream_t stream) {
		cudaMemsetAsync(contents.interactionsOwned, 0xFF,
			sizeof(InteractionToken) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound, stream);
		//cudaMemset(contents.scIdsQueryNonowned, 0xFF, sizeof(int) * TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound);
		CubWrappers::FillN(contents.scIdsQueryNonowned,
			TaskBuilderControlContents::maxTasksPerSc * nSuperclustersUpperbound, INT_MAX, stream);
	}

	~TaskBuilderControl() {
		cudaFree(contents.particlesBondedToParticle);
		cudaFree(contents.pclustersBondedToPcluster);
		cudaFree(contents.superclusterPositionSpheres);
		cudaFree(contents.nInteractionsOwned);
		cudaFree(contents.interactionsOwned);
		cudaFree(contents.nInteractionsNonowned);
		cudaFree(contents.scIdsQueryNonowned);
		cudaFree(contents.nResults);

		cudaFree(contents.nResultsPrefixsum);
		cudaFree(contents.nQueryBuffersPrefixsum);
	}
};
