#pragma once

#include "KernelConstants.cuh"

class NoBoundaryCondition {
public:
	__device__ __host__ void static applyBC(NodeIndex& origo) {}

	__device__ constexpr static NodeIndex applyBC(const NodeIndex& nodeindex, int nodesPerDim) { return nodeindex; }

	__device__ __host__ static void applyHyperpos(const NodeIndex& static_index, NodeIndex& movable_index) {}

	__device__ __host__ static inline void applyHyperposNM(const Float3& static_particle, Float3& movable_particle) {}

	__device__ __host__ static NodeIndex applyHyperpos_Return(const NodeIndex& static_index, const NodeIndex& movable_index) { return movable_index; }
};

class PeriodicBoundaryCondition {
public:
	__device__ constexpr static void applyBC(NodeIndex& origo) {
		origo.x += DeviceConstants::boxSize.blocksPerDim * ((origo.x < 0) - (origo.x >= DeviceConstants::boxSize.blocksPerDim));
		origo.y += DeviceConstants::boxSize.blocksPerDim * ((origo.y < 0) - (origo.y >= DeviceConstants::boxSize.blocksPerDim));
		origo.z += DeviceConstants::boxSize.blocksPerDim * ((origo.z < 0) - (origo.z >= DeviceConstants::boxSize.blocksPerDim));

	}

	__device__ constexpr static NodeIndex applyBC(const NodeIndex& nodeindex, int nodesPerDim) {
		NodeIndex output = nodeindex;
		output.x += nodesPerDim * ((nodeindex.x < 0) - (nodeindex.x >= nodesPerDim));
		output.y += nodesPerDim * ((nodeindex.y < 0) - (nodeindex.y >= nodesPerDim));
		output.z += nodesPerDim * ((nodeindex.z < 0) - (nodeindex.z >= nodesPerDim));
		return output;
	}

	__device__ constexpr static void applyHyperpos(const NodeIndex& static_index, NodeIndex& movable_index) {
		const NodeIndex difference = static_index - movable_index;
		movable_index.x += DeviceConstants::boxSize.blocksPerDim * (difference.x > (DeviceConstants::boxSize.blocksPerDimHalf));		// Dont need to +1 to account of uneven, this is correct (im pretty sure)
		movable_index.x -= DeviceConstants::boxSize.blocksPerDim * (difference.x < -(DeviceConstants::boxSize.blocksPerDimHalf));
        movable_index.y += DeviceConstants::boxSize.blocksPerDim * (difference.y > (DeviceConstants::boxSize.blocksPerDimHalf));
		movable_index.y -= DeviceConstants::boxSize.blocksPerDim * (difference.y < -(DeviceConstants::boxSize.blocksPerDimHalf));
		movable_index.z += DeviceConstants::boxSize.blocksPerDim * (difference.z > (DeviceConstants::boxSize.blocksPerDimHalf));
		movable_index.z -= DeviceConstants::boxSize.blocksPerDim * (difference.z < -(DeviceConstants::boxSize.blocksPerDimHalf));
	}

	__device__ constexpr static NodeIndex applyHyperpos_Return(const NodeIndex& static_index, const NodeIndex& movable_index) {
		NodeIndex hyperIndex = movable_index;
        const int halfBox = DeviceConstants::boxSize.blocksPerDimHalf;
		const NodeIndex difference = static_index - movable_index;

		hyperIndex.x += DeviceConstants::boxSize.blocksPerDim * ((difference.x > halfBox) - (difference.x < -halfBox));
		hyperIndex.y += DeviceConstants::boxSize.blocksPerDim * ((difference.y > halfBox) - (difference.y < -halfBox));
		hyperIndex.z += DeviceConstants::boxSize.blocksPerDim * ((difference.z > halfBox) - (difference.z < -halfBox));

		return hyperIndex;
	}

	__device__ constexpr static inline void applyHyperposNM(const Float3& static_particle, Float3& movable_particle) {
		const float boxlenhalf_nm = DeviceConstants::boxSize.boxSizeNM_f / 2.f;

		movable_particle.x += DeviceConstants::boxSize.boxSizeNM_f * ((static_particle.x - movable_particle.x) > boxlenhalf_nm);
		movable_particle.x -= DeviceConstants::boxSize.boxSizeNM_f * ((static_particle.x - movable_particle.x) < -boxlenhalf_nm);
		movable_particle.y += DeviceConstants::boxSize.boxSizeNM_f * ((static_particle.y - movable_particle.y) > boxlenhalf_nm);
		movable_particle.y -= DeviceConstants::boxSize.boxSizeNM_f * ((static_particle.y - movable_particle.y) < -boxlenhalf_nm);
		movable_particle.z += DeviceConstants::boxSize.boxSizeNM_f * ((static_particle.z - movable_particle.z) > boxlenhalf_nm);
		movable_particle.z -= DeviceConstants::boxSize.boxSizeNM_f * ((static_particle.z - movable_particle.z) < -boxlenhalf_nm);
	}

	__device__ constexpr static void applyBCNM(Float3& current_position) {	// Only changes position if position is outside of box;		
		current_position.x += DeviceConstants::boxSize.boxSizeNM_f * (current_position.x < 0.f);
		current_position.x -= DeviceConstants::boxSize.boxSizeNM_f * (current_position.x > DeviceConstants::boxSize.boxSizeNM_f);
		current_position.y += DeviceConstants::boxSize.boxSizeNM_f * (current_position.y < 0.f);
		current_position.y -= DeviceConstants::boxSize.boxSizeNM_f * (current_position.y > DeviceConstants::boxSize.boxSizeNM_f);
		current_position.z += DeviceConstants::boxSize.boxSizeNM_f * (current_position.z < 0.f);
		current_position.z -= DeviceConstants::boxSize.boxSizeNM_f * (current_position.z > DeviceConstants::boxSize.boxSizeNM_f);
	}
};