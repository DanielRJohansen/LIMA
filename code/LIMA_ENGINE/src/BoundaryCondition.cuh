#pragma once

#include "KernelConstants.cuh"

class NoBoundaryCondition {
public:
	__device__ __host__ void static applyBC(NodeIndex& origo) {}

	__device__ constexpr static NodeIndex applyBC(const NodeIndex& nodeindex, const Int3& nodesPerDim) { return nodeindex; }

	__device__ __host__ static void applyHyperpos(const NodeIndex& static_index, NodeIndex& movable_index) {}

	__device__ __host__ static inline void applyHyperposNM(const Float3& static_particle, Float3& movable_particle) {}

	__device__ __host__ static NodeIndex applyHyperpos_Return(const NodeIndex& static_index, const NodeIndex& movable_index) { return movable_index; }
};

class PeriodicBoundaryCondition {
public:
	__device__ constexpr static void applyBC(NodeIndex& origo) {
		origo.x += DeviceConstants::boxSize.blocksPerDim.x * ((origo.x < 0) - (origo.x >= DeviceConstants::boxSize.blocksPerDim.x));
		origo.y += DeviceConstants::boxSize.blocksPerDim.y * ((origo.y < 0) - (origo.y >= DeviceConstants::boxSize.blocksPerDim.y));
		origo.z += DeviceConstants::boxSize.blocksPerDim.z * ((origo.z < 0) - (origo.z >= DeviceConstants::boxSize.blocksPerDim.z));

	}

	__device__ constexpr static NodeIndex applyBC(const NodeIndex& nodeindex, const Int3& gridDim) {
		NodeIndex output = nodeindex;
		output.x += gridDim.x * ((nodeindex.x < 0) - (nodeindex.x >= gridDim.x));
		output.y += gridDim.y * ((nodeindex.y < 0) - (nodeindex.y >= gridDim.y));
		output.z += gridDim.z * ((nodeindex.z < 0) - (nodeindex.z >= gridDim.z));
		return output;
	}

	__device__ constexpr static void applyHyperpos(const NodeIndex& static_index, NodeIndex& movable_index) {
		const NodeIndex difference = static_index - movable_index;
		movable_index.x += DeviceConstants::boxSize.blocksPerDim.x * (difference.x > (DeviceConstants::boxSize.blocksPerDimHalf.x));		// Dont need to +1 to account of uneven, this is correct (im pretty sure)
		movable_index.x -= DeviceConstants::boxSize.blocksPerDim.x * (difference.x < -(DeviceConstants::boxSize.blocksPerDimHalf.x));
        movable_index.y += DeviceConstants::boxSize.blocksPerDim.y * (difference.y > (DeviceConstants::boxSize.blocksPerDimHalf.y));
		movable_index.y -= DeviceConstants::boxSize.blocksPerDim.y * (difference.y < -(DeviceConstants::boxSize.blocksPerDimHalf.y));
		movable_index.z += DeviceConstants::boxSize.blocksPerDim.z * (difference.z > (DeviceConstants::boxSize.blocksPerDimHalf.z));
		movable_index.z -= DeviceConstants::boxSize.blocksPerDim.z * (difference.z < -(DeviceConstants::boxSize.blocksPerDimHalf.z));
	}

	__device__ constexpr static NodeIndex applyHyperpos_Return(const NodeIndex& static_index, const NodeIndex& movable_index) {
		NodeIndex hyperIndex = movable_index;
        const Int3 halfBox = DeviceConstants::boxSize.blocksPerDimHalf;
		const NodeIndex difference = static_index - movable_index;

		hyperIndex.x += DeviceConstants::boxSize.blocksPerDim.x * ((difference.x > halfBox.x) - (difference.x < -halfBox.x));
		hyperIndex.y += DeviceConstants::boxSize.blocksPerDim.y * ((difference.y > halfBox.y) - (difference.y < -halfBox.y));
		hyperIndex.z += DeviceConstants::boxSize.blocksPerDim.z * ((difference.z > halfBox.z) - (difference.z < -halfBox.z));

		return hyperIndex;
	}

	__device__ constexpr static inline void applyHyperposNM(const Float3& static_particle, Float3& movable_particle) {
		const Float3 boxlenhalf_nm = DeviceConstants::boxSize.boxSizeNM_f * 0.5f;

		movable_particle.x += DeviceConstants::boxSize.boxSizeNM_f.x * ((static_particle.x - movable_particle.x) > boxlenhalf_nm.x);
		movable_particle.x -= DeviceConstants::boxSize.boxSizeNM_f.x * ((static_particle.x - movable_particle.x) < -boxlenhalf_nm.x);
		movable_particle.y += DeviceConstants::boxSize.boxSizeNM_f.y * ((static_particle.y - movable_particle.y) > boxlenhalf_nm.y);
		movable_particle.y -= DeviceConstants::boxSize.boxSizeNM_f.y * ((static_particle.y - movable_particle.y) < -boxlenhalf_nm.y);
		movable_particle.z += DeviceConstants::boxSize.boxSizeNM_f.z * ((static_particle.z - movable_particle.z) > boxlenhalf_nm.z);
		movable_particle.z -= DeviceConstants::boxSize.boxSizeNM_f.z * ((static_particle.z - movable_particle.z) < -boxlenhalf_nm.z);
	}

	__device__ constexpr static void applyBCNM(Float3& current_position) {	// Only changes position if position is outside of box;		
		current_position.x += DeviceConstants::boxSize.boxSizeNM_f.x * (current_position.x < 0.f);
		current_position.x -= DeviceConstants::boxSize.boxSizeNM_f.x * (current_position.x > DeviceConstants::boxSize.boxSizeNM_f.x);
		current_position.y += DeviceConstants::boxSize.boxSizeNM_f.y * (current_position.y < 0.f);
		current_position.y -= DeviceConstants::boxSize.boxSizeNM_f.y * (current_position.y > DeviceConstants::boxSize.boxSizeNM_f.y);
		current_position.z += DeviceConstants::boxSize.boxSizeNM_f.z * (current_position.z < 0.f);
		current_position.z -= DeviceConstants::boxSize.boxSizeNM_f.z * (current_position.z > DeviceConstants::boxSize.boxSizeNM_f.z);
	}
};