#pragma once

#include "KernelConstants.cuh"

class NoBoundaryCondition {
public:
	__device__ __host__ void static applyBC(NodeIndex& origo) {}

	__device__ constexpr static NodeIndex applyBC(const NodeIndex& nodeindex, const Int3& nodesPerDim) { return nodeindex; }
	static void constexpr applyBCNM(const Float3&) {}
	__device__ __host__ static void applyHyperpos(const NodeIndex& static_index, NodeIndex& movable_index) {}

	__device__ __host__ static inline void applyHyperposNM(const Float3& static_particle, Float3& movable_particle) {}
	
	__device__ static inline void ApplyBC(Float3& currentPosition, const Float3& boxSize, const Float3& boxSizeInv) {}
	__device__ constexpr static inline void ApplyHyperpos(const Float3& staticParticle, Float3& movableParticle, const Float3& boxSize, const Float3& boxSizeInv) {}
	__device__ static inline void ApplyHyperpos(const Float3& staticParticle, float& x, float& y, float& z, const Float3& boxSize, const Float3& boxSizeInv) {}
	/*__device__ static inline Float3 GetHyperposTranslation(const Float3& staticParticle, const Float3& movableParticle, const Float3& boxSize, const Float3& boxSizeInv) {
		return Float3{ 0, 0, 0 };
	}*/
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

	__device__ constexpr static inline void applyHyperposNM(const Float3& static_particle, Float3& movable_particle) {
		const Float3 boxlenhalf_nm = DeviceConstants::boxSize.boxSizeNM_f * 0.5f;

		movable_particle.x += DeviceConstants::boxSize.boxSizeNM_f.x * ((static_particle.x - movable_particle.x) > boxlenhalf_nm.x);
		movable_particle.x -= DeviceConstants::boxSize.boxSizeNM_f.x * ((static_particle.x - movable_particle.x) < -boxlenhalf_nm.x);
		movable_particle.y += DeviceConstants::boxSize.boxSizeNM_f.y * ((static_particle.y - movable_particle.y) > boxlenhalf_nm.y);
		movable_particle.y -= DeviceConstants::boxSize.boxSizeNM_f.y * ((static_particle.y - movable_particle.y) < -boxlenhalf_nm.y);
		movable_particle.z += DeviceConstants::boxSize.boxSizeNM_f.z * ((static_particle.z - movable_particle.z) > boxlenhalf_nm.z);
		movable_particle.z -= DeviceConstants::boxSize.boxSizeNM_f.z * ((static_particle.z - movable_particle.z) < -boxlenhalf_nm.z);
	}

	__device__  static inline void ApplyHyperpos(const Float3& staticParticle, Float3& movableParticle, const Float3& boxSize, const Float3& boxSizeInv) {
		const Float3 delta = staticParticle - movableParticle;

		movableParticle.x += boxSize.x * float(__float2int_rn(delta.x * boxSizeInv.x));
		movableParticle.y += boxSize.y * float(__float2int_rn(delta.y * boxSizeInv.y));
		movableParticle.z += boxSize.z * float(__float2int_rn(delta.z * boxSizeInv.z));
	}

	__device__ static inline void ApplyHyperpos(const Float3& staticParticle, float& x, float& y, float& z, const Float3& boxSize, const Float3& boxSizeInv) {
		const Float3 delta = staticParticle - Float3(x, y, z);

		x += boxSize.x * float(__float2int_rn(delta.x * boxSizeInv.x));
		y += boxSize.y * float(__float2int_rn(delta.y * boxSizeInv.y));
		z += boxSize.z * float(__float2int_rn(delta.z * boxSizeInv.z));
	}

	//__device__ static inline Float3 GetHyperposTranslation(const Float3& staticParticle, const Float3& movableParticle, const Float3& boxSize, const Float3& boxSizeInv) {
	//	const Float3 delta = staticParticle - movableParticle;
	//	return Float3(
	//		boxSize.x * float(__float2int_rn(delta.x * boxSizeInv.x)),
	//		boxSize.y * float(__float2int_rn(delta.y * boxSizeInv.y)),
	//		boxSize.z * float(__float2int_rn(delta.z * boxSizeInv.z))
	//	);
	//}

	__device__ static inline void ApplyBC(Float3& currentPosition, const Float3& boxSize, const Float3& boxSizeInv) {	
		currentPosition.x -= boxSize.x * floorf(currentPosition.x * boxSizeInv.x);
		currentPosition.y -= boxSize.y * floorf(currentPosition.y * boxSizeInv.y);
		currentPosition.z -= boxSize.z * floorf(currentPosition.z * boxSizeInv.z);
	}

// TODO: CHECK ALL THESE! Most are wrong, its CRITICAL we do >= not just >!!!
	__device__ constexpr static void applyBCNM(Float3& current_position) {	// Only changes position if position is outside of box;		
		current_position.x += DeviceConstants::boxSize.boxSizeNM_f.x * (current_position.x < 0.f);
		current_position.x -= DeviceConstants::boxSize.boxSizeNM_f.x * (current_position.x >= DeviceConstants::boxSize.boxSizeNM_f.x);
		current_position.y += DeviceConstants::boxSize.boxSizeNM_f.y * (current_position.y < 0.f);
		current_position.y -= DeviceConstants::boxSize.boxSizeNM_f.y * (current_position.y >= DeviceConstants::boxSize.boxSizeNM_f.y);
		current_position.z += DeviceConstants::boxSize.boxSizeNM_f.z * (current_position.z < 0.f);
		current_position.z -= DeviceConstants::boxSize.boxSizeNM_f.z * (current_position.z >= DeviceConstants::boxSize.boxSizeNM_f.z);
	}
};