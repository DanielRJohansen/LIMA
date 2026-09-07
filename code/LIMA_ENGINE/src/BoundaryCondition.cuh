#pragma once

class NoBoundaryCondition {
public:
	__device__ __host__ void static applyBC(NodeIndex& origo, const Int3& gridDim) {}

	__device__ constexpr static NodeIndex applyBC(const NodeIndex& nodeindex, const Int3& nodesPerDim) { return nodeindex; }
	static void constexpr applyBCNM(Float3& currentPosition, const Float3& boxSize) {}
	__device__ __host__ static void applyHyperpos(const NodeIndex& staticIndex, NodeIndex& movableIndex, const Int3& gridDim) {}

	__device__ __host__ static inline void applyHyperposNM(const Float3& staticParticle, Float3& movableParticle, const Float3& boxSize) {}
	
	__device__ static inline void ApplyBC(Float3& currentPosition, const Float3& boxSize, const Float3& boxSizeInv) {}
	__device__ constexpr static inline void ApplyHyperpos(const Float3& staticParticle, Float3& movableParticle, const Float3& boxSize, const Float3& boxSizeInv) {}
	__device__ static inline void ApplyHyperpos(const Float3& staticParticle, float& x, float& y, float& z, const Float3& boxSize, const Float3& boxSizeInv) {}
	/*__device__ static inline Float3 GetHyperposTranslation(const Float3& staticParticle, const Float3& movableParticle, const Float3& boxSize, const Float3& boxSizeInv) {
		return Float3{ 0, 0, 0 };
	}*/
};

class PeriodicBoundaryCondition {
public:
	__device__ constexpr static void applyBC(NodeIndex& origo, const Int3& gridDim) {
		origo.x += gridDim.x * ((origo.x < 0) - (origo.x >= gridDim.x));
		origo.y += gridDim.y * ((origo.y < 0) - (origo.y >= gridDim.y));
		origo.z += gridDim.z * ((origo.z < 0) - (origo.z >= gridDim.z));
	}

	__device__ constexpr static NodeIndex applyBC(const NodeIndex& nodeindex, const Int3& gridDim) {
		NodeIndex output = nodeindex;
		output.x += gridDim.x * ((nodeindex.x < 0) - (nodeindex.x >= gridDim.x));
		output.y += gridDim.y * ((nodeindex.y < 0) - (nodeindex.y >= gridDim.y));
		output.z += gridDim.z * ((nodeindex.z < 0) - (nodeindex.z >= gridDim.z));
		return output;
	}

	__device__ constexpr static void applyHyperpos(const NodeIndex& staticIndex, NodeIndex& movableIndex, const Int3& gridDim) {
		const NodeIndex difference = staticIndex - movableIndex;
		const Int3 gridDimHalf = gridDim / 2;
		movableIndex.x += gridDim.x * (difference.x > gridDimHalf.x);		// Dont need to +1 to account of uneven, this is correct (im pretty sure)
		movableIndex.x -= gridDim.x * (difference.x < -gridDimHalf.x);
        movableIndex.y += gridDim.y * (difference.y > gridDimHalf.y);
		movableIndex.y -= gridDim.y * (difference.y < -gridDimHalf.y);
		movableIndex.z += gridDim.z * (difference.z > gridDimHalf.z);
		movableIndex.z -= gridDim.z * (difference.z < -gridDimHalf.z);
	}

	__device__ constexpr static inline void applyHyperposNM(const Float3& staticParticle, Float3& movableParticle, const Float3& boxSize) {
		const Float3 boxSizeHalf = boxSize * 0.5f;

		movableParticle.x += boxSize.x * ((staticParticle.x - movableParticle.x) > boxSizeHalf.x);
		movableParticle.x -= boxSize.x * ((staticParticle.x - movableParticle.x) < -boxSizeHalf.x);
		movableParticle.y += boxSize.y * ((staticParticle.y - movableParticle.y) > boxSizeHalf.y);
		movableParticle.y -= boxSize.y * ((staticParticle.y - movableParticle.y) < -boxSizeHalf.y);
		movableParticle.z += boxSize.z * ((staticParticle.z - movableParticle.z) > boxSizeHalf.z);
		movableParticle.z -= boxSize.z * ((staticParticle.z - movableParticle.z) < -boxSizeHalf.z);
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
	__device__ constexpr static void applyBCNM(Float3& currentPosition, const Float3& boxSize) {	// Only changes position if position is outside of box;
		currentPosition.x += boxSize.x * (currentPosition.x < 0.f);
		currentPosition.x -= boxSize.x * (currentPosition.x >= boxSize.x);
		currentPosition.y += boxSize.y * (currentPosition.y < 0.f);
		currentPosition.y -= boxSize.y * (currentPosition.y >= boxSize.y);
		currentPosition.z += boxSize.z * (currentPosition.z < 0.f);
		currentPosition.z -= boxSize.z * (currentPosition.z >= boxSize.z);
	}
};
