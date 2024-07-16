#pragma once

#include "KernelConstants.cuh"

class NoBoundaryCondition {
public:
	__device__ __host__ void static applyBC(NodeIndex& origo) {}


	__device__ __host__ static void applyHyperpos(const NodeIndex& static_index, NodeIndex& movable_index) {}

	__device__ __host__ static inline void applyHyperposNM(const Float3& static_particle, Float3& movable_particle) {}

	__device__ __host__ static NodeIndex applyHyperpos_Return(const NodeIndex& static_index, const NodeIndex& movable_index) { return movable_index; }
};

class PeriodicBoundaryCondition {
public:
	__device__ __host__ void static applyBC(NodeIndex& origo) {
		origo.x += boxSize_device.blocksPerDim * (origo.x < 0);
		origo.x -= boxSize_device.blocksPerDim * (origo.x >= boxSize_device.blocksPerDim);
		origo.y += boxSize_device.blocksPerDim * (origo.y < 0);
		origo.y -= boxSize_device.blocksPerDim * (origo.y >= boxSize_device.blocksPerDim);
		origo.z += boxSize_device.blocksPerDim * (origo.z < 0);
		origo.z -= boxSize_device.blocksPerDim * (origo.z >= boxSize_device.blocksPerDim);
	}

	__device__ __host__ static void applyHyperpos(const NodeIndex& static_index, NodeIndex& movable_index) {
		const NodeIndex difference = static_index - movable_index;		
		movable_index.x += boxSize_device.blocksPerDim * (difference.x > (boxSize_device.blocksPerDim / 2));		// Dont need to +1 to account of uneven, this is correct (im pretty sure)
		movable_index.x -= boxSize_device.blocksPerDim * (difference.x < -(boxSize_device.blocksPerDim / 2));
		movable_index.y += boxSize_device.blocksPerDim * (difference.y > (boxSize_device.blocksPerDim / 2));
		movable_index.y -= boxSize_device.blocksPerDim * (difference.y < -(boxSize_device.blocksPerDim / 2));
		movable_index.z += boxSize_device.blocksPerDim * (difference.z > (boxSize_device.blocksPerDim / 2));
		movable_index.z -= boxSize_device.blocksPerDim * (difference.z < -(boxSize_device.blocksPerDim / 2));
	}

	__device__ __host__ static NodeIndex applyHyperpos_Return(const NodeIndex& static_index, const NodeIndex& movable_index) {
		NodeIndex hyperIndex = movable_index;

		const NodeIndex difference = static_index - movable_index;
		hyperIndex.x += boxSize_device.blocksPerDim * (difference.x > (boxSize_device.blocksPerDim / 2));		// Dont need to +1 to account of uneven, this is correct (im pretty sure)
		hyperIndex.x -= boxSize_device.blocksPerDim * (difference.x < -(boxSize_device.blocksPerDim / 2));
		hyperIndex.y += boxSize_device.blocksPerDim * (difference.y > (boxSize_device.blocksPerDim / 2));
		hyperIndex.y -= boxSize_device.blocksPerDim * (difference.y < -(boxSize_device.blocksPerDim / 2));
		hyperIndex.z += boxSize_device.blocksPerDim * (difference.z > (boxSize_device.blocksPerDim / 2));
		hyperIndex.z -= boxSize_device.blocksPerDim * (difference.z < -(boxSize_device.blocksPerDim / 2));

		return hyperIndex;
	}

	__device__ __host__ static inline void applyHyperposNM(const Float3& static_particle, Float3& movable_particle) {		
		const float boxlenhalf_nm = boxSize_device.boxSizeNM_f / 2.f;

		for (int i = 0; i < 3; i++) {
			movable_particle[i] += boxSize_device.boxSizeNM_f * ((static_particle[i] - movable_particle[i]) > boxlenhalf_nm);
			movable_particle[i] -= boxSize_device.boxSizeNM_f * ((static_particle[i] - movable_particle[i]) < -boxlenhalf_nm);
		}
	}

	__device__ __host__ static void applyBCNM(Float3& current_position) {	// Only changes position if position is outside of box;		
		current_position.x += boxSize_device.boxSizeNM_f * (current_position.x < 0.f);
		current_position.x -= boxSize_device.boxSizeNM_f * (current_position.x > boxSize_device.boxSizeNM_f);
		current_position.y += boxSize_device.boxSizeNM_f * (current_position.y < 0.f);
		current_position.y -= boxSize_device.boxSizeNM_f * (current_position.y > boxSize_device.boxSizeNM_f);
		current_position.z += boxSize_device.boxSizeNM_f * (current_position.z < 0.f);
		current_position.z -= boxSize_device.boxSizeNM_f * (current_position.z > boxSize_device.boxSizeNM_f);
	}
};

// This macro functions calls the correct kernel, depending on the BoundaryCondition select
#define LAUNCH_GENERIC_KERNEL(kernel, grid_size, block_size, bc_select, kernel_arg) \
	do { \
		switch (bc_select) { \
			case NoBC: \
			kernel<NoBoundaryCondition> <<<grid_size, block_size>>> (kernel_arg); \
			break; \
	\
			case PBC: \
			kernel<PeriodicBoundaryCondition> <<<grid_size, block_size>>> (kernel_arg); \
			break; \
	\
			default: \
			throw std::runtime_error("Unsupported boundary condition in LAUNCH_GENERIC_KERNEL"); \
		} \
	} \
	while(0);	// Dont question the do-while, its safer for some reason i barely understand

// This macro functions calls the correct kernel, depending on the BoundaryCondition select and the doEM flag
#define LAUNCH_GENERIC_KERNEL_2(kernel, grid_size, block_size, bc_select, doEM, kernel_arg) \
	do { \
		switch (bc_select) { \
			case NoBC: \
				if (doEM) { \
					kernel<NoBoundaryCondition, true> <<<grid_size, block_size>>> (kernel_arg); \
				} else { \
					kernel<NoBoundaryCondition, false> <<<grid_size, block_size>>> (kernel_arg); \
				} \
				break; \
			\
			case PBC: \
				if (doEM) { \
					kernel<PeriodicBoundaryCondition, true> <<<grid_size, block_size>>> (kernel_arg); \
				} else { \
					kernel<PeriodicBoundaryCondition, false> <<<grid_size, block_size>>> (kernel_arg); \
				} \
				break; \
			\
			default: \
				throw std::runtime_error("Unsupported boundary condition in LAUNCH_GENERIC_KERNEL"); \
		} \
	} \
	while(0);	// Dont question the do-while, its safer for some reason i barely understand

// Macro that calls the correct function, depending on the boundarycondition select
#define CALL_FUNCTION_WITH_BC(func, bc_select, ...) \
    do { \
        switch (bc_select) { \
            case NoBC: \
                func<NoBoundaryCondition>(__VA_ARGS__); \
                break; \
            case PBC: \
                func<PeriodicBoundaryCondition>(__VA_ARGS__); \
                break; \
            default: \
                throw std::runtime_error("Unsupported boundary condition in CALL_FUNCTION_WITH_BC"); \
        } \
    } while (0)

//// This is just pure ChatGPT magic now...
//#define CALL_FUNCTION_WITH_BC_RET(func, bc_select, ...) \
//    ([&]() -> decltype(auto) { \
//        switch (bc_select) { \
//            case NoBC: \
//                return func<NoBoundaryCondition>(__VA_ARGS__); \
//            case PBC: \
//                return func<PeriodicBoundaryCondition>(__VA_ARGS__); \
//            default: \
//                throw std::runtime_error("Unsupported boundary condition in CALL_FUNCTION_WITH_BC_RET"); \
//        } \
//    }())