#pragma once

#include <UserConstants.h>

class NoBoundaryCondition {
public:
	__host__ void static applyBC(NodeIndex& origo) {}

	__device__ __host__ static void applyBC(LimaPosition& position) {}

	__device__ __host__ static void applyHyperpos(const NodeIndex& static_index, NodeIndex& movable_index) {}

	__device__ __host__ static inline void applyHyperposNM(const Float3* static_particle, Float3* movable_particle) {}
};

class PeriodicBoundaryCondition {
public:
	__device__ __host__ void static applyBC(NodeIndex& origo) {
		origo.x += BOXGRID_N_NODES * (origo.x < 0);
		origo.x -= BOXGRID_N_NODES * (origo.x >= BOXGRID_N_NODES);
		origo.y += BOXGRID_N_NODES * (origo.y < 0);
		origo.y -= BOXGRID_N_NODES * (origo.y >= BOXGRID_N_NODES);
		origo.z += BOXGRID_N_NODES * (origo.z < 0);
		origo.z -= BOXGRID_N_NODES * (origo.z >= BOXGRID_N_NODES);
	}

	__device__ __host__ static void applyBC(LimaPosition& position) {
		// Offset position so we grab onto the correct node - NOT REALLY SURE ABOUT THIS...
		const int64_t offset = BOXGRID_NODE_LEN_i / 2; // + 1;
		position.x += BOX_LEN_i * (position.x + offset < 0);
		position.x -= BOX_LEN_i * (position.x + offset >= BOX_LEN_i);
		position.y += BOX_LEN_i * (position.y + offset < 0);
		position.y -= BOX_LEN_i * (position.y + offset >= BOX_LEN_i);
		position.z += BOX_LEN_i * (position.z + offset < 0);
		position.z -= BOX_LEN_i * (position.z + offset >= BOX_LEN_i);
	}

	__device__ __host__ static void applyHyperpos(const NodeIndex& static_index, NodeIndex& movable_index) {
		const NodeIndex difference = static_index - movable_index;
		movable_index.x += BOXGRID_N_NODES * (difference.x > (BOXGRID_N_NODES / 2));		// Dont need to +1 to account of uneven, this is correct (im pretty sure)
		movable_index.x -= BOXGRID_N_NODES * (difference.x < -(BOXGRID_N_NODES / 2));
		movable_index.y += BOXGRID_N_NODES * (difference.y > (BOXGRID_N_NODES / 2));
		movable_index.y -= BOXGRID_N_NODES * (difference.y < -(BOXGRID_N_NODES / 2));
		movable_index.z += BOXGRID_N_NODES * (difference.z > (BOXGRID_N_NODES / 2));
		movable_index.z -= BOXGRID_N_NODES * (difference.z < -(BOXGRID_N_NODES / 2));
	}

	__device__ __host__ static inline void applyHyperposNM(const Float3* static_particle, Float3* movable_particle) {
		const float boxlen_nm = BOX_LEN / NANO_TO_LIMA;
		const float boxlenhalf_nm = boxlen_nm / 2.f;

		for (int i = 0; i < 3; i++) {
			*movable_particle->placeAt(i) += boxlen_nm * ((static_particle->at(i) - movable_particle->at(i)) > boxlenhalf_nm);
			*movable_particle->placeAt(i) -= boxlen_nm * ((static_particle->at(i) - movable_particle->at(i)) < -boxlenhalf_nm);
		}
	}

	__device__ __host__ static void applyBCNM(Float3& current_position) {	// Only changes position if position is outside of box;
		current_position.x += BOX_LEN_NM * (current_position.x < 0.f);
		current_position.x -= BOX_LEN_NM * (current_position.x > BOX_LEN_NM);
		current_position.y += BOX_LEN_NM * (current_position.y < 0.f);
		current_position.y -= BOX_LEN_NM * (current_position.y > BOX_LEN_NM);
		current_position.z += BOX_LEN_NM * (current_position.z < 0.f);
		current_position.z -= BOX_LEN_NM * (current_position.z > BOX_LEN_NM);
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