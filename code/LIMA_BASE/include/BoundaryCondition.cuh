#pragma once


class NoBoundaryCondition {
public:
	__host__ void static applyBC(NodeIndex& origo) {}

	__host__ static void applyHyperpos(const LimaPosition& static_position, LimaPosition& movable_position) {}

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

	__device__ __host__ static void applyPC(NodeIndex& origo) {
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

	__host__ static void applyHyperpos(const LimaPosition& static_position, LimaPosition& movable_position) {
		const LimaPosition difference = static_position - movable_position;
		movable_position.x += BOX_LEN_i * (difference.x > BOX_LEN_i / 2);
		movable_position.x -= BOX_LEN_i * (difference.x < -BOX_LEN_i / 2);
		movable_position.y += BOX_LEN_i * (difference.y > BOX_LEN_i / 2);
		movable_position.y -= BOX_LEN_i * (difference.y < -BOX_LEN_i / 2);
		movable_position.z += BOX_LEN_i * (difference.z > BOX_LEN_i / 2);
		movable_position.z -= BOX_LEN_i * (difference.z < -BOX_LEN_i / 2);
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