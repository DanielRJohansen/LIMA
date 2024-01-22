#include "BoundaryConditionPublic.h"

class NoBoundaryCondition {
public:
	static void applyBC(NodeIndex& origo) {}

	static void applyHyperpos(const LimaPosition& static_position, LimaPosition& movable_position) {}

	static void applyBC(LimaPosition& position) {}

	static void applyBC(Float3& position, float boxlen_nm) {}

};

class PeriodicBoundaryCondition {
public:
	static void applyBC(NodeIndex& origo, int boxgrid_n_nodes) {
		origo.x += boxgrid_n_nodes * (origo.x < 0);
		origo.x -= boxgrid_n_nodes * (origo.x >= boxgrid_n_nodes);
		origo.y += boxgrid_n_nodes * (origo.y < 0);
		origo.y -= boxgrid_n_nodes * (origo.y >= boxgrid_n_nodes);
		origo.z += boxgrid_n_nodes * (origo.z < 0);
		origo.z -= boxgrid_n_nodes * (origo.z >= boxgrid_n_nodes);
	}

	static void applyHyperpos(const LimaPosition& static_position, LimaPosition& movable_position, float boxlen_nm) {
		const int boxlen_i = static_cast<int>(boxlen_nm * NANO_TO_LIMA);
		const LimaPosition difference = static_position - movable_position;
		movable_position.x += boxlen_i * (difference.x > boxlen_i / 2);
		movable_position.x -= boxlen_i * (difference.x < -boxlen_i / 2);
		movable_position.y += boxlen_i * (difference.y > boxlen_i / 2);
		movable_position.y -= boxlen_i * (difference.y < -boxlen_i / 2);
		movable_position.z += boxlen_i * (difference.z > boxlen_i / 2);
		movable_position.z -= boxlen_i * (difference.z < -boxlen_i / 2);
	}

	static void applyBC(LimaPosition& position, float boxlen_nm) {
		// Offset position so we grab onto the correct node - NOT REALLY SURE ABOUT THIS...
		const int boxlen_i = static_cast<int>(boxlen_nm * NANO_TO_LIMA);
		const int64_t offset = BOXGRID_NODE_LEN_i / 2; // + 1;
		position.x += boxlen_i * (position.x + offset < 0);
		position.x -= boxlen_i * (position.x + offset >= boxlen_i);
		position.y += boxlen_i * (position.y + offset < 0);
		position.y -= boxlen_i * (position.y + offset >= boxlen_i);
		position.z += boxlen_i * (position.z + offset < 0);
		position.z -= boxlen_i * (position.z + offset >= boxlen_i);
	}

	static void applyBC(Float3& position, float boxlen_nm) {
		position.x += boxlen_nm * (position.x < 0.f);
		position.x -= boxlen_nm * (position.x > boxlen_nm);
		position.y += boxlen_nm * (position.y < 0.f);
		position.y -= boxlen_nm * (position.y > boxlen_nm);
		position.z += boxlen_nm * (position.z < 0.f);
		position.z -= boxlen_nm * (position.z > boxlen_nm);
	}

	static void applyHyperposNM(const Float3* static_particle, Float3* movable_particle, float boxlen_nm) {
		const float boxlenhalf_nm = boxlen_nm / 2.f;

		for (int i = 0; i < 3; i++) {
			*movable_particle->placeAt(i) += boxlen_nm * ((static_particle->at(i) - movable_particle->at(i)) > boxlenhalf_nm);
			*movable_particle->placeAt(i) -= boxlen_nm * ((static_particle->at(i) - movable_particle->at(i)) < -boxlenhalf_nm);
		}
	}
};



void BoundaryConditionPublic::applyBC(NodeIndex& nodeindex, float boxlen_nm, BoundaryConditionSelect bc) {
	const int boxgrid_n_nodes = static_cast<int>(std::roundf(boxlen_nm * NANO_TO_LIMA / BOXGRID_NODE_LEN));
	switch (bc) {
	case NoBC: {
		NoBoundaryCondition::applyBC(nodeindex);
		break;
	}
	case PBC: {
		PeriodicBoundaryCondition::applyBC(nodeindex, boxgrid_n_nodes);
		break;
	}
	}
}
void BoundaryConditionPublic::applyBC(LimaPosition& position, float boxlen_nm, BoundaryConditionSelect bc) {
	const int boxgrid_n_nodes = static_cast<int>(std::roundf(boxlen_nm * NANO_TO_LIMA / BOXGRID_NODE_LEN));
	switch (bc) {
	case NoBC: {
		NoBoundaryCondition::applyBC(position);
		break;
	}
	case PBC: {
		PeriodicBoundaryCondition::applyBC(position, boxlen_nm);
		break;
	}
	}
}

void BoundaryConditionPublic::applyBCNM(Float3& pos_nm, float boxlen_nm, BoundaryConditionSelect bc) {
	switch (bc) {
	case NoBC: {
		NoBoundaryCondition::applyBC(pos_nm, boxlen_nm);
		break;
	}
	case PBC: {
		PeriodicBoundaryCondition::applyBC(pos_nm, boxlen_nm);
		break;
	}
	}
}

void BoundaryConditionPublic::applyHyperpos(const LimaPosition& static_position, LimaPosition& movable_position, float boxlen_nm, BoundaryConditionSelect bc) {
	switch (bc) {
	case NoBC: {
		NoBoundaryCondition::applyHyperpos(static_position, movable_position);
		break;
	}
	case PBC: {
		PeriodicBoundaryCondition::applyHyperpos(static_position, movable_position, boxlen_nm);
		break;
	}
	}
}

void BoundaryConditionPublic::applyHyperposNM(const Float3* static_position, Float3* movable_position, float boxlen_nm, BoundaryConditionSelect bc) {
	switch (bc) {
	case NoBC: {
		break;
	}
	case PBC: {
		PeriodicBoundaryCondition::applyHyperposNM(static_position, movable_position, boxlen_nm);
		break;
	}
	}
}