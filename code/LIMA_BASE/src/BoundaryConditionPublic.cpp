#include "BoundaryConditionPublic.h"

class NoBoundaryCondition {
public:
	static void constexpr applyBC(NodeIndex& origo) {}

	static void constexpr applyBC(Float3& position, float boxlen_nm) {}

};

class PeriodicBoundaryCondition {
public:
	static void constexpr applyBC(NodeIndex& origo, int boxgrid_n_nodes) {
		origo.x += boxgrid_n_nodes * (origo.x < 0);
		origo.x -= boxgrid_n_nodes * (origo.x >= boxgrid_n_nodes);
		origo.y += boxgrid_n_nodes * (origo.y < 0);
		origo.y -= boxgrid_n_nodes * (origo.y >= boxgrid_n_nodes);
		origo.z += boxgrid_n_nodes * (origo.z < 0);
		origo.z -= boxgrid_n_nodes * (origo.z >= boxgrid_n_nodes);
	}

	static void constexpr applyBC(Float3& position, float boxlen_nm) {
		position.x += boxlen_nm * (position.x < 0.f);
		position.x -= boxlen_nm * (position.x > boxlen_nm);
		position.y += boxlen_nm * (position.y < 0.f);
		position.y -= boxlen_nm * (position.y > boxlen_nm);
		position.z += boxlen_nm * (position.z < 0.f);
		position.z -= boxlen_nm * (position.z > boxlen_nm);
	}

	static void constexpr applyHyperposNM(const Float3& static_particle, Float3& movable_particle, float boxlen_nm) {
		const float boxlenhalf_nm = boxlen_nm / 2.f;

		movable_particle.x += boxlen_nm * ((static_particle.x - movable_particle.x) > boxlenhalf_nm);
		movable_particle.x -= boxlen_nm * ((static_particle.x - movable_particle.x) < -boxlenhalf_nm);
		movable_particle.y += boxlen_nm * ((static_particle.y - movable_particle.y) > boxlenhalf_nm);
		movable_particle.y -= boxlen_nm * ((static_particle.y - movable_particle.y) < -boxlenhalf_nm);
		movable_particle.z += boxlen_nm * ((static_particle.z - movable_particle.z) > boxlenhalf_nm);
		movable_particle.z -= boxlen_nm * ((static_particle.z - movable_particle.z) < -boxlenhalf_nm);
	}

	static void constexpr applyHyperpos(const NodeIndex& staticNodeindex, NodeIndex& movableNodeindex, int boxlen_nm) {
		const int boxlenHalfNM = boxlen_nm / 2;
		
		movableNodeindex.x += BoxGrid::NodesPerDim(boxlen_nm) * ((staticNodeindex.x - movableNodeindex.x) > boxlenHalfNM);
		movableNodeindex.x -= BoxGrid::NodesPerDim(boxlen_nm) * ((staticNodeindex.x - movableNodeindex.x) < -boxlenHalfNM);
		movableNodeindex.y += BoxGrid::NodesPerDim(boxlen_nm) * ((staticNodeindex.y - movableNodeindex.y) > boxlenHalfNM);
		movableNodeindex.y -= BoxGrid::NodesPerDim(boxlen_nm) * ((staticNodeindex.y - movableNodeindex.y) < -boxlenHalfNM);
		movableNodeindex.z += BoxGrid::NodesPerDim(boxlen_nm) * ((staticNodeindex.z - movableNodeindex.z) > boxlenHalfNM);
		movableNodeindex.z -= BoxGrid::NodesPerDim(boxlen_nm) * ((staticNodeindex.z - movableNodeindex.z) < -boxlenHalfNM);
		
	}
};



void BoundaryConditionPublic::applyBC(NodeIndex& nodeindex, int boxlenNM, BoundaryConditionSelect bc) {
	switch (bc) {
	case None: {
		NoBoundaryCondition::applyBC(nodeindex);
		break;
	}
	case PBC: {
		PeriodicBoundaryCondition::applyBC(nodeindex, BoxGrid::NodesPerDim(boxlenNM));
		break;
	}
	}
}
void BoundaryConditionPublic::applyBC(NodeIndex& nodeindex, int nNodesPerDim) {
	PeriodicBoundaryCondition::applyBC(nodeindex, nNodesPerDim);
}

void BoundaryConditionPublic::applyBCNM(Float3& pos_nm, float boxlen_nm, BoundaryConditionSelect bc) {
	switch (bc) {
	case None: {
		NoBoundaryCondition::applyBC(pos_nm, boxlen_nm);
		break;
	}
	case PBC: {
		PeriodicBoundaryCondition::applyBC(pos_nm, boxlen_nm);
		break;
	}
	}
}

void BoundaryConditionPublic::applyHyperposNM(const Float3& static_position, Float3& movable_position, float boxlen_nm, BoundaryConditionSelect bc) {
	switch (bc) {
	case None: {
		break;
	}
	case PBC: {
		PeriodicBoundaryCondition::applyHyperposNM(static_position, movable_position, boxlen_nm);
		break;
	}
	}
}

void BoundaryConditionPublic::applyHyperpos(const NodeIndex& staticNodeindex, NodeIndex& movableNodeindex, int boxlen_nm, BoundaryConditionSelect bc) {
	switch (bc) {
	case None: {
		break;
	}
	case PBC: {
		PeriodicBoundaryCondition::applyHyperpos(staticNodeindex, movableNodeindex, boxlen_nm);
		break;
	}
	}
}