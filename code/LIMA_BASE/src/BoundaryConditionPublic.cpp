#include "BoundaryConditionPublic.h"

class NoBoundaryCondition {
public:
	static void constexpr applyBC(NodeIndex& origo) {}

	static void constexpr applyBC(Float3& position, const Float3& boxlen_nm) {}
	static void constexpr applyBCNM(const Float3&) {}

};

class PeriodicBoundaryCondition {
public:
	static void constexpr applyBC(NodeIndex& origo, Int3 gridDim) {
		origo.x += gridDim.x * (origo.x < 0);
		origo.x -= gridDim.x * (origo.x >= gridDim.x);
		origo.y += gridDim.y * (origo.y < 0);
		origo.y -= gridDim.y * (origo.y >= gridDim.y);
		origo.z += gridDim.z * (origo.z < 0);
		origo.z -= gridDim.z * (origo.z >= gridDim.z)	;
	}

	static void constexpr applyBC(Float3& position, const Float3& boxlen_nm) {
		position.x += boxlen_nm.x * (position.x < 0.f);
		position.x -= boxlen_nm.x * (position.x > boxlen_nm.x);
		position.y += boxlen_nm.y * (position.y < 0.f);
		position.y -= boxlen_nm.y * (position.y > boxlen_nm.y);
		position.z += boxlen_nm.z * (position.z < 0.f);
		position.z -= boxlen_nm.z * (position.z > boxlen_nm.z);
	}

	static void constexpr applyHyperposNM(const Float3& static_particle, Float3& movable_particle, const Float3& boxlen_nm) {
		const Float3 boxlenhalf_nm = boxlen_nm / 2.f;

		movable_particle.x += boxlen_nm.x * ((static_particle.x - movable_particle.x) > boxlenhalf_nm.x);
		movable_particle.x -= boxlen_nm.x * ((static_particle.x - movable_particle.x) < -boxlenhalf_nm.x);
		movable_particle.y += boxlen_nm.y * ((static_particle.y - movable_particle.y) > boxlenhalf_nm.y);
		movable_particle.y -= boxlen_nm.y * ((static_particle.y - movable_particle.y) < -boxlenhalf_nm.y);
		movable_particle.z += boxlen_nm.z * ((static_particle.z - movable_particle.z) > boxlenhalf_nm.z);
		movable_particle.z -= boxlen_nm.z * ((static_particle.z - movable_particle.z) < -boxlenhalf_nm.z);
	}

	static void constexpr applyHyperpos(const NodeIndex& staticNodeindex, NodeIndex& movableNodeindex, const Int3& boxlen_nm) {
		const Int3 boxlenHalfNM = boxlen_nm / 2;
		
		movableNodeindex.x += BoxGrid::NodesPerDim(boxlen_nm.x) * ((staticNodeindex.x - movableNodeindex.x) > boxlenHalfNM.x);
		movableNodeindex.x -= BoxGrid::NodesPerDim(boxlen_nm.x) * ((staticNodeindex.x - movableNodeindex.x) < -boxlenHalfNM.x);
		movableNodeindex.y += BoxGrid::NodesPerDim(boxlen_nm.y) * ((staticNodeindex.y - movableNodeindex.y) > boxlenHalfNM.y);
		movableNodeindex.y -= BoxGrid::NodesPerDim(boxlen_nm.y) * ((staticNodeindex.y - movableNodeindex.y) < -boxlenHalfNM.y);
		movableNodeindex.z += BoxGrid::NodesPerDim(boxlen_nm.z) * ((staticNodeindex.z - movableNodeindex.z) > boxlenHalfNM.z);
		movableNodeindex.z -= BoxGrid::NodesPerDim(boxlen_nm.z) * ((staticNodeindex.z - movableNodeindex.z) < -boxlenHalfNM.z);
		
	}
};



void BoundaryConditionPublic::applyBC(NodeIndex& nodeindex, const Int3& boxlenNM, BoundaryConditionSelect bc) {
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
void BoundaryConditionPublic::applyBC(NodeIndex& nodeindex, const Int3& nNodesPerDim) {
	PeriodicBoundaryCondition::applyBC(nodeindex, nNodesPerDim);
}

void BoundaryConditionPublic::applyBCNM(Float3& pos_nm, const Float3& boxlen_nm, BoundaryConditionSelect bc) {
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

void BoundaryConditionPublic::applyHyperposNM(const Float3& static_position, Float3& movable_position, const Float3& boxlen_nm, BoundaryConditionSelect bc) {
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

void BoundaryConditionPublic::applyHyperpos(const NodeIndex& staticNodeindex, NodeIndex& movableNodeindex, const Int3& boxlen_nm, BoundaryConditionSelect bc) {
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