/// This file provides BC for non-engine entities. These functions should perfectly mimic their compiler-optimized equals in BoundaryCondition.cuh, 
/// but be independent of any compile-time constants for boxlen or similar.
/// Be careful of making any changes in this file, without making the simular changes to the similar function in BoundaryCondition.cuh
#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "BoxGrid.cuh"

namespace BoundaryConditionPublic {
	void applyBC(NodeIndex& nodeindex, const Int3& boxlenNM, BoundaryConditionSelect);
	void applyBC(NodeIndex& nodeindex, const Int3& nNodesPerDim);	// TODO temp, this is not the best way forward	
	void applyBCNM(Float3& pos_nm, const Float3& boxlenNM, BoundaryConditionSelect);

	void applyHyperpos(const NodeIndex& staticNodeindex, NodeIndex& movableNodeindex, const Int3& boxlen_nm, BoundaryConditionSelect);
	void applyHyperposNM(const Float3& static_particle, Float3& movable_particle, const Float3& boxlen_nm, BoundaryConditionSelect);
}