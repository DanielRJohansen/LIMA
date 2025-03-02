/// This file provides BC for non-engine entities. These functions should perfectly mimic their compiler-optimized equals in BoundaryCondition.cuh, 
/// but be independent of any compile-time constants for boxlen or similar.
/// Be careful of making any changes in this file, without making the simular changes to the similar function in BoundaryCondition.cuh
#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "BoxGrid.cuh"

namespace BoundaryConditionPublic {
	void applyBC(NodeIndex& nodeindex, int boxlenNM, BoundaryConditionSelect);
	void applyBC(NodeIndex& nodeindex, int nNodesPerDim);	// TODO temp, this is not the best way forward	
	void applyBCNM(Float3& pos_nm, float boxlenNM, BoundaryConditionSelect);

	void applyHyperpos(const NodeIndex& staticNodeindex, NodeIndex& movableNodeindex, int boxlen_nm, BoundaryConditionSelect);
	void applyHyperposNM(const Float3& static_particle, Float3& movable_particle, float boxlen_nm, BoundaryConditionSelect);
}