/// This file provides BC for non-engine entities. These functions should perfectly mimic their compiler-optimized equals in BoundaryCondition.cuh, 
/// but be independent of any compile-time constants for boxlen or similar.
/// Be careful of making any changes in this file, without making the simular changes to the similar function in BoundaryCondition.cuh
#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"


namespace BoundaryConditionPublic {
	void applyBC(NodeIndex& nodeindex, float boxlen_nm, BoundaryConditionSelect);
	void applyBC(LimaPosition& position, float boxlen_nm, BoundaryConditionSelect);
	void applyBCNM(Float3& pos_nm, float boxlen_nm, BoundaryConditionSelect);

	void applyHyperpos(const LimaPosition& static_position, LimaPosition& movable_position, float boxlen_nm, BoundaryConditionSelect bc);
	void applyHyperposNM(const Float3* static_particle, Float3* movable_particle, float boxlen_nm, BoundaryConditionSelect);
}