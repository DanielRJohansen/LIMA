#pragma once

#include<iostream>
//#include <cmath>

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"

#include "EngineUtilsWarnings.cuh"
#include "BoundaryConditionPublic.h"

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"




namespace LIMAPOSITIONSYSTEM {
	// -------------------------------------------------------- PBC and HyperPos -------------------------------------------------------- //




	// -------------------------------------------------------- LimaPosition Conversion -------------------------------------------------------- //


	
	template<typename T>
	T floor_div(T num, T denom) {
		static_assert(std::is_integral<T>::value, "floor_div requires integer types");
		return (num - (num < 0 ? denom - 1 : 0)) / denom;
	}

	// Converts to nodeindex, applies PBC
	__host__ static NodeIndex absolutePositionToNodeIndex(const PositionHighRes& position, BoundaryConditionSelect bc, float boxlen_nm) {
		const int64_t offset = BOXGRID_NODE_LEN_i / 2;
		/*NodeIndex nodeindex{
			static_cast<int>((position.x + offset) / BOXGRID_NODE_LEN_i),
			static_cast<int>((position.y + offset) / BOXGRID_NODE_LEN_i),
			static_cast<int>((position.z + offset) / BOXGRID_NODE_LEN_i)
		};*/
		NodeIndex nodeindex{
			(int)floor_div(position.x + offset, static_cast<int64_t>(BOXGRID_NODE_LEN_i)),
			(int)floor_div(position.y + offset, static_cast<int64_t>(BOXGRID_NODE_LEN_i)),
			(int)floor_div(position.z + offset, static_cast<int64_t>(BOXGRID_NODE_LEN_i))
		};
		BoundaryConditionPublic::applyBC(nodeindex, boxlen_nm, bc);

		return nodeindex;
	}


	/// <summary>
	/// Converts a nodeindex to a relative position in [lm]. ONLY safe to call with relatively small node indexes. 
	/// If the index has proponents larger than what a coord can represent, then :((( 
	/// </summary>
	/// <returns>Coord in [lm]</returns>
	__device__ __host__ static Coord nodeIndexToCoord(const NodeIndex& node_index) { 
		return Coord{ node_index.x, node_index.y, node_index.z } * BOXGRID_NODE_LEN_i; 
	}
	

	// Returns absolute position of nodeindex [nm]
	__device__ __host__ static Float3 nodeIndexToAbsolutePosition(const NodeIndex& node_index) {
		const float nodelen_nm = BOXGRID_NODE_LEN / NANO_TO_LIMA;
		return Float3{ 
			static_cast<float>(node_index.x) * nodelen_nm,
			static_cast<float>(node_index.y) * nodelen_nm,
			static_cast<float>(node_index.z) * nodelen_nm
		};
	}

	//template <typename BoundaryCondition>
	static Coord getRelativeCoord(const PositionHighRes& absolute_position, const NodeIndex& nodeindex, const int max_node_diff, float boxlen_nm, BoundaryConditionSelect bc) {
		// Subtract nodeindex from abs position to get relative position
		PositionHighRes hyperpos = absolute_position;
		//applyHyperpos<PeriodicBoundaryCondition>(createLimaPosition(nodeindex), hyperpos);
		BoundaryConditionPublic::applyHyperpos(PositionHighRes{ nodeindex,  BOXGRID_NODE_LEN_i}, hyperpos, boxlen_nm, bc);

		const PositionHighRes relpos = hyperpos - PositionHighRes{ nodeindex, BOXGRID_NODE_LEN_i };

		if (relpos.largestMagnitudeElement() > BOXGRID_NODE_LEN_i * static_cast<int64_t>(max_node_diff)) {
			throw std::runtime_error("Tried to place a position that was not correcly assigned a node");
		}

		return Coord{ static_cast<int32_t>(relpos.x), static_cast<int32_t>(relpos.y), static_cast<int32_t>(relpos.z) };
	}

	// relpos in LM
	__device__ __host__ static Float3 relposToAbsolutePosition(const Coord& relpos) {
		return relpos.toFloat3() / NANO_TO_LIMA;
	}


	__host__ static std::tuple<NodeIndex, Coord> absolutePositionPlacement(const PositionHighRes& position, float boxlen_nm, BoundaryConditionSelect bc) {
		const NodeIndex nodeindex = absolutePositionToNodeIndex(position, bc, boxlen_nm);	// TEMP
		const Coord relpos = getRelativeCoord(position, nodeindex, 1, boxlen_nm, bc);
		return std::make_tuple(nodeindex, relpos);
	}

	__device__ __host__ static Float3 getAbsolutePositionNM(const NodeIndex& nodeindex, const Coord& coord) {
		return nodeIndexToAbsolutePosition(nodeindex) + relposToAbsolutePosition(coord);
	}

	/// <summary>
	/// Transfer external coordinates to internal multi-range LIMA coordinates
	/// </summary>
	/// <param name="state">Absolute positions of particles as float [nm]</param>
	/// <param name="key_particle_index">Index of centermost particle of compound</param>
	static CompoundCoords positionCompound(const std::vector<PositionHighRes>& positions,  int key_particle_index, float boxlen_nm, BoundaryConditionSelect bc) {
		CompoundCoords compoundcoords{};

		// WARNING: It may become a problem that state and state_prev does not share an origo. That should be fixed..
		//compoundcoords.origo = absolutePositionToNodeIndex<BoundaryCondition>(positions[key_particle_index]);
		compoundcoords.origo = absolutePositionToNodeIndex(positions[key_particle_index], bc, boxlen_nm);	//TEMP

		for (int i = 0; i < positions.size(); i++) {
			// Allow some leeway, as different particles in compound may fit different gridnodes
			compoundcoords.rel_positions[i] = getRelativeCoord(positions[i], compoundcoords.origo, 3, boxlen_nm, bc);	
		}
		if (compoundcoords.rel_positions[key_particle_index].maxElement() > BOXGRID_NODE_LEN_i) {
			compoundcoords.rel_positions[key_particle_index].print('k');
			//printf("%d\n", compoundcoords.rel_positions[key_particle_index].maxElement());
			throw std::runtime_error("Failed to place compound correctly.");
		}
		return compoundcoords;
	}



	__device__ __host__ static bool canRepresentRelativeDist(const Coord& origo_a, const Coord& origo_b) {
		const auto diff = origo_a - origo_b;
		return std::abs(diff.x) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.y) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.z) < MAX_REPRESENTABLE_DIFF_NM;
	}

	// Get hyper index of "other"
	template <typename BoundaryCondition>
	__device__ __host__ static NodeIndex getHyperNodeIndex(const NodeIndex& self, const NodeIndex& other) {
		NodeIndex temp = other;
		BoundaryCondition::applyHyperpos(self, temp);
		//applyHyperpos<BoundaryCondition>(self, temp);
		return temp;
	}


	// The following two functions MUST ALWAYS be used together
	// Shift refers to the wanted difference in the relative positions, thus origo must move -shift.
	// Keep in mind that origo is in nm and rel_pos are in lm
	// ONLY CALL FROM THREAD 0
	__device__ static Coord shiftOrigo(CompoundCoords& coords, const int keyparticle_index) {
		//const Coord shift_nm = coords.rel_positions[keyparticle_index] / static_cast<int32_t>(NANO_TO_LIMA);	// OPTIM. If LIMA wasn't 100 femto, but rather a power of 2, we could do this much better!
		//const NodeIndex shift_node = coordToNodeIndex(coords.rel_positions[keyparticle_index]);

		const NodeIndex shift = NodeIndex{
			coords.rel_positions[keyparticle_index].x / (BOXGRID_NODE_LEN_i/2),	// /2 so we switch to new index once we are halfway there
			coords.rel_positions[keyparticle_index].y / (BOXGRID_NODE_LEN_i/2),
			coords.rel_positions[keyparticle_index].z / (BOXGRID_NODE_LEN_i/2)
		};
		EngineUtilsWarnings::verifyCompoundOrigoshiftDuringIntegrationIsValid(shift, coords.rel_positions[keyparticle_index]);
		coords.origo += shift;
		return -nodeIndexToCoord(shift);
	}

	__device__ static NodeIndex getOnehotDirection(const Coord relpos, const int32_t threshold) {
		const int32_t magnitude_x = std::abs(relpos.x);
		const int32_t magnitude_y = std::abs(relpos.y);
		const int32_t magnitude_z = std::abs(relpos.z);

		const bool out_of_bounds = 
			(relpos.x < -threshold || relpos.x >= threshold) ||
			(relpos.y < -threshold || relpos.y >= threshold) ||
			(relpos.z < -threshold || relpos.z >= threshold);

		const bool is_x_largest = (magnitude_x >= magnitude_y) && (magnitude_x >= magnitude_z);
		const bool is_y_largest = (magnitude_y >= magnitude_x) && (magnitude_y >= magnitude_z);
		const bool is_z_largest = (magnitude_z >= magnitude_x) && (magnitude_z >= magnitude_y);

		const int x_component = out_of_bounds * is_x_largest * (relpos.x < 0 ? -1 : 1);
		const int y_component = out_of_bounds * is_y_largest * (relpos.y < 0 ? -1 : 1);
		const int z_component = out_of_bounds * is_z_largest * (relpos.z < 0 ? -1 : 1);

		return NodeIndex{ x_component, y_component, z_component };
	}

	// Since coord is rel to 0,0,0 of a block, we need to offset the positions so they are scattered around the origo instead of above it
	// We also need a threshold of half a blocklen, otherwise we should not transfer, and return{0,0,0}
	__device__ static NodeIndex getTransferDirection(const Coord& relpos) {
		const int32_t blocklen_half = BOXGRID_NODE_LEN_i / 2;

		EngineUtilsWarnings::verifyValidRelpos(relpos);

		return getOnehotDirection(relpos, blocklen_half);
	}

	/// <summary>NListManager
	/// Shifts the position 1/2 blocklen so we can find the appropriate origo.
	/// Applies PBC to the solvent
	/// </summary>
	/// <param name="position">Absolute position of solvent [nm] </param>
	__host__ static SolventCoord createSolventcoordFromAbsolutePosition(const PositionHighRes& position, float boxlen_nm, BoundaryConditionSelect bc) {	// Find a way to do this without the the BC
		PositionHighRes hyperpos = position;
		//applyBC<BoundaryCondition>(hyperpos);
		BoundaryConditionPublic::applyBC(hyperpos, boxlen_nm, bc);


		//const auto [nodeindex, relpos] = LIMAPOSITIONSYSTEM::absolutePositionPlacement(position);
		NodeIndex nodeindex; Coord relpos;
		//std::tie(nodeindex, relpos) = LIMAPOSITIONSYSTEM::absolutePositionPlacement<BoundaryCondition>(hyperpos);
		std::tie(nodeindex, relpos) = absolutePositionPlacement(hyperpos, boxlen_nm, bc);

		SolventCoord solventcoord{ nodeindex, relpos };
		//applyBC<BoundaryCondition>(solventcoord);
		BoundaryConditionPublic::applyBC(solventcoord.origo, boxlen_nm, bc);
		return solventcoord;
	}

	template <typename BoundaryCondition>
	__host__ static float calcHyperDist(const NodeIndex& left, const NodeIndex& right) {
		const NodeIndex right_hyper = getHyperNodeIndex<BoundaryCondition>(left, right);
		const NodeIndex diff = right_hyper - left;
		return nodeIndexToAbsolutePosition(diff).len();
	}

	template <typename BoundaryCondition>
	__device__ __host__ static float calcHyperDistNM(const Float3* const p1, const Float3* const p2) {
		Float3 temp = *p2;	
		BoundaryCondition::applyHyperposNM(p1, &temp);
		return (*p1 - temp).len();
	}

	__host__ static Float3 GetPosition(const CompoundcoordsCircularQueue& coords, int step, int compoundIndex, int particleIndex) {
		return getAbsolutePositionNM(coords.getCoordArray(step, compoundIndex).origo, coords.getCoordArray(step, compoundIndex).rel_positions[particleIndex]);
	}


};


// gcc is being a bitch with threadIdx and blockIdx in .cuh files that are also included by .c++ files.
// This workaround is to have these functions as static class fucntinos instead of namespace, which avoid the issue somehow. fuck its annoying tho
class LIMAPOSITIONSYSTEM_HACK{
public:
	__device__ static void getRelativePositions(const Coord* coords, Float3* positions, const unsigned int n_elements) {
		if (threadIdx.x < n_elements)
			positions[threadIdx.x] = coords[threadIdx.x].toFloat3();
	}

	__device__ static void shiftRelPos(CompoundCoords& coords, const Coord& shift_lm) {
		coords.rel_positions[threadIdx.x] += shift_lm;
	}

	template <typename BoundaryCondition>
	__device__ static void applyBC(CompoundCoords& coords) {
		if (threadIdx.x != 0) { return; }
		BoundaryCondition::applyBC(coords.origo);
	}


	// This function is only used in bridge, and can be made alot smarter with that context. TODO
	// Calculate the shift in [lm] for all relpos belonging to right, so they will share origo with left
	template <typename BoundaryCondition>
	__device__ static Coord getRelativeShiftBetweenCoordarrays(CompoundcoordsCircularQueue* coordarray_circular_queue, int step, int compound_index_left, int compound_index_right) {
		NodeIndex& nodeindex_left = coordarray_circular_queue->getCoordarrayRef(step, compound_index_left)->origo;
		NodeIndex& nodeindex_right = coordarray_circular_queue->getCoordarrayRef(step, compound_index_right)->origo;

		const NodeIndex hypernodeindex_right = LIMAPOSITIONSYSTEM::getHyperNodeIndex<BoundaryCondition>(nodeindex_left, nodeindex_right);
		const NodeIndex nodeshift_right_to_left = nodeindex_left - hypernodeindex_right;

		EngineUtilsWarnings::verifyNodeIndexShiftIsSafe(nodeshift_right_to_left);

		// Calculate necessary shift in relative position for all particles of right, so they share origo with left
		//return (coord_origo_left - hyperorigo_right) * static_cast<uint32_t>(NANO_TO_LIMA);	// This fucks up when the diff is > ~20
		return LIMAPOSITIONSYSTEM::nodeIndexToCoord(nodeshift_right_to_left * -1);
	}

		// Calculate the necessary shift in LM of all elements of FROM, assuming the origo has been shifted to TO
	/*__device__ static Coord getRelShiftFromOrigoShift(const Coord& origo_from, const Coord& origo_to) {
		return (origo_from - origo_to) * static_cast<int32_t>(NANO_TO_LIMA);
	}*/
	__device__ static Coord getRelShiftFromOrigoShift(const NodeIndex& from, const NodeIndex& to) {
		EngineUtilsWarnings::verifyOrigoShiftIsValid(from, to);

		const NodeIndex origo_shift = from - to;
		return LIMAPOSITIONSYSTEM::nodeIndexToCoord(origo_shift);
	}
};
