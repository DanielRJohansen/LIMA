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
#include "BoxGrid.cuh"



namespace LIMAPOSITIONSYSTEM {
	// -------------------------------------------------------- PBC and HyperPos -------------------------------------------------------- //




	// -------------------------------------------------------- LimaPosition Conversion -------------------------------------------------------- //


	
	template<typename T>
	inline T floor_div(T num, T denom) {
		static_assert(std::is_integral<T>::value, "floor_div requires integer types");
		return (num - (num < 0 ? denom - 1 : 0)) / denom;
	}


	/// <summary>
	/// Use with care, will overflow if posLM is > 20 nm. Does NOT apply boundary condition
	/// </summary>
	/// <param name="posLM"></param>
	/// <returns></returns>
	__device__ inline NodeIndex PositionToNodeIndex(const Float3& posLM) {
		NodeIndex nodeindex{
			static_cast<int>(round((posLM.x) / static_cast<float>(BoxGrid::blocksizeLM))),
			static_cast<int>(round((posLM.y) / static_cast<float>(BoxGrid::blocksizeLM))),
			static_cast<int>(round((posLM.z) / static_cast<float>(BoxGrid::blocksizeLM)))
		};

		return nodeindex;
	}

	__device__ __host__ inline NodeIndex PositionToNodeIndexNM(const Float3& posNM) {
		NodeIndex nodeindex{
			static_cast<int>(round((posNM.x) / static_cast<float>(BoxGrid::blocksizeNM))),
			static_cast<int>(round((posNM.y) / static_cast<float>(BoxGrid::blocksizeNM))),
			static_cast<int>(round((posNM.z) / static_cast<float>(BoxGrid::blocksizeNM)))
		};

		return nodeindex;
	}

	/// <summary>
	/// Converts a nodeindex to a relative position in [lm]. ONLY safe to call with relatively small node indexes. 
	/// If the index has proponents larger than what a coord can represent, then :((( 
	/// </summary>
	/// <returns>Coord in [lm]</returns>
	__device__ __host__ static Coord nodeIndexToCoord(const NodeIndex& node_index) { 
		return Coord{ node_index.x, node_index.y, node_index.z } * static_cast<int32_t>(BoxGrid::blocksizeLM);	// TODO: Unsafe
	}
	

	// Returns absolute position of nodeindex [nm]
	__device__ __host__ static Float3 nodeIndexToAbsolutePosition(const NodeIndex& node_index) {
		const float nodelen_nm = static_cast<float>(BoxGrid::blocksizeNM);
		return Float3{ 
			static_cast<float>(node_index.x) * nodelen_nm,
			static_cast<float>(node_index.y) * nodelen_nm,
			static_cast<float>(node_index.z) * nodelen_nm
		};
	}

	static Coord getRelativeCoord(const Float3& absPosNM, const NodeIndex& nodeindex, const int max_node_diff, float boxlen_nm, BoundaryConditionSelect bc) {
		// Subtract nodeindex from abs position to get relative position
		Float3 hyperPos = absPosNM;
		const Float3 nodePos = nodeIndexToAbsolutePosition(nodeindex);
		BoundaryConditionPublic::applyHyperposNM(nodePos, hyperPos, boxlen_nm, bc);

		const Float3 relpos = hyperPos - nodePos;

		if (relpos.LargestMagnitudeElement() > static_cast<float>(max_node_diff)*BoxGrid::blocksizeNM) {
			/*auto absPos = absolute_position.toFloat3();
			auto hPos = hyperPos.toFloat3();
			auto nPos = nodePos.toFloat3();*/
			throw std::runtime_error("Tried to place a position that was not correcly assigned a node.");
			// Pos: " + absPos.toString() + " hyperpos : " + hPos.toString() + " nodePos : " + nPos.toString());
			//+ "% f % f % f Hyperpos % f % f % f node % f % f % f");
		}

		return Coord{ 
			static_cast<int32_t>(relpos.x * NANO_TO_LIMA), 
			static_cast<int32_t>(relpos.y * NANO_TO_LIMA), 
			static_cast<int32_t>(relpos.z * NANO_TO_LIMA) 
		};
	}

	// relpos in LM
	__device__ __host__ static Float3 relposToAbsolutePosition(const Coord& relpos) {
		return relpos.toFloat3() / NANO_TO_LIMA;
	}

	__host__ static std::tuple<NodeIndex, Coord> absolutePositionPlacement(const Float3& position, float boxlen_nm, BoundaryConditionSelect bc) {
		NodeIndex nodeindex = PositionToNodeIndexNM(position);	// TEMP
		BoundaryConditionPublic::applyBC(nodeindex, boxlen_nm, bc);

		const Coord relpos = getRelativeCoord(position, nodeindex, 1, boxlen_nm, bc);
		return std::make_tuple(nodeindex, relpos);
	}

	__device__ __host__ static Float3 GetAbsolutePositionNM(const NodeIndex& nodeindex, const Coord& coord) {
		return nodeIndexToAbsolutePosition(nodeindex) + relposToAbsolutePosition(coord);
	}
	__device__ __host__ static Float3 GetAbsolutePositionNM(const NodeIndex& nodeindex, const Float3& relposNM) {
		return nodeIndexToAbsolutePosition(nodeindex) + relposNM;
	}


	/// <summary>
	/// Transfer external coordinates to internal multi-range LIMA coordinates
	/// </summary>
	/// <param name="state">Absolute positions of particles as float [nm]</param>
	/// <param name="key_particle_index">Index of centermost particle of compound</param>
	static CompoundCoords positionCompound(const std::vector<Float3>& positions,  int key_particle_index, float boxlen_nm, BoundaryConditionSelect bc) {
		CompoundCoords compoundcoords{};

		compoundcoords.origo = PositionToNodeIndexNM(positions[key_particle_index]);
		BoundaryConditionPublic::applyBC(compoundcoords.origo, boxlen_nm, bc);

		for (int i = 0; i < positions.size(); i++) {
			// Allow some leeway, as different particles in compound may fit different gridnodes
			compoundcoords.rel_positions[i] = getRelativeCoord(positions[i], compoundcoords.origo, 3, boxlen_nm, bc);

		}
		if (compoundcoords.rel_positions[key_particle_index].maxElement() > static_cast<int32_t>(BoxGrid::blocksizeLM)) {
			compoundcoords.rel_positions[key_particle_index].print('k');
			//printf("%d\n", compoundcoords.rel_positions[key_particle_index].maxElement());
			throw std::runtime_error("Failed to place compound correctly.");
		}
		return compoundcoords;
	}



	//__device__ __host__ static bool canRepresentRelativeDist(const Coord& origo_a, const Coord& origo_b) {
	//	const auto diff = origo_a - origo_b;
	//	return std::abs(diff.x) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.y) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.z) < MAX_REPRESENTABLE_DIFF_NM;
	//}

	//// Get hyper index of "other"
	//template <typename BoundaryCondition>
	//__device__ __host__ static NodeIndex getHyperNodeIndex(const NodeIndex& self, const NodeIndex& other) {
	//	NodeIndex temp = other;
	//	BoundaryCondition::applyHyperpos(self, temp);
	//	//applyHyperpos<BoundaryCondition>(self, temp);
	//	return temp;
	//}


	// The following two functions MUST ALWAYS be used together
	// Shift refers to the wanted difference in the relative positions, thus origo must move -shift.
	// Keep in mind that origo is in nm and rel_pos are in lm
	// ONLY CALL FROM THREAD 0
	__device__ static Coord shiftOrigo(CompoundCoords& coords, const int keyparticle_index) {
		//const Coord shift_nm = coords.rel_positions[keyparticle_index] / static_cast<int32_t>(NANO_TO_LIMA);	// OPTIM. If LIMA wasn't 100 femto, but rather a power of 2, we could do this much better!
		//const NodeIndex shift_node = coordToNodeIndex(coords.rel_positions[keyparticle_index]);

		const NodeIndex shift = NodeIndex{
			coords.rel_positions[keyparticle_index].x / (static_cast<int>(BoxGrid::blocksizeLM) / 2),	// /2 so we switch to new index once we are halfway there
			coords.rel_positions[keyparticle_index].y / (static_cast<int>(BoxGrid::blocksizeLM) / 2),
			coords.rel_positions[keyparticle_index].z / (static_cast<int>(BoxGrid::blocksizeLM) / 2)
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
		const int32_t blocklen_half = BoxGrid::blocksizeLM / 2;

		EngineUtilsWarnings::verifyValidRelpos(relpos);

		return getOnehotDirection(relpos, blocklen_half);
	}

	/// <summary>NListManager
	/// Shifts the position 1/2 blocklen so we can find the appropriate origo.
	/// Applies PBC to the solvent
	/// </summary>
	/// <param name="position">Absolute position of solvent [nm] </param>
	__host__ static SolventCoord createSolventcoordFromAbsolutePosition(Float3 position, float boxlen_nm, BoundaryConditionSelect bc) {	// Find a way to do this without the the BC
		BoundaryConditionPublic::applyBCNM(position, boxlen_nm, bc);


		NodeIndex nodeindex; Coord relpos;
		std::tie(nodeindex, relpos) = absolutePositionPlacement(position, boxlen_nm, bc);

		SolventCoord solventcoord{ nodeindex, relpos };
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
	__device__ __host__ static float calcHyperDistNM(const Float3& p1, const Float3& p2) {
		Float3 temp = p2;	
		BoundaryCondition::applyHyperposNM(p1, temp);
		return (p1 - temp).len();
	}

	__host__ static Float3 GetPosition(const CompoundcoordsCircularQueue_Host& coords, int64_t step, int compoundIndex, int particleIndex) {
		return GetAbsolutePositionNM(coords.getCoordArray(step, compoundIndex).origo, coords.getCoordArray(step, compoundIndex).rel_positions[particleIndex]);
	}


	__device__ inline void LoadCompoundPositionAsLm(const CompoundCoords* const coordsGlobal, Float3& origoOut, Float3* relposOut, int nActiveParticles) {
		if (threadIdx.x == 0) {
			origoOut = nodeIndexToAbsolutePosition(coordsGlobal->origo);
		}
		relposOut[threadIdx.x] = coordsGlobal->rel_positions[threadIdx.x].toFloat3();
	}

	__device__ inline Float3 LoadRelposLmAndOrigo(const CompoundCoords* const coordsGlobal, Float3& origoOut) {
		if (threadIdx.x == 0) {
			origoOut = nodeIndexToAbsolutePosition(coordsGlobal->origo);
		}
		return coordsGlobal->rel_positions[threadIdx.x].toFloat3();
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
	__device__ static Coord getRelativeShiftBetweenCoordarrays(CompoundCoords* coordarray_circular_queue, int64_t step, int compound_index_left, int compound_index_right) {
		NodeIndex& nodeindex_left = CompoundcoordsCircularQueueUtils::getCoordarrayRef(coordarray_circular_queue, step, compound_index_left)->origo;
		NodeIndex& nodeindex_right = CompoundcoordsCircularQueueUtils::getCoordarrayRef(coordarray_circular_queue, step, compound_index_right)->origo;

		const NodeIndex hypernodeindex_right = BoundaryCondition::applyHyperpos_Return(nodeindex_left, nodeindex_right);
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
