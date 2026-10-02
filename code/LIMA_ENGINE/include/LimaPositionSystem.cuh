#pragma once

#include<iostream>
//#include <cmath>

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"

#include "BoundaryConditionPublic.h"

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "BoxGrid.cuh"



namespace LIMAPOSITIONSYSTEM {

	// -------------------------------------------------------- LimaPosition Conversion -------------------------------------------------------- //


	__device__ __host__ inline NodeIndex PositionToNodeIndexNM(const Float3& posNM) {
		NodeIndex nodeindex{
			static_cast<int>(round((posNM.x) / static_cast<float>(BoxGrid::blocksizeNM))),
			static_cast<int>(round((posNM.y) / static_cast<float>(BoxGrid::blocksizeNM))),
			static_cast<int>(round((posNM.z) / static_cast<float>(BoxGrid::blocksizeNM)))
		};

		return nodeindex;
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

	__device__ __host__ static Float3 GetAbsolutePositionNM(const NodeIndex& nodeindex, const Float3& relposNM) {
		return nodeIndexToAbsolutePosition(nodeindex) + relposNM;
	}


	/// <summary>
	/// Transfer external coordinates to internal multi-range LIMA coordinates
	/// </summary>
	/// <param name="state">Absolute positions of particles as float [nm]</param>
	/// <param name="key_particle_index">Index of centermost particle of compound</param>
	//static CompoundCoords positionCompound(const std::vector<Float3>& positions,  int key_particle_index, Int3 boxlen_nm, BoundaryConditionSelect bc) {
	//	CompoundCoords compoundcoords{};

	//	compoundcoords.origo = PositionToNodeIndexNM(positions[key_particle_index]);
	//	BoundaryConditionPublic::applyBC(compoundcoords.origo, boxlen_nm, bc);

	//	for (int i = 0; i < positions.size(); i++) {
	//		// Allow some leeway, as different particles in compound may fit different gridnodes
	//		compoundcoords.rel_positions[i] = getRelativeCoord(positions[i], compoundcoords.origo, 3, Float3::FromInt3(boxlen_nm), bc);

	//	}
	//	return compoundcoords;
	//}



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




	// Returns a one-hot vector of the largest magnitude axis, IF the abs of that axis is above threshold
	__device__ static NodeIndex GetTransferDirection(const Float3& pos, float threshold = 0.5f) { // optim consider making the threshold a template param
		const float ax = fabsf(pos.x);
		const float ay = fabsf(pos.y);
		const float az = fabsf(pos.z);

		// masks for which axis wins (ties resolved deterministically)
		const int mx = (ax >= ay) & (ax >= az) & (ax >= threshold);
		const int my = (ay > ax) & (ay >= az) & (ay >= threshold);
		const int mz = (az > ax) & (az > ay) & (az >= threshold);

		// sign without branching
		const int sx = (pos.x > 0.f) - (pos.x < 0.f);
		const int sy = (pos.y > 0.f) - (pos.y < 0.f);
		const int sz = (pos.z > 0.f) - (pos.z < 0.f);

		return {
			mx * sx,
			my * sy,
			mz * sz
		};
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

    template <typename BoundaryCondition>
    __device__ __host__ static float calcHyperDistSquaredNM(const Float3& p1, const Float3& p2) {
        Float3 temp = p2;
        BoundaryCondition::applyHyperposNM(p1, temp);
        return (p1 - temp).lenSquared();
    }

	//__host__ static Float3 GetPosition(const CompoundcoordsCircularQueue_Host& coords, int64_t step, int compoundIndex, int particleIndex) {
	//	return GetAbsolutePositionNM(coords.getCoordArray(step, compoundIndex).origo, coords.getCoordArray(step, compoundIndex).rel_positions[particleIndex]);
	//}



};


// gcc is being a bitch with threadIdx and blockIdx in .cuh files that are also included by .c++ files.
// This workaround is to have these functions as static class fucntinos instead of namespace, which avoid the issue somehow. fuck its annoying tho
class LIMAPOSITIONSYSTEM_HACK{
public:


	/*__device__ static Coord GetRelShiftFromOrigoShift_Coord(const NodeIndex& from, const NodeIndex& to) {
		EngineUtilsWarnings::verifyOrigoShiftIsValid(from, to);

		const NodeIndex origo_shift = from - to;
		return Coord{ origo_shift };
	}*/
    __device__ static constexpr Float3 GetRelShiftFromOrigoShift_Float3(const NodeIndex& from, const NodeIndex& to) {
        return NodeIndex{from-to}.toFloat3();
    }

};
