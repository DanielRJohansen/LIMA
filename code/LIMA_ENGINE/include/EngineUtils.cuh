#pragma once

#include<iostream>
#include "LIMA_BASE/include/LimaTypes.cuh"
#include "LIMA_BASE/include/Constants.cuh"

#include "LIMA_BASE/include/Simulation.cuh"
#include "LIMA_BASE/include/Forcefield.cuh"
#include <cmath>

//#include <cuda.h>
//#include <device_launch_parameters.h>
//#include <cuda_runtime_api.h>

namespace ForceCalc {
//	-
};

//#include "cuda/std/cmath"
//#include "cuda/std//utility"

namespace CPPD {
	__device__ __host__ constexpr int32_t ceil(float num) {
		return (static_cast<float>(static_cast<int32_t>(num)) == num)
			? static_cast<int32_t>(num)
			: static_cast<int32_t>(num) + ((num > 0) ? 1 : 0);
	}

	template <typename T>
	__device__ __host__ static T max(const T l, const T r) {
		return r > l ? r : l;
	}

	template <typename T>
	__device__ __host__ static T min(const T l, const T r) {
		return r < l ? r : l;
	}

	__device__ __host__ static int32_t abs(const int32_t val) {
		return val < 0 ? -val : val;
	}
}




namespace LIMAPOSITIONSYSTEM {
	// -------------------------------------------------------- PBC and HyperPos -------------------------------------------------------- //

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

	__device__ __host__ static void applyPBC(NodeIndex& origo) {
		origo.x += BOXGRID_N_NODES * (origo.x < 0);
		origo.x -= BOXGRID_N_NODES * (origo.x >= BOXGRID_N_NODES);
		origo.y += BOXGRID_N_NODES * (origo.y < 0);
		origo.y -= BOXGRID_N_NODES * (origo.y >= BOXGRID_N_NODES);
		origo.z += BOXGRID_N_NODES * (origo.z < 0);
		origo.z -= BOXGRID_N_NODES * (origo.z >= BOXGRID_N_NODES);
	}

	__device__ __host__ static void applyPBC(SolventCoord& coord) { applyPBC(coord.origo); }

	__device__ static void applyPBC(CompoundCoords& coords) {
		if (threadIdx.x != 0) { return; }
		applyPBC(coords.origo);
	}

	__device__ __host__ static void applyPBC(LimaPosition& position) {
		// Offset position so we grab onto the correct node - NOT REALLY SURE ABOUT THIS...
		const int64_t offset = BOXGRID_NODE_LEN_i / 2; // + 1;
		position.x += BOX_LEN_i * (position.x + offset < 0);
		position.x -= BOX_LEN_i * (position.x + offset >= BOX_LEN_i);
		position.y += BOX_LEN_i * (position.y + offset < 0);
		position.y -= BOX_LEN_i * (position.y + offset >= BOX_LEN_i);
		position.z += BOX_LEN_i * (position.z + offset < 0);
		position.z -= BOX_LEN_i * (position.z + offset >= BOX_LEN_i);
	}

	// -------------------------------------------------------- LimaPosition Conversion -------------------------------------------------------- //

	__host__ static LimaPosition createLimaPosition(const NodeIndex& nodeindex) {
		return LimaPosition{
			nodeindex.x * BOXGRID_NODE_LEN_i,
			nodeindex.y * BOXGRID_NODE_LEN_i,
			nodeindex.z * BOXGRID_NODE_LEN_i
		};
	}
	__host__ static LimaPosition createLimaPosition(const Float3& pos_nm) {
		const Float3 pos_lm = pos_nm * NANO_TO_LIMA;
		return LimaPosition{ static_cast<int64_t>(pos_lm.x), static_cast<int64_t>(pos_lm.y), static_cast<int64_t>(pos_lm.z) };
	}
	
	//// Safe to call with any Coord
	//__device__ __host__ static NodeIndex coordToNodeIndex(const Coord& coord) { 
	//	return NodeIndex{ 
	//		coord.x / BOXGRID_NODE_LEN_i, 
	//		coord.y / BOXGRID_NODE_LEN_i , 
	//		coord.z / BOXGRID_NODE_LEN_i 
	//	}; 
	//}

	// Converts to nodeindex, applies PBC
	__host__ static NodeIndex absolutePositionToNodeIndex(const LimaPosition& position) {
		int offset = BOXGRID_NODE_LEN_i / 2;
		NodeIndex nodeindex{
			static_cast<int>((position.x + offset) / BOXGRID_NODE_LEN_i),
			static_cast<int>((position.y + offset) / BOXGRID_NODE_LEN_i),
			static_cast<int>((position.z + offset) / BOXGRID_NODE_LEN_i)
		};
		applyPBC(nodeindex);
		return nodeindex;
	}

	// Used only for neighborlist. Might be temp. Input-position in nm!
	__host__ static NodeIndex absolutePositionToNodeIndex(const Float3& position) {
		float factor = NANO_TO_LIMA / BOXGRID_NODE_LEN;
		NodeIndex nodeindex{
			static_cast<int32_t>(std::roundf(position.x * factor)),
			static_cast<int32_t>(std::roundf(position.y * factor)),
			static_cast<int32_t>(std::roundf(position.z * factor))
		};
		applyPBC(nodeindex);
		return nodeindex;
	}

	/// <summary>
	/// Converts a nodeindex to a relative position in [lm]. ONLY safe to call with relatively small node indexes. 
	/// If the index has proponents larger than what a coord can represent, then :((( 
	/// </summary>
	/// <returns>Coord in [lm]</returns>
	__device__ __host__ static Coord nodeIndexToCoord(const NodeIndex& node_index) { return Coord{ node_index.x, node_index.y, node_index.z } * BOXGRID_NODE_LEN_i; }
	

	// Returns absolute position of nodeindex [nm]
	__device__ __host__ static Float3 nodeIndexToAbsolutePosition(const NodeIndex& node_index) {
		const float nodelen_nm = BOXGRID_NODE_LEN / NANO_TO_LIMA;
		return Float3{ 
			static_cast<float>(node_index.x) * nodelen_nm,
			static_cast<float>(node_index.y) * nodelen_nm,
			static_cast<float>(node_index.z) * nodelen_nm
		};
	}

	__host__ static Coord getRelativeCoord(const LimaPosition& absolute_position, const NodeIndex& nodeindex, const int max_node_diff=1) {
		//float epsilon = 1.1;
		//if (position.largestMagnitudeElement() / epsilon > BOXGRID_NODE_LEN / NANO_TO_LIMA) { // Check if position is somewhat correctly placed
		//	throw "Tried to place a position that was not correcly assigned a node"; 
		//}

		// Subtract nodeindex from abs position to get relative position
		LimaPosition hyperpos = absolute_position;
		LIMAPOSITIONSYSTEM::applyHyperpos(createLimaPosition(nodeindex), hyperpos);
		const LimaPosition relpos = hyperpos - createLimaPosition(nodeindex);
		auto absposf = absolute_position.toFloat3();
		auto abs_hyperposf = hyperpos.toFloat3();
		auto origof = createLimaPosition(nodeindex).toFloat3();
		auto relposf = relpos.toFloat3();
		auto p = createLimaPosition(nodeindex);

		if (relpos.largestMagnitudeElement() > BOXGRID_NODE_LEN_i * max_node_diff) {
			throw "Tried to place a position that was not correcly assigned a node";
		}

		return Coord{ static_cast<int32_t>(relpos.x), static_cast<int32_t>(relpos.y), static_cast<int32_t>(relpos.z) };
	}

	// relpos in LM
	__device__ __host__ static Float3 relposToAbsolutePosition(const Coord& relpos) {
		return relpos.toFloat3() / NANO_TO_LIMA;
	}

	__host__ static std::tuple<NodeIndex, Coord> absolutePositionPlacement(const LimaPosition& position) {
		const NodeIndex nodeindex = absolutePositionToNodeIndex(position);
		const Coord relpos = getRelativeCoord(position, nodeindex);
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
	/// <returns></returns>
	static CompoundCoords positionCompound(const std::vector<LimaPosition>& positions,  int key_particle_index=0) {
		CompoundCoords compoundcoords{};

		// WARNING: It may become a problem that state and state_prev does not share an origo. That should be fixed..
		compoundcoords.origo = absolutePositionToNodeIndex(positions[key_particle_index]);

		for (int i = 0; i < positions.size(); i++) {
			// Allow some leeway, as different particles in compound may fit different gridnodes
			compoundcoords.rel_positions[i] = getRelativeCoord(positions[i], compoundcoords.origo, 2);	
		}
		return compoundcoords;
	}

	// Transforms the relative coords to relative float positions.
	__device__ static void getRelativePositions(const Coord* coords, Float3* positions, const unsigned int n_elements) {
		if (threadIdx.x < n_elements)
			positions[threadIdx.x] = coords[threadIdx.x].toFloat3();
	}

	__device__ __host__ static bool canRepresentRelativeDist(const Coord& origo_a, const Coord& origo_b) {
		const auto diff = origo_a - origo_b;
		return std::abs(diff.x) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.y) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.z) < MAX_REPRESENTABLE_DIFF_NM;
	}




	
	// Calculate the necessary shift in LM of all elements of FROM, assuming the origo has been shifted to TO
	/*__device__ static Coord getRelShiftFromOrigoShift(const Coord& origo_from, const Coord& origo_to) {
		return (origo_from - origo_to) * static_cast<int32_t>(NANO_TO_LIMA);
	}*/
	__device__ static Coord getRelShiftFromOrigoShift(const NodeIndex& from, const NodeIndex& to) {
		return nodeIndexToCoord(from - to);
	}

	// Get hyper index of "other"
	__device__ __host__ static NodeIndex getHyperNodeIndex(const NodeIndex& self, const NodeIndex& other) {
		NodeIndex temp = other;
		applyHyperpos(self, temp);
		return temp;
	}


	// The following two functions MUST ALWAYS be used together
	// Shift refers to the wanted difference in the relative positions, thus origo must move -shift.
	// Keep in mind that origo is in nm and rel_pos are in lm
	// ONLY CALL FROM THREAD 0
	__device__ static Coord shiftOrigo(CompoundCoords& coords, const int keyparticle_index=0) {
		//const Coord shift_nm = coords.rel_positions[keyparticle_index] / static_cast<int32_t>(NANO_TO_LIMA);	// OPTIM. If LIMA wasn't 100 femto, but rather a power of 2, we could do this much better!
		//const NodeIndex shift_node = coordToNodeIndex(coords.rel_positions[keyparticle_index]);

		const NodeIndex shift_node = NodeIndex{
			coords.rel_positions[keyparticle_index].x / (BOXGRID_NODE_LEN_i/2),	// /2 so we switch to new index once we are halfway there
			coords.rel_positions[keyparticle_index].y / (BOXGRID_NODE_LEN_i/2),
			coords.rel_positions[keyparticle_index].z / (BOXGRID_NODE_LEN_i/2)
		};

		//coords.origo += shift_nm;
		coords.origo += shift_node;
		return -nodeIndexToCoord(shift_node);
		//return -shift_nm * static_cast<int32_t>(NANO_TO_LIMA);
	}
	__device__ static void shiftRelPos(CompoundCoords& coords, const Coord& shift_lm) {
		coords.rel_positions[threadIdx.x] += shift_lm;
	}

	// This function is only used in bridge, and can be made alot smarter with that context. TODO
	// Calculate the shift in [lm] for all relpos belonging to right, so they will share origo with left
	__device__ static Coord getRelativeShiftBetweenCoordarrays(CompoundCoords* coordarray_circular_queue, int step, int compound_index_left, int compound_index_right) {
		NodeIndex& nodeindex_left = CoordArrayQueueHelpers::getCoordarrayRef(coordarray_circular_queue, step, compound_index_left)->origo;
		NodeIndex& nodeindex_right = CoordArrayQueueHelpers::getCoordarrayRef(coordarray_circular_queue, step, compound_index_right)->origo;


		const NodeIndex hypernodeindex_right = LIMAPOSITIONSYSTEM::getHyperNodeIndex(nodeindex_left, nodeindex_right);
		const NodeIndex nodeshift_right_to_left = nodeindex_left - hypernodeindex_right;

#ifdef LIMASAFEMODE
		if (nodeshift_right_to_left.manhattanLen() > MAX_SAFE_SHIFT) {
			printf("Shifting compound further than what is safe!");
		}
#endif

		// Calculate necessary shift in relative position for all particles of right, so they share origo with left
		//return (coord_origo_left - hyperorigo_right) * static_cast<uint32_t>(NANO_TO_LIMA);	// This fucks up when the diff is > ~20
		return nodeIndexToCoord(nodeshift_right_to_left * -1);
	}






	// ReCenter origo of solvent, and the relative pos around said origo
	__device__ static void updateSolventcoord(SolventCoord& coord) {
		//const int shift_at = static_cast<int32_t>(NANO_TO_LIMA) / 2;
		//Coord shift_nm = coord.rel_position / static_cast<int32_t>(shift_at);	// OPTIM. If LIMA wasn't 100 femto, but rather a power of 2, we could do this much better! 
		//const NodeIndex node_shift = coord.rel_position / (BOXGRID_NODE_LEN_i / 2)
		//const NodeIndex node_shift = LIMAPOSITIONSYSTEM::coordToNodeIndex(coord.rel_position * 2);	// If the coord is more than halfway to the next node, we shift it.

		const NodeIndex node_shift{
			coord.rel_position.x / (BOXGRID_NODE_LEN_i / 2),	// If the coord is more than halfway to the next node, we shift it.
			coord.rel_position.y / (BOXGRID_NODE_LEN_i / 2),
			coord.rel_position.z / (BOXGRID_NODE_LEN_i / 2)
		};

		SolventCoord tmp = coord;
		coord.origo += node_shift;
		//coord.rel_position -= shift_nm * static_cast<int32_t>(NANO_TO_LIMA);
		coord.rel_position -= LIMAPOSITIONSYSTEM::nodeIndexToCoord(node_shift);
		
		if (blockIdx.x + threadIdx.x == 0 && node_shift.x != 0) {
			tmp.origo.print('o');
			tmp.rel_position.print('r');
			node_shift.print('s');
			coord.origo.print('O');
			coord.rel_position.print('R');
		}
	}

	// Get the relpos_prev, if the solvent was in the same solventblock last step
	__device__ static Coord getRelposPrev(SolventBlockGrid* solventblockgrid_circularqueue, const int solventblock_id, const int current_step) {
		const int step_prev = current_step == 0 ? STEPS_PER_SOLVENTBLOCKTRANSFER-1 : current_step - 1;	// Unnecessary if we use transfermodule's prevpos for first step!!!!!!!!! TODOTODO
		auto blockPtr = CoordArrayQueueHelpers::getSolventBlockPtr(solventblockgrid_circularqueue, step_prev, solventblock_id);
		return blockPtr->rel_pos[threadIdx.x];
	}








	//__device__ static Coord getOnehotDirectionEdgeblock(const Coord relpos, const Coord thresholdInward, const Coord thresholdOutward)

	__device__ static NodeIndex getOnehotDirection(const Coord relpos, const int32_t threshold) {
		const int32_t magnitude_x = std::abs(relpos.x);
		const int32_t magnitude_y = std::abs(relpos.y);
		const int32_t magnitude_z = std::abs(relpos.z);

		//if (magnitude_x <= threshold && magnitude_y <= threshold && magnitude_z <= threshold) { return Coord{ 0,0,0 }; }	// This sadly means that there are multiple correct positions for a particle :(

		// Including at bottomleftfront, excluding at toprightback
		if (   relpos.x < -threshold || relpos.x >= threshold
			|| relpos.y < -threshold || relpos.y >= threshold
			|| relpos.z < -threshold || relpos.z >= threshold)
		{
			// Determine which magnitude is the largest
			if (magnitude_x >= magnitude_y && magnitude_x >= magnitude_z) {
				// The x component has the largest magnitude
				return NodeIndex{ relpos.x < 0 ? -1 : 1, 0, 0 };
			}
			else if (magnitude_y >= magnitude_z) { // The y component has the largest magnitude			
				return NodeIndex{ 0, relpos.y < 0 ? -1 : 1, 0 };
			}
			else { // The z component has the largest magnitude		
				return NodeIndex{ 0, 0, relpos.z < 0 ? -1 : 1 };
			}
		}
		else {
			return NodeIndex{};
		}		
	}

	// Since coord is rel to 0,0,0 of a block, we need to offset the positions so they are scattered around the origo instead of above it
	// We also need a threshold of half a blocklen, otherwise we should not transfer, and return{0,0,0}
	__device__ static NodeIndex getTransferDirection(const Coord relpos) {
		const int32_t blocklen_half = BOXGRID_NODE_LEN_i / 2;
		const Coord rel_blockcenter{ blocklen_half };
		if (relpos.x < INT32_MIN + blocklen_half || relpos.y < INT32_MIN + blocklen_half || relpos.z < INT32_MIN + blocklen_half) {
			printf("\nWe have underflow!\n");
			relpos.print('R');
		}
		if (relpos.x > INT32_MAX - blocklen_half || relpos.y > INT32_MAX - blocklen_half || relpos.z > INT32_MAX - blocklen_half) {
			printf("\nWe have overflow!\n");
			relpos.print('R');
		}
		//return LIMAPOSITIONSYSTEM::getOnehotDirection(relpos + rel_blockcenter, blocklen_half);
		return LIMAPOSITIONSYSTEM::getOnehotDirection(relpos, blocklen_half);
	}

	/// <summary>
	/// Shifts the position 1/2 blocklen so we can find the appropriate origo.
	/// Applies PBC to the solvent
	/// </summary>
	/// <param name="position">Absolute position of solvent [nm] </param>
	__host__ static SolventCoord createSolventcoordFromAbsolutePosition(const LimaPosition& position) {
		LimaPosition hyperpos = position;
		applyPBC(hyperpos);

		//const auto [nodeindex, relpos] = LIMAPOSITIONSYSTEM::absolutePositionPlacement(position);
		NodeIndex nodeindex; Coord relpos;
		std::tie(nodeindex, relpos) = LIMAPOSITIONSYSTEM::absolutePositionPlacement(hyperpos);

		SolventCoord solventcoord{ nodeindex, relpos };
		applyPBC(solventcoord);
		return solventcoord;
	}

	__host__ static float calcHyperDist(const NodeIndex& left, const NodeIndex& right) {
		const NodeIndex right_hyper = getHyperNodeIndex(left, right);
		const NodeIndex diff = right_hyper - left;
		return nodeIndexToAbsolutePosition(diff).len();
	}
};








namespace EngineUtils {

	// -------------------------------------------------------- Functions in NM - doesnt fit well here! -------------------------------------------------------- //
	// Assumes positions in NM
	__device__ __host__ static inline void applyHyperposNM(const Float3* static_particle, Float3* movable_particle) {
		const float boxlen_nm = BOX_LEN / NANO_TO_LIMA;
		const float boxlenhalf_nm = boxlen_nm / 2.f;

		for (int i = 0; i < 3; i++) {
			*movable_particle->placeAt(i) += boxlen_nm * ((static_particle->at(i) - movable_particle->at(i)) > boxlenhalf_nm);
			*movable_particle->placeAt(i) -= boxlen_nm * ((static_particle->at(i) - movable_particle->at(i)) < -boxlenhalf_nm);
		}
	}

	// Assumes positions in NM
	__device__ __host__ static float calcHyperDistNM(const Float3* p1, const Float3* p2) {
		Float3 temp = *p2;
		applyHyperposNM(p1, &temp);
		return (*p1 - temp).len();
	}

	__device__ __host__ static float calcKineticEnergy(const float velocity, const float mass) {
		return 0.5f * mass * velocity * velocity;
	}

	// Positions in [nm] and time in [ns]
	__device__ __host__ static float calcKineticEnergy(const Float3* pos1, const Float3* pos2, const float mass, const float elapsed_time) {
		auto dist = calcHyperDistNM(pos1, pos2);
		if (dist > 1.f) {
			printf("DIST: %f\n", dist);
			pos1->print('1');
			pos2->print('2');
		}

		const float velocity = dist / elapsed_time;
		return calcKineticEnergy(velocity, mass);
		//const float kinE = 0.5f * mass * vel * vel;
		//return kinE;
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------- //









	// LimaPosition in [nm]
	__device__ __host__ static void applyPBCNM(Float3* current_position) {	// Only changes position if position is outside of box;
		for (int dim = 0; dim < 3; dim++) {
			*current_position->placeAt(dim) += BOX_LEN * (current_position->at(dim) < 0.f);
			*current_position->placeAt(dim) -= BOX_LEN * (current_position->at(dim) > BOX_LEN);
		}
	}



	//static float calcSpeedOfParticle(const float mass /*[kg]*/, const float temperature /*[K]*/) { // 
	//	const float R = 8.3144f;								// Ideal gas constants - J/(Kelvin*mol)
	//	const float v_rms = static_cast<float>(sqrt(3.f * R * temperature / mass));
	//	return v_rms;	// [m/s]
	//}
	
	//static float tempToVelocity(float temperature /*[K]*/, float mass /*[kg/mol]*/) {
	//	const float kNA = BOLTZMANNCONSTANT * AVOGADROSNUMBER;			// k * mol
	//	const float v_rms = std::sqrt(3 * kNA * temperature / mass);	// Calculate root-mean-square velocity
	//	return v_rms * std::sqrt(2);									// Convert to most probable velocity
	//}

	// Calculate mean speed of particles. [K], [kg/mol]
	//sqrt((8RT)/(piM))
	static float tempToVelocity(double temperature /*[K]*/, double mass /*[kg/mol]*/) {
		return sqrtf(static_cast<float>(3.0 * BOLTZMANNCONSTANT * temperature / (mass/AVOGADROSNUMBER)));
	}
	// http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/kintem.html

	static float kineticEnergyToTemperature(long double kineticEnergy /*[J]*/, int numParticles) {
		const double temperature = kineticEnergy * (2.0 / 3.0) / (BOLTZMANNCONSTANT * numParticles);
		return static_cast<float>(temperature);
	}

	// For solvents, compound_id = n_compounds and particle_id = solvent_index
	__device__ static uint32_t getLoggingIndexOfParticle(uint32_t step, uint32_t total_particles_upperbound, uint32_t compound_id, uint32_t particle_id_local) {
		const uint32_t steps_since_transfer = (step % STEPS_PER_LOGTRANSFER);
		const uint32_t step_offset = steps_since_transfer * total_particles_upperbound;
		const uint32_t compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
		return step_offset + compound_offset + particle_id_local;
	}

	// TODO: Marked for deletion
	__host__ static size_t getAlltimeIndexOfParticle(uint64_t step, uint32_t total_particles_upperbound, uint32_t compound_id, uint32_t particle_id_local) {
		const uint32_t step_offset = static_cast<uint32_t>(step) * total_particles_upperbound;
		const uint32_t compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
		return step_offset + compound_offset + particle_id_local;
	}

	__device__ int static getNewBlockId(const NodeIndex& transfer_direction, const NodeIndex& origo) {
		NodeIndex new_nodeindex = transfer_direction + origo;
		LIMAPOSITIONSYSTEM::applyPBC(new_nodeindex);
		return SolventBlockGrid::get1dIndex(new_nodeindex);
	}

	__device__ static SolventBlockTransfermodule* getTransfermoduleTargetPtr(SolventBlockTransfermodule* transfermodule_array, int blockId, const NodeIndex& transfer_direction) {
		NodeIndex new_nodeindex = SolventBlockGrid::get3dIndex(blockId) + transfer_direction;
		LIMAPOSITIONSYSTEM::applyPBC(new_nodeindex);
		auto index = SolventBlockGrid::get1dIndex(new_nodeindex);
		return &transfermodule_array[index];
	}



	// returns an int between -21 and 21
	__device__ static int genPseudoRandomNum(int& seed) {
		const unsigned int a = 1664525;
		const unsigned int c = 1013904223;
		seed = a * seed + c;
		return seed / 100000000;
	}

	// This function assumes that coord_tsub1 is already hyperpositioned to coord.
	//__device__ static Coord integratePosition(const Coord& coord, const Coord& coord_tsub1, const Float3* force, const float mass, const float dt, const float thermostat_scalar) {
	//	// All threads will follow the same path, so branch is no problem
	//	const Coord delta_pos = thermostat_scalar == 1.f
	//		? coord - coord_tsub1
	//		: Coord{ (coord - coord_tsub1).toFloat3() * thermostat_scalar};

	//	return coord + delta_pos + Coord{ *force * dt * dt / mass };
	//}

	// returns pos_tadd1
	__device__ static Coord integratePositionVVS(const Coord& pos, const Float3& vel, const Float3& force, const double mass, const double dt) {
		
		/*if (threadIdx.x == 0) {
			Float3 f = vel * dt + force * (0.5 / mass * dt * dt);
			Coord c = Coord{ vel * dt + force * (0.5 / mass * dt * dt) };
			f.print('f');
			c.print('c');
			printf("dt: %f\n", dt);
		}*/

		//const Coord pos_tadd1 = pos + Coord{ vel * dt + force * (0.5 / mass * dt * dt) };						// braking/stable version
		//const Coord pos_tadd1 = pos + Coord{ (vel * dt + force * (0.5 / mass * dt * dt)).round()};				// precise version

		const Coord pos_tadd1 = pos + Coord{ (vel * dt + force * (0.5 / mass * dt * dt)).round() };				// precise version

		return pos_tadd1;
	}
	__device__ static Float3 integrateVelocityVVS(const Float3& vel_tsub1, const Float3& force_tsub1, const Float3& force, const double dt, const double mass) {
		/*Float3 f = (force + force_tsub1) * dt * 0.5f / mass;
		Coord c{ (force + force_tsub1) * dt * 0.5f / mass };
		const Coord vel = vel_tsub1 + Coord{ (force + force_tsub1) * dt * 0.5f / mass };

		if (threadIdx.x == 0) {
			f.print('f');
			c.print('c');
			printf("dt: %f\n", dt);
		}
		*/
		const Float3 vel = vel_tsub1 + (force + force_tsub1) * (dt * 0.5 / mass);
		return vel;
	}

	__device__ inline void LogCompoundData(Compound& compound, Box* box, CompoundCoords& compound_coords, float* potE_sum, Float3& force, Float3& force_LJ_sol, uint32_t step, DatabuffersDevice* databuffers) {
		const uint32_t index = EngineUtils::getLoggingIndexOfParticle(step, box->boxparams.total_particles_upperbound, blockIdx.x, threadIdx.x);
		databuffers->traj_buffer[index] = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compound_coords.origo, compound_coords.rel_positions[threadIdx.x]); //LIMAPOSITIONSYSTEM::getGlobalPosition(compound_coords);
		databuffers->potE_buffer[index] = *potE_sum;
		databuffers->vel_buffer[index] = compound.vels_prev[threadIdx.x];
	}

	__device__ inline void LogSolventData(Box* box, const float& potE, const SolventBlock& solventblock, bool solvent_active, const Float3& force, const Float3 velocity, uint32_t step, DatabuffersDevice* databuffers) {
		if (solvent_active) {
			const uint32_t index = EngineUtils::getLoggingIndexOfParticle(step, box->boxparams.total_particles_upperbound, box->boxparams.n_compounds, solventblock.ids[threadIdx.x]);
			//box->traj_buffer[index] = SolventBlockHelpers::extractAbsolutePositionLM(solventblock);
			if (solventblock.ids[threadIdx.x] == 0) {
				//LIMAPOSITIONSYSTEM::getAbsolutePositionNM(solventblock.origo, solventblock.rel_pos[threadIdx.x]).print('0');
				//force.print('f');
			}
			databuffers->traj_buffer[index] = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(solventblock.origo, solventblock.rel_pos[threadIdx.x]);
			databuffers->potE_buffer[index] = potE;


			if (solventblock.ids[threadIdx.x] > 13000) printf("\nhiewr: %u\n", solventblock.ids[threadIdx.x]);

#ifdef USEDEBUGF3
			const auto debug_index = (box->step * box->total_particles_upperbound + box->n_compounds * MAX_COMPOUND_PARTICLES + solventblock.ids[threadIdx.x]) * DEBUGDATAF3_NVARS;
			//box->debugdataf3[debug_index] = Float3(solventblock.ids[threadIdx.x] * 10 + 1.f, solventblock.ids[threadIdx.x] * 10 + 2.f, solventblock.ids[threadIdx.x] * 10 + 3.f);
			box->debugdataf3[debug_index] = force;
			box->debugdataf3[debug_index + 1] = velocity;
			box->debugdataf3[debug_index + 2] = SolventBlockHelpers::extractAbsolutePositionLM(solventblock) / NANO_TO_LIMA;
#endif
		}
	}



	// Slow function, never for device
	__host__ static float calcDistance(const NodeIndex& o1, const Coord& relpos1, const NodeIndex& o2, const Coord& relpos2) {
		auto pos1 = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(o1, relpos1);
		auto pos2 = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(o2, relpos2);
		return (pos1 - pos2).len();
	}
};


namespace LIMAKERNELDEBUG {
	__device__ void static transferOut(STransferQueue* queue_global, const STransferQueue& queue_local, const NodeIndex& transferdir_queue, const int queue_index) {
		if (queue_global->rel_positions[threadIdx.x].x < -2 * static_cast<int32_t>(NANO_TO_LIMA) || queue_global->rel_positions[threadIdx.x].x > 2 * static_cast<int32_t>(NANO_TO_LIMA)
			|| queue_global->rel_positions[threadIdx.x].y < -2 * static_cast<int32_t>(NANO_TO_LIMA) || queue_global->rel_positions[threadIdx.x].y > 2 * static_cast<int32_t>(NANO_TO_LIMA)
			|| queue_global->rel_positions[threadIdx.x].z < -2 * static_cast<int32_t>(NANO_TO_LIMA) || queue_global->rel_positions[threadIdx.x].z > 2 * static_cast<int32_t>(NANO_TO_LIMA)
			) {
			printf("\n");
			transferdir_queue.print('t');
			queue_local.rel_positions[threadIdx.x].print('q');
			queue_global->rel_positions[threadIdx.x].print('Q');
		}

		if (threadIdx.x == 0) {
			if (queue_global->n_elements != 0) {
				printf("\nN elements was: %d in queue %d\n", queue_global->n_elements, queue_index);
				transferdir_queue.print('d');
			}

			queue_global->n_elements = queue_local.n_elements;
			if (queue_local.n_elements > 15) {
				printf("\nTransferring %d elements\n", queue_local.n_elements);
			}
		}
	}


	__device__ void static compoundIntegration(const Coord& relpos_prev, const Coord& relpos_next, const Float3& force, bool& critical_error_encountered) {
		const auto dif = (relpos_next - relpos_prev);
		const int32_t max_diff = BOXGRID_NODE_LEN_i / 20;
		if (std::abs(dif.x) > max_diff || std::abs(dif.y) > max_diff || std::abs(dif.z) > max_diff || force.len() > 3.f) {
			printf("\nParticle %d in compound %d is moving too fast\n", threadIdx.x, blockIdx.x);
			dif.printS('D');
			force.print('F');
			relpos_next.printS('R');
			critical_error_encountered = true;
		}
	}

	__device__ void static solventIntegration(const Coord& relpos_prev, const Coord& relpos_next, const Float3& force, bool& critical_error_encountered, int id) {
		const auto dif = (relpos_next - relpos_prev);
		const int32_t max_diff = BOXGRID_NODE_LEN_i / 20;
		if (std::abs(dif.x) > max_diff || std::abs(dif.y) > max_diff || std::abs(dif.z) > max_diff) {
			printf("\nSolvent %d moving too fast\n", id);
			dif.printS('D');
			force.print('F');
			critical_error_encountered = true;
		}
	}

	// 
	//__device__ bool static forceTooLarge(const Float3& force, const float threshold=0.5) { // 

	//}

};

namespace DEBUGUTILS {

	/// <summary>
	/// Puts the nearest solvent of each solvent in the out vector.
	/// </summary>
	void findAllNearestSolventSolvent(SolventBlockGrid* solventblockgrid, size_t n_solvents, std::vector<float>& out);
}



// LIMA algorithm Library
namespace LAL {
	__device__ constexpr int getBlellochTablesize(int n) {
		float nf = static_cast<float>(n);
		return CPPD::ceil(nf * log2f(nf) * 2.f);
	}
}