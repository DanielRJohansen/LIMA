#pragma once

#include<iostream>
#include "LimaTypes.cuh"
#include "Constants.cuh"

#include "Simulation.cuh"
#include "Forcefield.cuh"
#include <cmath>

#include <cuda.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>

namespace ForceCalc {
//	-
};

#include "cuda/std/cmath"
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

	__device__ __host__ static int32_t abs(const int32_t val) {
		return val < 0 ? -val : val;
	}
}




namespace LIMAPOSITIONSYSTEM {
	
	/// <summary>
	/// Transfer external coordinates to internal multi-range LIMA coordinates
	/// </summary>
	/// <param name="state">Absolute positions of particles as float [nm]</param>
	/// <param name="key_particle_index">Index of centermost particle of compound</param>
	/// <returns></returns>
	static CompoundCoords positionCompound(CompoundState& state,  int key_particle_index=0) {
		CompoundCoords compoundcoords{};

		// WARNING: It may become a problem that state and state_prev does not share an origo. That should be fixed..
		compoundcoords.origo = Coord{ state.positions[key_particle_index] };
		//compoundcoords.origo = Coord(0);	// Temp, use the one above in future

		double default_norm_dist = 1.0;	// By default, 2 nm has a relative distance of 1.0 (int float) and 2^29 (uint32_t

		for (int i = 0; i < state.n_particles; i++) {
			double x = (static_cast<double>(state.positions[i].x) - static_cast<double>(compoundcoords.origo.x)) / default_norm_dist;
			double y = (static_cast<double>(state.positions[i].y) - static_cast<double>(compoundcoords.origo.y)) / default_norm_dist;
			double z = (static_cast<double>(state.positions[i].z) - static_cast<double>(compoundcoords.origo.z)) / default_norm_dist;

			Float3 rel_pos_nm{ x, y, z };
			Coord rel_coord{ rel_pos_nm * NANO_TO_LIMA };
			compoundcoords.rel_positions[i] = rel_coord;


		}
		return compoundcoords;
	}

	// Returns position in LimaMetres
	__device__ static Float3 getGlobalPosition(const CompoundCoords& coords) {
		return coords.origo.toFloat3() * NANO_TO_LIMA + coords.rel_positions[threadIdx.x].toFloat3();
	}

	// Transforms the relative coords to relative float positions.
	__device__ static void getRelativePositions(const Coord* coords, Float3* positions, const int n_elements) {
		if (threadIdx.x < n_elements)
			positions[threadIdx.x] = coords[threadIdx.x].toFloat3();
	}

	__device__ __host__ static bool canRepresentRelativeDist(const Coord& origo_a, const Coord& origo_b) {
		const auto diff = origo_a - origo_b;
		return std::abs(diff.x) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.y) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.z) < MAX_REPRESENTABLE_DIFF_NM;
	}




	
	// Calculate the necessary shift in LM of all elements of FROM, assuming the origo has been shifted to TO
	__device__ static Coord getRelShiftFromOrigoShift(const Coord& origo_from, const Coord& origo_to) {
		return (origo_from - origo_to) * static_cast<int32_t>(NANO_TO_LIMA);
	}

	__device__ static void applyHyperpos(const Coord& static_coord, Coord& movable_coord) {
		movable_coord.x += static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.x - movable_coord.x) > static_cast<int32_t>( BOX_LEN_HALF_NM));
		movable_coord.x -= static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.x - movable_coord.x) < static_cast<int32_t>(-BOX_LEN_HALF_NM));

		movable_coord.y += static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.y - movable_coord.y) > static_cast<int32_t>(BOX_LEN_HALF_NM));
		movable_coord.y -= static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.y - movable_coord.y) < static_cast<int32_t>(-BOX_LEN_HALF_NM));

		movable_coord.z += static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.z - movable_coord.z) > static_cast<int32_t>(BOX_LEN_HALF_NM));
		movable_coord.z -= static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.z - movable_coord.z) < static_cast<int32_t>(-BOX_LEN_HALF_NM));
	}



	//static void alignCoordinates(CompoundCoords& lhs, CompoundCoords& rhs) {
	//	applyHyperpos(lhs, rhs);

	//	Coord common_origo = {
	//		std::min(lhs.origo.x, rhs.origo.x) + ((std::max(lhs.origo.x, rhs.origo.x) - std::min(lhs.origo.x, rhs.origo.x)) / 2),
	//		std::min(lhs.origo.y, rhs.origo.y) + ((std::max(lhs.origo.y, rhs.origo.y) - std::min(lhs.origo.y, rhs.origo.y)) / 2),
	//		std::min(lhs.origo.z, rhs.origo.z) + ((std::max(lhs.origo.z, rhs.origo.z) - std::min(lhs.origo.z, rhs.origo.z)) / 2)
	//	};
	//}

	// Gets hyperorigo of other
	__device__ static Coord getHyperOrigo(const Coord& self, const Coord& other) {
		Coord temp = other;
		applyHyperpos(self, temp);
		return temp;
	}

	// The following two functions MUST ALWAYS be used together
	// Shift refers to the wanted difference in the relative positions, thus origo must move -shift.
	// Keep in mind that origo is in nm and rel_pos are in lm
	// ONLY CALL FROM THREAD 0
	__device__ static Coord shiftOrigo(CompoundCoords& coords, const int keyparticle_index=0) {
		const Coord shift_nm = coords.rel_positions[keyparticle_index] / static_cast<int32_t>(NANO_TO_LIMA);	// OPTIM. If LIMA wasn't 100 femto, but rather a power of 2, we could do this much better!
		coords.origo += shift_nm;
		return -shift_nm * static_cast<int32_t>(NANO_TO_LIMA);
	}
	__device__ static void shiftRelPos(CompoundCoords& coords, const Coord& shift_lm) {
		coords.rel_positions[threadIdx.x] += shift_lm;
	}

	__device__ static Coord getRelativeShiftBetweenCoordarrays(CompoundCoords* coordarray_circular_queue, int step, int compound_index_left, int compound_index_right) {
		Coord& coord_origo_left = CoordArrayQueueHelpers::getCoordarrayPtr(coordarray_circular_queue, step, compound_index_left)->origo;
		Coord& coord_origo_right = CoordArrayQueueHelpers::getCoordarrayPtr(coordarray_circular_queue, step, compound_index_right)->origo;

		Coord hyperorigo_right = LIMAPOSITIONSYSTEM::getHyperOrigo(coord_origo_left, coord_origo_right);

		// Calculate necessary shift in relative position for all particles in other, so they share origo with left
		return (coord_origo_left - hyperorigo_right) * static_cast<uint32_t>(NANO_TO_LIMA);	// This fucks up when the diff is > ~20
	}

	//// Calculates the relative position of movable_solvent, relative to another staticcoords's origo.
	//// Returns false if it is not possible to represent the position as coord. In that case, we should avoid
	//// following computations..
	//__device__ static Coord getRelativeHyperposition(const SolventCoord& static_solvent, const SolventCoord& movable_solvent) {
	//	const Coord hyperorigo_other = LIMAPOSITIONSYSTEM::getHyperOrigo(static_solvent.origo, movable_solvent.origo);

	//	// calc Relative Position Shift from the origo-shift
	//	//const Coord relPosShiftOfMovable = LIMAPOSITIONSYSTEM::getRelShift(static_solvent.origo, hyperorigo_other);
	//	// 
	//	//const Coord relPosShiftOfMovable = LIMAPOSITIONSYSTEM::getRelShiftFromOrigoShift(static_solvent.origo, hyperorigo_movable);
	//	const Coord relPosShiftOfMovable = LIMAPOSITIONSYSTEM::getRelShiftFromOrigoShift(hyperorigo_other, static_solvent.origo);
	//	return movable_solvent.rel_position + relPosShiftOfMovable;
	//}

	__device__ __host__ static void applyPBC(Coord& origo) {
		origo.x += BOX_LEN_NM_INT * (origo.x < 0);
		origo.x -= BOX_LEN_NM_INT * (origo.x >= BOX_LEN_NM_INT);
		origo.y += BOX_LEN_NM_INT * (origo.y < 0);
		origo.y -= BOX_LEN_NM_INT * (origo.y >= BOX_LEN_NM_INT);
		origo.z += BOX_LEN_NM_INT * (origo.z < 0);
		origo.z -= BOX_LEN_NM_INT * (origo.z >= BOX_LEN_NM_INT);
	}

	__device__ __host__ static void applyPBC(SolventCoord& coord) { applyPBC(coord.origo); }

	__device__ static void applyPBC(CompoundCoords& coords) {
		if (threadIdx.x != 0) { return; }
		applyPBC(coords.origo);		
	}





	// ReCenter origo of solvent, and the relative pos around said origo
	__device__ static void updateSolventcoord(SolventCoord& coord) {
		const int shift_at = static_cast<int32_t>(NANO_TO_LIMA) / 2;
		Coord shift_nm = coord.rel_position / static_cast<int32_t>(shift_at);	// OPTIM. If LIMA wasn't 100 femto, but rather a power of 2, we could do this much better! 

		SolventCoord tmp = coord;
		coord.origo += shift_nm;
		coord.rel_position -= shift_nm * static_cast<int32_t>(NANO_TO_LIMA);
		
		if (blockIdx.x + threadIdx.x == 0 && shift_nm.x != 0) {
			tmp.origo.print('o');
			tmp.rel_position.print('r');
			shift_nm.print('s');
			coord.origo.print('O');
			coord.rel_position.print('R');
		}
	}

	// Get the relpos_prev, if the solvent was in the same solventblock last step
	__device__ static Coord getRelposPrev(SolventBlockGrid* solventblockgrid_circularqueue, const int solventblock_id, const int current_step) {
		const int step_prev = current_step == 0 ? 0 : current_step - 1;	// Unnecessary if we use transfermodule's prevpos for first step!!!!!!!!! TODOTODO
		auto blockPtr = CoordArrayQueueHelpers::getSolventBlockPtr(solventblockgrid_circularqueue, step_prev, solventblock_id);
		return blockPtr->rel_pos[threadIdx.x];
	}








	//__device__ static Coord getOnehotDirectionEdgeblock(const Coord relpos, const Coord thresholdInward, const Coord thresholdOutward)

	__device__ static Coord getOnehotDirection(const Coord relpos, const int32_t threshold) {
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
				return Coord{ relpos.x < 0 ? -1 : 1, 0, 0 };
			}
			else if (magnitude_y >= magnitude_z) { // The y component has the largest magnitude			
				return Coord{ 0, relpos.y < 0 ? -1 : 1, 0 };
			}
			else { // The z component has the largest magnitude		
				return Coord{ 0, 0, relpos.z < 0 ? -1 : 1 };
			}
		}
		else {
			return Coord{ 0 };
		}		
	}

	// Since coord is rel to 0,0,0 of a block, we need to offset the positions so they are scattered around the origo instead of above it
	// We also need a threshold of half a blocklen, otherwise we should not transfer, and return{0,0,0}
	__device__ static Coord getTransferDirection(const Coord relpos) {
		//const int32_t blocklen_half = static_cast<int32_t>(NANO_TO_LIMA) / 2;
		const int32_t blocklen_half = static_cast<int32_t>(SolventBlockGrid::node_len) / 2;
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
	/// Shifts the position 1/2 blocklen so we can find the appropriate origo by 
	/// </summary>
	/// <param name="position">Absolute position of solvent [nm] </param>
	__host__ static SolventCoord createSolventcoordFromAbsolutePosition(const Float3& position) {
		//const Float3 blockcenter_relative{ SolventBlockGrid::node_len / 2.f / NANO_TO_LIMA };	// [nm]
		//const Float3 position_adjusted = position + blockcenter_relative;						// [nm]

		//const Float3 origo_f = position_adjusted.piecewiseRound();								// [nm]
		const Float3 origo_f = position.piecewiseRound();
		const Float3 relpos_f = (position - origo_f) * NANO_TO_LIMA;							// [lm}

		SolventCoord solventcoord{ Coord{origo_f}, Coord{relpos_f } };
		applyPBC(solventcoord);
		return solventcoord;
	}





	//__device__ static bool willSolventRe

	// Get the relpos_prev, if the solvent was NOT in the same solventblock last step
	//__device__ static Coord getRelposPrevAfterTransfer(SolventBlockGrid* solventblockgrid_circularqueue, const int solventblock_id, const int step) {


	//__device__ static void applyPBC(Compound* compound);
};








namespace EngineUtils {
	__device__ __host__ static inline void applyHyperpos(const Float3* static_particle, Float3* movable_particle) {
		for (int i = 0; i < 3; i++) {
			//*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
			//*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
			*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);
		}
	}

	__device__ __host__ static float calcHyperDist(const Float3* p1, const Float3* p2) {
		Float3 temp = *p2;
		applyHyperpos(p1, &temp);
		return (*p1 - temp).len();
	}

	__device__ __host__ static float calcHyperDistDenormalized(const Float3* p1, const Float3* p2) {
		Float3 temp = *p2;
		applyHyperpos(p1, &temp);
		return (*p1 - temp).len() * NORMALIZER;
	}

	__device__ __host__ static float calcKineticEnergy(const Float3* pos1, const Float3* pos2, const float mass, const float elapsed_time) {
		const float vel = calcHyperDist(pos1, pos2) / elapsed_time * NORMALIZER;
		const float kinE = 0.5f * mass * vel * vel;
		return kinE;
	}

	__device__ __host__ static void applyPBC(Float3* current_position) {	// Only changes position if position is outside of box;
		for (int dim = 0; dim < 3; dim++) {
			/**current_position->placeAt(dim) += BOX_LEN * (current_position->at(dim) < 0.f);
			*current_position->placeAt(dim) -= BOX_LEN * (current_position->at(dim) > BOX_LEN);*/
			*current_position->placeAt(dim) += BOX_LEN * (current_position->at(dim) < 0.f);
			*current_position->placeAt(dim) -= BOX_LEN * (current_position->at(dim) > BOX_LEN);
		}
	}

	static void __host__ genericErrorCheck(const char* text) {
		cudaDeviceSynchronize();
		cudaError_t cuda_status = cudaGetLastError();
		if (cuda_status != cudaSuccess) {
			std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
			fprintf(stderr, text);
			exit(1);
		}
	}

	static float calcSpeedOfParticle(const float mass /*[kg]*/, const float temperature /*[K]*/) { // 
		const float R = 8.3144f;								// Ideal gas constants - J/(Kelvin*mol)
		const float v_rms = static_cast<float>(sqrt(3.f * R * temperature / mass));
		return v_rms;	// [m/s]
	}

	// For solvents, compound_id = n_compounds and particle_id = solvent_index
	__device__ static uint32_t getLoggingIndexOfParticle(uint32_t step, uint32_t total_particles_upperbound, uint32_t compound_id, uint32_t particle_id_local) {
		const uint32_t steps_since_transfer = (step % STEPS_PER_LOGTRANSFER);
		const uint32_t step_offset = steps_since_transfer * total_particles_upperbound;
		const uint32_t compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
		return step_offset + compound_offset + particle_id_local;
	}

	__host__ static size_t getAlltimeIndexOfParticle(uint64_t step, uint32_t total_particles_upperbound, uint32_t compound_id, uint32_t particle_id_local) {
		const uint32_t step_offset = static_cast<uint32_t>(step) * total_particles_upperbound;
		const uint32_t compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
		return step_offset + compound_offset + particle_id_local;
	}

	__device__ int static getNewBlockId(const Coord& transfer_direction, const Coord& origo) {
		auto newCoord3d = transfer_direction + origo;
		LIMAPOSITIONSYSTEM::applyPBC(newCoord3d);
		return SolventBlockGrid::get1dIndex(newCoord3d);
	}

	__device__ static SolventBlockTransfermodule* getTransfermoduleTargetPtr(SolventBlockTransfermodule* transfermodule_array, int blockId, const Coord& transfer_direction) {
		Coord new3dIndex = SolventBlockGrid::get3dIndex(blockId) + transfer_direction;
		LIMAPOSITIONSYSTEM::applyPBC(new3dIndex);
		auto index = SolventBlockGrid::get1dIndex(new3dIndex);
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
	__device__ static Coord integratePosition(const Coord& coord, const Coord& coord_tsub1, const Float3* force, const float mass, const double dt, const float thermostat_scalar) {

		return coord + (coord - coord_tsub1) * thermostat_scalar + Coord{ *force * dt * dt / mass };
		//return coord + Coord{ 1000000, 0, 0 };
		if (threadIdx.x == 0 && blockIdx.x == 0) {
			//force->print('f');
			//printf("dt %f mass %f\n", dt, mass);
			//(*force*dt*dt / mass).print('F');
			//uint32_t diff = coord.x - x - dx;
			//printf("x  %d  dx %d  force %.10f ddx %d    x_ %d   dif %d\n", x, dx, force->x, ddx, dx + ddx, diff);
		}
	}


};


namespace LIMADEBUG {
	__device__ void static transferOut(STransferQueue* queue_global, const STransferQueue& queue_local, const Coord& transferdir_queue, const int queue_index) {
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

};





// LIMA algorithm Library
namespace LAL {
	__device__ constexpr int getBlellochTablesize(int n) {
		float nf = static_cast<float>(n);
		return CPPD::ceil(nf * log2f(nf) * 2.f);
	}
}