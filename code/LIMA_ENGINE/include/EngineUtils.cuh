#pragma once

#include<iostream>
#include "LimaTypes.cuh"
#include "Constants.cuh"

#include "Simulation.cuh"
#include "Forcefield.cuh"
#include <math.h>

#include <cuda.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>

namespace ForceCalc {
//	-
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
			*current_position->placeAt(dim) += BOX_LEN *(current_position->at(dim) < 0.f);
			*current_position->placeAt(dim) -= BOX_LEN * (current_position->at(dim) > BOX_LEN);
		}
	}

	static void __host__ genericErrorCheck(const char* text) {
		cudaDeviceSynchronize();
		cudaError_t cuda_status = cudaGetLastError();
		if (cuda_status != cudaSuccess) {
			std::cout << "\nCuda error code: " << cuda_status << std::endl;
			fprintf(stderr, text);
			exit(1);
		}
	}

	//__device__ __host__ static void applyABC(Coord* pos) {

	//}

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


};







namespace LIMAPOSITIONSYSTEM {
	//static Coord coordFromAbsPos(Float3 pos)
	static CompoundCoords positionCompound(CompoundState& state,  int key_particle_index=0) {
		CompoundCoords compoundcoords{};
		//Float3& key_pos = state.positions[key_particle_index];
		//compoundcoords.origo = Float3(key_pos);	
		compoundcoords.origo = Float3(0);	// Temp, use the one above in future

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

	__device__ static Float3 getGlobalPositionNM(CompoundCoords& coords) {
		return coords.origo.toFloat3() + coords.rel_positions[threadIdx.x].toFloat3() / NANO_TO_LIMA;
	}

	__device__ static Float3 getGlobalPositionFM(CompoundCoords& coords) {
		return coords.origo.toFloat3() * NANO_TO_FEMTO + coords.rel_positions[threadIdx.x].toFloat3() * LIMA_TO_FEMTO;
	}

	// Returns position in LimaMetres
	__device__ static Float3 getGlobalPosition(const CompoundCoords& coords) {
		return coords.origo.toFloat3() * NANO_TO_LIMA + coords.rel_positions[threadIdx.x].toFloat3();
	}

	// Returns positions in LimaMetres
	__device__ static void getGlobalPositions(CompoundCoords& coords, CompoundState& state) {
		state.positions[threadIdx.x] = getGlobalPosition(coords);
	}

	__device__ __host__ static bool canRepresentRelativeDist(const Coord& origo_a, const Coord& origo_b) {
		const auto diff = origo_a - origo_b;
		return std::abs(diff.x) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.y) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.z) < MAX_REPRESENTABLE_DIFF_NM;
	}

	__device__ static void getRelativePositions(Coord* coords, Float3* positions) {
		positions[threadIdx.x] = coords[threadIdx.x].toFloat3();
	}

	__device__ static void getRelativePositions(CompoundCoords& coords, CompoundState& state) {
		//state.positions[threadIdx.x] = coords.rel_positions[threadIdx.x].toFloat3();
		//or 
		getRelativePositions(coords.rel_positions, state.positions);
	}

	//static void applyHyperpos(CompoundCoords& lhs, CompoundCoords& rhs) {
	//	// TODO: IMplement
	//}
	
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
		//Coord shift_nm = -coords.rel_positions[keyparticle_index] / static_cast<uint32_t>(NANO_TO_LIMA);	// OPTIM. If LIMA wasn't 100 femto, but rather a power of 2, we could do this much better!

		Coord shift_nm = coords.rel_positions[keyparticle_index] / static_cast<int32_t>(NANO_TO_LIMA);	// OPTIM. If LIMA wasn't 100 femto, but rather a power of 2, we could do this much better!
		
		if (threadIdx.x == 0 && blockIdx.x == 1) {
			//printf("x %d shift.x %d  shift_out%d\n", coords.rel_positions[keyparticle_index].x, shift_nm.x, shift_nm.x * static_cast<int32_t>(NANO_TO_LIMA));
			//printf("%d %d %d\n", shift_nm.x, shift_nm.y, shift_nm.z);
		}
		
		coords.origo += shift_nm;
		return -shift_nm * static_cast<int32_t>(NANO_TO_LIMA);
	}


	__device__ static void shiftRelPos(CompoundCoords& coords, const Coord& shift_lm) {
		coords.rel_positions[threadIdx.x] += shift_lm;
	}


	static void moveCoordinate(Coord& coord, Float3 delta /*[nm]*/) {	// TODO: some checks to make this safe
		Coord delta_c{ delta / LIMA_TO_FEMTO };
		coord += delta_c;
	}

	__device__ static void applyPBC(CompoundCoords& coords) {
		if (threadIdx.x != 0) { return; }
		coords.origo.x += static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.x < 0);
		coords.origo.x -= static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.x > static_cast<int32_t>(BOX_LEN_NM));
		coords.origo.y += static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.y < 0);
		coords.origo.y -= static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.y > static_cast<int32_t>(BOX_LEN_NM));
		coords.origo.z += static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.z < 0);
		coords.origo.z -= static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.z > static_cast<int32_t>(BOX_LEN_NM));
	}

	// DANGEROUS: Only wraps around once! Min is assumed to be 0
	__device__ static int getIntegarWraparound(int num, const int max) {
		num += max * (num < 0);
		num -= max * (num > max);
		return num;
	}

	__device__ static Coord getRelativeShiftBetweenCoordarrays(CompoundCoords* coordarray_circular_queue, int step, int compound_index_left, int compound_index_right) {
		Coord& coord_origo_left = CoordArrayQueueHelpers::getCoordarrayPtr(coordarray_circular_queue, step, compound_index_left)->origo;
		Coord& coord_origo_right = CoordArrayQueueHelpers::getCoordarrayPtr(coordarray_circular_queue, step, compound_index_right)->origo;

		Coord hyperorigo_right = LIMAPOSITIONSYSTEM::getHyperOrigo(coord_origo_left, coord_origo_right);

		// Calculate necessary shift in relative position for all particles in other, so they share origo with left
		return (coord_origo_left - hyperorigo_right) * static_cast<uint32_t>(NANO_TO_LIMA);	// This fucks up when the diff is > ~20
	}

	// Calculates the relative position of movable_solvent, relative to another staticcoords's origo.
	// Returns false if it is not possible to represent the position as coord. In that case, we should avoid
	// following computations..
	__device__ static Coord getRelativeHyperposition(const SolventCoord& static_solvent, const SolventCoord& movable_solvent) {
		const Coord hyperorigo_other = LIMAPOSITIONSYSTEM::getHyperOrigo(static_solvent.origo, movable_solvent.origo);

		// calc Relative Position Shift from the origo-shift
		//const Coord relPosShiftOfMovable = LIMAPOSITIONSYSTEM::getRelShift(static_solvent.origo, hyperorigo_other);
		// 
		//const Coord relPosShiftOfMovable = LIMAPOSITIONSYSTEM::getRelShiftFromOrigoShift(static_solvent.origo, hyperorigo_movable);
		const Coord relPosShiftOfMovable = LIMAPOSITIONSYSTEM::getRelShiftFromOrigoShift(hyperorigo_other, static_solvent.origo);
		return movable_solvent.rel_position + relPosShiftOfMovable;
	}

	__device__ __host__ static void applyPBC(SolventCoord& coord) {	// Only changes position if position is outside of box;
		Coord& pos = coord.origo;
		pos.x += BOX_LEN_NM * (pos.x < 0);
		pos.x -= BOX_LEN_NM * (pos.x >= BOX_LEN_NM);
		pos.y += BOX_LEN_NM * (pos.y < 0);
		pos.y -= BOX_LEN_NM * (pos.y >= BOX_LEN_NM);
		pos.z += BOX_LEN_NM * (pos.z < 0);
		pos.z -= BOX_LEN_NM * (pos.z >= BOX_LEN_NM);
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
	__device__ static Coord getRelposPrev(SolventBlockGrid* solventblockgrid_circularqueue, const int solventblock_id, const int step) {
		const int step_prev = step == 0 ? 0 : step - 1;
		auto blockPtr = CoordArrayQueueHelpers::getSolventBlockPtr(solventblockgrid_circularqueue, step_prev, solventblock_id);
		return blockPtr->rel_pos[threadIdx.x];
	}

	// Get the relpos_prev, if the solvent was NOT in the same solventblock last step
	//__device__ static Coord getRelposPrevAfterTransfer(SolventBlockGrid* solventblockgrid_circularqueue, const int solventblock_id, const int step) {


	//__device__ static void applyPBC(Compound* compound);
};

namespace CPPD {
	constexpr int32_t ceil(float num) {
		return (static_cast<float>(static_cast<int32_t>(num)) == num)
			? static_cast<int32_t>(num)
			: static_cast<int32_t>(num) + ((num > 0) ? 1 : 0);
	}
}