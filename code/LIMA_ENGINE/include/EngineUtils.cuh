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

	__device__ static void applyPBC(Float3* current_position) {	// Only changes position if position is outside of box;
		for (int dim = 0; dim < 3; dim++) {
			/**current_position->placeAt(dim) += BOX_LEN * (current_position->at(dim) < 0.f);
			*current_position->placeAt(dim) -= BOX_LEN * (current_position->at(dim) > BOX_LEN);*/
			*current_position->placeAt(dim) += BOX_LEN_RELATIVE * (current_position->at(dim) < 0.f);
			*current_position->placeAt(dim) -= BOX_LEN_RELATIVE * (current_position->at(dim) > BOX_LEN_RELATIVE);
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

	__device__ static void getRelativePositions(CompoundCoords& coords, CompoundState& state) {
		state.positions[threadIdx.x] = coords.rel_positions[threadIdx.x].toFloat3();
	}

	static void applyHyperpos(CompoundCoords& lhs, CompoundCoords& rhs) {
		// TODO: IMplement
	}
	
	// For coordinates of OTHER, we find the value in LM that each coord must be shifted, to be aligned with coordinates of self
	__device__ static Coord getRelShift(const Coord& origo_self, const Coord& origo_other) {
		return (origo_self - origo_other) * static_cast<int32_t>(NANO_TO_LIMA);
	}

	__device__ static void applyHyperpos(const Coord& static_coord, Coord& movable_coord) {
		movable_coord.x += static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.x - movable_coord.x) > static_cast<int32_t>( BOX_LEN_HALF_NM));
		movable_coord.x -= static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.x - movable_coord.x) < static_cast<int32_t>(-BOX_LEN_HALF_NM));

		movable_coord.y += static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.y - movable_coord.y) > static_cast<int32_t>(BOX_LEN_HALF_NM));
		movable_coord.y -= static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.y - movable_coord.y) < static_cast<int32_t>(-BOX_LEN_HALF_NM));

		movable_coord.z += static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.z - movable_coord.z) > static_cast<int32_t>(BOX_LEN_HALF_NM));
		movable_coord.z -= static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.z - movable_coord.z) < static_cast<int32_t>(-BOX_LEN_HALF_NM));
	}


	static void alignCoordinates(CompoundCoords& lhs, CompoundCoords& rhs) {
		applyHyperpos(lhs, rhs);

		Coord common_origo = {
			std::min(lhs.origo.x, rhs.origo.x) + ((std::max(lhs.origo.x, rhs.origo.x) - std::min(lhs.origo.x, rhs.origo.x)) / 2),
			std::min(lhs.origo.y, rhs.origo.y) + ((std::max(lhs.origo.y, rhs.origo.y) - std::min(lhs.origo.y, rhs.origo.y)) / 2),
			std::min(lhs.origo.z, rhs.origo.z) + ((std::max(lhs.origo.z, rhs.origo.z) - std::min(lhs.origo.z, rhs.origo.z)) / 2)
		};
	}

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

	//__device__ static void applyPBC(Compound* compound);
};