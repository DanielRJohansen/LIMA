#pragma once

#include<iostream>
#include "LimaTypes.cuh"
#include "Constants.cuh"

#include "Simulation.cuh"
#include "Forcefield.cuh"
#include <math.h>

namespace ForceCalc {
//	-
};

namespace EngineUtils {
	__device__ __host__ static inline void applyHyperpos(const Float3* static_particle, Float3* movable_particle) {
		for (int i = 0; i < 3; i++) {
			//*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
			//*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			*movable_particle->placeAt(i) += BOX_LEN_RELATIVE * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_RELATIVE_HALF);
			*movable_particle->placeAt(i) -= BOX_LEN_RELATIVE * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_RELATIVE_HALF);
		}
	}

	__device__ __host__ static float calcHyperDist(const Float3* p1, const Float3* p2) {
		Float3 temp = *p2;
		//applyHyperpos(p1, &temp);
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
		//double mean_velocity = M / (2 * k_B * T);				// This is bullshit. Only for sol mass
		const float v_rms = static_cast<float>(sqrt(3.f * R * temperature / mass));
		return v_rms;
	}



};







namespace LIMAPOSITIONSYSTEM {
	//static Coord coordFromAbsPos(Float3 pos)
	static CompoundCoords positionCompound(CompoundState& state,  int key_particle_index=0) {
		CompoundCoords compoundcoords{};
		Float3& key_pos = state.positions[key_particle_index];
		compoundcoords.origo = Float3(key_pos);

		double default_norm_dist = 1.0;	// By default, 2 nm has a relative distance of 1.0 (int float) and 2^29 (uint32_t

		for (int i = 0; i < state.n_particles; i++) {
			double x = (static_cast<double>(state.positions[i].x) - static_cast<double>(compoundcoords.origo.x)) / default_norm_dist;
			double y = (static_cast<double>(state.positions[i].y) - static_cast<double>(compoundcoords.origo.y)) / default_norm_dist;
			double z = (static_cast<double>(state.positions[i].z) - static_cast<double>(compoundcoords.origo.z)) / default_norm_dist;

			Float3 rel_pos{ x, y, z };
			Coord rel_coord{ rel_pos * static_cast<float>(1 << 29) };
			compoundcoords.rel_positions[i] = rel_coord;

			if (rel_pos.len() > 1.f) {
				throw "Compound spans too large a distance";
			}
		}
		return compoundcoords;
	}

	__device__ static Float3 getGlobalPositionNM(CompoundCoords& coords) {
		return coords.origo.toFloat3() / 1e+6f + coords.rel_positions[threadIdx.x].toFloat3() / static_cast<float>(1 << 29);
	}

	__device__ static Float3 getGlobalPositionFM(CompoundCoords& coords) {
		return coords.origo.toFloat3() * NANO_TO_FEMTO*0.f + coords.rel_positions[threadIdx.x].toFloat3() * LIMASCALE_TO_FEMTO;
	}

	__device__ static void getGlobalPositionsFM(CompoundCoords& coords, CompoundState& state) {
		//state.positions[threadIdx.x] = coords.origo.toFloat3() / 1e+6f + coords.rel_positions[threadIdx.x].toFloat3() / static_cast<float>(1 << 29);
		state.positions[threadIdx.x] = getGlobalPositionFM(coords);
	}

	static void applyHyperpos(CompoundCoords& lhs, CompoundCoords& rhs) {
		// TODO: IMplement
	}

	static void alignCoordinates(CompoundCoords& lhs, CompoundCoords& rhs) {
		applyHyperpos(lhs, rhs);

		Coord common_origo = {
			std::min(lhs.origo.x, rhs.origo.x) + ((std::max(lhs.origo.x, rhs.origo.x) - std::min(lhs.origo.x, rhs.origo.x)) / 2),
			std::min(lhs.origo.y, rhs.origo.y) + ((std::max(lhs.origo.y, rhs.origo.y) - std::min(lhs.origo.y, rhs.origo.y)) / 2),
			std::min(lhs.origo.z, rhs.origo.z) + ((std::max(lhs.origo.z, rhs.origo.z) - std::min(lhs.origo.z, rhs.origo.z)) / 2)
		};

	}

	static void moveCoordinate(Coord& coord, Float3 delta /*[nm]*/) {	// TODO: some checks to make this safe
		Coord delta_c{ delta / LIMASCALE_TO_FEMTO };
		coord += delta_c;
	}


	//__device__ static void applyPBC(Compound* compound);
};