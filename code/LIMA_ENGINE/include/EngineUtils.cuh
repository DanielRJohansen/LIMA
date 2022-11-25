#pragma once

#include<iostream>
#include "LimaTypes.cuh"
#include "Constants.cuh"

#include "Simulation.cuh"
#include "Forcefield.cuh"

namespace ForceCalc {
//	-
};

namespace EngineUtils {
	__device__ __host__ static inline void applyHyperpos(const Float3* static_particle, Float3* movable_particle) {
		//#pragma unroll
		for (int i = 0; i < 3; i++) {
			*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
			*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		}
	}

	__device__ __host__ static inline float calcHyperDist(const Float3* p1, const Float3* p2) {
		Float3 temp = *p2;
		applyHyperpos(p1, &temp);
		return (*p1 - temp).len();
	}

	__device__ __host__ static float calcKineticEnergy(const Float3* pos1, const Float3* pos2, const float mass, const float elapsed_time) {
		const float vel = calcHyperDist(pos1, pos2) / elapsed_time;
		const float kinE = 0.5f * mass * vel * vel;
		return kinE;
	}

	__device__ static void applyPBC(Float3* current_position) {	// Only changes position if position is outside of box;
		for (int dim = 0; dim < 3; dim++) {
			*current_position->placeAt(dim) += BOX_LEN * (current_position->at(dim) < 0.f);
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

	__device__ __host__ static void applyABC(Coord* pos) {

	}

	static float calcSpeedOfParticle(const float mass /*[kg]*/, const float temperature /*[K]*/) { // 
		const float R = 8.3144;								// Ideal gas constants - J/(Kelvin*mol)
		//double mean_velocity = M / (2 * k_B * T);				// This is bullshit. Only for sol mass
		const float v_rms = static_cast<float>(sqrt(3 * R * temperature / mass));
		return v_rms;
	}



};