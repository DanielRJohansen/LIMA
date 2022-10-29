#pragma once

#include<iostream>
#include "LimaTypes.cuh"
#include "Constants.cuh"
struct EngineUtils {
	static void __device__ __host__ applyHyperpos(Float3* static_particle, Float3* movable_particle) {
		//#pragma unroll
		for (int i = 0; i < 3; i++) {
			*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
			*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		}
	}

	static float __device__ __host__ calcKineticEnergy(Float3* pos1, Float3* pos2, float mass, double dt) {	// pos1/2 MUST be 2 steps apart!!!!
		EngineUtils::applyHyperpos(pos1, pos2);

		if ((*pos1 - *pos2).len() > 1) {
			//printf("KinE Dist over 1 nm!\n");
			//pos1->print('1');
			//pos2->print('2');
		}


		float vel = (*pos1 - *pos2).len() * (float)(0.5f / dt);
		float kinE = 0.5f * mass * vel * vel;
		return kinE;
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
};