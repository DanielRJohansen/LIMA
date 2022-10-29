#pragma once

#include<iostream>
#include "LimaTypes.cuh"
#include "Constants.cuh"

#include "Simulation.cuh"
#include "Forcefield.cuh"



namespace EngineUtils {
	static inline void __device__ __host__ applyHyperpos(Float3* static_particle, Float3* movable_particle) {
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

	static Float3 __host__ getBoxTemperature(Simulation* simulation, ForceField& forcefield_host) {

		const uint64_t step = simulation->getStep() - 1;

		long double temp_sum = 0;				// [k]

		long double kinE_sum = 0;
		float biggest_contribution = 0;


		const int step_offset_a = step * simulation->total_particles_upperbound;
		const int step_offset_b = (step - 2) * simulation->total_particles_upperbound;
		const int solvent_offset = MAX_COMPOUND_PARTICLES * simulation->n_compounds;


		for (int c = 0; c < simulation->n_compounds; c++) {
			int compound_offset = c * MAX_COMPOUND_PARTICLES;
			for (int i = 0; i < simulation->compounds_host[c].n_particles; i++) {	// i gotta move this somewhere else....

				Float3 posa = simulation->traj_buffer[i + compound_offset + step_offset_a];
				Float3 posb = simulation->traj_buffer[i + compound_offset + step_offset_b];
				float kinE = EngineUtils::calcKineticEnergy(&posa, &posb, forcefield_host.particle_parameters[simulation->compounds_host[c].atom_types[i]].mass, simulation->dt);			// Doesnt work, use forcefield_host!!
				//float kinE = EngineUtils::calcKineticEnergy(&posa, &posb, forcefield_device.particle_parameters[simulation->box->compounds[c].atom_types[i]].mass, simulation->dt);			// Doesnt work, use forcefield_host!!
				//printf("kinE %f\n", kinE);
				//printf("mass %f\n", forcefield_host.particle_parameters[simulation->compounds_host[c].atom_types[i]].mass);
				biggest_contribution = std::max(biggest_contribution, kinE);

				kinE_sum += kinE;
				//particles_total++;
			}
		}
		//printf("\nKin e %Lf\n", kinE_sum);
		//kinE_sum = 0;
		for (int i = 0; i < simulation->n_solvents; i++) {
			Float3 posa = simulation->traj_buffer[i + solvent_offset + step_offset_a];
			Float3 posb = simulation->traj_buffer[i + solvent_offset + step_offset_b];
			float kinE = EngineUtils::calcKineticEnergy(&posa, &posb, SOLVENT_MASS, simulation->dt);
			biggest_contribution = max(biggest_contribution, kinE);
			kinE_sum += kinE;
		}
		//double avg_kinE = kinE_sum / (long double)particles_total;
		double avg_kinE = kinE_sum / (long double)simulation->total_particles;
		float temperature = avg_kinE * 2.f / (3.f * 8.3145);
		//printf("\nTemp: %f\n", temperature);
		return Float3(temperature, biggest_contribution, avg_kinE);
	}
};