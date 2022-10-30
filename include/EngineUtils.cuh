#pragma once

#include<iostream>
#include "LimaTypes.cuh"
#include "Constants.cuh"

#include "Simulation.cuh"
#include "Forcefield.cuh"



namespace EngineUtils {
	__device__ __host__ static inline void applyHyperpos(Float3* static_particle, Float3* movable_particle) {
		//#pragma unroll
		for (int i = 0; i < 3; i++) {
			*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
			*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		}
	}

	__device__ __host__ static float calcKineticEnergy(Float3* pos1, Float3* pos2, float mass, double dt) {	// pos1/2 MUST be 2 steps apart!!!!
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

	static Float3 __host__ getBoxTemperature(Simulation* simulation, ForceField& forcefield_host) {
		const uint64_t step = simulation->getStep() - 1;

		long double kinE_sum = 0;	// [k]
		float biggest_contribution = 0;


		const uint64_t step_offset_a = step * simulation->total_particles_upperbound;
		const uint64_t step_offset_b = (step - 2) * simulation->total_particles_upperbound;
		const uint64_t solvent_offset = MAX_COMPOUND_PARTICLES * simulation->n_compounds;


		for (int c = 0; c < simulation->n_compounds; c++) {
			uint64_t compound_offset = c * MAX_COMPOUND_PARTICLES;
			for (uint64_t i = 0; i < simulation->compounds_host[c].n_particles; i++) {	// i gotta move this somewhere else....

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
			float kinE = EngineUtils::calcKineticEnergy(&posa, &posb, forcefield_host.particle_parameters[0].mass, simulation->dt);
			biggest_contribution = max(biggest_contribution, kinE);
			kinE_sum += static_cast<float>(kinE);
		}
		//double avg_kinE = kinE_sum / (long double)particles_total;
		float avg_kinE = static_cast<float>(kinE_sum / static_cast<long double>(simulation->total_particles));
		float temperature = avg_kinE * 2.f / (3.f * 8.3145f);
		//printf("\nTemp: %f\n", temperature);
		return Float3(temperature, biggest_contribution, avg_kinE);
	}
};