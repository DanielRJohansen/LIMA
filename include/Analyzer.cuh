#pragma once

#include <iostream>
#include <chrono>
#include <cstring>
#include <string>
#include "Simulation.cuh"
#include "Engine.cuh"


//#include "Forcefield.cuh"



class Analyzer {
public:
	Analyzer() {}

	struct AnalyzedPackage {
		AnalyzedPackage(Float3* e_ptr, int e_cnt, float* t_ptr, int t_cnt) : n_energy_values(e_cnt), n_temperature_values(t_cnt) {
			energy_data = e_ptr;
			temperature_data = t_ptr;
		}
		Float3* energy_data; // potE, kinE, totalE
		int n_energy_values;

		float* temperature_data;
		int n_temperature_values;
		~AnalyzedPackage() {
			delete[] energy_data;
			delete[] temperature_data;
		}
	};

	AnalyzedPackage analyzeEnergy(Simulation* simulation); // Prints a file of doubles: [step, molecule, atom, coordinate_dim]

	Float3* analyzeSolvateEnergy(Simulation* simulation, uint64_t n_steps);
	Float3* analyzeCompoundEnergy(Simulation* simulation, uint64_t n_steps);



private:
	Engine engine;



	Float3* traj_buffer_device;
	float* potE_buffer_device;
};
