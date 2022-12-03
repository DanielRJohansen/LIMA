#pragma once

#include "LimaTypes.cuh"
#include <iostream>
#include <chrono>
#include <cstring>
#include <string>
#include "Simulation.cuh"
#include "Engine.cuh"
#include <vector>

//#include "Forcefield.cuh"



class Analyzer {
public:
	Analyzer() {}

	struct AnalyzedPackage {
		AnalyzedPackage() = default;
		AnalyzedPackage(std::vector<Float3>& avg_energy, float* t_ptr, int t_cnt) {
			//energy_data.resize(e_cnt);
			//memcpy(energy_data.data(), e_ptr, e_cnt);
			energy_data = avg_energy;
			auto e_cnt = energy_data.size();

			temperature_data.resize(t_cnt);
			memcpy(temperature_data.data(), t_ptr, t_cnt);

			pot_energy.resize(e_cnt);
			kin_energy.resize(e_cnt);
			total_energy.resize(e_cnt);
			for (int i = 0; i < e_cnt; i++) {
				pot_energy[i] = energy_data[i].x;
				kin_energy[i] = energy_data[i].y;
				total_energy[i] = energy_data[i].z;
			}
		}
		/*: n_energy_values(e_cnt), n_temperature_values(t_cnt) {
			energy_data = e_ptr;
			temperature_data = t_ptr;
		}*/
		std::vector<float> pot_energy;
		std::vector<float> kin_energy;
		std::vector<float> total_energy;
		std::vector<Float3> energy_data; // potE, kinE, totalE
		int n_energy_values = 0;

		std::vector<float> temperature_data;
		int n_temperature_values = 0;
		~AnalyzedPackage() {
			//delete[] energy_data;
			//delete[] temperature_data;
		}
	};

	AnalyzedPackage analyzeEnergy(Simulation* simulation); // Prints a file of doubles: [step, molecule, atom, coordinate_dim]

	Float3* analyzeSolvateEnergy(Simulation* simulation, uint64_t n_steps);
	std::vector<Float3> analyzeCompoundEnergy(Simulation* simulation, uint64_t n_steps);

	void moveAndPadData(Simulation* sim, uint64_t steps_in_kernel, uint64_t step_offset);

	static void printEnergy(AnalyzedPackage* package);
	static float getStdDevNorm(std::vector<float>& vec);

private:
	Engine engine;



	Float3* traj_buffer_device = nullptr;
	float* potE_buffer_device = nullptr;
};
