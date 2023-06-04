#pragma once

#include "LIMA_BASE/include/LimaTypes.cuh"
#include <iostream>
#include <chrono>
#include <cstring>
#include <string>
#include "LIMA_BASE/include/Simulation.cuh"
#include "LIMA_ENGINE/include/Engine.cuh"
#include <vector>

//#include "Forcefield.cuh"



class Analyzer {
public:
	Analyzer() {}

	struct AnalyzedPackage {
		AnalyzedPackage() = default;
		AnalyzedPackage(std::vector<Float3>& avg_energy, float* t_ptr, int t_cnt) {	// TODO: t_ptr should be a vector
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

		std::vector<float> pot_energy;
		std::vector<float> kin_energy;
		std::vector<float> total_energy;
		std::vector<Float3> energy_data; // potE, kinE, totalE
		int n_energy_values = 0;

		std::vector<float> temperature_data;
		int n_temperature_values = 0;
	};

	AnalyzedPackage analyzeEnergy(Simulation* simulation); // Prints a file of doubles: [step, molecule, atom, coordinate_dim]

	Float3* analyzeSolvateEnergy(Simulation* simulation, uint64_t n_steps);
	std::vector<Float3> analyzeCompoundEnergy(Simulation* simulation, uint64_t n_steps);

	void moveAndPadData(Simulation* sim, uint64_t steps_in_kernel, uint64_t step_offset);

	static void printEnergy(AnalyzedPackage* package);
	static float getVarianceCoefficient(std::vector<float>& vec);

private:
	Engine engine;



	Float3* traj_buffer_device = nullptr;
	float* potE_buffer_device = nullptr;
};
