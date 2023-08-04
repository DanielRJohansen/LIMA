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
	Analyzer(std::unique_ptr<LimaLogger> logger) : m_logger(std::move(logger)) {}

	struct AnalyzedPackage {
		AnalyzedPackage() = default;
		AnalyzedPackage(std::vector<Float3>& avg_energy, std::vector<float> temperature);

		std::vector<float> pot_energy;
		std::vector<float> kin_energy;
		std::vector<float> total_energy;
		std::vector<Float3> energy_data; // potE, kinE, totalE

		float energy_gradient;
		float variance_coefficient;

		float mean_energy{};

		std::vector<float> temperature_data;
	};

	AnalyzedPackage analyzeEnergy(Simulation* simulation); // Prints a file of doubles: [step, molecule, atom, coordinate_dim]

	static void printEnergy(AnalyzedPackage* package);
	//static float getVarianceCoefficient(const std::vector<float>& total_energy);


	// Temp dev function
	static void findAndDumpPiecewiseEnergies(const Simulation& sim, const std::string& workdir);

	
	



private:
	std::vector<Float3> analyzeSolvateEnergy(Simulation* simulation, uint64_t n_steps);
	std::vector<Float3> analyzeCompoundEnergy(Simulation* simulation, uint64_t n_steps);

	float* potE_buffer_device = nullptr;
	float* vel_buffer_device = nullptr;

	std::unique_ptr<LimaLogger> m_logger;
};
