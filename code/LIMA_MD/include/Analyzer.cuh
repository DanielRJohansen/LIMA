#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "Utilities.h"
#include <filesystem>
#include <vector>
#include <string>


namespace SimAnalysis {
	struct AnalyzedPackage {
		AnalyzedPackage() = default;
		AnalyzedPackage(std::vector<Float3>& avg_energy, std::vector<float> temperature);

		void Print() const;

		std::vector<float> pot_energy;
		std::vector<float> kin_energy;
		std::vector<float> total_energy;
		std::vector<Float3> energy_data; // potE, kinE, totalE

		float energy_gradient = 0.f;
		float variance_coefficient = 0.f;

		float mean_energy{};

		std::vector<float> temperature_data;
	};

	AnalyzedPackage analyzeEnergy(Simulation* simulation); // Prints a file of doubles: [step, molecule, atom, coordinate_dim]


	void PlotPotentialEnergyDistribution(const Simulation& sim, const std::filesystem::path& dir, const std::vector<int>& stepsToPlot);

	int CountOscillations(std::vector<float>& data);
}