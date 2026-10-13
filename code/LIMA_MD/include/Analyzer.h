#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "Utilities.h"
#include "MDFiles.h"
#include <filesystem>
#include <functional>
#include <vector>
#include <string>

struct SimulationResult;

namespace SimAnalysis {
	struct AnalyzedPackage {
		AnalyzedPackage() = default;
		AnalyzedPackage(std::vector<Float3> avgEnergy, std::vector<float> temperature);

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
	void AnalyzeEnergy(SimulationResult& result);
	// Potential, kinetic and total energy [kJ/mol] and temperature [K] of each logged step, as an .xvg file like gmx energy writes
	void WriteEnergySeries(const Simulation& simulation, const std::filesystem::path& outputPath);
	void DensityProfile(const Simulation& simulation, const std::filesystem::path& outputPath, bool show = false);
	struct DensityProfileGroup {
		std::string composition;
		float temperature = 0.f;
		std::vector<std::filesystem::path> profiles;
	};
	void CompareDensityProfiles(const std::vector<DensityProfileGroup>& groups, const std::filesystem::path& outputPath, bool show = true);
	std::vector<Float3> GetForces(const Simulation& simulation, int64_t step);

	struct SelfDiffusion {
		double coefficient = 0.0;	// [m^2/s]
		double temperature = 0.0;	// [K] mean over the analyzed part of the run
	};
	// Self-diffusion of the selected atoms from their mean squared displacement, MSD(t) = 6Dt, ignoring the first
	// skipPs of the trajectory. No finite-size correction, so D is somewhat lower than in an infinite system
	SelfDiffusion SelfDiffusionCoefficient(const Simulation& simulation,
		const std::function<bool(const GroRecord&)>& selection, double skipPs);


	void PlotPotentialEnergyDistribution(const Simulation& sim, const std::filesystem::path& dir, const std::vector<int>& stepsToPlot);

	int CountOscillations(std::vector<float>& data);
}
