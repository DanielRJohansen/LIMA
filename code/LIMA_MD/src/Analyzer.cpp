#include "Analyzer.h"

#include "PhysicsUtils.cuh"
#include "Printer.h"
#include "Statistics.h"

#include <algorithm>
#include <numeric>

std::vector<Float3> SimAnalysis::GetForces(const Simulation& simulation, int64_t step) {
	int atomCount = 0;
	for (const auto& metadata : simulation.box->persistentClustersMetadata)
		for (const int globalId : metadata.particleIdsGlobal)
			atomCount = (std::max)(atomCount, globalId + 1);

	std::vector<Float3> forces(atomCount); // [kJ/mol/nm]
	for (int clusterId = 0; clusterId < simulation.box->persistentClusters.size(); clusterId++) {
		for (int particleId = 0; particleId < PersistentCluster::maxParticles; particleId++) {
			const int globalId = simulation.box->persistentClustersMetadata[clusterId].particleIdsGlobal[particleId];
			if (globalId >= 0)
				forces[globalId] = simulation.forceBuffer->GetDatapointAtStep(clusterId, particleId, step) / KILO;
		}
	}
	return forces;
}

SimAnalysis::AnalyzedPackage SimAnalysis::analyzeEnergy(Simulation* simulation) {
	const auto& metadata = simulation->box->persistentClustersMetadata;
	if (metadata.empty())
		return {};

	const int64_t nEntries = LIMALOGSYSTEM::getMostRecentDataentryIndex(
		simulation->getStep(), simulation->simParams.data_logging_interval);
	if (nEntries < 2)
		return {};

	struct Particle {
		size_t bufferOffset;
		float mass;
	};
	std::vector<Particle> particles;
	particles.reserve(metadata.size() * PersistentCluster::maxParticles);
	for (size_t clusterId = 0; clusterId < metadata.size(); clusterId++) {
		for (size_t particleId = 0; particleId < PersistentCluster::maxParticles; particleId++) {
			if (metadata[clusterId].particleIdsGlobal[particleId] >= 0) {
				particles.push_back({
					clusterId * PersistentCluster::maxParticles + particleId,
					metadata[clusterId].mass[particleId] });
			}
		}
	}
	if (particles.empty())
		return {};

	const size_t valuesPerEntry = metadata.size() * PersistentCluster::maxParticles;
	const auto& potentialEnergies = simulation->potE_buffer->data();
	const auto& velocities = simulation->vel_buffer->data();
	std::vector<Float3> averageEnergies(nEntries);
	for (int64_t entry = 0; entry < nEntries; entry++) {
		const size_t entryOffset = entry * valuesPerEntry;
		double potentialSum = 0.;
		double kineticSum = 0.;
		double totalSum = 0.;
		for (const auto& particle : particles) {
			const size_t index = entryOffset + particle.bufferOffset;
			const float potentialEnergy = potentialEnergies[index];
			const float kineticEnergy = PhysicsUtils::calcKineticEnergy(velocities[index], particle.mass);
			potentialSum += potentialEnergy;
			kineticSum += kineticEnergy;
			totalSum += potentialEnergy + kineticEnergy;
		}
		const double particleCount = static_cast<double>(particles.size());
		averageEnergies[entry] = Float3{
			potentialSum / particleCount,
			kineticSum / particleCount,
			totalSum / particleCount };
	}

	return AnalyzedPackage(std::move(averageEnergies), simulation->temperature_buffer);
}

namespace {
	float GetVarianceCoefficient(const std::vector<float>& values) {
		if (values.empty())
			return 0.f;
		const float standardDeviation = Statistics::StdDev(values);
		const float mean = Statistics::Mean(values);
		return standardDeviation == 0.f && mean == 0.f ? 0.f : standardDeviation / std::abs(mean);
	}

	void PrintRow(const std::string& title, const std::vector<float>& values) {
		if (values.empty())
			return;
		LIMA_Printer::printTableRow(title, {
			*std::ranges::min_element(values), *std::ranges::max_element(values),
			Statistics::StdDev(values), (values.back() - values.front()) / values.front() });
	}

	float CalculateSlopeLinearRegression(const std::vector<float>& values, float mean) {
		const size_t count = values.size();
		float sumX = 0.f;
		float sumY = 0.f;
		float sumXY = 0.f;
		float sumXX = 0.f;
		for (size_t index = 0; index < count; index++) {
			sumX += index;
			sumY += values[index];
			sumXY += index * values[index];
			sumXX += index * index;
		}
		const float slope = (count * sumXY - sumX * sumY) / (count * sumXX - sumX * sumX);
		return slope / mean;
	}
}

void SimAnalysis::AnalyzedPackage::Print() const {
	LIMA_Printer::printTableRow({ "", "min", "max", "Std. deviation", "Change 0->n" });
	PrintRow("potE", pot_energy);
	PrintRow("kinE", kin_energy);
	PrintRow("totalE", total_energy);
}

SimAnalysis::AnalyzedPackage::AnalyzedPackage(
	std::vector<Float3> avgEnergy, std::vector<float> temperature)
	: energy_data(std::move(avgEnergy)), temperature_data(std::move(temperature))
{
	pot_energy.reserve(energy_data.size());
	kin_energy.reserve(energy_data.size());
	total_energy.reserve(energy_data.size());
	for (const Float3& energy : energy_data) {
		pot_energy.push_back(energy.x);
		kin_energy.push_back(energy.y);
		total_energy.push_back(energy.z);
	}

	mean_energy = Statistics::Mean(total_energy);
	energy_gradient = CalculateSlopeLinearRegression(total_energy, mean_energy);
	variance_coefficient = GetVarianceCoefficient(total_energy);
}

void SimAnalysis::PlotPotentialEnergyDistribution(
	const Simulation&, const std::filesystem::path&, const std::vector<int>&) {}

int SimAnalysis::CountOscillations(std::vector<float>& data) {
	if (data.size() < 2)
		return 0;
	const float mean = std::accumulate(data.begin(), data.end(), 0.f) / data.size();
	int count = 0;
	bool aboveMean = data.front() > mean;
	for (size_t index = 1; index < data.size(); index++) {
		const bool currentAboveMean = data[index] > mean;
		if (!aboveMean && currentAboveMean)
			count++;
		aboveMean = currentAboveMean;
	}
	return count;
}
