#include "Analyzer.h"

#include "Environment.h"
#include "PhysicsUtils.cuh"
#include "Printer.h"
#include "Statistics.h"
#include "BoxImageBuilder.h"
#include "Filehandling.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <numeric>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

namespace {
	void LaunchDetachedPython(const std::filesystem::path& script,
		const std::filesystem::path& input, bool show) {
#ifdef _WIN32
		std::wstring command = L"python \"" + script.wstring() + L"\" --comparison \""
			+ input.wstring() + L"\"" + (show ? L" --show" : L"");
		STARTUPINFOW startupInfo{};
		startupInfo.cb = sizeof(startupInfo);
		PROCESS_INFORMATION processInfo{};
		if (!CreateProcessW(nullptr, command.data(), nullptr, nullptr, FALSE,
			DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP, nullptr, nullptr, &startupInfo, &processInfo))
			throw std::runtime_error("Failed to launch density-profile comparison");
		CloseHandle(processInfo.hThread);
		CloseHandle(processInfo.hProcess);
#else
		std::string command = std::format("python \"{}\" --comparison \"{}\"{} >/dev/null 2>&1 &",
			script.string(), input.string(), show ? " --show" : "");
		if (std::system(command.c_str()) != 0)
			throw std::runtime_error("Failed to launch density-profile comparison");
#endif
	}
}

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

void SimAnalysis::AnalyzeEnergy(SimulationResult& result) {
	result.analysis = analyzeEnergy(result.simulation.get());
}

void SimAnalysis::DensityProfile(const Simulation& simulation, const std::filesystem::path& outputPath, bool show) {
	if (!simulation.traj_buffer || !simulation.boxImage)
		throw std::invalid_argument("Density profile requires trajectory data and atom metadata");
	const int loggingInterval = simulation.simParams.data_logging_interval;
	if (loggingInterval <= 0)
		throw std::invalid_argument("Density profile requires trajectory logging");
	const int64_t nFrames = LIMALOGSYSTEM::getMostRecentDataentryIndex(simulation.getStep(), loggingInterval);
	if (nFrames <= 0)
		throw std::invalid_argument("Density profile requires at least one logged trajectory frame");

	const Float3 boxSize = simulation.box->boxparams.BoxSizeFloat();
	if (boxSize.x <= 0.f || boxSize.y <= 0.f || boxSize.z <= 0.f)
		throw std::invalid_argument("Density profile requires a non-empty box");

	constexpr int nBins = 80;
	std::array<std::vector<float>, 3> densities;
	for (auto& density : densities) density.assign(nBins, 0.f);
	const float binWidth = boxSize.z / nBins;
	const auto Group = [&simulation](const PersistentClusterMeta& metadata, int particleId) -> std::optional<size_t> {
		if (metadata.isSolvent) return 0;
		const int globalId = metadata.particleIdsGlobal[particleId];
		if (globalId < 0 || globalId >= simulation.boxImage->grofile.atoms.size()) return std::nullopt;
		const auto& atom = simulation.boxImage->grofile.atoms[globalId];
		if (atom.atomName == "P") return 1;
		if (atom.atomName.View().starts_with('C')) return 2;
		return std::nullopt;
	};

	for (int64_t frame = 0; frame < nFrames; ++frame) {
		for (size_t clusterId = 0; clusterId < simulation.box->persistentClustersMetadata.size(); ++clusterId) {
			const auto& metadata = simulation.box->persistentClustersMetadata[clusterId];
			for (int particleId = 0; particleId < PersistentCluster::maxParticles; ++particleId) {
				const auto group = Group(metadata, particleId);
				if (!group) continue;
				const Float3 position = simulation.traj_buffer->GetDatapoint(static_cast<int>(clusterId), particleId, frame);
				const float z = position.z - std::floor(position.z / boxSize.z) * boxSize.z;
				const int bin = std::clamp(static_cast<int>(z / binWidth), 0, nBins - 1);
				densities[*group][bin] += 1.f;
			}
		}
	}

	const float scale = 1.f / (nFrames * binWidth * boxSize.x * boxSize.y);
	for (auto& density : densities)
		for (float& value : density)
			value *= scale;

	if (!outputPath.parent_path().empty())
		std::filesystem::create_directories(outputPath.parent_path());
	std::ofstream file(outputPath);
	if (!file.is_open())
		throw std::runtime_error(std::format("Failed to write density profile {}", outputPath.string()));
	file << "z_nm,water_density,head_density,tail_density\n";
	for (int bin = 0; bin < nBins; ++bin)
		file << (bin + .5f) * binWidth << ',' << densities[0][bin] << ',' << densities[1][bin] << ',' << densities[2][bin] << '\n';
	file.close();

	const auto script = FileUtils::GetLimaDir() / "dev" / "PyTools" / "DensityProfile.py";
	std::string command = std::format("python \"{}\" \"{}\"", script.string(), outputPath.string());
	if (show) command += " --show";
	if (std::system(command.c_str()) != 0)
		throw std::runtime_error("Matplotlib failed to render density profile");
}

void SimAnalysis::CompareDensityProfiles(const std::vector<DensityProfileGroup>& groups, const std::filesystem::path& outputPath, bool show) {
	if (groups.empty())
		throw std::invalid_argument("Density profile comparison requires at least one group");
	if (!outputPath.parent_path().empty())
		std::filesystem::create_directories(outputPath.parent_path());
	std::ofstream output(outputPath);
	if (!output.is_open())
		throw std::runtime_error(std::format("Failed to write density profile comparison {}", outputPath.string()));
	output << "composition,temperature,z_nm,water_density,head_density,tail_density\n";

	for (const auto& group : groups) {
		if (group.profiles.empty())
			throw std::invalid_argument("Density profile comparison group has no profiles");
		std::vector<std::array<float, 4>> averages;
		for (const auto& profilePath : group.profiles) {
			std::ifstream input(profilePath);
			if (!input.is_open())
				throw std::runtime_error(std::format("Failed to read density profile {}", profilePath.string()));
			std::string line;
			std::getline(input, line);
			for (size_t bin = 0; std::getline(input, line); ++bin) {
				std::array<float, 4> values{};
				std::istringstream row(line);
				char comma;
				if (!(row >> values[0] >> comma >> values[1] >> comma >> values[2] >> comma >> values[3]))
					throw std::runtime_error(std::format("Invalid density profile row in {}", profilePath.string()));
				if (averages.size() <= bin) averages.emplace_back();
				for (size_t column = 0; column < values.size(); ++column)
					averages[bin][column] += values[column];
			}
		}
		for (auto& average : averages) {
			for (size_t column = 0; column < average.size(); ++column)
				average[column] /= static_cast<float>(group.profiles.size());
			output << group.composition << ',' << group.temperature;
			for (const float value : average) output << ',' << value;
			output << '\n';
		}
	}
	output.close();

	const auto script = FileUtils::GetLimaDir() / "dev" / "PyTools" / "DensityProfile.py";
	LaunchDetachedPython(script, outputPath, show);
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
