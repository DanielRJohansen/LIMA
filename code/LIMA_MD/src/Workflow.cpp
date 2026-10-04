#include "Workflow.h"
#include "Format.h"

#include "PhysicsUtils.cuh"

#include <random>

namespace {
	void InitializeVelocities(Simulation& simulation, float temperature, uint32_t seed) {
		if (temperature <= 0.f)
			throw std::invalid_argument("Initial temperature must be positive");

		std::mt19937 generator(seed);
		std::normal_distribution<float> normal;
		Float3 totalMomentum{};
		float totalMass = 0.f;
		for (size_t clusterId = 0; clusterId < simulation.box->persistentClustersMetadata.size(); ++clusterId) {
			const PersistentClusterMeta& metadata = simulation.box->persistentClustersMetadata[clusterId];
			PersistentclusterInterimState& state = simulation.box->pclusterInterimStates[clusterId];
			for (int particleId = 0; particleId < PersistentCluster::maxParticles; ++particleId) {
				const float mass = metadata.mass[particleId];
				if (metadata.particleIdsGlobal[particleId] < 0 || mass <= 0.f) continue;
				const float componentDeviation = PhysicsUtils::tempToVelocity(temperature, mass) / std::sqrt(3.f);
				state.vels_prev[particleId] = Float3{
					normal(generator), normal(generator), normal(generator) } * componentDeviation;
				totalMomentum += state.vels_prev[particleId] * mass;
				totalMass += mass;
			}
		}

		if (totalMass <= 0.f)
			throw std::runtime_error("Cannot initialize velocities for an empty system");
		const Float3 centerOfMassVelocity = totalMomentum / totalMass;
		double kineticEnergy = 0.;
		for (size_t clusterId = 0; clusterId < simulation.box->persistentClustersMetadata.size(); ++clusterId) {
			const PersistentClusterMeta& metadata = simulation.box->persistentClustersMetadata[clusterId];
			PersistentclusterInterimState& state = simulation.box->pclusterInterimStates[clusterId];
			for (int particleId = 0; particleId < PersistentCluster::maxParticles; ++particleId) {
				const float mass = metadata.mass[particleId];
				if (metadata.particleIdsGlobal[particleId] < 0 || mass <= 0.f) continue;
				state.vels_prev[particleId] -= centerOfMassVelocity;
				kineticEnergy += PhysicsUtils::calcKineticEnergy(
					state.vels_prev[particleId].len(), mass);
			}
		}
		const float actualTemperature = PhysicsUtils::kineticEnergyToTemperature(
			kineticEnergy, simulation.box->boxparams.degreesOfFreedom);
		if (actualTemperature <= 0.f)
			throw std::runtime_error("Failed to generate nonzero initial velocities");
		const float scale = std::sqrt(temperature / actualTemperature);
		for (PersistentclusterInterimState& state : simulation.box->pclusterInterimStates)
			for (Float3& velocity : state.vels_prev) velocity *= scale;
		simulation.temperature_buffer.push_back(temperature);
	}
}

std::vector<Programs::WorkflowInput> Programs::MakeMembraneInputs(
	const std::vector<Lipids::Selection>& compositions, const std::vector<int>& seeds,
	Float3 boxSize, MembraneGeometry::Figure geometry, bool solvate) {
	std::vector<WorkflowInput> inputs;
	inputs.reserve(compositions.size() * seeds.size());
	for (const Lipids::Selection& composition : compositions) {
		const std::string compositionName = Lipids::NameSelection(composition);
		for (const int seed : seeds) {
			const std::string name = compositionName + Lima::Format("_seed{}", seed);
			inputs.push_back({
				.name = name,
				.tags = { { "composition", compositionName }, { "seed", std::to_string(seed) } },
				.MakeJob = [composition, boxSize, geometry, seed, solvate](
					fs::path runDir, SimParams params, EnvMode mode) {
					return MakeMembraneJob(std::move(runDir), composition, boxSize, geometry,
						seed, std::move(params), mode, solvate);
				}
			});
		}
	}
	return inputs;
}

Programs::SimulationWorkflow::SimulationWorkflow(fs::path workDir, EnvMode mode)
	: workDir(std::move(workDir)), mode(mode) {}

void Programs::SimulationWorkflow::AddInputs(std::vector<WorkflowInput> newInputs) {
	inputs.insert(inputs.end(), std::make_move_iterator(newInputs.begin()),
		std::make_move_iterator(newInputs.end()));
}

void Programs::SimulationWorkflow::AddStage(WorkflowStage stage) {
	stages.push_back(std::move(stage));
}

void Programs::SimulationWorkflow::CompareDensityProfiles(std::string compositionTag,
	std::string temperatureTag, fs::path output) {
	densityProfileComparison.emplace(
		std::move(compositionTag), std::move(temperatureTag), std::move(output));
}

void Programs::SimulationWorkflow::Run(Environment& environment) {
	if (inputs.empty()) throw std::invalid_argument("A simulation workflow requires at least one input");
	if (stages.empty()) throw std::invalid_argument("A simulation workflow requires at least one stage");

	struct PendingRun {
		std::string name;
		std::map<std::string, std::string> tags;
		fs::path runDir;
		SimulationHandle handle;
	};
	std::vector<PendingRun> pending;
	const WorkflowStage& firstStage = stages.front();
	if (!firstStage.variants.empty())
		throw std::invalid_argument("The first workflow stage cannot have variants");
	pending.reserve(inputs.size());
	for (const WorkflowInput& input : inputs) {
		const fs::path runDir = workDir / input.name;
		SimulationJob job = input.MakeJob(runDir, firstStage.params, mode);
		job.name = input.name;
		job.outputs = firstStage.outputs;
		pending.push_back({ input.name, input.tags, runDir, environment.Submit(std::move(job)) });
	}

	for (size_t stageId = 1; stageId < stages.size(); ++stageId) {
		const WorkflowStage& stage = stages[stageId];
		std::vector<WorkflowVariant> variants = stage.variants;
		if (variants.empty()) variants.push_back({ .name = stage.name });
		std::vector<PendingRun> next;
		next.reserve(pending.size() * variants.size());
		for (PendingRun& parent : pending) {
			SimulationResult completed = parent.handle.Get();
			MolecularSystem parentSystem = completed.FinalSystem();
			for (size_t variantId = 0; variantId < variants.size(); ++variantId) {
				const WorkflowVariant& variant = variants[variantId];
				SimParams params = stage.params;
				if (variant.Configure) variant.Configure(params);
				MolecularSystem system = variantId + 1 == variants.size()
					? std::move(parentSystem) : parentSystem;
				const std::string runName = variant.name.empty() ? stage.name : variant.name;
				const fs::path runDir = parent.runDir / runName;
				SimulationJob job = MakeSimulationJob(runDir, std::move(system), std::move(params), mode);
				job.name = parent.name + " " + runName;
				job.outputs = stage.outputs;
				if (stage.initializeAtReferenceTemperature) {
					const uint32_t velocitySeed = static_cast<uint32_t>(std::hash<std::string>{}(job.name));
					job.configureSimulation = [velocitySeed](Simulation& simulation) {
						InitializeVelocities(simulation, simulation.simParams.ref_t, velocitySeed);
					};
				}
				std::map<std::string, std::string> tags = parent.tags;
				for (const auto& [key, value] : variant.tags)
					tags.insert_or_assign(key, value);
				next.push_back({ job.name, std::move(tags), runDir,
					environment.Submit(std::move(job)) });
			}
		}
		pending = std::move(next);
	}

	for (PendingRun& run : pending) run.handle.Get();
	if (!densityProfileComparison) return;
	if (!stages.back().outputs.contains(OutputSelect::DensityProfile))
		throw std::invalid_argument("Density-profile comparison requires density-profile output");

	std::map<std::pair<std::string, std::string>, std::vector<fs::path>> groupedProfiles;
	for (const PendingRun& run : pending) {
		const auto composition = run.tags.find(densityProfileComparison->compositionTag);
		const auto temperature = run.tags.find(densityProfileComparison->temperatureTag);
		if (composition == run.tags.end() || temperature == run.tags.end())
			throw std::invalid_argument("Density-profile comparison references a missing workflow tag");
		groupedProfiles[{ composition->second, temperature->second }].push_back(
			run.runDir / "density_profile.csv");
	}
	std::vector<SimAnalysis::DensityProfileGroup> groups;
	groups.reserve(groupedProfiles.size());
	for (auto& [key, profiles] : groupedProfiles)
		groups.push_back({ key.first, std::stof(key.second), std::move(profiles) });
	SimAnalysis::CompareDensityProfiles(groups, workDir / densityProfileComparison->output);
}
