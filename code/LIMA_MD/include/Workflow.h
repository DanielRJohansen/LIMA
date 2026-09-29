#pragma once

#include "Programs.h"

#include <map>

namespace Programs {
	struct WorkflowInput {
		std::string name;
		std::map<std::string, std::string> tags;
		std::function<SimulationJob(fs::path, SimParams, EnvMode)> MakeJob;
	};

	struct WorkflowVariant {
		std::string name;
		std::map<std::string, std::string> tags;
		std::function<void(SimParams&)> Configure;
	};

	struct WorkflowStage {
		std::string name;
		SimParams params;
		std::vector<WorkflowVariant> variants;
		std::set<OutputSelect> outputs;
		bool initializeAtReferenceTemperature = false;
	};

	std::vector<WorkflowInput> MakeMembraneInputs(const std::vector<Lipids::Selection>& compositions,
		const std::vector<int>& seeds, Float3 boxSize, MembraneGeometry::Figure geometry,
		bool solvate = false);

	class SimulationWorkflow {
	public:
		SimulationWorkflow(fs::path workDir, EnvMode mode = EnvMode::Headless);
		void AddInputs(std::vector<WorkflowInput> newInputs);
		void AddStage(WorkflowStage stage);
		void CompareDensityProfiles(std::string compositionTag, std::string temperatureTag,
			fs::path output = "density_profile_comparison.csv");
		void Run(Environment& environment = Environment::Get());

	private:
		struct DensityProfileComparison {
			std::string compositionTag;
			std::string temperatureTag;
			fs::path output;
		};

		fs::path workDir;
		EnvMode mode;
		std::vector<WorkflowInput> inputs;
		std::vector<WorkflowStage> stages;
		std::optional<DensityProfileComparison> densityProfileComparison;
	};
}
