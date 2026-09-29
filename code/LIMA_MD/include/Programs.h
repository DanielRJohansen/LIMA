#pragma once

#include "SimulationBuilder.h"
#include "Environment.h"

#include <map>
#include <string_view>

class MoleculeHullCollection;

namespace Programs {
	enum class WaterModel { Tip3p, Tip4p, Tips3p, Tip5p, Spc, Spce };

	struct GmxConversionResult {
		GroFile grofile;
		TopologyFile topology;
		// Ordered like topology.GetSystem().molecules. Each file contains the
		// position restraints for the corresponding molecule.
		std::vector<GenericItpFile> positionRestraints;
	};

	WaterModel ParseWaterModel(std::string_view name);

	void GetForcefieldParams(const GroFile&, const TopologyFile&, const fs::path& workdir);

	MoleculeHullCollection MakeLipidVesicle(GroFile&, TopologyFile&, Lipids::Selection, float vesicleRadius, 
		Float3 vesicleCenter, std::optional<int> numLipids=std::nullopt);

	void MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize, bool renderProgress);

	void StaticbodyEnergyMinimize(GroFile&, const TopologyFile&, bool render);

	SimulationJob MakeMembraneJob(fs::path workDir, Lipids::Selection composition,
		Float3 boxSize, MembraneGeometry::Figure geometry, int seed,
		SimParams params, EnvMode mode = EnvMode::Headless);

	SimulationJob MakeSimulationJob(fs::path workDir, MolecularSystem system,
		SimParams params, EnvMode mode = EnvMode::Headless);

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
	};

	std::vector<WorkflowInput> MakeMembraneInputs(const std::vector<Lipids::Selection>& compositions,
		const std::vector<int>& seeds, Float3 boxSize, MembraneGeometry::Figure geometry);

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

	/// Build in-memory CHARMM27 coordinates, topology, and position restraints
	/// from a protein PDB or mmCIF structure.
	GmxConversionResult ToGmx(const fs::path& structureFile, WaterModel waterModel = WaterModel::Tip3p);
}
