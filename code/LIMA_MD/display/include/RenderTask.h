#pragma once

#include "Backbone.h"
#include "LipidModel.h"
#include "ParticleCategory.h"
#include "MDFiles.h"
#include "MoleculeHull.cuh"
#include "Simulation.cuh"

#include <array>
#include <memory>
#include <set>
#include <variant>
#include <vector>

using SimulationId = int;

struct BoxImage;

namespace Rendering {	
	struct MoleculeInfo {
		std::string name;
		size_t number = 0;
		size_t typeCount = 0;
		std::vector<int> atomIds;
	};

	std::vector<MoleculeInfo> GetMoleculeInfo(const BoxImage& boxImage);
	std::optional<MoleculeInfo> GetMoleculeInfo(const BoxImage& boxImage, int atomId);
	// Empty when the simulation has no box image
	std::vector<MoleculeInfo> GetMoleculeInfo(const Simulation& simulation);
	// One molecule per residue, for coordinate-only rendering
	std::vector<MoleculeInfo> GetMoleculeInfo(const GroFile& grofile);
	struct NoTask {};

	struct FreeTask{};

	enum class LipidPart : uint8_t { None, Body, Head };

	struct AtomRenderData {
		char atomLetter = ' ';
		float charge = 0.f;
		int groupId = -1;
		bool isSolvent = false;
		// Set from the lipid models when the task is prepared for display
		LipidPart lipidPart = LipidPart::None;
		LipidHeadgroup lipidHeadgroup = LipidHeadgroup::Anionic;
	};

	struct AtomRenderTask {	
		std::vector<Float3> positions;
		std::vector<AtomRenderData> atoms;
		std::vector<int> packedPositionIndices;
		Float3 boxSize{};
		SimStatus simStatus;
		BackboneChains backboneChains;
		LipidModels lipidModels;
		ParticleCategories particleCategories;
		std::set<int> highlightedAtoms;
		std::vector<MoleculeInfo> molecules;
		bool showSolvents = true;
		std::optional<ColoringMethod> representation; // Initial representation, otherwise the display keeps its current

		AtomRenderTask(const GroFile& grofile, bool showSolvents = true);
		AtomRenderTask(
			const std::vector<PersistentCluster>& pclusters,
			const std::vector<PersistentClusterMeta>& pcMeta,
			const BoxParams& boxparams,
			SimStatus simStatus = {},
			BackboneChains backboneChains = {},
			std::vector<MoleculeInfo> molecules = {},
			LipidModels lipidModels = {},
			ParticleCategories particleCategories = {});
	};

	// Covalent bonds guessed from distances in the task's current positions, for the licorice representation.
	// Bonds across the periodic boundary are left out, as a stick through the whole box would be wrong
	std::vector<std::array<int, 2>> FindBonds(const AtomRenderTask& task);

	// One entry of the display's right-click menu
	struct SelectionOption {
		std::string label;
		std::string count;	// Shown right-aligned
		MoleculeInfo selection;	// Empty atomIds clears the selection
	};

	// What can be selected starting from atomId: the atom, its element, its molecule, its kind of molecule, and
	// every category in the system. atomId -1 (empty space) gives only the system-wide options
	std::vector<SelectionOption> SelectionOptions(const AtomRenderTask& task, int atomId, bool hasSelection);

	struct SimulationTaskUpdate {
		const Float3* const positions = nullptr;
		const float* const forceMagnitudes = nullptr;
		SimStatus simStatus;
	};

	struct MoleculehullTask {	
		const MoleculeHullCollection& molCollection;
		Float3 boxSize{};
	};

	using Task = std::variant<NoTask, FreeTask, std::unique_ptr<AtomRenderTask>, std::unique_ptr<SimulationTaskUpdate>, std::unique_ptr<MoleculehullTask>>;
}
