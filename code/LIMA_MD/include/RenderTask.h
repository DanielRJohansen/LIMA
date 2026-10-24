#pragma once

#include "Backbone.h"
#include "LipidModel.h"
#include "MDFiles.h"
#include "MoleculeHull.cuh"
#include "Simulation.cuh"

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
			LipidModels lipidModels = {});
	};

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
