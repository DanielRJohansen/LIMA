#pragma once

#include "Backbone.h"
#include "MDFiles.h"
#include "MoleculeHull.cuh"
#include "Simulation.cuh"

#include <memory>
#include <set>
#include <variant>
#include <vector>

namespace Rendering {
	struct NoTask {};

	struct AtomRenderData {
		char atomLetter = ' ';
		float charge = 0.f;
		int groupId = -1;
		bool isSolvent = false;
	};

	struct AtomRenderTask {
		std::vector<Float3> positions;
		std::vector<AtomRenderData> atoms;
		std::vector<int> packedPositionIndices;
		Float3 boxSize{};
		SimStatus simStatus;
		BackboneChains backboneChains;
		std::set<int> highlightedAtoms;
		bool showSolvents = true;

		AtomRenderTask(const GroFile& grofile, bool showSolvents = true);
		AtomRenderTask(
			const std::vector<PersistentCluster>& pclusters,
			const std::vector<PersistentClusterMeta>& pcMeta,
			const BoxParams& boxparams,
			SimStatus simStatus = {},
			BackboneChains backboneChains = {});
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

	using Task = std::variant<NoTask, std::unique_ptr<AtomRenderTask>, std::unique_ptr<SimulationTaskUpdate>, std::unique_ptr<MoleculehullTask>>;
}
