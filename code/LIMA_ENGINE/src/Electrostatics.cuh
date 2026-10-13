#pragma once

#include "ForceAccumulator.cuh"
#include <memory>
#include <vector>

struct EngineSimulationData;
struct SimParams;
struct CapacityOverflow;

namespace Electrostatics {
	struct ShortRangeParameters {
		bool enabled = false;
		float ewaldKappa = 0.f;
	};

	// Both backends consume the current supercluster layout and write the same MD/EM outputs.
	// Backend allocations and numerical details never escape this interface.
	class Controller {
	public:
		virtual ~Controller() = default;
		virtual void CalcCharges(SuperCluster* scData, SuperClusterMeta* scMeta, int nSuperclusters,
			ForceEnergy* forceEnergy, ForceAccumulator forceAcc) = 0;
		virtual const CapacityOverflow* Overflow() const { return nullptr; }
	};

	ShortRangeParameters GetShortRangeParameters(const SimParams& params);
	std::unique_ptr<Controller> Create(const std::vector<EngineSimulationData>& simulations,
		const SimParams& params, cudaStream_t& stream);
}
