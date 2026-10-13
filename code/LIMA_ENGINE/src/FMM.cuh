#pragma once

#include "Electrostatics.cuh"

namespace FMM {
	class Controller final : public Electrostatics::Controller {
		struct Data;
		std::unique_ptr<Data> data;
	public:
		Controller(const std::vector<EngineSimulationData>& simulations, const SimParams& params, cudaStream_t stream);
		~Controller() override;
		void CalcCharges(SuperCluster* scData, SuperClusterMeta* scMeta, int nSuperclusters,
			ForceEnergy* forceEnergy, ForceAccumulator forceAcc) override;
	};
}
