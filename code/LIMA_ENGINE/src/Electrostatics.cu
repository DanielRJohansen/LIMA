#include "Electrostatics.cuh"
#include "PME.cuh"
#include "FMM.cuh"

Electrostatics::ShortRangeParameters Electrostatics::GetShortRangeParameters(const SimParams& params) {
	switch (params.longrange_method.value) {
	case LongRangeElectrostaticsMethod::PME:
		// Preserve PME's existing real-space policy even without a reciprocal
		// controller. Fresh disabled-ES boxes have their charges cleared by BoxBuilder;
		// reused boxes retain the historical screened short-range behavior.
		return PME::Controller::ShortRange(params.cutoff_nm);
	case LongRangeElectrostaticsMethod::FMM:
		return {};
	}
	throw std::invalid_argument("Invalid longrange_method");
}

std::unique_ptr<Electrostatics::Controller> Electrostatics::Create(
	const std::vector<EngineSimulationData>& simulations, const SimParams& params, cudaStream_t& stream) {
	switch (params.longrange_method.value) {
	case LongRangeElectrostaticsMethod::PME:
		return std::make_unique<PME::Controller>(simulations, params.cutoff_nm, stream);
	case LongRangeElectrostaticsMethod::FMM:
		return std::make_unique<FMM::Controller>(simulations, params, stream);
	}
	throw std::invalid_argument("Invalid longrange_method");
}
