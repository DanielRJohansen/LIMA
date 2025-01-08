//#pragma once // Only allowed to be included by engine.cu
#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "SimulationDevice.cuh"

#include <cufft.h>
#include <memory>

namespace ChargeBlock {	struct ChargeblockBuffers; }
// TODO: Do i need to account for e0, vacuum/spaceial permitivity here? Probably....



namespace PME {
	const int gridpointsPerNm = 10;
	constexpr float gridpointsPerNm_f = static_cast<float>(gridpointsPerNm);
	constexpr float invCellVolume = static_cast<float>(gridpointsPerNm * gridpointsPerNm * gridpointsPerNm);


	class Controller {
		int gridpointsPerDim = -1;
		size_t nGridpointsRealspace = 0;
		int nGridpointsReciprocalspace = -1;
		const float ewaldKappa;
		float boxlenNm{};
		const int nChargeblocks;

		// Always applied constant per particle
		float selfenergyCorrection;

		// FFT
		float* realspaceGrid;
		cufftComplex* fourierspaceGrid;
		float* greensFunctionScalars;

		// Chargeblocks 
		std::unique_ptr<ChargeBlock::ChargeblockBuffers> chargeblockBuffers;

		cufftHandle planForward;
		cufftHandle planInverse;

		// For system with a net charge, we apply to correction to each realspaceGridnode
		//LAL::optional<float> backgroundchargeCorrection;

		void CalcEnergyCorrection(const Box& box);

	public:

		Controller(int boxLenNm, const Box& box, float cutoffNM);
		~Controller();

		void CalcCharges(const BoxConfig& config, const BoxState& state, int nCompounds, ForceEnergy* const forceEnergy);

	private:
		//Just for debugging
		void PlotPotentialSlices();

	};
}