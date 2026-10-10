#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "SimulationData.h"
#include "CapacityOverflow.cuh"
#include "ForceAccumulator.cuh"
namespace ChargeBlock { struct ChargeblockBuffers; }

#include <cufft.h>
#include <memory>

// TODO: Do i need to account for e0, vacuum/spaceial permitivity here? Probably....


namespace PME {
	const int gridpointsPerNm = 10;
	constexpr float gridpointsPerNm_f = static_cast<float>(gridpointsPerNm);
	constexpr float invCellVolume = static_cast<float>(gridpointsPerNm * gridpointsPerNm * gridpointsPerNm);


	class Controller {
		Int3 gridpointsPerDim{};
		size_t nGridpointsRealspace = 0;
		int nGridpointsReciprocalspace = -1;	// TODO: Does this also need to be size_t?
		const float ewaldKappa;
		Float3 boxlenNm{};
		const int nChargeblocks;

		// Always applied constant per particle
		CudaBuffer<float> selfenergyCorrections;
		CudaBuffer<int> simulationSlots;
		std::vector<int> activeSimulationIds;
		int batchCount = 0;

		// FFT
		float* realspaceGrid = nullptr;
		cufftComplex* fourierspaceGrid = nullptr;
		float* greensFunctionScalars;

		// Chargeblocks 
		std::unique_ptr<ChargeBlock::ChargeblockBuffers> chargeblockBuffers;

		// The Greens multiply could be folded into the inverse FFT's load with a cuFFT LTO callback (cufftXtSetJITCallback),
		// saving ApplyGreensFunctionKernel (~60 us/step on STMV). Not done: it JIT-links at plan creation, and complicates the
		// static cuFFT linking on Linux
		cufftHandle planForward = 0;
		cufftHandle planInverse = 0;

		cudaStream_t& stream;
		// For system with a net charge, we apply to correction to each realspaceGridnode
		//LAL::optional<float> backgroundchargeCorrection;

		static float CalcEnergyCorrection(const Box& box, float ewaldKappa);

	public:

		Controller(const std::vector<EngineSimulationData>& simulations, float cutoffNM, cudaStream_t& stream);
		void SetActiveSimulations(const std::vector<EngineSimulationData>& simulations);
		~Controller();

		// MD adds the forces to forceAcc, EM stores them in forceEnergy (pcluster layout)
		void CalcCharges(SuperCluster* scData, SuperClusterMeta* scMeta, int nSuperclusters, ForceEnergy* forceEnergy, ForceAccumulator forceAcc);

		// Charge-block capacity status, or nullptr while no simulation is active
		const CapacityOverflow* Overflow() const;

	private:
		//Just for debugging
		void PlotPotentialSlices();

	};
}
