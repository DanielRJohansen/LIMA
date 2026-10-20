#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "SimulationData.h"
#include "CapacityOverflow.cuh"
#include "ForceAccumulator.cuh"
namespace ChargeBlock { struct ChargeblockBuffers; }

#include "platform/GpuFft.h"
#include <memory>

namespace PME {
	// 0.125 nm grid spacing. With 4th order B-splines this gives ~0.3% RMS reciprocal force error on a solvated protein,
	// measured against an exact Ewald sum (0.14% at 0.1 nm)
	const int gridpointsPerNm = 8;
	constexpr float gridpointsPerNm_f = static_cast<float>(gridpointsPerNm);
	constexpr float invCellVolume = static_cast<float>(gridpointsPerNm * gridpointsPerNm * gridpointsPerNm);


	class Controller {
		Int3 gridpointsPerDim{};
		size_t nGridpointsRealspace = 0;
		int nGridpointsReciprocalspace = -1;
		const float ewaldKappa;
		Float3 boxlenNm{};
		Int3 blocksPerDim{};	// The box length rounded up, independent of the engine's cells
		Float3 gridScale{};		// Per axis, from nm to the grid's coordinates, where the chargeblocks are 1 by 1 by 1
		const int nChargeblocks;

		CudaBuffer<int> simulationSlots;
		std::vector<int> activeSimulationIds;
		int batchCount = 0;

		// FFT
		float* realspaceGrid = nullptr;
		cufftComplex* fourierspaceGrid = nullptr;
		float* greensFunctionScalars;
		float* reciprocalVirialScalars;

		// Chargeblocks 
		std::unique_ptr<ChargeBlock::ChargeblockBuffers> chargeblockBuffers;

		// The Greens multiply could be folded into the inverse FFT's load with a cuFFT LTO callback (cufftXtSetJITCallback),
		// saving ApplyGreensFunctionKernel (~45 us/step on STMV). Not done: it JIT-links at plan creation, and complicates the
		// static cuFFT linking on Linux
		cufftHandle planForward = 0;
		cufftHandle planInverse = 0;

		cudaStream_t& stream;
		// For system with a net charge, we apply to correction to each realspaceGridnode
		//LAL::optional<float> backgroundchargeCorrection;

	public:

		Controller(const std::vector<EngineSimulationData>& simulations, float cutoffNM, cudaStream_t& stream);
		void SetActiveSimulations(const std::vector<EngineSimulationData>& simulations);
		~Controller();

		// MD adds the forces to forceAcc, EM stores them in forceEnergy (pcluster layout)
		void CalcCharges(SuperCluster* scData, SuperClusterMeta* scMeta, int nSuperclusters, ForceEnergy* forceEnergy, ForceAccumulator forceAcc,
			long long* virialPartials = nullptr);
		int VirialPartialCount() const { return (nGridpointsReciprocalspace * batchCount + 63) / 64; }
		bool RequiresRebuild(Float3 boxLength) const;
		void UpdateBoxLength(Float3 boxLength);

		// Charge-block capacity status, or nullptr while no simulation is active
		const CapacityOverflow* Overflow() const;

	private:
		//Just for debugging
		void PlotPotentialSlices();

	};
}
