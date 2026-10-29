#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "SimulationData.h"
#include "CapacityOverflow.cuh"
#include "ForceAccumulator.cuh"
namespace ChargeBlock { struct ChargeblockBuffers; }

#include "platform/GpuFft.h"
#include "PmeFft.h"
#include <memory>

namespace PME {
	// Gridpoints per chargeblock and axis. A chargeblock is 1 nm in the grid's scaled coordinates, so this is per scaled nm.
	// The real spacing follows from how many chargeblocks tile the box, see ChargeblocksPerDim
	const int gridpointsPerNm = 8;
	constexpr float gridpointsPerNm_f = static_cast<float>(gridpointsPerNm);
	constexpr float invCellVolume = static_cast<float>(gridpointsPerNm * gridpointsPerNm * gridpointsPerNm);


	class Controller {
		Int3 gridpointsPerDim{};
		size_t nGridpointsRealspace = 0;
		int nGridpointsReciprocalspace = -1;
		const float ewaldKappa;
		const float gridSpacing;	// The largest grid spacing [nm]
		Float3 boxlenNm{};
		Int3 blocksPerDim{};	// See ChargeblocksPerDim, independent of the engine's cells
		Float3 gridScale{};		// Per axis, from nm to the grid's coordinates, where the chargeblocks are 1 by 1 by 1
		const int nChargeblocks;

		CudaBuffer<int> simulationSlots;
		std::vector<int> activeSimulationIds;
		int batchCount = 0;

		// FFT
		float* realspaceGrid = nullptr;
		int* fixedPointGrid = nullptr;		// The charges spread to the grid, in fixed point so their sums are deterministic, see SpreadChargesKernel
		cufftComplex* fourierspaceGrid = nullptr;
		float* greensFunctionScalars = nullptr;	// Without the z pass
		float* reciprocalVirialScalars = nullptr;

		// Chargeblocks 
		std::unique_ptr<ChargeBlock::ChargeblockBuffers> chargeblockBuffers;

		// With the z pass the plans transform the planes along x and y, else the whole grid. Without it, the Greens multiply could be
		// folded into the inverse FFT's load with a cuFFT LTO callback (cufftXtSetJITCallback), saving ApplyGreensFunctionKernel.
		// Not done: it JIT-links at plan creation, and complicates the static cuFFT linking on Linux
		cufftHandle planForward = 0;
		cufftHandle planInverse = 0;

		// PmeFft's z pass, see PmeFft.cuh, for grids whose z axis it takes
		bool zPass = false;
		PmeFft::Axis zAxis{};
		float2* zTwiddles = nullptr;
		float2* greensAxisTables = nullptr;
		void UploadGreensAxisTables();

		cudaStream_t& stream;
		// For system with a net charge, we apply to correction to each realspaceGridnode
		//LAL::optional<float> backgroundchargeCorrection;

	public:

		// fourierSpacing: the largest grid spacing [nm], 0 for the default, see SimParams::fourier_spacing
		Controller(const std::vector<EngineSimulationData>& simulations, float cutoffNM, float fourierSpacing, cudaStream_t& stream);
		void SetActiveSimulations(const std::vector<EngineSimulationData>& simulations);
		~Controller();

		// MD adds the forces to forceAcc, EM stores them in forceEnergy (pcluster layout)
		void CalcCharges(SuperCluster* scData, SuperClusterMeta* scMeta, int nSuperclusters, ForceEnergy* forceEnergy, ForceAccumulator forceAcc,
			long long* virialPartials = nullptr);
		int VirialPartialCount() const;
		bool RequiresRebuild(Float3 boxLength) const;
		void UpdateBoxLength(Float3 boxLength);

		// Charge-block capacity status, or nullptr while no simulation is active
		const CapacityOverflow* Overflow() const;

	private:
		//Just for debugging
		void PlotPotentialSlices();

	};
}
