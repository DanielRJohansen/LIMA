//#pragma once // Only allowed to be included by engine.cu
#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "BoundaryCondition.cuh"
#include "BoxGrid.cuh"
#include "SimulationDevice.cuh"
#include "PhysicsUtils.cuh"

#include "Constants.h"
#include "Filehandling.h"
#include "Utilities.h"

#include <cuda_runtime.h>
#include <cufft.h>


struct ChargeBlock;

// TODO: Do i need to account for e0, vacuum/spaceial permitivity here? Probably....

namespace PME {
	const int gridpointsPerNm = 10;
	constexpr float gridpointsPerNm_f = static_cast<float>(gridpointsPerNm);

	class Controller {
		// FFT
		float* realspaceGrid;
		cufftComplex* fourierspaceGrid;
		float* greensFunctionScalars;

		// 
		ChargeBlock* chargeBlocks;

		int gridpointsPerDim = -1;
		size_t nGridpointsRealspace = 0;
		int nGridpointsReciprocalspace = -1;
		const float ewaldKappa;
		float boxlenNm{};

		// Always applied constant per particle
		float selfenergyCorrection;

		cufftHandle planForward;
		cufftHandle planInverse;

		// For system with a net charge, we apply to correction to each realspaceGridnode
		//LAL::optional<float> backgroundchargeCorrection;

		void ForwardFFT() {
			cufftResult result = cufftExecR2C(planForward, realspaceGrid, fourierspaceGrid);
			if (result != CUFFT_SUCCESS) {
				fprintf(stderr, "cufftExecR2C failed with error code %d\n", result);
			}
		}

		void InverseFFT() {
			cufftResult result = cufftExecC2R(planInverse, fourierspaceGrid, realspaceGrid);
			if (result != CUFFT_SUCCESS) {
				fprintf(stderr, "cufftExecC2R failed with error code %d\n", result);
			}
		}

		void CalcEnergyCorrection(const Box& box) {
			double chargeSquaredSum = 0;
			double chargeSum = 0;
			for (const Compound& compound : box.compounds) {
				for (int pid = 0; pid < compound.n_particles; ++pid) {
					const double charge = compound.atom_charges[pid];		// [kC/mol]
					chargeSquaredSum += charge * charge;
					chargeSum += charge;
				}
			}

			selfenergyCorrection = static_cast<float>(-ewaldKappa / sqrt(PI) * chargeSquaredSum * PhysicsUtils::modifiedCoulombConstant_Potential);

			/*if (std::abs(chargeSum) > 1e-6) {
				backgroundchargeCorrection = chargeSum / static_cast<double>
			}*/
			// Only relevant for systems with a net charge
			//const float volume = static_cast<float>(box.boxparams.boxSize) * static_cast<float>(box.boxparams.boxSize) * static_cast<float>(box.boxparams.boxSize);	
			//const float backgroundEnergyCorrection = static_cast<float>(-PI * chargeSum * chargeSum / (kappa * kappa * volume) * PhysicsUtils::modifiedCoulombConstant_Potential);
			//return selfEnergyCorrection + backgroundEnergyCorrection;		// [J/mol]
		}

	public:

		Controller(int boxLenNm, const Box& box, float cutoffNM);
		~Controller();

		void CalcCharges(const BoxConfig& config, const BoxState& state, int nCompounds, ForceEnergy* const forceEnergy);

	private:
		//Just for debugging
		void PlotPotentialSlices() {

			std::vector<float> gridHost;
			GenericCopyToHost(realspaceGrid, gridHost, nGridpointsRealspace);

			int centerSlice = 100;
			int numSlices = 1;
			int spacing = 10;

			std::vector<float> combinedData;
			std::vector<int> sliceIndices;

			for (int i = -numSlices; i <= numSlices; ++i) {
				int sliceIndex = centerSlice + i * spacing;
				if (sliceIndex < 0 || sliceIndex >= gridpointsPerDim) {
					throw std::runtime_error("error");
				}

				int firstIndex = gridpointsPerDim * gridpointsPerDim * sliceIndex;
				int lastIndex = firstIndex + gridpointsPerDim * gridpointsPerDim;
				combinedData.insert(combinedData.end(), gridHost.begin() + firstIndex, gridHost.begin() + lastIndex);
				sliceIndices.push_back(sliceIndex);
			}

			// Save all slices to one file
			FileUtils::WriteVectorToBinaryFile("C:/Users/Daniel/git_repo/LIMA_data/Pool/PmePot_AllSlices.bin", combinedData);

			// Call the Python script to plot the slices
			std::string pyscriptPath = (FileUtils::GetLimaDir() / "dev" / "PyTools" / "Plot2dVec.py").string();
			std::string command = "python " + pyscriptPath + " " + std::to_string(numSlices * 2 + 1) + " " + std::to_string(gridpointsPerDim) + " " + std::to_string(centerSlice) + " " + std::to_string(spacing);
			std::system(command.c_str());
		}

	};
}