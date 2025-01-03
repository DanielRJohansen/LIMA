//#pragma once // Only allowed to be included by engine.cu

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"
#include "BoundaryCondition.cuh"

#include "Filehandling.h"
#include "PhysicsUtils.cuh"

#include <cuda_runtime.h>
#include <cufft.h>
#include "Utilities.h"


// TODO: Do i need to account for e0, vacuum/spaceial permitivity here? Probably....

namespace PME {
	const int gridpointsPerNm = 20;
	constexpr float gridpointsPerNm_f = static_cast<float>(gridpointsPerNm);

	__device__ __host__ int GetGridIndexRealspace(const Int3& gridIndex, int gridpointsPerDim) {
		return gridIndex.x + gridIndex.y * gridpointsPerDim + gridIndex.z * gridpointsPerDim * gridpointsPerDim;
	}
	__device__ __host__ int GetGridIndexReciprocalspace(const Int3& gridIndex, int gridpointsPerDim, int nGridpointsHalfdim) {
		return gridIndex.x + gridIndex.y * nGridpointsHalfdim + gridIndex.z * nGridpointsHalfdim * gridpointsPerDim;
	}
	__device__ NodeIndex Get3dIndexRealspace(int index1d, int gridpointsPerDim) {
		int z = index1d / (gridpointsPerDim * gridpointsPerDim);
		index1d -= z * gridpointsPerDim * gridpointsPerDim;
		int y = index1d / gridpointsPerDim;
		index1d -= y * gridpointsPerDim;
		int x = index1d;
		return NodeIndex{ x, y, z };
	}
	__device__ NodeIndex Get3dIndexReciprocalspace(int index1d, int gridpointsPerDim, int nGridpointsHalfdim) {
		int z = index1d / (nGridpointsHalfdim * gridpointsPerDim);
		index1d -= z * nGridpointsHalfdim * gridpointsPerDim;
		int y = index1d / nGridpointsHalfdim;
		index1d -= y * nGridpointsHalfdim;
		int x = index1d;
		return NodeIndex{ x, y, z };
	}



	__global__ void DistributeChargesToGridKernel(const BoxConfig config, const BoxState state, float* chargeGrid, int gridpointsPerDim) {
		if (threadIdx.x >= config.compounds[blockIdx.x].n_particles)
			return;

		const float charge = config.compoundsAtomCharges[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];	// [kC/mol]
		if (charge == 0.f)
			return;

		const NodeIndex origo = state.compoundOrigos[blockIdx.x];
		const Float3 relpos = state.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
		Float3 absPos = relpos + origo.toFloat3();
		PeriodicBoundaryCondition::applyBCNM(absPos);

		// Map absolute position to fractional grid coordinates
		Float3 gridPos = absPos * gridpointsPerNm_f;
		int ix = static_cast<int>(floorf(gridPos.x));
		int iy = static_cast<int>(floorf(gridPos.y));
		int iz = static_cast<int>(floorf(gridPos.z));

		float fx = gridPos.x - static_cast<float>(ix);
		float fy = gridPos.y - static_cast<float>(iy);
		float fz = gridPos.z - static_cast<float>(iz);
		
		// Cubic B-spline weights along x
		float wx[4];
		wx[0] = (1.f - fx) * (1.f - fx) * (1.f - fx) / 6.f;
		wx[1] = (4.f - 6.f * fx * fx + 3.f * fx * fx * fx) / 6.f;
		wx[2] = (1.f + 3.f * fx + 3.f * fx * fx - 3.f * fx * fx * fx) / 6.f;
		wx[3] = (fx * fx * fx) / 6.f;

		// Cubic B-spline weights along y
		float wy[4];
		wy[0] = (1.f - fy) * (1.f - fy) * (1.f - fy) / 6.f;
		wy[1] = (4.f - 6.f * fy * fy + 3.f * fy * fy * fy) / 6.f;
		wy[2] = (1.f + 3.f * fy + 3.f * fy * fy - 3.f * fy * fy * fy) / 6.f;
		wy[3] = (fy * fy * fy) / 6.f;

		// Cubic B-spline weights along z
		float wz[4];
		wz[0] = (1.f - fz) * (1.f - fz) * (1.f - fz) / 6.f;
		wz[1] = (4.f - 6.f * fz * fz + 3.f * fz * fz * fz) / 6.f;
		wz[2] = (1.f + 3.f * fz + 3.f * fz * fz - 3.f * fz * fz * fz) / 6.f;
		wz[3] = (fz * fz * fz) / 6.f;


		// Distribute the charge to the surrounding 4x4x4 cells
		for (int dx = 0; dx < 4; dx++) {
			int X = -1 + ix + dx;
			float wxCur = wx[dx];
			for (int dy = 0; dy < 4; dy++) {
				int Y = -1 + iy + dy;
				float wxyCur = wxCur * wy[dy];
				for (int dz = 0; dz < 4; dz++) {
					int Z = -1 + iz + dz;
					//printf("Weight %f\n", wxyCur * wz[dz]);

					const NodeIndex index3d = PeriodicBoundaryCondition::applyBC(NodeIndex{ X,Y,Z }, gridpointsPerDim);
					const int index1D = GetGridIndexRealspace(index3d, gridpointsPerDim);
					chargeGrid[index1D] += charge * wxyCur * wz[dz]; // TODO This must be applyBC
				}
			}
		}
	}

	__global__ void InterpolateForcesAndPotentialKernel(
		const BoxConfig config, 
		const BoxState state, 
		const float* potentialGrid, 
		int gridpointsPerDim, 
		ForceEnergy* const forceEnergies, 
		float selfenergyCorrection			// [J/mol]
	)
	{
		if (threadIdx.x >= config.compounds[blockIdx.x].n_particles) {
			return;
		}

		const float charge = config.compoundsAtomCharges[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];	// [kC/mol]
		if (charge == 0.f)
			return;

		const NodeIndex origo = state.compoundOrigos[blockIdx.x];
		const Float3 relpos = state.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
		Float3 absPos = relpos + origo.toFloat3();
		PeriodicBoundaryCondition::applyBCNM(absPos);

		const Float3 gridPos = absPos * gridpointsPerNm_f;
		int ix = static_cast<int>(floorf(gridPos.x));
		int iy = static_cast<int>(floorf(gridPos.y));
		int iz = static_cast<int>(floorf(gridPos.z));

		float fx = gridPos.x - static_cast<float>(ix);
		float fy = gridPos.y - static_cast<float>(iy);
		float fz = gridPos.z - static_cast<float>(iz);

		float wx[4];
		wx[0] = (1.f - fx) * (1.f - fx) * (1.f - fx) / 6.f;
		wx[1] = (4.f - 6.f * fx * fx + 3.f * fx * fx * fx) / 6.f;
		wx[2] = (1.f + 3.f * fx + 3.f * fx * fx - 3.f * fx * fx * fx) / 6.f;
		wx[3] = (fx * fx * fx) / 6.f;

		float wy[4];
		wy[0] = (1.f - fy) * (1.f - fy) * (1.f - fy) / 6.f;
		wy[1] = (4.f - 6.f * fy * fy + 3.f * fy * fy * fy) / 6.f;
		wy[2] = (1.f + 3.f * fy + 3.f * fy * fy - 3.f * fy * fy * fy) / 6.f;
		wy[3] = (fy * fy * fy) / 6.f;

		float wz[4];
		wz[0] = (1.f - fz) * (1.f - fz) * (1.f - fz) / 6.f;
		wz[1] = (4.f - 6.f * fz * fz + 3.f * fz * fz * fz) / 6.f;
		wz[2] = (1.f + 3.f * fz + 3.f * fz * fz - 3.f * fz * fz * fz) / 6.f;
		wz[3] = (fz * fz * fz) / 6.f;

		Float3 force{};			// [J/mol/nm]
		float potential{};		// [J/mol]

		for (int dx = 0; dx < 4; dx++) {
			int X = ix - 1 + dx;
			float wxCur = wx[dx];
			for (int dy = 0; dy < 4; dy++) {
				int Y = iy - 1 + dy;
				float wxyCur = wxCur * wy[dy];
				for (int dz = 0; dz < 4; dz++) {
					int Z = iz - 1 + dz;
					float wxyzCur = wxyCur * wz[dz];

					const NodeIndex node = PeriodicBoundaryCondition::applyBC(NodeIndex{ X, Y, Z }, gridpointsPerDim);
					const int gridIndex = GetGridIndexRealspace(node, gridpointsPerDim);

					float phi = potentialGrid[gridIndex];

					NodeIndex plusX  = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x + 1, node.y,     node.z }, gridpointsPerDim);
					NodeIndex minusX = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x - 1, node.y,     node.z }, gridpointsPerDim);
					NodeIndex plusY  = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y + 1, node.z }, gridpointsPerDim);
					NodeIndex minusY = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y - 1, node.z }, gridpointsPerDim);
					NodeIndex plusZ  = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y,     node.z + 1 }, gridpointsPerDim);
					NodeIndex minusZ = PeriodicBoundaryCondition::applyBC(NodeIndex{ node.x,     node.y,     node.z - 1 }, gridpointsPerDim);

					float phi_plusX = potentialGrid[GetGridIndexRealspace(plusX, gridpointsPerDim)];// / (gridpointsPerDim * gridpointsPerDim * gridpointsPerDim);
					float phi_minusX = potentialGrid[GetGridIndexRealspace(minusX, gridpointsPerDim)];// / (gridpointsPerDim * gridpointsPerDim * gridpointsPerDim);
					float phi_plusY = potentialGrid[GetGridIndexRealspace(plusY, gridpointsPerDim)];// / (gridpointsPerDim * gridpointsPerDim * gridpointsPerDim);
					float phi_minusY = potentialGrid[GetGridIndexRealspace(minusY, gridpointsPerDim)];// / (gridpointsPerDim * gridpointsPerDim * gridpointsPerDim);
					float phi_plusZ = potentialGrid[GetGridIndexRealspace(plusZ, gridpointsPerDim)];// / (gridpointsPerDim * gridpointsPerDim * gridpointsPerDim);
					float phi_minusZ = potentialGrid[GetGridIndexRealspace(minusZ, gridpointsPerDim)];// / (gridpointsPerDim * gridpointsPerDim * gridpointsPerDim);

					float E_x = -(phi_plusX - phi_minusX) * (gridpointsPerNm / 2.0f);
					float E_y = -(phi_plusY - phi_minusY) * (gridpointsPerNm / 2.0f);
					float E_z = -(phi_plusZ - phi_minusZ) * (gridpointsPerNm / 2.0f);

					force += Float3{ E_x, E_y, E_z } *wxyzCur;
					potential += phi * wxyzCur;
				}
			}
		}

		// Now add self charge to calculations
		force *= charge;
		potential *= charge;

		// Ewald self-energy correction
		potential += selfenergyCorrection;

		potential *= 0.5f; // Potential is halved because we computing for both this and the other particle's

		forceEnergies[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = ForceEnergy{ force, potential};
	}



	__global__ void PrecomputeGreensFunctionKernel(float* d_greensFunction, int gridpointsPerDim, 
		double boxLen,		// [nm]
		double ewaldKappa	// [nm^-1]
	) {
		const int halfNodes = gridpointsPerDim / 2;
		int nGridpointsHalfdim = gridpointsPerDim / 2 + 1;
		const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
		/*if (index1D > nGridpointsHalfdim * gridpointsPerDim * gridpointsPerDim)
			return;*/
		NodeIndex freqIndex = Get3dIndexReciprocalspace(index1D, gridpointsPerDim, nGridpointsHalfdim);

		if (freqIndex.x >= nGridpointsHalfdim || freqIndex.y >= gridpointsPerDim || freqIndex.z >= gridpointsPerDim)
			return;

		// Remap frequencies to negative for indices > N/2
		int kxIndex = freqIndex.x;
		int kyShiftedIndex = (freqIndex.y <= halfNodes) ? freqIndex.y : freqIndex.y - gridpointsPerDim;
		int kzShiftedIndex = (freqIndex.z <= halfNodes) ? freqIndex.z : freqIndex.z - gridpointsPerDim;

		// Ewald kappa fixed
		double volume = boxLen * boxLen * boxLen;				// [nm^3]
		double delta = boxLen / (double)gridpointsPerDim;		// [nm]

		// Physical wavevectors
		double kx = (2.0 * PI * (double)kxIndex) / boxLen;
		double ky = (2.0 * PI * (double)kyShiftedIndex) / boxLen;
		double kz = (2.0 * PI * (double)kzShiftedIndex) / boxLen;

		double kSquared = kx * kx + ky * ky + kz * kz;

		double currentGreensValue = 0.0f;

		// Compute B-spline structure factor (4th order)
		double kHalfX = kx * (delta * 0.5);
		double kHalfY = ky * (delta * 0.5);
		double kHalfZ = kz * (delta * 0.5);

		const double epsilon = 1e-8;

		auto splineFactor = [epsilon](double kh) {
			if (fabs(kh) < epsilon) return 1.0;
			double ratio = sin(kh) / kh;
			return pow(ratio, 4);			
			};

		double Sx = splineFactor(kHalfX);
		double Sy = splineFactor(kHalfY);
		double Sz = splineFactor(kHalfZ);
		
		double splineCorrection = (Sx * Sy * Sz);
		splineCorrection = splineCorrection * splineCorrection; // squared for forward+back interpolation

		//splineCorrection = 1;

		if (kSquared > epsilon) {
			currentGreensValue = (4.0 * PI / (volume * kSquared)) 
				* exp(-kSquared / (4.0 * ewaldKappa * ewaldKappa)) 
				* splineCorrection
				* PhysicsUtilsDevice::modifiedCoulombConstant_Force
				;
		}

		//if (freqIndex.y == halfNodes || freqIndex.z == halfNodes) {
		//	currentGreensValue = 0.0;
		//}

		d_greensFunction[index1D] = static_cast<float>(currentGreensValue);
	}


	__global__ void ApplyGreensFunctionKernel(
		cufftComplex* d_reciprocalFreqData, 
		const float* d_greensFunctionArray, 
		int gridpointsPerDim, 
		float boxLength							// [nm]
	)
	{
		int nGridpointsHalfdim = gridpointsPerDim / 2 + 1;// TODO: add comment here
		const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
		if (index1D > nGridpointsHalfdim * gridpointsPerDim * gridpointsPerDim)
			return;


		// Retrieve the precomputed Green’s function value at this frequency point.
		const float currentGreenValue = d_greensFunctionArray[index1D];

		// The frequency-domain data is complex, and we multiply by the real Green’s function to apply the long-range correction.
		cufftComplex currentFreqValue = d_reciprocalFreqData[index1D];

		// Scaling the complex frequency data by the Green’s function imposes the correct electrostatic filtering.
		currentFreqValue.x *= currentGreenValue;
		currentFreqValue.y *= currentGreenValue;

		d_reciprocalFreqData[index1D] = currentFreqValue;
	}


	class Controller {

		float* chargeGrid;
		cufftComplex* potentialGrid;
		float* greensFunctionScalars;

		int gridpointsPerDim=-1;
		size_t nGridpointsRealspace = -1;
		int nGridpointsReciprocalspace = -1;
		const float ewaldKappa;
		float boxlenNm{};

		// Always applied constant per particle
		float selfenergyCorrection;

		cufftHandle planForward;
		cufftHandle planInverse;

		// For system with a net charge, we apply to correction to each chargegridnode
		//LAL::optional<float> backgroundchargeCorrection;

		void ForwardFFT() {
			cufftResult result = cufftExecR2C(planForward, chargeGrid, potentialGrid);
			if (result != CUFFT_SUCCESS) {
				fprintf(stderr, "cufftExecR2C failed with error code %d\n", result);
			}
		}

		void InverseFFT() {		
			cufftResult result = cufftExecC2R(planInverse, potentialGrid, chargeGrid);
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

		Controller(int boxLenNm, const Box& box, float cutoffNM) : boxlenNm(boxLenNm), ewaldKappa(PhysicsUtils::CalcEwaldkappa(cutoffNM))
		{
			gridpointsPerDim = boxLenNm * gridpointsPerNm;
			nGridpointsRealspace = gridpointsPerDim * gridpointsPerDim * gridpointsPerDim;
			nGridpointsReciprocalspace = gridpointsPerDim * gridpointsPerDim * (gridpointsPerDim / 2 + 1);

			if (nGridpointsRealspace > INT32_MAX)
				throw std::runtime_error("Ewald grid too large to index with integers");

			const size_t byteSize = nGridpointsRealspace * sizeof(float) + nGridpointsReciprocalspace * sizeof(float) * 3; // 3= 2 from cufftComplex, 1 from greensfunction
			if (byteSize > 10'000'000'000)
				throw std::runtime_error("Ewald grid too large");

			CalcEnergyCorrection(box);

			cufftPlan3d(&planForward, gridpointsPerDim, gridpointsPerDim, gridpointsPerDim, CUFFT_R2C);
			cufftPlan3d(&planInverse, gridpointsPerDim, gridpointsPerDim, gridpointsPerDim, CUFFT_C2R);

			cudaMalloc(&chargeGrid, nGridpointsRealspace * sizeof(float));
			cudaMalloc(&potentialGrid, nGridpointsReciprocalspace * sizeof(cufftComplex));
			cudaMalloc(&greensFunctionScalars, nGridpointsReciprocalspace * sizeof(float));
			cudaDeviceSynchronize();

			const int nBlocks = (nGridpointsReciprocalspace + 63) / 64;
			PrecomputeGreensFunctionKernel<<<nBlocks, 64 >>>(greensFunctionScalars, gridpointsPerDim, boxLenNm, ewaldKappa);
			LIMA_UTILS::genericErrorCheck("PrecomputeGreensFunctionKernel failed!");
		}
		~Controller() {
			cudaDeviceSynchronize();

			cufftDestroy(planForward);
			cufftDestroy(planInverse);

			cudaFree(chargeGrid);
			cudaFree(potentialGrid);
			cudaFree(greensFunctionScalars);
		}

		void CalcCharges(const BoxConfig& config, const BoxState& state, int nCompounds, ForceEnergy* const forceEnergy) {
			if (nCompounds == 0)
				return;

			cudaMemset(chargeGrid, 0, nGridpointsRealspace * sizeof(float));// maybe do this async, after the last kernel?
			cudaMemset(potentialGrid, 0, nGridpointsReciprocalspace * sizeof(cufftComplex));
			cudaDeviceSynchronize();
			DistributeChargesToGridKernel<<<nCompounds, MAX_COMPOUND_PARTICLES>> > (config, state, chargeGrid, gridpointsPerDim);
			LIMA_UTILS::genericErrorCheckNoSync("DistributeChargesToGridKernel failed!");
			cudaDeviceSynchronize();

			ForwardFFT();
			cudaDeviceSynchronize();

			ApplyGreensFunctionKernel<<<(nGridpointsReciprocalspace + 63) / 64, 64 >> > (potentialGrid, greensFunctionScalars, gridpointsPerDim, boxlenNm);
			LIMA_UTILS::genericErrorCheckNoSync("ApplyGreensFunctionKernel failed!");

			InverseFFT();
			cudaDeviceSynchronize();

			//PlotPotentialSlices();

			InterpolateForcesAndPotentialKernel << <nCompounds, MAX_COMPOUND_PARTICLES >> > (config, state, chargeGrid, gridpointsPerDim, forceEnergy, selfenergyCorrection);
			LIMA_UTILS::genericErrorCheckNoSync("InterpolateForcesAndPotentialKernel failed!");

			cudaDeviceSynchronize();
		}
		
	private:
		//Just for debugging
		void PlotPotentialSlices() {

			std::vector<float> gridHost;
			GenericCopyToHost(chargeGrid, gridHost, nGridpointsRealspace);

			int centerSlice = 200;
			int numSlices = 1;
			int spacing = 20;

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