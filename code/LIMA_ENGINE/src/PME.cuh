

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"
#include "BoundaryCondition.cuh"

#include "Filehandling.h"
#include "PhysicsUtils.cuh"

#include <cuda_runtime.h>
#include <cufft.h>
//#include <optional>
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

		const double epsilon = 1e-24;

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

		splineCorrection = 1;

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
		int nGridpointsRealspace = -1;
		int nGridpointsReciprocalspace = -1;
		const float ewaldKappa;
		float boxlenNm{};

		// Always applied constant per particle
		float selfenergyCorrection;

		// For system with a net charge, we apply to correction to each chargegridnode
		//LAL::optional<float> backgroundchargeCorrection;

		void ForwardFFT() {
			cufftHandle plan; // TODO move this to member var. Also check for error when creating
			cufftPlan3d(&plan, gridpointsPerDim, gridpointsPerDim, gridpointsPerDim, CUFFT_R2C);
			cufftResult result = cufftExecR2C(plan, chargeGrid, potentialGrid);
			if (result != CUFFT_SUCCESS) {
				fprintf(stderr, "cufftExecR2C failed with error code %d\n", result);
			}
			cufftDestroy(plan);
		}

		void InverseFFT() {
			cufftHandle plan;
			cufftPlan3d(&plan, gridpointsPerDim, gridpointsPerDim, gridpointsPerDim, CUFFT_C2R);

			cufftResult result = cufftExecC2R(plan, potentialGrid, chargeGrid);
			if (result != CUFFT_SUCCESS) {
				fprintf(stderr, "cufftExecC2R failed with error code %d\n", result);
			}
			cufftDestroy(plan);
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


			cudaMalloc(&chargeGrid, nGridpointsRealspace * sizeof(float));
			cudaMalloc(&potentialGrid, nGridpointsReciprocalspace * sizeof(cufftComplex));
			cudaMalloc(&greensFunctionScalars, nGridpointsReciprocalspace * sizeof(float));
			cudaDeviceSynchronize();

			const int nBlocks = (nGridpointsReciprocalspace + 63) / 64;
			PrecomputeGreensFunctionKernel<<<nBlocks, 64 >>>(greensFunctionScalars, gridpointsPerDim, boxLenNm, ewaldKappa);
			LIMA_UTILS::genericErrorCheck("PrecomputeGreensFunctionKernel failed!");
		}
		~Controller() {
			cudaFree(chargeGrid);
			cudaFree(potentialGrid);
			cudaFree(greensFunctionScalars);
		}

		void CalcCharges(const BoxConfig& config, const BoxState& state, int nCompounds, ForceEnergy* const forceEnergy) {
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


/// Code to plot the potential grid
//

//
////const int firstIndex = gridpointsPerDim * gridpointsPerDim * 200;
////const int lastIndex = firstIndex + gridpointsPerDim * gridpointsPerDim;
////std::vector<float> pots;
////GenericCopyToHost(chargeGrid, pots, nGridpointsRealspace);
////std::vector<float> data = std::vector<float>(pots.begin() + firstIndex, pots.begin() + lastIndex);
////FileUtils::WriteVectorToBinaryFile(R"(C:\Users\Daniel\git_repo\LIMA_data\Pool\PmePot.bin)", data);
//
//{
//	//int i = 0;
//	//const int y = 92;
//	//for (int z = 0; z < gridpointsPerDim; z++) {
//	//	for (int x = 0; x < gridpointsPerDim; x++) {
//	//		const int index1D = GetGridIndexRealspace(Int3{ x,y,z }, gridpointsPerDim);
//	//		data[i] = pots[index1D];
//	//	}
//	//}
//	//FileUtils::WriteVectorToBinaryFile(R"(C:\Users\Daniel\git_repo\LIMA_data\Pool\PmePot.bin)", data);
//}





//	__global__ void InterpolateForcesAndPotentialKernel1(const BoxConfig config, const BoxState state, const float* potentialGrid, int gridpointsPerDim, ForceEnergy* const forceEnergies) 
//	{
//		if (threadIdx.x >= config.compounds[blockIdx.x].n_particles) {
//			return;
//		}
//
//		const float charge = config.compoundsAtomCharges[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
//		if (charge == 0.f)
//			return;
//
//		const NodeIndex origo = state.compoundOrigos[blockIdx.x];
//		const Float3 relpos = state.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
//		Float3 absPos = relpos + origo.toFloat3();
//		PeriodicBoundaryCondition::applyBCNM(absPos);
//
//		const Float3 gridnodePosition = (absPos * gridpointsPerNm_f);
//		const NodeIndex gridnode000 = PeriodicBoundaryCondition::applyBC({ static_cast<int>(floor(gridnodePosition.x)), static_cast<int>(floor(gridnodePosition.y)), static_cast<int>(floor(gridnodePosition.z)) }, gridpointsPerDim);
//
//
//		//if (blockIdx.x != 1)
//		//	return;
//
////		printf("gridnode000 %d %d %d\n", gridnode000.x, gridnode000.y, gridnode000.z);
//
//		const Float3 fraction = gridnodePosition - gridnode000.toFloat3();
//
//		// Identify the eight surrounding grid points
//		NodeIndex gridIndices[8];
//		gridIndices[0] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 0, 0 }, gridpointsPerDim); // c000
//		gridIndices[1] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 0, 0 }, gridpointsPerDim); // c100
//		gridIndices[2] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 1, 0 }, gridpointsPerDim); // c010
//		gridIndices[3] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 1, 0 }, gridpointsPerDim); // c110
//		gridIndices[4] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 0, 1 }, gridpointsPerDim); // c001
//		gridIndices[5] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 0, 1 }, gridpointsPerDim); // c101
//		gridIndices[6] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 1, 1 }, gridpointsPerDim); // c011
//		gridIndices[7] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 1, 1 }, gridpointsPerDim); // c111
//
//		float weights[8];
//		weights[0] = (1.0f - fraction.x) * (1.0f - fraction.y) * (1.0f - fraction.z);
//		weights[1] = fraction.x * (1.0f - fraction.y) * (1.0f - fraction.z);
//		weights[2] = (1.0f - fraction.x) * fraction.y * (1.0f - fraction.z);
//		weights[3] = fraction.x * fraction.y * (1.0f - fraction.z);
//		weights[4] = (1.0f - fraction.x) * (1.0f - fraction.y) * fraction.z;
//		weights[5] = fraction.x * (1.0f - fraction.y) * fraction.z;
//		weights[6] = (1.0f - fraction.x) * fraction.y * fraction.z;
//		weights[7] = fraction.x * fraction.y * fraction.z;
//
//		Float3 force{};
//		float potential{};
//
//		// Loop over the eight surrounding grid nodes to interpolate forces and potential
//		//#pragma unroll
//		for (int i = 0; i < 8; i++) {
//			// Compute the linear index for the current grid node using zyx-major ordering
//			//int gridIndex = getGridIndex(gridIndices[i], gridpointsPerDim); // WHY: Converts 3D grid indices to 1D index for memory access
//			const int gridIndex = GetGridIndexRealspace(gridIndices[i], gridpointsPerDim);
//
//
//			// Retrieve the potential at the current grid node
//			float phi = potentialGrid[gridIndex]; // WHY: Potential at the grid node is needed for both force and potential interpolation
//
//			// Compute neighboring grid indices for finite differences to estimate the electric field (gradient of potential)
//			NodeIndex plusX = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x + 1, gridIndices[i].y, gridIndices[i].z }, gridpointsPerDim);
//			NodeIndex minusX = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x - 1, gridIndices[i].y, gridIndices[i].z }, gridpointsPerDim);
//			NodeIndex plusY = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x, gridIndices[i].y + 1, gridIndices[i].z }, gridpointsPerDim);
//			NodeIndex minusY = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x, gridIndices[i].y - 1, gridIndices[i].z }, gridpointsPerDim);
//			NodeIndex plusZ = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x, gridIndices[i].y, gridIndices[i].z + 1 }, gridpointsPerDim);
//			NodeIndex minusZ = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x, gridIndices[i].y, gridIndices[i].z - 1 }, gridpointsPerDim);
//
//			// Retrieve potentials at grid nodes neighboring the current query
//			float phi_plusX = potentialGrid[GetGridIndexRealspace(plusX, gridpointsPerDim)];
//			float phi_minusX = potentialGrid[GetGridIndexRealspace(minusX, gridpointsPerDim)];
//			float phi_plusY = potentialGrid[GetGridIndexRealspace(plusY, gridpointsPerDim)];
//			float phi_minusY = potentialGrid[GetGridIndexRealspace(minusY, gridpointsPerDim)];
//			float phi_plusZ = potentialGrid[GetGridIndexRealspace(plusZ, gridpointsPerDim)];
//			float phi_minusZ = potentialGrid[GetGridIndexRealspace(minusZ, gridpointsPerDim)];
//
//			// Compute electric field components using central finite differences (E = -grad(phi))
//			float E_x = -(phi_plusX - phi_minusX) * (gridpointsPerNm / 2.0f);
//			float E_y = -(phi_plusY - phi_minusY) * (gridpointsPerNm / 2.0f);
//			float E_z = -(phi_plusZ - phi_minusZ) * (gridpointsPerNm / 2.0f);
//
//			// Accumulate the weighted electric field contributions to compute the particle's force
//			force += Float3{ E_x, E_y, E_z } * weights[i];		
//
//			if (blockIdx.x == 0)
//				printf("weights %f phi %f\n", weights[i], phi);
//
//			//if (i == 0)
//			//	printf("E %f %f %f\n", E_x, E_y, E_z);
//
//			// Accumulate the weighted potential contributions to compute the particle's potential
//			potential += phi * weights[i];
//		}
//
//		// Store the computed force and potential in the output arrays
//		forceEnergies[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = ForceEnergy{ force * charge, potential * charge };
//	}