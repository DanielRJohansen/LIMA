

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"
#include "BoundaryCondition.cuh"

#include "PhysicsUtils.cuh"

#include <cuda_runtime.h>
#include <cufft.h>

















namespace PME {
	const int gridpointsPerNm = 20;
	constexpr float gridpointsPerNm_f = 20.f;

	__device__ int getGridIndex(const Int3& gridIndex, int gridpointsPerDim) {
		return gridIndex.x + gridIndex.y * gridpointsPerDim + gridIndex.z * gridpointsPerDim * gridpointsPerDim;
	}

	__global__ void DistributeChargesToGridKernel(const BoxConfig config, const BoxState state, float* chargeGrid, int gridpointsPerDim) {
		if (threadIdx.x >= config.compounds[blockIdx.x].n_particles) {
			return;
		}

		const float charge = config.compoundsAtomCharges[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
		if (charge == 0.f)
			return;

		const NodeIndex origo = state.compoundOrigos[blockIdx.x];
		const Float3 relpos = state.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
		Float3 absPos = relpos + origo.toFloat3();
		PeriodicBoundaryCondition::applyBCNM(absPos);

		const Float3 gridnodePosition = (absPos * gridpointsPerNm_f);
		const NodeIndex gridnode000{ static_cast<int>(floor(gridnodePosition.x)), static_cast<int>(floor(gridnodePosition.y)), static_cast<int>(floor(gridnodePosition.z)) };

		const Float3 fraction = gridnodePosition - gridnode000.toFloat3();
		//fraction.print('f');


		// Identify the eight surrounding grid points
		NodeIndex gridIndices[8];
		gridIndices[0] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 0, 0 }, gridpointsPerDim); // c000
		gridIndices[1] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 0, 0 }, gridpointsPerDim); // c100
		gridIndices[2] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 1, 0 }, gridpointsPerDim); // c010
		gridIndices[3] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 1, 0 }, gridpointsPerDim); // c110
		gridIndices[4] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 0, 1 }, gridpointsPerDim); // c001
		gridIndices[5] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 0, 1 }, gridpointsPerDim); // c101
		gridIndices[6] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 1, 1 }, gridpointsPerDim); // c011
		gridIndices[7] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 1, 1 }, gridpointsPerDim); // c111

		float weights[8];
		weights[0] = (1.0f - fraction.x) * (1.0f - fraction.y) * (1.0f - fraction.z);
		weights[1] = fraction.x * (1.0f - fraction.y) * (1.0f - fraction.z);
		weights[2] = (1.0f - fraction.x) * fraction.y * (1.0f - fraction.z);
		weights[3] = fraction.x * fraction.y * (1.0f - fraction.z);
		weights[4] = (1.0f - fraction.x) * (1.0f - fraction.y) * fraction.z;
		weights[5] = fraction.x * (1.0f - fraction.y) * fraction.z;
		weights[6] = (1.0f - fraction.x) * fraction.y * fraction.z;
		weights[7] = fraction.x * fraction.y * fraction.z;



		#pragma unroll
		for (int i = 0; i < 8; i++) {
			const int gridIndex = getGridIndex(gridIndices[i], gridpointsPerDim);
			/*if (charge > 0.f)
				printf("weighted charge %f\n", charge * weights[i]);*/
			chargeGrid[gridIndex] += charge * weights[i];
		}
	}

	__global__ void InterpolateForcesAndPotentialKernel(const BoxConfig config, const BoxState state, const float* potentialGrid, int gridpointsPerDim, ForceEnergy* const forceEnergies) 
	{
		if (threadIdx.x >= config.compounds[blockIdx.x].n_particles) {
			return;
		}

		const float charge = config.compoundsAtomCharges[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
		if (charge == 0.f)
			return;

		const NodeIndex origo = state.compoundOrigos[blockIdx.x];
		const Float3 relpos = state.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
		Float3 absPos = relpos + origo.toFloat3();
		PeriodicBoundaryCondition::applyBCNM(absPos);

		const Float3 gridnodePosition = (absPos * gridpointsPerNm_f);
		const NodeIndex gridnode000 = PeriodicBoundaryCondition::applyBC({ static_cast<int>(floor(gridnodePosition.x)), static_cast<int>(floor(gridnodePosition.y)), static_cast<int>(floor(gridnodePosition.z)) }, gridpointsPerDim);

		const Float3 fraction = gridnodePosition - gridnode000.toFloat3();

		// Identify the eight surrounding grid points
		NodeIndex gridIndices[8];
		gridIndices[0] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 0, 0 }, gridpointsPerDim); // c000
		gridIndices[1] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 0, 0 }, gridpointsPerDim); // c100
		gridIndices[2] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 1, 0 }, gridpointsPerDim); // c010
		gridIndices[3] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 1, 0 }, gridpointsPerDim); // c110
		gridIndices[4] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 0, 1 }, gridpointsPerDim); // c001
		gridIndices[5] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 0, 1 }, gridpointsPerDim); // c101
		gridIndices[6] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 0, 1, 1 }, gridpointsPerDim); // c011
		gridIndices[7] = PeriodicBoundaryCondition::applyBC(gridnode000 + NodeIndex{ 1, 1, 1 }, gridpointsPerDim); // c111

		float weights[8];
		weights[0] = (1.0f - fraction.x) * (1.0f - fraction.y) * (1.0f - fraction.z);
		weights[1] = fraction.x * (1.0f - fraction.y) * (1.0f - fraction.z);
		weights[2] = (1.0f - fraction.x) * fraction.y * (1.0f - fraction.z);
		weights[3] = fraction.x * fraction.y * (1.0f - fraction.z);
		weights[4] = (1.0f - fraction.x) * (1.0f - fraction.y) * fraction.z;
		weights[5] = fraction.x * (1.0f - fraction.y) * fraction.z;
		weights[6] = (1.0f - fraction.x) * fraction.y * fraction.z;
		weights[7] = fraction.x * fraction.y * fraction.z;

		Float3 force{};
		float potential{};

		// Loop over the eight surrounding grid nodes to interpolate forces and potential
		//#pragma unroll
		for (int i = 0; i < 8; i++) {
			// Compute the linear index for the current grid node using zyx-major ordering
			int gridIndex = getGridIndex(gridIndices[i], gridpointsPerDim); // WHY: Converts 3D grid indices to 1D index for memory access

			// Retrieve the potential at the current grid node
			float phi = potentialGrid[gridIndex]; // WHY: Potential at the grid node is needed for both force and potential interpolation

			// Compute neighboring grid indices for finite differences to estimate the electric field (gradient of potential)
			NodeIndex plusX = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x + 1, gridIndices[i].y, gridIndices[i].z }, gridpointsPerDim);
			NodeIndex minusX = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x - 1, gridIndices[i].y, gridIndices[i].z }, gridpointsPerDim);
			NodeIndex plusY = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x, gridIndices[i].y + 1, gridIndices[i].z }, gridpointsPerDim);
			NodeIndex minusY = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x, gridIndices[i].y - 1, gridIndices[i].z }, gridpointsPerDim);
			NodeIndex plusZ = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x, gridIndices[i].y, gridIndices[i].z + 1 }, gridpointsPerDim);
			NodeIndex minusZ = PeriodicBoundaryCondition::applyBC(NodeIndex{ gridIndices[i].x, gridIndices[i].y, gridIndices[i].z - 1 }, gridpointsPerDim);

			// Retrieve potentials at grid nodes neighboring the current query
			float phi_plusX = potentialGrid[getGridIndex(plusX, gridpointsPerDim)];
			float phi_minusX = potentialGrid[getGridIndex(minusX, gridpointsPerDim)];
			float phi_plusY = potentialGrid[getGridIndex(plusY, gridpointsPerDim)];
			float phi_minusY = potentialGrid[getGridIndex(minusY, gridpointsPerDim)];
			float phi_plusZ = potentialGrid[getGridIndex(plusZ, gridpointsPerDim)];
			float phi_minusZ = potentialGrid[getGridIndex(minusZ, gridpointsPerDim)];

			// Compute electric field components using central finite differences (E = -grad(phi))
			float E_x = -(phi_plusX - phi_minusX) * (gridpointsPerNm / 2.0f);
			float E_y = -(phi_plusY - phi_minusY) * (gridpointsPerNm / 2.0f);
			float E_z = -(phi_plusZ - phi_minusZ) * (gridpointsPerNm / 2.0f);

			// Accumulate the weighted electric field contributions to compute the particle's force
			force += Float3{ E_x, E_y, E_z } * weights[i];		

			// Accumulate the weighted potential contributions to compute the particle's potential
			potential += phi * weights[i];
		}

		// Store the computed force and potential in the output arrays
		forceEnergies[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = ForceEnergy{ force * charge, potential * charge };
	}

	__global__ void PrecomputeGreensFunctionKernel(float* d_greensFunction, int gridnodesPerDim, float boxLen, float ewaldKappa) {
		int halfNodes = gridnodesPerDim / 2;
		int reciprocalComplexCountZ = halfNodes + 1; // Only need N/2+1 points in Z due to Hermitian symmetry.

		// Each dimension (X, Y, Z) corresponds to a frequency component in reciprocal space.
		const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
		NodeIndex freqIndex = BoxGrid::Get3dIndexWithNNodes(index1D, gridnodesPerDim);

		if (freqIndex.x >= gridnodesPerDim || freqIndex.y >= gridnodesPerDim || freqIndex.z >= reciprocalComplexCountZ) {
			return;
		}

		// Compute the wavevector indices.Frequencies above N / 2 represent negative frequencies.
		// This re-mapping ensures that kx and ky run from negative to positive frequencies symmetrically.
		int kxShiftedIndex = (freqIndex.x <= halfNodes) ? freqIndex.x : freqIndex.x - gridnodesPerDim;
		int kyShiftedIndex = (freqIndex.y <= halfNodes) ? freqIndex.y : freqIndex.y - gridnodesPerDim;
		int kzIndex = freqIndex.z;

		// Convert these index-based frequencies to actual wavevector values (in nm^-1).
		// Multiplying by 2pi and dividing by the box length maps the discrete index to a physical wavevector.
		float kx = (2.0f * PI * (float)kxShiftedIndex) / boxLen;
		float ky = (2.0f * PI * (float)kyShiftedIndex) / boxLen;
		float kz = (2.0f * PI * (float)kzIndex) / boxLen;

		// Compute k^2.This determines how strong the damping and division by k^2 in the Green’s function will be.
		float kSquared = kx * kx + ky * ky + kz * kz;

		// Initialize G(k) to zero. If k=0, G(k) should remain zero to avoid division by zero.
		float currentGreensValue = 0.0f;
		if (kSquared > 1e-14f) { // Avoid division by zero at k=0
			//double exponentialFactor = exp(-kSquared / (4.0f * ewaldKappa * ewaldKappa)); 
			//printf("ef %f k2 %f\n", exponentialFactor, kSquared);
			//currentGreensValue = (4.0f * PI * exponentialFactor) / kSquared;
			currentGreensValue = exp(-kSquared / (4.0f * ewaldKappa * ewaldKappa)) / kSquared;
		}

		// Yes this is correct here even tho zdim is only half, because it is still the major
		d_greensFunction[index1D] = currentGreensValue;
	}

	__global__ void ApplyGreensFunctionKernel(cufftComplex* d_reciprocalFreqData, const float* d_greensFunctionArray, int realSpaceNodes, float boxLength)
	{
		int halfNodes = realSpaceNodes / 2;
		int reciprocalComplexCountZ = halfNodes + 1; // Only half + 1 needed due to Hermitian symmetry.

		const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
		NodeIndex freqIndex = BoxGrid::Get3dIndex(index1D, realSpaceNodes);

		// If the thread’s index falls outside the valid range, we do nothing.
		if (freqIndex.x >= realSpaceNodes || freqIndex.y >= realSpaceNodes || freqIndex.z >= reciprocalComplexCountZ) {
			return; // Avoid out-of-bounds memory access.
		}

		// Retrieve the precomputed Green’s function value at this frequency point.
		const float currentGreenValue = d_greensFunctionArray[index1D];

		// The frequency-domain data is complex, and we multiply by the real Green’s function to apply the long-range correction.
		cufftComplex currentFreqValue = d_reciprocalFreqData[index1D];

		/*if (currentFreqValue.x != 0.f || currentFreqValue.y != 0.f)
			printf("currentFreqValue %f %f greenvalue %f\n", currentFreqValue.x, currentFreqValue.y, currentGreenValue);*/

		// Scaling the complex frequency data by the Green’s function imposes the correct electrostatic filtering.
		currentFreqValue.x *= currentGreenValue;
		currentFreqValue.y *= currentGreenValue;

		d_reciprocalFreqData[index1D] = currentFreqValue;
	}

	__global__ void NormalizePotentialKernel(float* const d_potential, float totalElements) {
		int idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < totalElements) {
			d_potential[idx] = d_potential[idx] / totalElements;
			//d_potential[idx] *= PhysicsUtilsDevice::modifiedCoulombConstant_Potential;
		}
	}





	class Controller {

		float* chargeGrid;
		cufftComplex* potentialGrid;
		float* greensFunctionScalars;

		int gridpointsPerDim=-1;
		int nGridpointsRealspace = -1;
		int nGridpointsReciprocalspace = -1;
		float kappa = 1.f/1.2; // cutoff between LR and SR. TODO: should be the same as LR cutoff from params
		float boxlenNm{};

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
	


	public:

		Controller(int boxLenNm) :boxlenNm(boxLenNm){
			gridpointsPerDim = boxLenNm * gridpointsPerNm;
			nGridpointsRealspace = gridpointsPerDim * gridpointsPerDim * gridpointsPerDim;
			nGridpointsReciprocalspace = gridpointsPerDim * gridpointsPerDim * (gridpointsPerDim / 2 + 1);

			if (nGridpointsRealspace > INT32_MAX)
				throw std::runtime_error("Ewald grid too large to index with integers");

			const size_t byteSize = nGridpointsRealspace * sizeof(float)*3;
			if (byteSize > 8'000'000'000)
				throw std::runtime_error("Ewald grid too large");

			cudaMalloc(&chargeGrid, nGridpointsRealspace * sizeof(float));
			cudaMalloc(&potentialGrid, nGridpointsReciprocalspace * sizeof(cufftComplex));
			cudaMalloc(&greensFunctionScalars, nGridpointsRealspace * sizeof(float));
			cudaDeviceSynchronize();

			const int nBlocks = (nGridpointsReciprocalspace + 63) / 64;
			PrecomputeGreensFunctionKernel<<<nBlocks, 64 >>>(greensFunctionScalars, gridpointsPerDim, boxLenNm, kappa);
			LIMA_UTILS::genericErrorCheck("PrecomputeGreensFunctionKernel failed!");

		}
		~Controller() {
			cudaFree(chargeGrid);
			cudaFree(potentialGrid);
			cudaFree(greensFunctionScalars);
		}

		void CalcCharges(const BoxConfig& config, const BoxState& state, int nCompounds, ForceEnergy* const forceEnergy) {
			cudaMemset(chargeGrid, 0, gridpointsPerDim * gridpointsPerDim * gridpointsPerDim * sizeof(float));// maybe do this async, after the last kernel?
			cudaDeviceSynchronize();
			DistributeChargesToGridKernel<<<nCompounds, MAX_COMPOUND_PARTICLES>> > (config, state, chargeGrid, gridpointsPerDim);
			LIMA_UTILS::genericErrorCheckNoSync("DistributeChargesToGridKernel failed!");

			ForwardFFT();
			cudaDeviceSynchronize();

			ApplyGreensFunctionKernel<<<(nGridpointsReciprocalspace + 63) / 64, 64 >> > (potentialGrid, greensFunctionScalars, gridpointsPerDim, boxlenNm);
			LIMA_UTILS::genericErrorCheckNoSync("ApplyGreensFunctionKernel failed!");

			InverseFFT();
			cudaDeviceSynchronize();

			NormalizePotentialKernel<<<(nGridpointsRealspace + 63) / 64, 64 >> > (chargeGrid, static_cast<float>(nGridpointsRealspace));
			LIMA_UTILS::genericErrorCheckNoSync("NormalizePotentialKernel failed!");

			InterpolateForcesAndPotentialKernel << <nCompounds, MAX_COMPOUND_PARTICLES >> > (config, state, chargeGrid, gridpointsPerDim, forceEnergy);
			LIMA_UTILS::genericErrorCheckNoSync("InterpolateForcesAndPotentialKernel failed!");
		}
	};

}