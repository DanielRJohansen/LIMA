

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"
#include "BoundaryCondition.cuh"

#include "Filehandling.h"
#include "PhysicsUtils.cuh"

#include <cuda_runtime.h>
#include <cufft.h>




namespace PME {
	const int gridpointsPerNm = 20;
	constexpr float gridpointsPerNm_f = static_cast<float>(gridpointsPerNm);

	//__device__ __host__ int getGridIndex(const Int3& gridIndex, int gridpointsPerDim) {
	//	return gridIndex.x + gridIndex.y * gridpointsPerDim + gridIndex.z * gridpointsPerDim * gridpointsPerDim;
	//}


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

		const float charge = config.compoundsAtomCharges[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
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

		//printf("%f %f %f\n", fx, fy, fz);
		
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
			int X = (ix + dx + gridpointsPerDim) % gridpointsPerDim;
			float wxCur = wx[dx];
			for (int dy = 0; dy < 4; dy++) {
				int Y = (iy + dy + gridpointsPerDim) % gridpointsPerDim;
				float wxyCur = wxCur * wy[dy];
				for (int dz = 0; dz < 4; dz++) {
					int Z = (iz + dz + gridpointsPerDim) % gridpointsPerDim;
					const int index1D = GetGridIndexRealspace(Int3{ X,Y,Z }, gridpointsPerDim);
					chargeGrid[index1D] = charge * wxyCur * wz[dz]; // TODO This must be applyBC

					//atomicAdd(&chargeGrid[getGridIndex(Int3{ X,Y,Z }, gridpointsPerDim)], charge * wxyCur * wz[dz]);
				}
			}
		}


	}


	__global__ void InterpolateForcesAndPotentialKernel1(const BoxConfig config,
		const BoxState state,
		const float* __restrict__ potentialGrid,
		int gridpointsPerDim,
		ForceEnergy* __restrict__ forceEnergies)
	{
		int particleIdx = threadIdx.x;
		int compoundIdx = blockIdx.x;
		if (particleIdx >= config.compounds[compoundIdx].n_particles)
			return;

		int outIdx = compoundIdx * MAX_COMPOUND_PARTICLES + particleIdx;

		// Compute the particle's absolute position
		NodeIndex compoundOrigin = state.compoundOrigos[compoundIdx];
		Float3 relativePos = state.compoundsRelposNm[compoundIdx * MAX_COMPOUND_PARTICLES + particleIdx];
		Float3 absolutePos = relativePos + compoundOrigin.toFloat3();
		PeriodicBoundaryCondition::applyBCNM(absolutePos);

		// Map to fractional grid coordinates
		Float3 gridPos = absolutePos * gridpointsPerNm_f;
		int ix = (int)floorf(gridPos.x);
		int iy = (int)floorf(gridPos.y);
		int iz = (int)floorf(gridPos.z);

		float fx = gridPos.x - ix;
		float fy = gridPos.y - iy;
		float fz = gridPos.z - iz;

		// B-spline weights and their derivatives for x
		float w0x = ((1.f - fx) * (1.f - fx) * (1.f - fx)) / 6.f;
		float w1x = (4.f - 6.f * fx * fx + 3.f * fx * fx * fx) / 6.f;
		float w2x = (1.f + 3.f * fx + 3.f * fx * fx - 3.f * fx * fx * fx) / 6.f;
		float w3x = (fx * fx * fx) / 6.f;

		float dw0x = -((1.f - fx) * (1.f - fx)) / 2.f;                      // d/dx of w0x
		float dw1x = (-12.f * fx + 9.f * fx * fx) / 6.f;                   // d/dx of w1x
		float dw2x = (3.f + 6.f * fx - 9.f * fx * fx) / 6.f;               // d/dx of w2x
		float dw3x = (3.f * fx * fx) / 6.f;                                // d/dx of w3x

		float wx[4] = { w0x, w1x, w2x, w3x };
		float dwx[4] = { dw0x, dw1x, dw2x, dw3x };

		// B-spline weights and their derivatives for y
		float w0y = ((1.f - fy) * (1.f - fy) * (1.f - fy)) / 6.f;
		float w1y = (4.f - 6.f * fy * fy + 3.f * fy * fy * fy) / 6.f;
		float w2y = (1.f + 3.f * fy + 3.f * fy * fy - 3.f * fy * fy * fy) / 6.f;
		float w3y = (fy * fy * fy) / 6.f;

		float dw0y = -((1.f - fy) * (1.f - fy)) / 2.f;
		float dw1y = (-12.f * fy + 9.f * fy * fy) / 6.f;
		float dw2y = (3.f + 6.f * fy - 9.f * fy * fy) / 6.f;
		float dw3y = (3.f * fy * fy) / 6.f;

		float wy[4] = { w0y, w1y, w2y, w3y };
		float dwy[4] = { dw0y, dw1y, dw2y, dw3y };

		// B-spline weights and their derivatives for z
		float w0z = ((1.f - fz) * (1.f - fz) * (1.f - fz)) / 6.f;
		float w1z = (4.f - 6.f * fz * fz + 3.f * fz * fz * fz) / 6.f;
		float w2z = (1.f + 3.f * fz + 3.f * fz * fz - 3.f * fz * fz * fz) / 6.f;
		float w3z = (fz * fz * fz) / 6.f;

		float dw0z = -((1.f - fz) * (1.f - fz)) / 2.f;
		float dw1z = (-12.f * fz + 9.f * fz * fz) / 6.f;
		float dw2z = (3.f + 6.f * fz - 9.f * fz * fz) / 6.f;
		float dw3z = (3.f * fz * fz) / 6.f;

		float wz[4] = { w0z, w1z, w2z, w3z };
		float dwz[4] = { dw0z, dw1z, dw2z, dw3z };

		// Interpolate potential and its gradient
		float potential = 0.f;
		float dpx = 0.f;
		float dpy = 0.f;
		float dpz = 0.f;

		for (int dx = 0; dx < 4; dx++) {
			int gx = (ix + dx + gridpointsPerDim) % gridpointsPerDim;
			float wxCur = wx[dx];
			float dwxCur = dwx[dx];
			for (int dy = 0; dy < 4; dy++) {
				int gy = (iy + dy + gridpointsPerDim) % gridpointsPerDim;
				float wyCur = wy[dy];
				float dwyCur = dwy[dy];
				for (int dz = 0; dz < 4; dz++) {
					int gz = (iz + dz + gridpointsPerDim) % gridpointsPerDim;
					float val = potentialGrid[GetGridIndexRealspace(Int3{ gx, gy, gz }, gridpointsPerDim)];

					float wzCur = wz[dz];
					float dwzCur = dwz[dz];

					// Accumulate potential
					potential += val * wxCur * wyCur * wzCur;

					// Accumulate gradients
					dpx += val * dwxCur * wyCur * wzCur;
					dpy += val * wxCur * dwyCur * wzCur;
					dpz += val * wxCur * wyCur * dwzCur;
				}
			}
		}

		// Force = -grad(potential)
		Float3 force = Float3(-dpx, -dpy, -dpz);

		/*forceEnergies[outIdx].force = force;
		forceEnergies[outIdx].potE = potential;*/
		forceEnergies[outIdx] = ForceEnergy(force, potential);
		// TODO: Apply correct unit scaling or additional constants if needed
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


		/*if (blockIdx.x != 0)
			return;*/

//		printf("gridnode000 %d %d %d\n", gridnode000.x, gridnode000.y, gridnode000.z);

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
			//int gridIndex = getGridIndex(gridIndices[i], gridpointsPerDim); // WHY: Converts 3D grid indices to 1D index for memory access
			const int gridIndex = GetGridIndexRealspace(gridIndices[i], gridpointsPerDim);


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
			float phi_plusX = potentialGrid[GetGridIndexRealspace(plusX, gridpointsPerDim)];
			float phi_minusX = potentialGrid[GetGridIndexRealspace(minusX, gridpointsPerDim)];
			float phi_plusY = potentialGrid[GetGridIndexRealspace(plusY, gridpointsPerDim)];
			float phi_minusY = potentialGrid[GetGridIndexRealspace(minusY, gridpointsPerDim)];
			float phi_plusZ = potentialGrid[GetGridIndexRealspace(plusZ, gridpointsPerDim)];
			float phi_minusZ = potentialGrid[GetGridIndexRealspace(minusZ, gridpointsPerDim)];

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

	__global__ void PrecomputeGreensFunctionKernel(float* d_greensFunction, int gridpointsPerDim, 
		float boxLen,		// [nm]
		float ewaldKappa	// [nm^-1]
	) {

		const int halfNodes = gridpointsPerDim / 2;
		int nGridpointsHalfdim = gridpointsPerDim / 2 + 1;// TODO: add comment here
		const int index1D = blockIdx.x * blockDim.x + threadIdx.x;
		if (index1D > nGridpointsHalfdim * gridpointsPerDim * gridpointsPerDim)
			return;
		NodeIndex freqIndex = Get3dIndexReciprocalspace(index1D, gridpointsPerDim, nGridpointsHalfdim);

		// Remap frequencies to negative for indices > N/2
		/*int kxShiftedIndex = (freqIndex.x <= halfNodes) ? freqIndex.x : freqIndex.x - gridpointsPerDim;
		int kyShiftedIndex = (freqIndex.y <= halfNodes) ? freqIndex.y : freqIndex.y - gridpointsPerDim;
		int kzIndex = freqIndex.z;*/

		int kxIndex = freqIndex.x;
		int kyShiftedIndex = (freqIndex.y <= halfNodes) ? freqIndex.y : freqIndex.y - gridpointsPerDim;
		int kzShiftedIndex = (freqIndex.z <= halfNodes) ? freqIndex.z : freqIndex.z - gridpointsPerDim;

		// Ewald kappa fixed
		float volume = boxLen * boxLen * boxLen;			// [nm^3]
		float delta = boxLen / (float)gridpointsPerDim;		// [nm]

		// Physical wavevectors
		float kx = (2.0f * PI * (float)kxIndex) / boxLen;
		float ky = (2.0f * PI * (float)kyShiftedIndex) / boxLen;
		float kz = (2.0f * PI * (float)kzShiftedIndex) / boxLen;


		float kSquared = kx * kx + ky * ky + kz * kz;

		float currentGreensValue = 0.0f;

		// Compute B-spline structure factor (4th order)
		float kHalfX = kx * (delta * 0.5f);
		float kHalfY = ky * (delta * 0.5f);
		float kHalfZ = kz * (delta * 0.5f);

		auto splineFactor = [](float kh) {
			if (fabs(kh) < 1e-14f) return 1.0f;
			float ratio = __sinf(kh) / kh;
			return powf(ratio, 4);
			};

		float Sx = splineFactor(kHalfX);
		float Sy = splineFactor(kHalfY);
		float Sz = splineFactor(kHalfZ);
		
		float splineCorrection = (Sx * Sy * Sz);
		splineCorrection = splineCorrection * splineCorrection; // squared for forward+back interpolation

		//printf("spC %f\n", splineCorrection);

		if (kSquared > 1e-14f) {
			currentGreensValue = (4.0f * PI / (volume * kSquared)) 
				* expf(-kSquared / (4.0f * ewaldKappa * ewaldKappa)) 
				* splineCorrection
				* PhysicsUtilsDevice::modifiedCoulombConstant_Force
				;
		}

		d_greensFunction[index1D] = currentGreensValue;
	}


	__global__ void ApplyGreensFunctionKernel(cufftComplex* d_reciprocalFreqData, const float* d_greensFunctionArray, int gridpointsPerDim, float boxLength)
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

	__global__ void NormalizePotentialKernel(float* const d_potential, float totalElements, int nGridnodesPerDim) {
		int idx = blockIdx.x * blockDim.x + threadIdx.x;
		NodeIndex idx3 = BoxGrid::Get3dIndexWithNNodes(idx, nGridnodesPerDim);
		if (idx < totalElements) {
			//d_potential[idx] = d_potential[idx] / totalElements;
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
			cudaMalloc(&greensFunctionScalars, nGridpointsReciprocalspace * sizeof(float));
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

			NormalizePotentialKernel<<<(nGridpointsRealspace + 63) / 64, 64 >> > (chargeGrid, static_cast<float>(nGridpointsRealspace), gridpointsPerDim);
			LIMA_UTILS::genericErrorCheckNoSync("NormalizePotentialKernel failed!");
			cudaDeviceSynchronize();

		
			{
				//std::vector<float> pots;
				//GenericCopyToHost(chargeGrid, pots, nGridpointsRealspace);

				//int centerSlice = 200;
				//int numSlices = 2;
				//int spacing = 20;

				//std::vector<float> combinedData;
				//std::vector<int> sliceIndices;

				//for (int i = -numSlices; i <= numSlices; ++i) {
				//	int sliceIndex = centerSlice + i * spacing;
				//	if (sliceIndex < 0 || sliceIndex >= gridpointsPerDim) {
				//		throw std::runtime_error("error");
				//	}

				//	int firstIndex = gridpointsPerDim * gridpointsPerDim * sliceIndex;
				//	int lastIndex = firstIndex + gridpointsPerDim * gridpointsPerDim;
				//	combinedData.insert(combinedData.end(), pots.begin() + firstIndex, pots.begin() + lastIndex);
				//	sliceIndices.push_back(sliceIndex);
				//}

				//// Save all slices to one file
				//FileUtils::WriteVectorToBinaryFile("C:/Users/Daniel/git_repo/LIMA_data/Pool/PmePot_AllSlices.bin", combinedData);

				//// Call the Python script to plot the slices
				//std::string pyscriptPath = (FileUtils::GetLimaDir() / "dev" / "PyTools" / "Plot2dVec.py").string();
				//std::string command = "python " + pyscriptPath + " " + std::to_string(numSlices * 2 + 1) + " " + std::to_string(gridpointsPerDim) + " " + std::to_string(centerSlice) + " " + std::to_string(spacing);
				//std::system(command.c_str());
			}

			//const int firstIndex = gridpointsPerDim * gridpointsPerDim * 200;
			//const int lastIndex = firstIndex + gridpointsPerDim * gridpointsPerDim;
			//std::vector<float> pots;
			//GenericCopyToHost(chargeGrid, pots, nGridpointsRealspace);
			//std::vector<float> data = std::vector<float>(pots.begin() + firstIndex, pots.begin() + lastIndex);
			//FileUtils::WriteVectorToBinaryFile(R"(C:\Users\Daniel\git_repo\LIMA_data\Pool\PmePot.bin)", data);

			{
				//int i = 0;
				//const int y = 92;
				//for (int z = 0; z < gridpointsPerDim; z++) {
				//	for (int x = 0; x < gridpointsPerDim; x++) {
				//		const int index1D = GetGridIndexRealspace(Int3{ x,y,z }, gridpointsPerDim);
				//		data[i] = pots[index1D];
				//	}
				//}
				//FileUtils::WriteVectorToBinaryFile(R"(C:\Users\Daniel\git_repo\LIMA_data\Pool\PmePot.bin)", data);
			}

			InterpolateForcesAndPotentialKernel << <nCompounds, MAX_COMPOUND_PARTICLES >> > (config, state, chargeGrid, gridpointsPerDim, forceEnergy);
			LIMA_UTILS::genericErrorCheckNoSync("InterpolateForcesAndPotentialKernel failed!");
		}
	};

}