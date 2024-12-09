

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





	__device__ __host__ float Spline(float x, int p) {
		if ((x <= 0.f) || x >= static_cast<float>(p))
			return 0.f;
		if (p==2)
			return 1.f - fabs(x - 1.f);
		return (x * spline(x, p - 1) + (p - x) * spline(x - 1.f, p - 1)) / static_cast<float>(p - 1);
	}

	__device__ __host__ float BCoeff(int p, int K, int k) {
		if ((p % 2 == 1) && (2 * abs(k) == K)) return 0.f;
		float c = 0.f;
		float s = 0.f;
		for (int q = 0; q <= p - 2; q++) {
			c += spline(q + 1.f, p) * cos((2.f * PI * k * q) / K);
			s += spline(q + 1.f, p) * sin((2.f * PI * k * q) / K);
		}
		return 1.f / (c*c + s*s);
	}

	__global__ void InitFFTGrid(int K, float* const D, int pmax, Float3 volume) {
		const int idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx > K * K * K)
			return;
		const NodeIndex k = BoxGrid::Get3dIndex(idx, K) - NodeIndex{ K - 1, K - 1, K - 1 };

		Int3 kp{
			(k.x + K) / K,
			(k.y + K) / K,
			(k.z + K) / K
		};
		const float m = (k.x / volume.x * k.x / volume.x) + (k.y / volume.y * k.y / volume.y) + (k.z / volume.z * k.z / volume.z);

		const int index = BoxGrid::Get3dIndex(idx, K);
		if (m > 0.f) {
			const float skal = 1.f;
			constfloat G = 1.f;
			D[index] = exp(-m*PI/G* PI / G)* skal /
				(m*PI**volume.x*volume.y*volume.z) *
				BCoeff(pmax, K, k.x)*
				BCoeff(pmax, K, k.y)*
				BCoeff(pmax, K, k.z);
		}
		else {
			D[index] = 0.f;
		}
	}

	__global__ void ApplyD((cufftComplex* Q, float* D) {
		const int idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx > K * K * K)
			return;
		
		Q[idx].x *= D[idx];
		Q[idx].y *= D[idx];
	}



	__global__ void InterpolateToMesh(const BoxConfig config, const BoxState state, float* chargeGrid, int gridpointsPerDim) {
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

	}

	void temp(float D*, int* K, int pmax, Float3 l) {
		NodeIndex m = 
	}




	class Controller {

		float* chargeGrid;
		cufftComplex* potentialGrid;
		float* greensFunctionScalars;

		int gridpointsPerDim = -1;
		int nGridpointsRealspace = -1;
		int nGridpointsReciprocalspace = -1;
		float kappa = 1.f / 1.2; // cutoff between LR and SR. TODO: should be the same as LR cutoff from params
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

		Controller(int boxLenNm) :boxlenNm(boxLenNm) {
			gridpointsPerDim = boxLenNm * gridpointsPerNm;
			nGridpointsRealspace = gridpointsPerDim * gridpointsPerDim * gridpointsPerDim;
			nGridpointsReciprocalspace = gridpointsPerDim * gridpointsPerDim * (gridpointsPerDim / 2 + 1);

			if (nGridpointsRealspace > INT32_MAX)
				throw std::runtime_error("Ewald grid too large to index with integers");

			const size_t byteSize = nGridpointsRealspace * sizeof(float) * 3;
			if (byteSize > 8'000'000'000)
				throw std::runtime_error("Ewald grid too large");

			cudaMalloc(&chargeGrid, nGridpointsRealspace * sizeof(float));
			cudaMalloc(&potentialGrid, nGridpointsReciprocalspace * sizeof(cufftComplex));
			cudaMalloc(&greensFunctionScalars, nGridpointsRealspace * sizeof(float));
			cudaDeviceSynchronize();

			const int nBlocks = (nGridpointsReciprocalspace + 63) / 64;
			PrecomputeGreensFunctionKernel << <nBlocks, 64 >> > (greensFunctionScalars, gridpointsPerDim, boxLenNm, kappa);
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
			DistributeChargesToGridKernel << <nCompounds, MAX_COMPOUND_PARTICLES >> > (config, state, chargeGrid, gridpointsPerDim);
			LIMA_UTILS::genericErrorCheckNoSync("DistributeChargesToGridKernel failed!");

			ForwardFFT();
			cudaDeviceSynchronize();

			ApplyGreensFunctionKernel << <(nGridpointsReciprocalspace + 63) / 64, 64 >> > (potentialGrid, greensFunctionScalars, gridpointsPerDim, boxlenNm);
			LIMA_UTILS::genericErrorCheckNoSync("ApplyGreensFunctionKernel failed!");

			InverseFFT();
			cudaDeviceSynchronize();

			NormalizePotentialKernel << <(nGridpointsRealspace + 63) / 64, 64 >> > (chargeGrid, static_cast<float>(nGridpointsRealspace));
			LIMA_UTILS::genericErrorCheckNoSync("NormalizePotentialKernel failed!");

			InterpolateForcesAndPotentialKernel << <nCompounds, MAX_COMPOUND_PARTICLES >> > (config, state, chargeGrid, gridpointsPerDim, forceEnergy);
			LIMA_UTILS::genericErrorCheckNoSync("InterpolateForcesAndPotentialKernel failed!");
		}
	};

}