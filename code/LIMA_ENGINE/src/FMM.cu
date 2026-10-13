#include "FMM.cuh"
#include "SimulationData.h"
#include "PhysicsUtils.cuh"
#include "Utilities.h"
#include <cub/device/device_radix_sort.cuh>
#include <mma.h>
#include <cuda_pipeline.h>
#include <cuda_fp16.h>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <type_traits>

// Tensor Chebyshev interpolation FMM (P2M, M2M, M2L, L2L, L2P and direct P2P).
// The 27 nearest root images are handled by the ordinary tree. A smooth periodic
// root operator represents all remaining images with conducting boundary conditions
// and a uniform neutralizing background. Ewald is used only to construct this small,
// position-independent operator; there is no mesh, FFT or reciprocal particle sum.
namespace FMM {
namespace {
	constexpr int maxOrder = 8;
	constexpr double pi = std::numbers::pi;
	struct Wave { double x, y, z, weight; };
	template<typename Real>
	struct CompensatedSum {
		Real value = 0., correction = 0.;
		__device__ void Add(Real contribution) {
			if constexpr (std::is_same_v<Real, float>) {
				// Preserve compensation even when the rest of the engine uses fast math.
				const float adjusted = __fsub_rn(contribution, correction);
				const float next = __fadd_rn(value, adjusted);
				correction = __fsub_rn(__fsub_rn(next, value), adjusted);
				value = next;
			}
			else value += contribution;
		}
	};

	struct Interpolation {
		int order, count;
		double nodes[maxOrder], invDenominator[maxOrder];
		template<int p, bool withDerivative, typename Real>
		__device__ void Basis(double x, Real* value, Real* derivative = nullptr) const {
			#pragma unroll
			for (int i = 0; i < p; ++i) {
				Real v = invDenominator[i], d = 0.;
				#pragma unroll
				for (int j = 0; j < p; ++j) if (i != j) {
					const Real distance = Real(x) - Real(nodes[j]);
					if constexpr (withDerivative) d = d * distance + v;
					v *= distance;
				}
				value[i] = v;
				if constexpr (withDerivative) derivative[i] = d;
			}
		}
		__device__ double3 Node(int i) const {
			return {nodes[i % order], nodes[(i / order) % order], nodes[i / (order * order)]};
		}
	};

	__host__ __device__ int Index(int x, int y, int z, int n) { return (z * n + y) * n + x; }
	__device__ int3 Coordinates(int i, int n) { return {i % n, (i / n) % n, i / (n * n)}; }
	// Every tree level has a power-of-two side length, including periodic images.
	__device__ int Wrap(int i, int n) { return i & (n - 1); }
	__device__ double Squared(double3 v) { return v.x * v.x + v.y * v.y + v.z * v.z; }

	__global__ void InitializeOperators(Interpolation interpolation, Float3 box, float* translations) {
		const int k = interpolation.count;
		const size_t id = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
		if (id >= size_t(343) * k * k) return;
		const int target = id % k, source = (id / k) % k, offset = id / (k * k);
		const double3 a = interpolation.Node(target), b = interpolation.Node(source);
		const int3 c = Coordinates(offset, 7);
		const double3 r{box.x * (c.x - 3 + (a.x - b.x) * .5),
			box.y * (c.y - 3 + (a.y - b.y) * .5), box.z * (c.z - 3 + (a.z - b.z) * .5)};
		translations[id] = (abs(c.x - 3) > 1 || abs(c.y - 3) > 1 || abs(c.z - 3) > 1)
			? static_cast<float>(1. / sqrt(Squared(r))) : 0.f;
	}

	__global__ void InitializePeriodicOperator(Interpolation interpolation, Float3 box, const Wave* waves,
		int nWaves, double alpha, float* root) {
		const int k = interpolation.count, id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= k * k) return;
		const double3 a = interpolation.Node(id % k), b = interpolation.Node(id / k);
		const double3 r{(a.x - b.x) * box.x * .5, (a.y - b.y) * box.y * .5, (a.z - b.z) * box.z * .5};
		double value = -pi / (alpha * alpha * double(box.x) * box.y * box.z);
		for (int z = -4; z <= 4; ++z) for (int y = -4; y <= 4; ++y) for (int x = -4; x <= 4; ++x) {
			const double3 d{r.x + x * box.x, r.y + y * box.y, r.z + z * box.z};
			const double radius = sqrt(Squared(d));
			if (abs(x) <= 1 && abs(y) <= 1 && abs(z) <= 1)
				value -= radius == 0. ? 2. * alpha / sqrt(pi) : erf(alpha * radius) / radius;
			else if (alpha * radius < 7.) value += erfc(alpha * radius) / radius;
		}
		for (int i = 0; i < nWaves; ++i) {
			const Wave w = waves[i];
			value += w.weight * cos(w.x * r.x + w.y * r.y + w.z * r.z);
		}
		root[id] = value;
	}

	__global__ void InitializeParticles(unsigned int* keys, int* values, int* slots, int count) {
		const int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= count) return;
		keys[i] = UINT_MAX;
		values[i] = i;
		slots[i] = -1;
	}

	__global__ void OperatorCovariance(int k, const float* matrices, double* covariance) {
		const int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= k * k) return;
		const int row = id / k, column = id % k;
		if (column > row) return;
		double sum = 0.;
		for (int j = 0; j < 343 * k; ++j) sum += double(matrices[size_t(j) * k + row]) * matrices[size_t(j) * k + column];
		covariance[row * k + column] = covariance[column * k + row] = sum;
	}

	__global__ void ProjectOperators(int k, int rank, const float* basis, const float* matrices, float* projected) {
		const size_t id = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
		if (id >= size_t(343) * k * rank) return;
		const int component = id % rank, source = (id / rank) % k, offset = id / (rank * k);
		double sum = 0.;
		for (int j = 0; j < k; ++j) sum += double(basis[j * rank + component]) * matrices[(size_t(offset) * k + source) * k + j];
		projected[id] = sum;
	}

	__global__ void CompressOperators(int k, int rank, const float* basis, const float* projected, float* compressed) {
		const size_t id = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
		if (id >= size_t(343) * rank * rank) return;
		const int target = id % rank, source = (id / rank) % rank, offset = id / (rank * rank);
		double sum = 0.;
		for (int j = 0; j < k; ++j) sum += double(projected[(size_t(offset) * k + j) * rank + target]) * basis[j * rank + source];
		compressed[id] = sum;
	}

	__global__ void PackTranslations(int rank, float scale, const float* matrix, __half* high, __half* low) {
		const int stride = rank == 96 ? 96 : ((rank + 63) / 64) * 64;
		const size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
		if (i >= size_t(343) * rank * stride) return;
		const int row = i % stride;
		const float value = row < rank ? matrix[(i / stride) * rank + row] * scale : 0.f;
		high[i] = __float2half_rn(value);
		low[i] = __float2half_rn((value - __half2float(high[i])) * 4096.f);
	}

	__global__ void PackFactors(int count, const float* values, __half* high, __half* low) {
		const int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= count) return;
		high[i] = __float2half_rn(values[i]);
		low[i] = __float2half_rn((values[i] - __half2float(high[i])) * 4096.f);
	}

	__global__ void MaximumMoment(int count, const float* moments, unsigned int* maximum) {
		__shared__ float warpMaxima[8];
		float value = 0.f;
		for (int i = blockIdx.x * 1024 + threadIdx.x; i < min(count, int(blockIdx.x + 1) * 1024); i += 256)
			value = fmaxf(value, fabsf(moments[size_t(blockIdx.y) * count + i]));
		for (int distance = 16; distance > 0; distance /= 2) value = fmaxf(value, __shfl_down_sync(0xffffffff, value, distance));
		if (threadIdx.x % 32 == 0) warpMaxima[threadIdx.x / 32] = value;
		__syncthreads();
		if (threadIdx.x < 32) {
			value = threadIdx.x < 8 ? warpMaxima[threadIdx.x] : 0.f;
			for (int distance = 16; distance > 0; distance /= 2) value = fmaxf(value, __shfl_down_sync(0xffffffff, value, distance));
			// Positive IEEE floats have the same order as their unsigned bits.
			if (threadIdx.x == 0) atomicMax(maximum + blockIdx.y, __float_as_uint(value));
		}
	}

	__global__ void PackMoments(int count, const float* moments, const unsigned int* maximum, __half* high, __half* low) {
		const int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= count) return;
		int exponent = 0;
		frexpf(__uint_as_float(maximum[blockIdx.y]), &exponent);
		const size_t index = size_t(blockIdx.y) * count + i;
		const float value = ldexpf(moments[index], -exponent);
		high[index] = __float2half_rn(value);
		low[index] = __float2half_rn((value - __half2float(high[index])) * 4096.f);
	}

	__global__ void TransformCoefficients(int inputCount, int outputCount, const float* basis, const float* input, float* output) {
		const int node = threadIdx.x;
		if (node >= outputCount) return;
		float a = 0.f, b = 0.f, c = 0.f, d = 0.f;
		int j = 0;
		for (; j + 3 < inputCount; j += 4) {
			a += basis[j * outputCount + node] * input[size_t(blockIdx.x) * inputCount + j];
			b += basis[(j + 1) * outputCount + node] * input[size_t(blockIdx.x) * inputCount + j + 1];
			c += basis[(j + 2) * outputCount + node] * input[size_t(blockIdx.x) * inputCount + j + 2];
			d += basis[(j + 3) * outputCount + node] * input[size_t(blockIdx.x) * inputCount + j + 3];
		}
		for (; j < inputCount; ++j) a += basis[j * outputCount + node] * input[size_t(blockIdx.x) * inputCount + j];
		const double sum = double(a) + b + c + d;
		const size_t index = size_t(blockIdx.x) * outputCount + node;
		const float value = sum;
		output[index] = value;
	}

	__global__ void PackBasis(int inputCount, int outputCount, const float* basis, __half* high, __half* low) {
		const int stride = outputCount == 96 ? 96 : ((outputCount + 63) / 64) * 64;
		const int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= ((inputCount + 15) / 16) * 16 * stride) return;
		const int row = i % stride, column = i / stride;
		const float value = row < outputCount && column < inputCount ? basis[column * outputCount + row] : 0.f;
		high[i] = __float2half_rn(value);
		low[i] = __float2half_rn((value - __half2float(high[i])) * 4096.f);
	}

	// Apply a shared basis to 32 cells together, reusing its coefficients through
	// tensor-core matrix products instead of independent matrix-vector products.
	template<int rows>
	__global__ void TransformCoefficientsTensor(int inputCount, int outputCount, int cells, const __half* basisHigh,
		const __half* basisLow, const float* input, float* output) {
#if __CUDA_ARCH__ >= 800
		using namespace nvcuda;
		constexpr int columns = 32, threads = rows * 2, aStride = rows + 16, bStride = 24;
		__shared__ __align__(32) __half aHigh[aStride * 16], aLow[aStride * 16];
		__shared__ __align__(32) __half bHigh[bStride * columns], bLow[bStride * columns];
		__shared__ __align__(32) float result[rows * columns];
		__shared__ float warpMaxima[8];
		__shared__ int exponent;
		const int tid = threadIdx.x, warp = tid / 32, firstCell = blockIdx.x * columns, firstRow = blockIdx.y * rows;
		const int matrixStride = ((outputCount + rows - 1) / rows) * rows;
		float maximum = 0.f;
		for (int i = tid; i < inputCount * columns; i += threads) {
			const int cell = firstCell + i / inputCount;
			if (cell < cells) maximum = fmaxf(maximum, fabsf(input[size_t(cell) * inputCount + i % inputCount]));
		}
		for (int distance = 16; distance > 0; distance /= 2) maximum = fmaxf(maximum, __shfl_down_sync(0xffffffff, maximum, distance));
		if (tid % 32 == 0) warpMaxima[warp] = maximum;
		__syncthreads();
		if (tid < 32) {
			maximum = tid < threads / 32 ? warpMaxima[tid] : 0.f;
			for (int distance = 16; distance > 0; distance /= 2) maximum = fmaxf(maximum, __shfl_down_sync(0xffffffff, maximum, distance));
			if (tid == 0) frexpf(maximum, &exponent);
		}
		__syncthreads();
		wmma::fragment<wmma::accumulator, 16, 16, 16, float> high[columns / 16], low[columns / 16];
		#pragma unroll
		for (int tile = 0; tile < columns / 16; ++tile) {
			wmma::fill_fragment(high[tile], 0.f); wmma::fill_fragment(low[tile], 0.f);
		}
		for (int begin = 0; begin < inputCount; begin += 16) {
			for (int vector = tid; vector < (rows / 8) * 16; vector += threads) {
				const int row = firstRow + (vector % (rows / 8)) * 8, column = begin + vector / (rows / 8);
				const size_t globalIndex = size_t(column) * matrixStride + row;
				const int sharedIndex = (vector / (rows / 8)) * aStride + (vector % (rows / 8)) * 8;
				__pipeline_memcpy_async(aHigh + sharedIndex, basisHigh + globalIndex, 16);
				__pipeline_memcpy_async(aLow + sharedIndex, basisLow + globalIndex, 16);
			}
			for (int i = tid; i < columns * 16; i += threads) {
				const int cell = firstCell + i / 16, component = begin + i % 16;
				const float value = cell < cells && component < inputCount ? ldexpf(input[size_t(cell) * inputCount + component], -exponent) : 0.f;
				const int sharedIndex = (i / 16) * bStride + i % 16;
				bHigh[sharedIndex] = __float2half_rn(value);
				bLow[sharedIndex] = __float2half_rn((value - __half2float(bHigh[sharedIndex])) * 4096.f);
			}
			__pipeline_commit(); __pipeline_wait_prior(0); __syncthreads();
			wmma::fragment<wmma::matrix_a, 16, 16, 16, __half, wmma::col_major> ah, al;
			wmma::fragment<wmma::matrix_b, 16, 16, 16, __half, wmma::col_major> bh, bl;
			wmma::load_matrix_sync(ah, aHigh + warp * 16, aStride);
			wmma::load_matrix_sync(al, aLow + warp * 16, aStride);
			#pragma unroll
			for (int tile = 0; tile < columns / 16; ++tile) {
				wmma::load_matrix_sync(bh, bHigh + tile * 16 * bStride, bStride);
				wmma::load_matrix_sync(bl, bLow + tile * 16 * bStride, bStride);
				wmma::mma_sync(high[tile], ah, bh, high[tile]);
				wmma::mma_sync(low[tile], ah, bl, low[tile]);
				wmma::mma_sync(low[tile], al, bh, low[tile]);
			}
			__syncthreads();
		}
		#pragma unroll
		for (int tile = 0; tile < columns / 16; ++tile) {
			for (int j = 0; j < high[tile].num_elements; ++j)
				high[tile].x[j] = ldexpf(__fadd_rn(high[tile].x[j], low[tile].x[j] * (1.f / 4096.f)), exponent);
			wmma::store_matrix_sync(result + warp * 16 * columns + tile * 16, high[tile], columns, wmma::mem_row_major);
		}
		__syncthreads();
		for (int i = tid; i < rows * columns; i += threads) {
			const int row = firstRow + i / columns, cell = firstCell + i % columns;
			if (row < outputCount && cell < cells) output[size_t(cell) * outputCount + row] = result[i];
		}
#endif
	}

	__global__ void GatherParticles(const SuperCluster* clusters, const SuperClusterMeta* meta, int count,
		Float3 box, int n, const int* simulationSlots, unsigned int* keys, float4* particles, int* slots) {
		const int slot = blockIdx.x * blockDim.x + threadIdx.x, sc = slot / SuperCluster::maxParticles;
		const int lane = slot % SuperCluster::maxParticles;
		if (sc >= count || !clusters[sc].Valid(lane)) return;
		const int simulation = simulationSlots[meta[sc].simulationId];
		if (simulation < 0) return;
		const int i = meta[sc]._pclusterIds[lane] * PersistentCluster::maxParticles + meta[sc].indexInPcluster[lane];
		float4 p = clusters[sc].posCharge[lane];
		p.x -= floorf(p.x / box.x) * box.x;
		p.y -= floorf(p.y / box.y) * box.y;
		p.z -= floorf(p.z / box.z) * box.z;
		const int x = min(n - 1, max(0, int(p.x * n / box.x)));
		const int y = min(n - 1, max(0, int(p.y * n / box.y)));
		const int z = min(n - 1, max(0, int(p.z * n / box.z)));
		keys[i] = simulation * n * n * n + Index(x, y, z, n);
		particles[i] = p;
		slots[i] = slot;
	}

	__global__ void FindLeafRanges(const unsigned int* keys, int count, int2* ranges) {
		const int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= count || keys[i] == UINT_MAX) return;
		const unsigned int key = keys[i];
		if (i == 0 || keys[i - 1] != key) ranges[key].x = i;
		if (i + 1 == count || keys[i + 1] != key) ranges[key].y = i + 1;
	}

	template<int p, typename Real>
	__global__ __launch_bounds__(512) void ParticleToMultipole(Interpolation interpolation, double3 leafScale, int n, const int2* ranges,
		const int* sorted, const float4* particles, float* moments) {
		const int node = threadIdx.x, cell = blockIdx.x;
		constexpr int k = p * p * p;
		// Each particle's three one-dimensional bases are shared by all tensor nodes.
		// Tiling also bounds shared storage independently of leaf occupancy.
		__shared__ Real basis[3][32][maxOrder];
		__shared__ float charges[32];
		const int3 c = Coordinates(cell % (n * n * n), n);
		const int2 range = ranges[cell];
		double sum = 0.;
		for (int begin = range.x; begin < range.y; begin += 32) {
			const int count = min(32, range.y - begin);
			if (node < count) {
				const float4 particle = particles[sorted[begin + node]];
				Real x[p], y[p], z[p];
				interpolation.Basis<p, false>(double(particle.x) * leafScale.x - (2 * c.x + 1), x);
				interpolation.Basis<p, false>(double(particle.y) * leafScale.y - (2 * c.y + 1), y);
				interpolation.Basis<p, false>(double(particle.z) * leafScale.z - (2 * c.z + 1), z);
				#pragma unroll
				for (int a = 0; a < p; ++a) {
					basis[0][node][a] = x[a]; basis[1][node][a] = y[a]; basis[2][node][a] = z[a];
				}
				charges[node] = particle.w;
			}
			__syncthreads();
			Real tile = 0.;
			if (node < k) for (int j = 0; j < count; ++j)
				tile += charges[j] * basis[0][j][node % p] * basis[1][j][(node / p) % p] * basis[2][j][node / (p * p)];
			sum += tile;
			__syncthreads();
		}
		if (node < k) moments[size_t(cell) * k + node] = sum;
	}

	template<int p>
	__global__ void ParticleToMultipoleTensor(Interpolation interpolation, double3 leafScale, int n, const int2* ranges,
		const int* sorted, const float4* particles, float* moments) {
#if __CUDA_ARCH__ >= 800
		using namespace nvcuda;
		constexpr int columns = ((p * p + 15) / 16) * 16, threads = (columns / 16) * 32, basisStride = p | 1;
		__shared__ float basis[3][32][basisStride], charges[32];
		__shared__ __align__(32) __half aHigh[16 * 32], aLow[16 * 32];
		__shared__ __align__(32) __half bHigh[32 * columns], bLow[32 * columns];
		__shared__ __align__(32) float output[16 * columns];
		__shared__ int exponent;
		const int tid = threadIdx.x, warp = tid / 32, cell = blockIdx.x;
		const int3 c = Coordinates(cell % (n * n * n), n);
		const int2 range = ranges[cell];
		using Accumulator = wmma::fragment<wmma::accumulator, 16, 16, 16, float>;
		Accumulator high, low;
		CompensatedSum<float> total[Accumulator::num_elements];
		for (int begin = range.x; begin < range.y; begin += 32) {
			const int count = min(32, range.y - begin);
			float maximum = 0.f;
			if (tid < count) {
				const float4 particle = particles[sorted[begin + tid]];
				float x[p], y[p], z[p];
				interpolation.Basis<p, false>(double(particle.x) * leafScale.x - (2 * c.x + 1), x);
				interpolation.Basis<p, false>(double(particle.y) * leafScale.y - (2 * c.y + 1), y);
				interpolation.Basis<p, false>(double(particle.z) * leafScale.z - (2 * c.z + 1), z);
				#pragma unroll
				for (int a = 0; a < p; ++a) {
					basis[0][tid][a] = x[a]; basis[1][tid][a] = y[a]; basis[2][tid][a] = z[a];
				}
				charges[tid] = particle.w;
				maximum = fabsf(particle.w);
			}
			if (tid < 32) {
				for (int distance = 16; distance > 0; distance /= 2) maximum = fmaxf(maximum, __shfl_down_sync(0xffffffff, maximum, distance));
				if (tid == 0) frexpf(maximum, &exponent);
			}
			__syncthreads();
			for (int i = tid; i < 16 * 32; i += threads) {
				const int particle = i / 16, a = i % 16;
				const float value = particle < count && a < p ? ldexpf(charges[particle], -exponent) * basis[0][particle][a] : 0.f;
				aHigh[i] = __float2half_rn(value);
				aLow[i] = __float2half_rn((value - __half2float(aHigh[i])) * 4096.f);
			}
			for (int i = tid; i < 32 * columns; i += threads) {
				const int particle = i % 32, node = i / 32;
				const float value = particle < count && node < p * p ? basis[1][particle][node % p] * basis[2][particle][node / p] : 0.f;
				bHigh[i] = __float2half_rn(value);
				bLow[i] = __float2half_rn((value - __half2float(bHigh[i])) * 4096.f);
			}
			__syncthreads();
			wmma::fill_fragment(high, 0.f); wmma::fill_fragment(low, 0.f);
			#pragma unroll
			for (int q = 0; q < 32; q += 16) {
				wmma::fragment<wmma::matrix_a, 16, 16, 16, __half, wmma::col_major> ah, al;
				wmma::fragment<wmma::matrix_b, 16, 16, 16, __half, wmma::col_major> bh, bl;
				wmma::load_matrix_sync(ah, aHigh + q * 16, 16);
				wmma::load_matrix_sync(al, aLow + q * 16, 16);
				wmma::load_matrix_sync(bh, bHigh + warp * 16 * 32 + q, 32);
				wmma::load_matrix_sync(bl, bLow + warp * 16 * 32 + q, 32);
				wmma::mma_sync(high, ah, bh, high);
				wmma::mma_sync(low, ah, bl, low);
				wmma::mma_sync(low, al, bh, low);
			}
			#pragma unroll
			for (int j = 0; j < Accumulator::num_elements; ++j)
				total[j].Add(ldexpf(__fadd_rn(high.x[j], low.x[j] * (1.f / 4096.f)), exponent));
			__syncthreads();
		}
		#pragma unroll
		for (int j = 0; j < Accumulator::num_elements; ++j) high.x[j] = total[j].value;
		wmma::store_matrix_sync(output + warp * 16, high, columns, wmma::mem_row_major);
		__syncthreads();
		for (int i = tid; i < p * p * p; i += threads)
			moments[size_t(cell) * (p * p * p) + i] = output[(i % p) * columns + i / p];
#endif
	}

	template<typename Real>
	__device__ Real TransferTensor(int p, int child, bool upward, const float* transfer, float* work, float* scratch) {
		const int node = threadIdx.x, k = p * p * p;
		const int x = node % p, y = (node / p) % p, z = node / (p * p);
		Real sum = 0.;
		if (node < k) for (int j = 0; j < p; ++j)
			sum += Real(transfer[(child & 1) * p * p + (upward ? j * p + x : x * p + j)]) * work[(z * p + y) * p + j];
		if (node < k) scratch[node] = sum;
		__syncthreads();
		sum = 0.;
		if (node < k) for (int j = 0; j < p; ++j)
			sum += Real(transfer[((child >> 1) & 1) * p * p + (upward ? j * p + y : y * p + j)]) * scratch[(z * p + j) * p + x];
		if (node < k) work[node] = sum;
		__syncthreads();
		sum = 0.;
		if (node < k) for (int j = 0; j < p; ++j)
			sum += Real(transfer[(child >> 2) * p * p + (upward ? j * p + z : z * p + j)]) * work[(j * p + y) * p + x];
		__syncthreads();
		return sum;
	}

	template<typename Real>
	__global__ void MultipoleToMultipole(int p, int n, const float* transfer, const float* children, float* parents) {
		const int node = threadIdx.x;
		const int k = p * p * p;
		__shared__ float work[512], scratch[512];
		const int simulation = blockIdx.x / (n * n * n);
		const int3 c = Coordinates(blockIdx.x % (n * n * n), n);
		double sum = 0.;
		for (int child = 0; child < 8; ++child) {
			const int index = simulation * 8 * n * n * n + Index(2 * c.x + (child & 1),
				2 * c.y + ((child >> 1) & 1), 2 * c.z + (child >> 2), 2 * n);
			if (node < k) work[node] = children[size_t(index) * k + node];
			__syncthreads();
			sum += TransferTensor<Real>(p, child, true, transfer, work, scratch);
		}
		if (node < k) parents[size_t(blockIdx.x) * k + node] = sum;
	}

	__global__ void RootToLocal(int k, const float* matrix, const float* moments, float* local) {
		const int node = threadIdx.x;
		if (node >= k) return;
		double sum = 0.;
		for (int j = 0; j < k; ++j) sum += double(matrix[j * k + node]) * moments[blockIdx.x * k + j];
		local[blockIdx.x * k + node] = sum;
	}

	template<int k>
	__global__ void MultipoleToLocal(int n, const float* matrix, const float* moments, float* local) {
		const int node = threadIdx.x;
		if (node >= k) return;
		const int cell = blockIdx.x % (n * n * n), simulation = blockIdx.x / (n * n * n);
		const int3 c = Coordinates(cell, n);
		double sum = 0.;
		// Children of the 27 neighboring parents, excluding the 27 near neighbors.
		for (int z = 2 * (c.z / 2 - 1); z <= 2 * (c.z / 2 + 1) + 1; ++z)
		for (int y = 2 * (c.y / 2 - 1); y <= 2 * (c.y / 2 + 1) + 1; ++y)
		for (int x = 2 * (c.x / 2 - 1); x <= 2 * (c.x / 2 + 1) + 1; ++x) {
			if (abs(c.x - x) <= 1 && abs(c.y - y) <= 1 && abs(c.z - z) <= 1) continue;
			const int offset = Index(c.x - x + 3, c.y - y + 3, c.z - z + 3, 7);
			const int source = simulation * n * n * n + Index(Wrap(x, n), Wrap(y, n), Wrap(z, n), n);
			// Independent dot-product lanes hide load/FMA latency and shorten each
			// float reduction. The 189 source-cell contributions accumulate in double.
			float a = 0.f, b = 0.f, d = 0.f, e = 0.f;
			int j = 0;
			for (; j + 3 < k; j += 4) {
				a += matrix[(size_t(offset) * k + j) * k + node] * moments[size_t(source) * k + j];
				b += matrix[(size_t(offset) * k + j + 1) * k + node] * moments[size_t(source) * k + j + 1];
				d += matrix[(size_t(offset) * k + j + 2) * k + node] * moments[size_t(source) * k + j + 2];
				e += matrix[(size_t(offset) * k + j + 3) * k + node] * moments[size_t(source) * k + j + 3];
			}
			for (; j < k; ++j) a += matrix[(size_t(offset) * k + j) * k + node] * moments[size_t(source) * k + j];
			sum += double(a) + b + d + e;
		}
		local[size_t(blockIdx.x) * k + node] = sum * n;
	}

	// Cells of the same child parity share each translation matrix.
	// Scaled high/low half pairs retain the low mantissa of float inputs. The
	// source-cell reduction is compensated, with explicit rounding intrinsics
	// so fast-math cannot reassociate away the compensation.
	template<int rows, bool recompressed>
	__global__ void MultipoleToLocalTensor(int inputCount, int n, float matrixScaleInv, const __half* matrixHigh, const __half* matrixLow,
		const __half* momentsHigh, const __half* momentsLow, const unsigned int* maximum, float* local,
		const int* factorRanks, const __half* leftHigh, const __half* leftLow, const __half* rightHigh, const __half* rightLow) {
#if __CUDA_ARCH__ >= 800
		using namespace nvcuda;
		const int k = recompressed ? 96 : inputCount;
		constexpr int columns = 32, threads = rows * 2, aStride = rows + 16, bStride = 48;
		constexpr int aSize = recompressed ? 40 * 96 : aStride * 32;
		__shared__ __align__(32) __half aHigh[aSize], aLow[aSize];
		__shared__ __align__(32) __half bHigh[bStride * columns], bLow[bStride * columns];
		__shared__ __align__(32) float output[rows * columns];
		__shared__ int sources[columns], targets[columns];
		const int tid = threadIdx.x, warp = tid / 32;
		const int parents = (n / 2) * (n / 2) * (n / 2), groups = (parents + columns - 1) / columns;
		const int parity = blockIdx.x % 8, group = (blockIdx.x / 8) % groups, simulation = blockIdx.x / (8 * groups);
		const int firstRow = blockIdx.y * rows;
		const int matrixStride = k == 96 ? 96 : ((k + 63) / 64) * 64;
		const int3 child{parity & 1, (parity >> 1) & 1, parity >> 2};
		int3 c{};
		if (tid < columns) {
			const int parent = group * columns + tid;
			const int3 p = Coordinates(parent % parents, n / 2);
			c = {2 * p.x + child.x, 2 * p.y + child.y, 2 * p.z + child.z};
			targets[tid] = parent < parents ? simulation * n * n * n + Index(c.x, c.y, c.z, n) : -1;
		}
		int exponent = 0;
		frexpf(__uint_as_float(maximum[simulation]), &exponent);
		const double scale = ldexp(double(matrixScaleInv) * n, exponent);
		using Accumulator = wmma::fragment<wmma::accumulator, 16, 16, 16, float>;
		Accumulator high[columns / 16], low[columns / 16];
		float total[columns / 16][Accumulator::num_elements]{}, correction[columns / 16][Accumulator::num_elements]{};
		for (int parentGroup = blockIdx.z; parentGroup < 27; parentGroup += gridDim.z)
		for (int sourceChild = 0; sourceChild < 8; ++sourceChild) {
			const int dx = 2 * (parentGroup % 3 - 1) + (sourceChild & 1) - child.x;
			const int dy = 2 * ((parentGroup / 3) % 3 - 1) + ((sourceChild >> 1) & 1) - child.y;
			const int dz = 2 * (parentGroup / 9 - 1) + (sourceChild >> 2) - child.z;
			if (abs(dx) <= 1 && abs(dy) <= 1 && abs(dz) <= 1) continue;
			const int offset = Index(3 - dx, 3 - dy, 3 - dz, 7);
			if (tid < columns) sources[tid] = simulation * n * n * n + Index(Wrap(c.x + dx, n), Wrap(c.y + dy, n), Wrap(c.z + dz, n), n);
			#pragma unroll
			for (int tile = 0; tile < columns / 16; ++tile) {
				wmma::fill_fragment(high[tile], 0.f); wmma::fill_fragment(low[tile], 0.f);
			}
			__syncthreads();
			const int factorRank = recompressed && factorRanks ? factorRanks[offset] : 0;
			if (factorRank) {
				const int factorRow = factorRank == 16 ? 0 : (warp & 1) * 16;
				const int factorColumn = factorRank == 16 ? warp * 16 : (warp >> 1) * 16;
				const int rowShift = factorRank == 16 ? 1 : 2, firstStride = factorRank + 8;
				// The final output scratch can hold the entire source tile until V
				// has been applied. Load all 96 inputs once, without inner barriers.
				auto* firstHigh = reinterpret_cast<__half*>(output);
				auto* firstLow = firstHigh + columns * 96;
				// Apply V to all source vectors first. Only the warps owning its
				// 16/32 rows compute; every thread still helps load the input tiles.
				for (int vector = tid; vector < (factorRank / 8) * k; vector += threads) {
					const int row = (vector & ((factorRank / 8) - 1)) * 8, column = vector >> rowShift;
					const size_t globalIndex = (size_t(offset) * k + column) * 32 + row;
					const int sharedIndex = column * firstStride + row;
					__pipeline_memcpy_async(aHigh + sharedIndex, rightHigh + globalIndex, 16);
					__pipeline_memcpy_async(aLow + sharedIndex, rightLow + globalIndex, 16);
				}
				for (int vector = tid; vector < columns * (k / 8); vector += threads) {
					const int column = vector / (k / 8), row = (vector % (k / 8)) * 8;
					const int sharedIndex = column * k + row;
					if (column < columns && targets[column] >= 0) {
						const size_t globalIndex = size_t(sources[column]) * k + row;
						__pipeline_memcpy_async(firstHigh + sharedIndex, momentsHigh + globalIndex, 16);
						__pipeline_memcpy_async(firstLow + sharedIndex, momentsLow + globalIndex, 16);
					}
					else if (column < columns) {
						*reinterpret_cast<int4*>(firstHigh + sharedIndex) = {};
						*reinterpret_cast<int4*>(firstLow + sharedIndex) = {};
					}
				}
				__pipeline_commit(); __pipeline_wait_prior(0); __syncthreads();
				if (factorColumn < columns) for (int q = 0; q < k; q += 16) {
					wmma::fragment<wmma::matrix_a, 16, 16, 16, __half, wmma::col_major> ah, al;
					wmma::fragment<wmma::matrix_b, 16, 16, 16, __half, wmma::col_major> bh, bl;
					wmma::load_matrix_sync(ah, aHigh + factorRow + q * firstStride, firstStride);
					wmma::load_matrix_sync(al, aLow + factorRow + q * firstStride, firstStride);
					wmma::load_matrix_sync(bh, firstHigh + factorColumn * k + q, k);
					wmma::load_matrix_sync(bl, firstLow + factorColumn * k + q, k);
					wmma::mma_sync(high[0], ah, bh, high[0]);
					wmma::mma_sync(low[0], ah, bl, low[0]);
					wmma::mma_sync(low[0], al, bh, low[0]);
				}
				__syncthreads();
				if (factorColumn < columns) {
					for (int j = 0; j < Accumulator::num_elements; ++j)
						high[0].x[j] = __fadd_rn(high[0].x[j], low[0].x[j] * (1.f / 4096.f));
					wmma::store_matrix_sync(output + factorRow * columns + factorColumn, high[0], columns, wmma::mem_row_major);
				}
				__syncthreads();
				for (int index = tid; index < factorRank * columns; index += threads) {
					// Factors and normalized moments are bounded by one, so the
					// intermediate is bounded by k=96. This power-of-two scale keeps
					// the second high/low product in the same normalized range.
					const float value = output[index] * (1.f / 128.f);
					const int sharedIndex = (index % columns) * bStride + index / columns;
					bHigh[sharedIndex] = __float2half_rn(value);
					bLow[sharedIndex] = __float2half_rn((value - __half2float(bHigh[sharedIndex])) * 4096.f);
				}
				for (int vector = tid; vector < (rows / 8) * factorRank; vector += threads) {
					const int row = firstRow + (vector % (rows / 8)) * 8, column = vector / (rows / 8);
					const size_t globalIndex = (size_t(offset) * 32 + column) * k + row;
					const int sharedIndex = column * aStride + (vector % (rows / 8)) * 8;
					__pipeline_memcpy_async(aHigh + sharedIndex, leftHigh + globalIndex, 16);
					__pipeline_memcpy_async(aLow + sharedIndex, leftLow + globalIndex, 16);
				}
				#pragma unroll
				for (int tile = 0; tile < columns / 16; ++tile) {
					wmma::fill_fragment(high[tile], 0.f); wmma::fill_fragment(low[tile], 0.f);
				}
				__pipeline_commit(); __pipeline_wait_prior(0); __syncthreads();
				for (int q = 0; q < factorRank; q += 16) {
					wmma::fragment<wmma::matrix_a, 16, 16, 16, __half, wmma::col_major> ah, al;
					wmma::fragment<wmma::matrix_b, 16, 16, 16, __half, wmma::col_major> bh, bl;
					wmma::load_matrix_sync(ah, aHigh + warp * 16 + q * aStride, aStride);
					wmma::load_matrix_sync(al, aLow + warp * 16 + q * aStride, aStride);
					#pragma unroll
					for (int tile = 0; tile < columns / 16; ++tile) {
						wmma::load_matrix_sync(bh, bHigh + tile * 16 * bStride + q, bStride);
						wmma::load_matrix_sync(bl, bLow + tile * 16 * bStride + q, bStride);
						wmma::mma_sync(high[tile], ah, bh, high[tile]);
						wmma::mma_sync(low[tile], ah, bl, low[tile]);
						wmma::mma_sync(low[tile], al, bh, low[tile]);
					}
				}
				__syncthreads();
			}
			else {
				for (int begin = 0; begin < k; begin += 32) {
					for (int vector = tid; vector < (rows / 8) * 32; vector += threads) {
						const int row = firstRow + (vector % (rows / 8)) * 8, column = begin + vector / (rows / 8);
						const size_t globalIndex = (size_t(offset) * k + column) * matrixStride + row;
						const int sharedIndex = (vector / (rows / 8)) * aStride + (vector % (rows / 8)) * 8;
						__pipeline_memcpy_async(aHigh + sharedIndex, matrixHigh + globalIndex, 16);
						__pipeline_memcpy_async(aLow + sharedIndex, matrixLow + globalIndex, 16);
					}
					const int column = tid / 4, row = begin + (tid % 4) * 8;
					const int sharedIndex = column * bStride + (tid % 4) * 8;
					if (column < columns && targets[column] >= 0) {
						const size_t globalIndex = size_t(sources[column]) * k + row;
						__pipeline_memcpy_async(bHigh + sharedIndex, momentsHigh + globalIndex, 16);
						__pipeline_memcpy_async(bLow + sharedIndex, momentsLow + globalIndex, 16);
					}
					else if (column < columns) {
						*reinterpret_cast<int4*>(bHigh + sharedIndex) = {};
						*reinterpret_cast<int4*>(bLow + sharedIndex) = {};
					}
					__pipeline_commit();
					__pipeline_wait_prior(0);
					__syncthreads();
					for (int q = 0; q < 32; q += 16) {
						wmma::fragment<wmma::matrix_a, 16, 16, 16, __half, wmma::col_major> ah, al;
						wmma::fragment<wmma::matrix_b, 16, 16, 16, __half, wmma::col_major> bh, bl;
						wmma::load_matrix_sync(ah, aHigh + warp * 16 + q * aStride, aStride);
						wmma::load_matrix_sync(al, aLow + warp * 16 + q * aStride, aStride);
						#pragma unroll
						for (int tile = 0; tile < columns / 16; ++tile) {
							wmma::load_matrix_sync(bh, bHigh + tile * 16 * bStride + q, bStride);
							wmma::load_matrix_sync(bl, bLow + tile * 16 * bStride + q, bStride);
							wmma::mma_sync(high[tile], ah, bh, high[tile]);
							wmma::mma_sync(low[tile], ah, bl, low[tile]);
							wmma::mma_sync(low[tile], al, bh, low[tile]);
						}
					}
					__syncthreads();
				}
			}
			#pragma unroll
			for (int tile = 0; tile < columns / 16; ++tile) for (int j = 0; j < Accumulator::num_elements; ++j) {
				const float contribution = __fadd_rn(high[tile].x[j], low[tile].x[j] * (1.f / 4096.f)) * (factorRank ? 128.f : 1.f);
				const float value = __fsub_rn(contribution, correction[tile][j]);
				const float next = __fadd_rn(total[tile][j], value);
				correction[tile][j] = __fsub_rn(__fsub_rn(next, total[tile][j]), value);
				total[tile][j] = next;
			}
		}
		#pragma unroll
		for (int tile = 0; tile < columns / 16; ++tile) {
			for (int j = 0; j < Accumulator::num_elements; ++j) high[tile].x[j] = total[tile][j] * scale;
			wmma::store_matrix_sync(output + warp * 16 * columns + tile * 16, high[tile], columns, wmma::mem_row_major);
		}
		__syncthreads();
		for (int index = tid; index < rows * columns; index += threads) {
			const int row = firstRow + index / columns, target = targets[index % columns];
			const size_t groupOffset = size_t(blockIdx.z) * (gridDim.x / (8 * groups)) * n * n * n * k;
			if (row < k && target >= 0) local[groupOffset + size_t(target) * k + row] = output[index];
		}
#endif
	}

	__global__ void ReduceTranslations(int count, int groups, float* local) {
		const int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= count) return;
		double sum = 0.;
		for (int group = 0; group < groups; ++group) sum += local[size_t(group) * count + i];
		local[i] = sum;
	}

	template<typename Real>
	__global__ void LocalToLocal(int p, int n, const float* transfer, const float* parents, float* children) {
		const int node = threadIdx.x;
		const int k = p * p * p;
		__shared__ float work[512], scratch[512];
		const int simulation = blockIdx.x / (n * n * n);
		const int3 c = Coordinates(blockIdx.x % (n * n * n), n);
		const int child = (c.x & 1) + 2 * (c.y & 1) + 4 * (c.z & 1);
		const int parent = simulation * (n / 2) * (n / 2) * (n / 2) + Index(c.x / 2, c.y / 2, c.z / 2, n / 2);
		if (node < k) work[node] = parents[size_t(parent) * k + node];
		__syncthreads();
		const double sum = TransferTensor<Real>(p, child, false, transfer, work, scratch);
		if (node < k) children[size_t(blockIdx.x) * k + node] += sum;
	}

	template<bool withPotE, typename Real, typename Vector>
	__device__ void AddPair(float3 delta, float charge, Real& potential, Vector& field) {
		const float r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
		if (r2 == 0.f || charge == 0.f) return;
		// Positions and charges are floats. Keep the long reduction in double,
		// but evaluate individual pairs using the GPU reciprocal-square-root unit.
		const float inverse = rsqrtf(r2), value = charge * inverse;
		const float scalar = value * inverse * inverse;
		if constexpr (withPotE) potential += value;
		field.x += scalar * delta.x;
		field.y += scalar * delta.y;
		field.z += scalar * delta.z;
	}

	template<bool withPotE>
	__global__ void ParticleToParticle(Float3 box, int n, const int2* ranges, const int* sorted,
		const float4* particles, double4* nearField) {
		__shared__ float4 sources[64];
		const int lane = threadIdx.x, cells = n * n * n;
		const int simulation = blockIdx.x / cells;
		const int3 c = Coordinates(blockIdx.x % cells, n);
		const int2 targets = ranges[blockIdx.x];
		for (int begin = targets.x; begin < targets.y; begin += 64) {
			const bool valid = begin + lane < targets.y;
			const int i = valid ? sorted[begin + lane] : 0;
			const float4 particle = valid ? particles[i] : float4{};
			const float3 position{particle.x - c.x * box.x / n, particle.y - c.y * box.y / n, particle.z - c.z * box.z / n};
			double potential = 0.;
			double3 field{};
			for (int z = c.z - 1; z <= c.z + 1; ++z)
			for (int y = c.y - 1; y <= c.y + 1; ++y)
			for (int x = c.x - 1; x <= c.x + 1; ++x) {
				const int wx = Wrap(x, n), wy = Wrap(y, n), wz = Wrap(z, n);
				const int2 range = ranges[simulation * cells + Index(wx, wy, wz, n)];
				// Cell-relative coordinates preserve boundary separations without
				// subtracting a full box (or doing double arithmetic) per pair.
				const float3 target{position.x - (x - c.x) * box.x / n,
					position.y - (y - c.y) * box.y / n, position.z - (z - c.z) * box.z / n};
				for (int source = range.x; source < range.y; source += 64) {
					const int count = min(64, range.y - source);
					if (lane < count) {
						float4 other = particles[sorted[source + lane]];
						other.x -= wx * box.x / n; other.y -= wy * box.y / n; other.z -= wz * box.z / n;
						sources[lane] = other;
					}
					__syncthreads();
					float tilePotential = 0.f;
					float3 tileField{};
					if (valid) for (int j = 0; j < count; ++j) {
						const float4 other = sources[j];
						const float3 delta{target.x - other.x, target.y - other.y, target.z - other.z};
						AddPair<withPotE>(delta, other.w, tilePotential, tileField);
					}
					if constexpr (withPotE) potential += tilePotential;
					field.x += tileField.x; field.y += tileField.y; field.z += tileField.z;
					__syncthreads();
				}
			}
			if (valid) nearField[i] = {field.x, field.y, field.z, potential};
		}
	}

	template<int p, typename Real, bool withPotE>
	__global__ void Evaluate(Interpolation interpolation, Float3 box, double3 leafScale, double3 boxInverse, int n, int count, const float4* particles,
		const int* slots, const unsigned int* keys, const int* sorted,
		const float* local, const double4* nearField, const int2* exclusionRanges, const int* exclusions, ForceEnergy* energies, ForceAccumulator forceAcc) {
		const int index = blockIdx.x * blockDim.x + threadIdx.x;
		if (index >= count) return;
		const int i = sorted[index];
		if (slots[i] < 0) return;
		const float4 particle = particles[i];
		const int key = keys[i], cells = n * n * n;
		const int3 c = Coordinates(key % cells, n);
		Real x[p], y[p], z[p], dx[p], dy[p], dz[p];
		interpolation.Basis<p, true>(double(particle.x) * leafScale.x - (2 * c.x + 1), x, dx);
		interpolation.Basis<p, true>(double(particle.y) * leafScale.y - (2 * c.y + 1), y, dy);
		interpolation.Basis<p, true>(double(particle.z) * leafScale.z - (2 * c.z + 1), z, dz);
		CompensatedSum<Real> farPotential, farX, farY, farZ;
		// Keep one coefficient plane live; fully unrolling all planes exhausts registers.
		#pragma unroll 1
		for (int d = 0; d < p; ++d) {
			Real planeValue = 0., planeX = 0., planeY = 0.;
			#pragma unroll
			for (int b = 0; b < p; ++b) {
				Real value = 0., derivative = 0.;
				#pragma unroll
				for (int a = 0; a < p; ++a) {
					const Real v = local[size_t(key) * (p * p * p) + (d * p + b) * p + a];
					value += v * x[a]; derivative += v * dx[a];
				}
				planeValue += value * y[b];
				planeX += derivative * y[b];
				planeY += value * dy[b];
			}
			if constexpr (withPotE) farPotential.Add(planeValue * z[d]);
			farX.Add(planeX * z[d]);
			farY.Add(planeY * z[d]);
			farZ.Add(planeValue * dz[d]);
		}
		const double4 near = nearField[i];
		double potential = withPotE ? near.w + farPotential.value : 0.;
		double3 field{near.x - double(farX.value) * leafScale.x,
			near.y - double(farY.value) * leafScale.y, near.z - double(farZ.value) * leafScale.z};
		// Remove only the nearest image of an excluded pair. Its periodic images remain.
		const int2 excluded = exclusionRanges[i];
		for (int index = excluded.x; index < excluded.y; ++index) {
			const float4 other = particles[exclusions[index]];
			double3 r{double(particle.x) - other.x, double(particle.y) - other.y, double(particle.z) - other.z};
			r.x -= round(r.x * boxInverse.x) * box.x;
			r.y -= round(r.y * boxInverse.y) * box.y;
			r.z -= round(r.z * boxInverse.z) * box.z;
			AddPair<withPotE>({float(r.x), float(r.y), float(r.z)}, -other.w, potential, field);
		}
		const double scale = double(particle.w) * PhysicsUtils::modifiedCoulombConstant;
		const ForceEnergy result{Float3{float(scale * field.x), float(scale * field.y), float(scale * field.z)}, float(.5 * scale * potential)};
		if (forceAcc.fx) forceAcc.Add<withPotE>(slots[i], result);
		else energies[i] = result;
	}
}

namespace {
	int TreeDepth(int particles, int target) {
		int depth = 1;
		for (int64_t cells = 8; cells * target < particles; cells *= 8) ++depth;
		return depth;
	}

	int TranslationGroups(int n) {
		// Geometry alone selects the partition, preserving standalone/batch results.
		return n <= 8 ? 27 : n == 16 ? 3 : 1;
	}
}

struct Operators {
	CudaBuffer<float> translations, root;
	CudaBuffer<float> transfers1D;
	CudaBuffer<float> compressed, project, expand;
	CudaBuffer<__half> translationsHigh, translationsLow;
	CudaBuffer<__half> projectHigh, projectLow, expandHigh, expandLow;
	CudaBuffer<__half> leftHigh, leftLow, rightHigh, rightLow;
	CudaBuffer<int> factorRanks;
	float matrixScaleInv = 1.f;
	int rank = 0;

	void Recompress(int order, cudaStream_t stream) {
		// Two padded factors only save arithmetic below half the shared-basis rank.
		// Retain the dense operator when complete-pivot cross approximation needs
		// more than 32 columns. Bound the full residual, rather than a pivot alone.
		if (rank != 96 || order > 6) return;
		constexpr int capacity = 32;
		const auto matrices = compressed.GetData();
		std::vector<float> left(size_t(343) * rank * capacity, 0.f), right(left.size(), 0.f);
		std::vector<int> ranks(343, 0);
		const double tolerance = std::pow(10., (-order - 6) * .5);
		for (int offset = 0; offset < 343; ++offset) {
			std::vector<double> residual(size_t(rank) * rank);
			double norm = 0.;
			for (size_t i = 0; i < residual.size(); ++i) {
				residual[i] = matrices[size_t(offset) * rank * rank + i] / matrixScaleInv;
				norm += residual[i] * residual[i];
			}
			if (norm == 0.) continue;
			const double threshold = tolerance * tolerance * norm;
			int count = 0;
			bool bounded = true;
			while (norm > threshold && count < capacity) {
				const auto pivot = std::max_element(residual.begin(), residual.end(),
					[](double a, double b) { return std::abs(a) < std::abs(b); });
				const int pivotRow = int(pivot - residual.begin()) % rank, pivotColumn = int(pivot - residual.begin()) / rank;
				const double inverse = 1. / std::sqrt(std::abs(*pivot)), sign = std::copysign(1., *pivot);
				std::vector<double> u(rank), v(rank);
				for (int i = 0; i < rank; ++i) {
					u[i] = residual[pivotColumn * rank + i] * inverse;
					v[i] = residual[i * rank + pivotRow] * inverse * sign;
					bounded = bounded && std::abs(u[i]) <= 1. && std::abs(v[i]) <= 1.;
					left[(size_t(offset) * capacity + count) * rank + i] = u[i];
					right[(size_t(offset) * rank + i) * capacity + count] = v[i];
				}
				norm = 0.;
				for (int column = 0; column < rank; ++column) for (int row = 0; row < rank; ++row) {
					double& value = residual[column * rank + row];
					value -= u[row] * v[column];
					norm += value * value;
				}
				++count;
			}
			if (norm <= threshold && bounded) ranks[offset] = ((count + 15) / 16) * 16;
		}
		factorRanks.SetData(ranks);
		const auto Pack = [&](const std::vector<float>& values, CudaBuffer<__half>& high, CudaBuffer<__half>& low) {
			CudaBuffer<float> buffer;
			buffer.SetData(values);
			high.Expand(values.size()); low.Expand(values.size());
			PackFactors<<<(values.size() + 127) / 128, 128, 0, stream>>>(int(values.size()), buffer.Get(), high.Get(), low.Get());
			LIMA_UTILS::genericErrorCheck(stream, "FMM individual operator compression");
		};
		Pack(left, leftHigh, leftLow); Pack(right, rightHigh, rightLow);
	}

	void Compress(int k, int order, Float3 box, cudaStream_t stream) {
		int exponent = 0;
		std::frexp(k / std::min({box.x, box.y, box.z}), &exponent);
		matrixScaleInv = std::ldexp(1.f, exponent);
		CudaBuffer<double> covariance;
		covariance.Expand(size_t(k) * k);
		OperatorCovariance<<<(k * k + 127) / 128, 128, 0, stream>>>(k, translations.Get(), covariance.Get());
		LIMA_UTILS::genericErrorCheck(stream, "FMM operator covariance");
		auto matrix = covariance.GetData();
		std::vector<double> vectors(size_t(k) * k, 0.);
		for (int i = 0; i < k; ++i) vectors[i * k + i] = 1.;
		// Cyclic Jacobi diagonalization of the small symmetric covariance matrix.
		// This initialization-only basis is shared by every level and simulation.
		for (int sweep = 0; sweep < 32; ++sweep) {
			double maximum = 0., diagonal = 0.;
			for (int i = 0; i < k; ++i) diagonal = std::max(diagonal, std::abs(matrix[i * k + i]));
			for (int p = 0; p < k; ++p) for (int q = p + 1; q < k; ++q) {
				const double off = matrix[p * k + q];
				maximum = std::max(maximum, std::abs(off));
				if (std::abs(off) <= diagonal * 1e-17) continue;
				const double tau = (matrix[q * k + q] - matrix[p * k + p]) / (2. * off);
				const double t = std::copysign(1., tau) / (std::abs(tau) + std::hypot(1., tau));
				const double c = 1. / std::hypot(1., t), s = t * c;
				matrix[p * k + p] -= t * off; matrix[q * k + q] += t * off;
				matrix[p * k + q] = matrix[q * k + p] = 0.;
				for (int j = 0; j < k; ++j) {
					if (j != p && j != q) {
						const double a = matrix[j * k + p], b = matrix[j * k + q];
						matrix[j * k + p] = matrix[p * k + j] = c * a - s * b;
						matrix[j * k + q] = matrix[q * k + j] = s * a + c * b;
					}
					const double a = vectors[j * k + p], b = vectors[j * k + q];
					vectors[j * k + p] = c * a - s * b;
					vectors[j * k + q] = s * a + c * b;
				}
			}
			if (maximum <= diagonal * 1e-16) break;
			if (sweep == 31) throw std::runtime_error("FMM compression eigensolver did not converge");
		}
		std::vector<int> indices(k);
		std::iota(indices.begin(), indices.end(), 0);
		std::stable_sort(indices.begin(), indices.end(), [&](int a, int b) { return matrix[a * k + a] > matrix[b * k + b]; });
		const double tolerance = std::pow(10., -order - 6) * matrix[indices[0] * k + indices[0]];
		while (rank < k && matrix[indices[rank] * k + indices[rank]] > tolerance) ++rank;
		// Keep tensor tiles aligned; zero-padding introduces no extra approximation.
		rank = ((rank + 31) / 32) * 32;
		std::vector<float> projection(size_t(k) * rank, 0.f), expansion(size_t(k) * rank, 0.f);
		for (int component = 0; component < std::min(k, rank); ++component)
			for (int node = 0; node < k; ++node) {
				const float value = vectors[node * k + indices[component]];
				projection[node * rank + component] = expansion[component * k + node] = value;
			}
		project.SetData(projection); expand.SetData(expansion);
		if (order <= 6) {
			const auto Pack = [&](int inputCount, int outputCount, const float* basis, CudaBuffer<__half>& high, CudaBuffer<__half>& low) {
				const int stride = outputCount == 96 ? 96 : ((outputCount + 63) / 64) * 64;
				const int count = ((inputCount + 15) / 16) * 16 * stride;
				high.Expand(count); low.Expand(count);
				PackBasis<<<(count + 127) / 128, 128, 0, stream>>>(inputCount, outputCount, basis, high.Get(), low.Get());
			};
			Pack(k, rank, project.Get(), projectHigh, projectLow);
			Pack(rank, k, expand.Get(), expandHigh, expandLow);
		}
		CudaBuffer<float> projected;
		projected.Expand(size_t(343) * k * rank);
		compressed.Expand(size_t(343) * rank * rank);
		ProjectOperators<<<(projected.Size() + 127) / 128, 128, 0, stream>>>(k, rank, project.Get(), translations.Get(), projected.Get());
		CompressOperators<<<(compressed.Size() + 127) / 128, 128, 0, stream>>>(k, rank, project.Get(), projected.Get(), compressed.Get());
		translationsHigh.Expand(size_t(343) * rank * (rank == 96 ? 96 : ((rank + 63) / 64) * 64));
		translationsLow.Expand(translationsHigh.Size());
		PackTranslations<<<(translationsHigh.Size() + 127) / 128, 128, 0, stream>>>(rank, 1.f / matrixScaleInv, compressed.Get(), translationsHigh.Get(), translationsLow.Get());
		LIMA_UTILS::genericErrorCheck(stream, "FMM operator compression");
		Recompress(order, stream);
	}
};

struct Tree {
	Interpolation interpolation{};
	Float3 box;
	cudaStream_t stream;
	int depth = 1, n = 2, batchCount = 0, particleSlots = 0;
	bool tensorTranslations = false;
	std::vector<size_t> offsets;
	std::shared_ptr<Operators> operators;
	CudaBuffer<float> moments, locals;
	CudaBuffer<float> compressedMoments, compressedLocals;
	CudaBuffer<__half> momentsHigh, momentsLow;
	CudaBuffer<unsigned int> maximumMoment;
	CudaBuffer<float4> particles;
	CudaBuffer<double4> nearField;
	CudaBuffer<unsigned int> keys, sortedKeys;
	CudaBuffer<int> values, sortedValues, slots, simulationSlots, exclusions;
	CudaBuffer<int2> ranges, exclusionRanges;
	CudaBuffer<char> sortStorage;
	size_t sortBytes = 0;

	void AllocateTranslationBuffers() {
		const size_t count = size_t(batchCount) * n * n * n * operators->rank;
		compressedMoments.Expand(count);
		momentsHigh.Expand(count); momentsLow.Expand(count);
		size_t partialCount = 0;
		for (int side = 2; side <= n; side *= 2)
			partialCount = std::max(partialCount, size_t(batchCount) * side * side * side * TranslationGroups(side) * operators->rank);
		compressedLocals.Expand(partialCount);
	}

	Tree(const std::vector<EngineSimulationData>& simulations, const SimParams& params, cudaStream_t inputStream, int selectedDepth,
		std::shared_ptr<Operators>& sharedOperators)
		: box(NodeIndex(simulations.front().simulation->box->boxparams.boxSize).toFloat3()), stream(inputStream) {
		if (params.fmm_order < 4 || params.fmm_order > maxOrder)
			throw std::invalid_argument("fmm_order must be between 4 and 8");
		if (params.fmm_leaf_particles <= 0) throw std::invalid_argument("fmm_leaf_particles must be positive");
		if (params.bc_select != PBC) throw std::invalid_argument("FMM currently requires PBC");
		interpolation.order = params.fmm_order;
		int device = 0;
		cudaGetDevice(&device);
		cudaDeviceProp properties{};
		cudaGetDeviceProperties(&properties, device);
		tensorTranslations = properties.major >= 8;
		const int p = interpolation.order, k = p * p * p;
		interpolation.count = k;
		for (int i = 0; i < p; ++i) interpolation.nodes[i] = std::cos(pi * (2 * i + 1) / (2 * p));
		for (int i = 0; i < p; ++i) {
			interpolation.invDenominator[i] = 1.;
			for (int j = 0; j < p; ++j) if (i != j)
				interpolation.invDenominator[i] /= interpolation.nodes[i] - interpolation.nodes[j];
		}
		int maxParticles = 0;
		std::vector<int> simSlots(simulations.size(), -1);
		for (size_t i = 0; i < simulations.size(); ++i) {
			const auto& simulation = simulations[i];
			particleSlots = std::max(particleSlots, (simulation.device.pclusters.offset + simulation.device.pclusters.count) * PersistentCluster::maxParticles);
			if (!simulation.device.active || TreeDepth(simulation.device.particles.count, params.fmm_leaf_particles) != selectedDepth) continue;
			simSlots[i] = batchCount++;
			maxParticles = std::max(maxParticles, simulation.device.particles.count);
		}
		if (!batchCount) return;
		maximumMoment.Expand(batchCount);
		while (int64_t(n) * n * n * params.fmm_leaf_particles.value < maxParticles) { n *= 2; ++depth; }
		if (int64_t(n) * n * n * batchCount > INT_MAX) throw std::invalid_argument("FMM tree exceeds index range");
		std::vector<int2> excludedRanges(particleSlots);
		std::vector<int> excluded;
		for (size_t simulationId = 0; simulationId < simulations.size(); ++simulationId) {
			const auto& simulation = simulations[simulationId];
			if (simSlots[simulationId] < 0) continue;
			const auto& b = *simulation.simulation->box;
			std::vector<int> ids(b.boxparams.totalParticles, -1);
			for (size_t pc = 0; pc < b.persistentClustersMetadata.size(); ++pc)
				for (int lane = 0; lane < b.persistentClustersMetadata[pc].nParticles; ++lane)
					ids.at(b.persistentClustersMetadata[pc].particleIdsGlobal[lane]) = (simulation.device.pclusters.offset + int(pc)) * PersistentCluster::maxParticles + lane;
			for (size_t id = 0; id < ids.size(); ++id) {
				const int begin = int(excluded.size());
				for (int j = 0; j < 32; ++j) {
					const int other = b.particlesBondedToParticle[id].Get(j);
					if (!ParticlesBondedToParticle::IsValue(other)) break;
					if (other != int(id)) excluded.push_back(ids.at(other));
				}
				excludedRanges.at(ids[id]) = {begin, int(excluded.size())};
			}
		}
		simulationSlots.SetData(simSlots);
		exclusionRanges.SetData(excludedRanges);
		exclusions.SetData(excluded);
		particles.Expand(particleSlots); keys.Expand(particleSlots); sortedKeys.Expand(particleSlots);
		nearField.Expand(particleSlots);
		values.Expand(particleSlots); sortedValues.Expand(particleSlots); slots.Expand(particleSlots);
		ranges.Expand(size_t(batchCount) * n * n * n);
		size_t total = 0;
		for (int level = 0, side = 1; level <= depth; ++level, side *= 2) {
			offsets.push_back(total);
			total += size_t(batchCount) * side * side * side * k;
		}
		moments.Expand(total); locals.Expand(total);

		cub::DeviceRadixSort::SortPairs(nullptr, sortBytes, keys.Get(), sortedKeys.Get(), values.Get(), sortedValues.Get(), particleSlots, 0, 32, stream);
		sortStorage.Expand(sortBytes);
		if (sharedOperators) {
			operators = sharedOperators;
			AllocateTranslationBuffers();
			return;
		}
		operators = std::make_shared<Operators>();
		auto& translations = operators->translations;
		auto& root = operators->root;
		translations.Expand(size_t(343) * k * k); root.Expand(size_t(k) * k);
		InitializeOperators<<<(translations.Size() + 127) / 128, 128, 0, stream>>>(interpolation, box, translations.Get());
		std::vector<float> transfers1D(2 * p * p);
		for (int child = 0; child < 2; ++child) for (int source = 0; source < p; ++source) for (int target = 0; target < p; ++target) {
			const double position = .5 * (interpolation.nodes[source] + (child ? 1. : -1.));
			double value = interpolation.invDenominator[target];
			for (int j = 0; j < p; ++j) if (j != target) value *= position - interpolation.nodes[j];
			transfers1D[(child * p + source) * p + target] = value;
		}
		operators->transfers1D.SetData(transfers1D);
		const double alpha = 2. / std::min({box.x, box.y, box.z});
		const double maxWave = 2. * alpha * std::sqrt(36.);
		const int nx = int(std::ceil(maxWave * box.x / (2. * pi)));
		const int ny = int(std::ceil(maxWave * box.y / (2. * pi)));
		const int nz = int(std::ceil(maxWave * box.z / (2. * pi)));
		std::vector<Wave> waves;
		for (int z = -nz; z <= nz; ++z) for (int y = -ny; y <= ny; ++y) for (int x = -nx; x <= nx; ++x) {
			if (x == 0 && y == 0 && z == 0) continue;
			const double kx = 2. * pi * x / box.x, ky = 2. * pi * y / box.y, kz = 2. * pi * z / box.z;
			const double squared = kx * kx + ky * ky + kz * kz;
			if (squared <= maxWave * maxWave)
				waves.push_back({kx, ky, kz, 4. * pi / (double(box.x) * box.y * box.z) * std::exp(-squared / (4. * alpha * alpha)) / squared});
		}
		CudaBuffer<Wave> waveBuffer;
		waveBuffer.SetData(waves);
		InitializePeriodicOperator<<<(k * k + 127) / 128, 128, 0, stream>>>(interpolation, box, waveBuffer.Get(), int(waves.size()), alpha, root.Get());
		LIMA_UTILS::genericErrorCheck(stream, "FMM initialization");
		operators->Compress(k, p, box, stream);
		AllocateTranslationBuffers();
		sharedOperators = operators;
	}
};

struct Controller::Data {
	std::vector<std::unique_ptr<Tree>> trees;
	Data(const std::vector<EngineSimulationData>& simulations, const SimParams& params, cudaStream_t stream) {
		if (params.fmm_leaf_particles <= 0) throw std::invalid_argument("fmm_leaf_particles must be positive");
		std::set<int> depths;
		for (const auto& simulation : simulations)
			if (simulation.device.active) depths.insert(TreeDepth(simulation.device.particles.count, params.fmm_leaf_particles));
		std::shared_ptr<Operators> operators;
		for (int depth : depths) trees.push_back(std::make_unique<Tree>(simulations, params, stream, depth, operators));
	}
};

Controller::Controller(const std::vector<EngineSimulationData>& simulations, const SimParams& params, cudaStream_t stream)
	: data(std::make_unique<Data>(simulations, params, stream)) {}
Controller::~Controller() = default;

template<int order, bool withPotE>
static void Calculate(Tree& d, SuperCluster* scData, SuperClusterMeta* scMeta, int nSuperclusters,
	ForceEnergy* forceEnergy, ForceAccumulator forceAcc) {
	if (!d.batchCount || !nSuperclusters) return;
	using Real = std::conditional_t<(order <= 6), float, double>;
	const double3 leafScale{2. * d.n / d.box.x, 2. * d.n / d.box.y, 2. * d.n / d.box.z};
	const double3 boxInverse{1. / d.box.x, 1. / d.box.y, 1. / d.box.z};
	const int k = d.interpolation.count, threads = ((k + 31) / 32) * 32;
	const int blocks = (d.particleSlots + 127) / 128;
	InitializeParticles<<<blocks, 128, 0, d.stream>>>(d.keys.Get(), d.values.Get(), d.slots.Get(), d.particleSlots);
	GatherParticles<<<(nSuperclusters * SuperCluster::maxParticles + 127) / 128, 128, 0, d.stream>>>(scData, scMeta,
		nSuperclusters, d.box, d.n, d.simulationSlots.Get(), d.keys.Get(), d.particles.Get(), d.slots.Get());
	cub::DeviceRadixSort::SortPairs(d.sortStorage.Get(), d.sortBytes, d.keys.Get(), d.sortedKeys.Get(),
		d.values.Get(), d.sortedValues.Get(), d.particleSlots, 0, 32, d.stream);
	cudaMemsetAsync(d.ranges.Get(), 0, d.ranges.Size() * sizeof(int2), d.stream);
	FindLeafRanges<<<blocks, 128, 0, d.stream>>>(d.sortedKeys.Get(), d.particleSlots, d.ranges.Get());
	const auto Deposit = [&]() {
		ParticleToMultipole<order, Real><<<d.batchCount * d.n * d.n * d.n, threads, 0, d.stream>>>(d.interpolation, leafScale, d.n,
			d.ranges.Get(), d.sortedValues.Get(), d.particles.Get(), d.moments.Get() + d.offsets.back());
	};
	if constexpr (order <= 6) {
		if (d.tensorTranslations) ParticleToMultipoleTensor<order><<<d.batchCount * d.n * d.n * d.n, ((order * order + 15) / 16) * 32, 0, d.stream>>>(
			d.interpolation, leafScale, d.n, d.ranges.Get(), d.sortedValues.Get(), d.particles.Get(), d.moments.Get() + d.offsets.back());
		else Deposit();
	}
	else Deposit();
	for (int level = d.depth - 1; level >= 0; --level) {
		const int n = 1 << level;
		MultipoleToMultipole<Real><<<d.batchCount * n * n * n, threads, 0, d.stream>>>(d.interpolation.order, n, d.operators->transfers1D.Get(),
			d.moments.Get() + d.offsets[level + 1], d.moments.Get() + d.offsets[level]);
	}
	RootToLocal<<<d.batchCount, threads, 0, d.stream>>>(k, d.operators->root.Get(), d.moments.Get(), d.locals.Get());
	for (int level = 1; level <= d.depth; ++level) {
		const int n = 1 << level;
		const auto Translate = [&]<int order>() {
			MultipoleToLocal<order * order * order><<<d.batchCount * n * n * n, threads, 0, d.stream>>>(n,
				d.operators->translations.Get(), d.moments.Get() + d.offsets[level], d.locals.Get() + d.offsets[level]);
		};
		if (d.tensorTranslations) {
			const int parents = (n / 2) * (n / 2) * (n / 2);
			const int rank = d.operators->rank;
			const int cells = d.batchCount * n * n * n;
			const auto Transform = [&](int inputCount, int outputCount, const float* basis, const __half* high, const __half* low,
				const float* input, float* output) {
				if constexpr (order <= 6) if (n >= 8) {
					const auto Tensor = [&]<int rows>() {
						TransformCoefficientsTensor<rows><<<dim3((cells + 31) / 32, (outputCount + rows - 1) / rows), rows * 2, 0, d.stream>>>(
							inputCount, outputCount, cells, high, low, input, output);
					};
					if (outputCount == 96) Tensor.operator()<96>();
					else Tensor.operator()<64>();
					return;
				}
				TransformCoefficients<<<cells, ((outputCount + 31) / 32) * 32, 0, d.stream>>>(inputCount, outputCount, basis, input, output);
			};
			Transform(k, rank, d.operators->project.Get(), d.operators->projectHigh.Get(), d.operators->projectLow.Get(),
				d.moments.Get() + d.offsets[level], d.compressedMoments.Get());
			const int countPerSimulation = n * n * n * rank;
			cudaMemsetAsync(d.maximumMoment.Get(), 0, d.batchCount * sizeof(unsigned int), d.stream);
			MaximumMoment<<<dim3((countPerSimulation + 1023) / 1024, d.batchCount), 256, 0, d.stream>>>(
				countPerSimulation, d.compressedMoments.Get(), d.maximumMoment.Get());
			PackMoments<<<dim3((countPerSimulation + 127) / 128, d.batchCount), 128, 0, d.stream>>>(
				countPerSimulation, d.compressedMoments.Get(), d.maximumMoment.Get(), d.momentsHigh.Get(), d.momentsLow.Get());
			const int groups = TranslationGroups(n);
			const auto TranslateTensor = [&]<int rows>() {
				MultipoleToLocalTensor<rows, rows == 96><<<dim3(d.batchCount * 8 * ((parents + 31) / 32), (rank + rows - 1) / rows, groups), rows * 2, 0, d.stream>>>(
					rank, n, d.operators->matrixScaleInv, d.operators->translationsHigh.Get(), d.operators->translationsLow.Get(),
					d.momentsHigh.Get(), d.momentsLow.Get(), d.maximumMoment.Get(), d.compressedLocals.Get(),
					d.operators->factorRanks.Get(), d.operators->leftHigh.Get(), d.operators->leftLow.Get(),
					d.operators->rightHigh.Get(), d.operators->rightLow.Get());
			};
			if (rank == 96) TranslateTensor.operator()<96>();
			else TranslateTensor.operator()<64>();
			if (groups > 1) {
				const int count = d.batchCount * n * n * n * rank;
				ReduceTranslations<<<(count + 127) / 128, 128, 0, d.stream>>>(count, groups, d.compressedLocals.Get());
			}
			Transform(rank, k, d.operators->expand.Get(), d.operators->expandHigh.Get(), d.operators->expandLow.Get(),
				d.compressedLocals.Get(), d.locals.Get() + d.offsets[level]);
		}
		else switch (d.interpolation.order) {
		case 4: Translate.operator()<4>(); break;
		case 5: Translate.operator()<5>(); break;
		case 6: Translate.operator()<6>(); break;
		case 7: Translate.operator()<7>(); break;
		case 8: Translate.operator()<8>(); break;
		}
		LocalToLocal<Real><<<d.batchCount * n * n * n, threads, 0, d.stream>>>(d.interpolation.order, n, d.operators->transfers1D.Get(),
			d.locals.Get() + d.offsets[level - 1], d.locals.Get() + d.offsets[level]);
	}
	ParticleToParticle<withPotE><<<d.batchCount * d.n * d.n * d.n, 64, 0, d.stream>>>(d.box, d.n, d.ranges.Get(),
		d.sortedValues.Get(), d.particles.Get(), d.nearField.Get());
	Evaluate<order, Real, withPotE><<<blocks, 128, 0, d.stream>>>(d.interpolation, d.box, leafScale, boxInverse, d.n, d.particleSlots, d.particles.Get(), d.slots.Get(),
		d.keys.Get(), d.sortedValues.Get(), d.locals.Get() + d.offsets.back(),
		d.nearField.Get(), d.exclusionRanges.Get(), d.exclusions.Get(), forceEnergy, forceAcc);
	LIMA_UTILS::genericErrorCheckNoSync("FMM evaluation");
}

void Controller::CalcCharges(SuperCluster* scData, SuperClusterMeta* scMeta, int nSuperclusters,
	ForceEnergy* forceEnergy, ForceAccumulator forceAcc) {
	for (auto& tree : data->trees) {
		const auto Run = [&]<int order>() {
			// EM writes persistent ForceEnergy values and always needs potential energy.
			if (!forceAcc.fx || forceAcc.potE) Calculate<order, true>(*tree, scData, scMeta, nSuperclusters, forceEnergy, forceAcc);
			else Calculate<order, false>(*tree, scData, scMeta, nSuperclusters, forceEnergy, forceAcc);
		};
		switch (tree->interpolation.order) {
		case 4: Run.operator()<4>(); break;
		case 5: Run.operator()<5>(); break;
		case 6: Run.operator()<6>(); break;
		case 7: Run.operator()<7>(); break;
		case 8: Run.operator()<8>(); break;
		}
	}
}

}
