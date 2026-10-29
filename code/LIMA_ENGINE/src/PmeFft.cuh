#pragma once

// PME's z pass: the forward FFT along z, the Greens function multiply (and the reciprocal virial), and the inverse FFT along z, in
// one kernel, ZGreensKernel. cuFFT transforms the planes along x and y around it. It replaces the z passes of cuFFT's 3D transforms
// and ApplyGreensFunctionKernel, so the grid is read and written 6 times per step rather than 8: on small grids the kernel launches
// cost, on large ones the passes (villin -2.1%, adh_dodec -1.8%, stmv neutral). Grids whose z axis is longer than maxLength, or has a
// prime factor above 5, keep cuFFT's 3D transforms.
//
// The 1D transforms are Stockham autosort FFTs with radices 4, 2, 3 and 5, in shared memory, on a pencil of consecutive x frequencies:
// element z of column c at data[z * width + c], so the global loads are coalesced. Twiddles come from a table computed in double on
// the host. The Greens function is separable but for its 1 / k^2, so it is computed from per axis tables rather than loaded per mode

#include "PmeFft.h"
#include "LimaTypes.cuh"
#include "Virial.cuh"

// The kernels, included by PME.cu only
namespace PmeFft {
	__device__ __forceinline__ float2 Add(float2 a, float2 b) { return float2{ a.x + b.x, a.y + b.y }; }
	__device__ __forceinline__ float2 Sub(float2 a, float2 b) { return float2{ a.x - b.x, a.y - b.y }; }
	__device__ __forceinline__ float2 Scale(float2 a, float s) { return float2{ a.x * s, a.y * s }; }
	__device__ __forceinline__ float2 Mul(float2 a, float2 b) { return float2{ a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x }; }
	__device__ __forceinline__ float2 Conj(float2 a) { return float2{ a.x, -a.y }; }
	// -i a in the forward transform, i a in the inverse
	template <bool inverse>
	__device__ __forceinline__ float2 RotateQuarter(float2 a) { return inverse ? float2{ -a.y, a.x } : float2{ a.y, -a.x }; }

	template <int R, bool inverse>
	__device__ __forceinline__ void Butterfly(float2 (&v)[R]) {
		if constexpr (R == 2) {
			const float2 a = v[0], b = v[1];
			v[0] = Add(a, b);
			v[1] = Sub(a, b);
		}
		else if constexpr (R == 4) {
			const float2 a = Add(v[0], v[2]), b = Sub(v[0], v[2]), c = Add(v[1], v[3]), d = RotateQuarter<inverse>(Sub(v[1], v[3]));
			v[0] = Add(a, c);
			v[1] = Add(b, d);
			v[2] = Sub(a, c);
			v[3] = Sub(b, d);
		}
		else if constexpr (R == 3) {
			const float2 t = Add(v[1], v[2]);
			const float2 m = Sub(v[0], Scale(t, 0.5f));
			const float2 s = RotateQuarter<inverse>(Scale(Sub(v[1], v[2]), 0.86602540378443865f));
			v[0] = Add(v[0], t);
			v[1] = Add(m, s);
			v[2] = Sub(m, s);
		}
		else if constexpr (R == 5) {
			constexpr float c1 = 0.30901699437494742f, c2 = -0.80901699437494742f, s1 = 0.95105651629515357f, s2 = 0.58778525229247313f;
			const float2 t1 = Add(v[1], v[4]), t2 = Add(v[2], v[3]), t3 = Sub(v[1], v[4]), t4 = Sub(v[2], v[3]);
			const float2 a1 = Add(v[0], Add(Scale(t1, c1), Scale(t2, c2)));
			const float2 a2 = Add(v[0], Add(Scale(t1, c2), Scale(t2, c1)));
			const float2 b1 = RotateQuarter<inverse>(Add(Scale(t3, s1), Scale(t4, s2)));
			const float2 b2 = RotateQuarter<inverse>(Sub(Scale(t3, s2), Scale(t4, s1)));
			v[0] = Add(v[0], Add(t1, t2));
			v[1] = Add(a1, b1);
			v[2] = Add(a2, b2);
			v[3] = Sub(a2, b2);
			v[4] = Sub(a1, b1);
		}
	}

	// One Stockham stage, after ns points of each transform are combined: butterflies j of the length / R, with k = j % ns
	template <int R, bool inverse>
	__device__ void Stage(const float2* in, float2* out, int count, int length, int ns, const float2* twiddles) {
		const int nButterflies = length / R;
		const int twiddleStep = length / (ns * R);
		for (int b = threadIdx.x; b < count * nButterflies; b += blockDim.x) {
			const int t = b % count, j = b / count, k = j % ns;
			float2 v[R];
#pragma unroll
			for (int r = 0; r < R; r++)
				v[r] = in[(j + r * nButterflies) * count + t];
#pragma unroll
			for (int r = 1; r < R; r++) {
				const float2 w = twiddles[k * r * twiddleStep];
				v[r] = Mul(v[r], inverse ? Conj(w) : w);
			}
			Butterfly<R, inverse>(v);
			const int outBase = (j - k) * R + k;
#pragma unroll
			for (int r = 0; r < R; r++)
				out[(outBase + r * ns) * count + t] = v[r];
		}
	}

	// count transforms of axis.length, unnormalized, from data, using work as the second buffer. Returns the buffer holding the result.
	// Synchronizes before reading data and after writing the result
	template <bool inverse>
	__device__ float2* Transform(float2* data, float2* work, int count, const Axis& axis, const float2* twiddles) {
		int ns = 1;
		for (int s = 0; s < axis.nStages; s++) {
			__syncthreads();
			switch (axis.radices[s]) {
			case 2: Stage<2, inverse>(data, work, count, axis.length, ns, twiddles); break;
			case 3: Stage<3, inverse>(data, work, count, axis.length, ns, twiddles); break;
			case 4: Stage<4, inverse>(data, work, count, axis.length, ns, twiddles); break;
			default: Stage<5, inverse>(data, work, count, axis.length, ns, twiddles); break;
			}
			ns *= axis.radices[s];
			float2* const swap = data;
			data = work;
			work = swap;
		}
		__syncthreads();
		return data;
	}

	// A pencil of width x frequencies of one y row per block, gridDim = (ceil(nxh / width), ny, batches), blockDim = threads. Shared:
	// nz twiddles and 2 buffers of nz * width, see SharedBytes. Between the forward and inverse transform each mode is multiplied by
	// its Greens function value: greensScale / k^2 times the product of its axis factors, where axisTables holds (k^2, factor) per x, y
	// and z frequency index, see Controller::UploadGreensAxisTables. As ApplyGreensFunctionKernel, which includes the 1 / N
	// normalization. With computeVirial the block's part of the reciprocal virial goes to virialPartials
	template <bool computeVirial>
	__global__ void __launch_bounds__(threads) ZGreensKernel(float2* const fourier, const Axis axis, const float2* const twiddlesZ, int nx, int ny,
		const float2* const axisTables, float greensScale, float virialFactor, float inverseTwoKappaSq, long long* const virialPartials) {
		extern __shared__ float2 shared[];
		const int nz = axis.length, nxh = nx / 2 + 1, width = PencilWidth(nz);
		float2* const twiddles = shared;
		float2* const data = twiddles + nz;
		float2* const work = data + nz * width;
		for (int i = threadIdx.x; i < nz; i += blockDim.x)
			twiddles[i] = twiddlesZ[i];
		const int firstX = blockIdx.x * width;
		const int y = blockIdx.y;
		float2* const grid = fourier + size_t(blockIdx.z) * nz * ny * nxh;
		for (int i = threadIdx.x; i < nz * width; i += blockDim.x) {
			const int z = i / width, c = i % width;
			data[i] = firstX + c < nxh ? grid[(size_t(z) * ny + y) * nxh + firstX + c] : float2{};
		}
		float2* const modes = Transform<false>(data, work, width, axis, twiddles);

		float virial = 0.f;
		for (int i = threadIdx.x; i < nz * width; i += blockDim.x) {
			const int z = i / width, kx = firstX + i % width;
			if (kx >= nxh)
				continue;
			const float2 tx = axisTables[kx], ty = axisTables[nxh + y], tz = axisTables[nxh + ny + z];
			const float kSq = tx.x + ty.x + tz.x;
			const float greens = kSq > 0.f ? greensScale / kSq * (tx.y * ty.y * tz.y) : 0.f;
			const float2 rho = modes[i];
			if constexpr (computeVirial) {
				// The modes not stored, kx beyond nx / 2, mirror those with 0 < kx < nx / 2
				const float weight = kx == 0 || (nx % 2 == 0 && kx == nx / 2) ? 1.f : 2.f;
				virial += 0.5f * weight * (rho.x * rho.x + rho.y * rho.y) * greens * (1.f - kSq * inverseTwoKappaSq) * virialFactor;
			}
			modes[i] = Scale(rho, greens);
		}
		if constexpr (computeVirial) {
			__shared__ float warpVirials[threads / 32];
			for (int offset = 16; offset > 0; offset >>= 1)
				virial += __shfl_down_sync(LIMA_FULL_WARP_MASK, virial, offset);
			if ((threadIdx.x & 31) == 0)
				warpVirials[threadIdx.x >> 5] = virial;
			__syncthreads();
			if (threadIdx.x == 0) {
				float sum = 0.f;
				for (int w = 0; w < threads / 32; w++)
					sum += warpVirials[w];
				virialPartials[(size_t(blockIdx.z) * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x] = Virial::ToFixed(sum);
			}
		}

		const float2* const result = Transform<true>(modes, modes == data ? work : data, width, axis, twiddles);
		for (int i = threadIdx.x; i < nz * width; i += blockDim.x) {
			const int z = i / width, c = i % width;
			if (firstX + c < nxh)
				grid[(size_t(z) * ny + y) * nxh + firstX + c] = result[i];
		}
	}
}
