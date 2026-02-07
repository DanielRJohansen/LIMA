#pragma once

#include <cuda_runtime.h>

namespace LAL {
	__device__ __host__ constexpr int32_t ceil(float num) {
		return (static_cast<float>(static_cast<int32_t>(num)) == num)
			? static_cast<int32_t>(num)
			: static_cast<int32_t>(num) + ((num > 0) ? 1 : 0);
	}

	template <typename T>
	__device__ __host__ static constexpr T max(const T l, const T r) {
		return r > l ? r : l;
	}

	template <typename T>
	__device__ __host__ static T min(const T l, const T r) {
		return r < l ? r : l;
	}

	template <typename T>
	constexpr static T abs(const T val) {
		return val < 0 ? -val : val;
	}

	__device__ inline bool IsEven(int a) {
		return !(a & 1);
	}

	__device__ inline float lerp(float v0, float v1, float t) {
		return fma(t, v1, fma(-t, v0, v0));
	}

	__device__ inline void SequentialPrefixSum(int* const data, int nElements) {
		if (threadIdx.x == 0) {
			for (int i = 1; i < nElements; i++) {
				data[i] += data[i - 1];
			}
		}
		__syncthreads();
		data[threadIdx.x] -= 1;
	}

	constexpr bool Fequal(float a, float b, float eps = 1e-6f) noexcept {
		return abs(a - b) <= eps;
	}

	// TODO These functions are NOT what their names elude they are, fix that
	// SLOW - Returns sum of actives before, thus must be -1 for 0-based index :)
	__device__ inline void doSequentialPrefixSum(uint8_t* onehot_remainers, int n_elements) {
		for (int i = 1; i < n_elements; i++) {
			if (threadIdx.x == i) {
				onehot_remainers[i] += onehot_remainers[i - 1];
				//KernelHelpersWarnings::verifyOnehotRemaindersIsValid(onehot_remainers, i);
			}
			__syncthreads();
		}
	}

	__device__ inline uint8_t computePrefixSum(const bool remain, uint8_t* utility_buffer, int n_elements) {
		utility_buffer[threadIdx.x] = static_cast<uint8_t>(remain);
		__syncthreads();

		doSequentialPrefixSum(utility_buffer, n_elements);
		//doBlellochPrefixSum

		const uint8_t solventindex_new = utility_buffer[threadIdx.x] - 1; // Underflow here doesn't matter, as the underflowing threads wont remain anyways :)
		return solventindex_new;
	}

	template <typename T>
	__device__ void ExclusiveScan(T* data, int nElements) {
		const int tid = threadIdx.x;

		// Up-sweep phase (reduce)
		for (int stride = 1; stride < nElements; stride *= 2) {
			int index = (tid + 1) * stride * 2 - 1;
			if (index < nElements) {
				data[index] += data[index - stride];
			}
			__syncthreads();
		}

		// Clear the last element for exclusive scan
		if (tid == 0) {
			data[nElements - 1] = 0;
		}
		__syncthreads();

		// Down-sweep phase
		for (int stride = nElements / 2; stride > 0; stride /= 2) {
			int index = (tid + 1) * stride * 2 - 1;
			if (index < nElements) {
				T temp = data[index - stride];
				data[index - stride] = data[index];
				data[index] += temp;
			}
			__syncthreads();
		}
	}


	template<typename T>
	__device__ inline void distributedSummation(T* arrayptr, int array_len) {				// Places the result at pos 0 of input_array
		T temp;			// This is a lazy soluation, but maybe it is also fast? Definitely simple..
		for (int i = 1; i < array_len; i *= 2) {	// Distributed averaging							// Make a generic and SAFER function for this, PLEASE OK??
			if ((threadIdx.x + i) < array_len) {
				temp = arrayptr[threadIdx.x] + arrayptr[threadIdx.x + i];
			}
			__syncthreads();
			arrayptr[threadIdx.x] = temp;
			__syncthreads();
		}
	}

	template <typename T, typename Accessor = decltype([](const T& a) { return a; }) >
	__device__ inline void Sort(T* data, int nElements, Accessor accessor = [](const T& a) { return a; }) {
		int tid = threadIdx.x;

		for (int k = 2; k <= nElements; k <<= 1) {
			for (int j = k >> 1; j > 0; j >>= 1) {
				int ixj = tid ^ j;
				if (ixj > tid) {
					if ((tid & k) == 0) {
						if (accessor(data[tid]) > accessor(data[ixj])) {
							// Swap data[tid] and data[ixj]
							T temp = data[tid];
							data[tid] = data[ixj];
							data[ixj] = temp;
						}
					}
					else {
						if (accessor(data[tid]) < accessor(data[ixj])) {
							// Swap data[tid] and data[ixj]
							T temp = data[tid];
							data[tid] = data[ixj];
							data[ixj] = temp;
						}
					}
				}
				__syncthreads(); // Synchronize to ensure all threads complete this step before moving on
			}
		}
	}

	// Assumes nValues is a power of two and <= blockDim.x
	template <int nValues>
	__device__ inline void Sort(int* keys, Float3* attachedData)
	{
		static_assert(nValues == 16 || nValues == 32 || nValues == 64 || nValues == 128);
		for (int k = 2; k <= nValues; k <<= 1) {
			for (int j = k >> 1; j > 0; j >>= 1) {
				int i = threadIdx.x;
				if (i < nValues) {
					int ixj = i ^ j;
					if (ixj > i) {
						bool ascending = ((i & k) == 0);
						int key_i = keys[i];
						int key_j = keys[ixj];

						if ((ascending && key_i > key_j) ||
							(!ascending && key_i < key_j)) {

							keys[i] = key_j;
							keys[ixj] = key_i;

							Float3 tmp = attachedData[i];
							attachedData[i] = attachedData[ixj];
							attachedData[ixj] = tmp;
						}
					}
				}
				__syncthreads();
			}
		}
	}


//	// Must always be called by blocks with blockdim=32,1,1
//	template <int nBins, int nValuesPerBin, typename T>
//	__device__ __forceinline__
//		void SortBins(T* keys, int* ids)
//	{
//		//static_assert(nValuesPerBin <= 32);
//		static_assert((nBins*nValuesPerBin) % 32 == 0, "Total must be a multiple of 32");
//		static_assert((nValuesPerBin & (nValuesPerBin - 1)) == 0, "nValuesPerBin must be power-of-two for bitonic.");
//
//		bool active = (threadIdx.x / nValuesPerBin) == (threadIdx.x + 1) / nValuesPerBin;
//		int nIterations = nValuesPerBin / 2 + 1;
//		int nBatches = (nBins * nValuesPerBin) / 32;
//
//#pragma unroll
//		for (int iter = 0; iter < nIterations; iter++) {
//#pragma unroll
//			for (int batchIndex = 0; batchIndex < nBatches; batchIndex++){			
//				int i = batchIndex * 32 + threadIdx.x;
//				int j = batchIndex * 32 + threadIdx.x + 1;
//
//				if (threadIdx.x % 2 == 0 && active) {					
//					if (keys[i] > keys[j]) {
//						// Swap keys
//						T tempKey = keys[i];
//						keys[i] = keys[j];
//						keys[j] = tempKey;
//						// Swap ids
//						T tempId = ids[i];
//						ids[i] = ids[j];
//						ids[j] = tempId;
//					}
//				}
//				__syncthreads();
//
//				if (threadIdx.x % 2 == 1 && active) {
//					if (keys[i] > keys[j]) {
//						// Swap keys
//						T tempKey = keys[i];
//						keys[i] = keys[j];
//						keys[j] = tempKey;
//						// Swap ids
//						T tempId = ids[i];
//						ids[i] = ids[j];
//						ids[j] = tempId;
//					}
//				}
//				__syncthreads();
//			}
//		}
//	}

	template <int nBins, int nValuesPerBin, typename T>
	__device__ __forceinline__
		void SortBins(T* keys, int* ids)
	{
		constexpr int totalValues = nBins * nValuesPerBin;
		static_assert(totalValues % 32 == 0, "Total must be a multiple of 32");
		static_assert((nValuesPerBin & (nValuesPerBin - 1)) == 0, "nValuesPerBin must be power-of-two");

		constexpr int nSegments = totalValues / 32;
		constexpr int binMask = nValuesPerBin - 1;

		int nIterations = nValuesPerBin / 2 + 1;

		for (int i = 0; i < nIterations; i++) {

			// even phase
			for (int seg = 0; seg < nSegments; ++seg) {
				int idx = threadIdx.x + seg * 32;
				int inBin = idx & binMask;

				if ((inBin & 1) == 0 && (inBin + 1) < nValuesPerBin) {
					int j = idx + 1;
					if (keys[idx] > keys[j]) {
						T  tempKey = keys[idx]; keys[idx] = keys[j]; keys[j] = tempKey;
						int tempId = ids[idx];  ids[idx] = ids[j];  ids[j] = tempId;
					}
				}
			}
			__syncthreads();

			// odd phase
			for (int seg = 0; seg < nSegments; ++seg) {
				int idx = threadIdx.x + seg * 32;
				int inBin = idx & binMask;

				if ((inBin & 1) == 1 && (inBin + 1) < nValuesPerBin) {
					int j = idx + 1;
					if (keys[idx] > keys[j]) {
						T  tempKey = keys[idx]; keys[idx] = keys[j]; keys[j] = tempKey;
						int tempId = ids[idx];  ids[idx] = ids[j];  ids[j] = tempId;
					}
				}
			}
			__syncthreads();
		}
	}

//	template <int nBins, int nValuesPerBin, typename T>
//	__device__ __forceinline__ void SortBins(T* keys, int* ids)
//	{
//		constexpr int kThreads = 32;
//		constexpr int totalValues = nBins * nValuesPerBin;
//
//		static_assert(totalValues % kThreads == 0, "Total must be a multiple of 32");
//		static_assert(nValuesPerBin >= 2, "nValuesPerBin must be >= 2");
//		static_assert((nValuesPerBin & (nValuesPerBin - 1)) == 0, "nValuesPerBin must be power-of-two");
//
//		constexpr int kItersPerThread = totalValues / kThreads;
//		constexpr int kMaskInBin = nValuesPerBin - 1;
//
//		const int lane = threadIdx.x;
//
//		auto SwapAt = [&](int a, int b) {
//			T  tk = keys[a]; keys[a] = keys[b]; keys[b] = tk;
//			int ti = ids[a];  ids[a] = ids[b];  ids[b] = ti;
//		};
//
//		// Odd-even transposition sort per bin.
//#pragma unroll
//		for (int pass = 0; pass < nValuesPerBin; ++pass)
//		{
//			// Even phase: (0,1)(2,3)...
//#pragma unroll
//			for (int it = 0; it < kItersPerThread; ++it)
//			{
//				const int idx = lane + it * kThreads;
//				const int inBin = idx & kMaskInBin;
//				if (((inBin & 1) == 0) && (inBin + 1 < nValuesPerBin))
//				{
//					const int j = idx + 1;
//					if (keys[idx] > keys[j]) SwapAt(idx, j);
//				}
//			}
//			__syncthreads();
//
//			// Odd phase: (1,2)(3,4)...
//#pragma unroll
//			for (int it = 0; it < kItersPerThread; ++it)
//			{
//				const int idx = lane + it * kThreads;
//				const int inBin = idx & kMaskInBin;
//				if (((inBin & 1) == 1) && (inBin + 1 < nValuesPerBin))
//				{
//					const int j = idx + 1;
//					if (keys[idx] > keys[j]) SwapAt(idx, j);
//				}
//			}
//			__syncthreads();
//		}
//	}


}