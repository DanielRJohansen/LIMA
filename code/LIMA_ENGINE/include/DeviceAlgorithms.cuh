#pragma once

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

	__device__ __host__ static int32_t abs(const int32_t val) {
		return val < 0 ? -val : val;
	}

	__device__ inline bool IsEven(int a) {
		return !(a & 1);
	}

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
}

