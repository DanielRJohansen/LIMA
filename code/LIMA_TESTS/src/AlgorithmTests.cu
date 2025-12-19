#include "AlgorithmTests.h"

//#include "DeviceAlgorithms.cuh"
//#include <random>



template <int nBins, int nValuesPerBin>
__global__ void TestSortKernel32(int* keys, int* ids)
{
	//static_assert(nBins * nValuesPerBin == 32);

	//// exactly one warp
	//Sort<nBins, nValuesPerBin>(keys, ids);
}



namespace KernelAlgorithms {

	LimaUnittestResult WarpSort32_Unittest(EnvMode envmode)
		constexpr int totalValues = 64;
		constexpr int iters = 20'000;

		auto runCase = [&](int nBins, int nValuesPerBin)
			{
				ASSERT(nBins * nValuesPerBin == 32, "invalid config");

				std::vector<int> hKeys(totalValues);
				std::vector<int> hIds(totalValues);
				std::vector<int> refKeys(totalValues);
				std::vector<int> refIds(totalValues);

				int* dKeys;
				int* dIds;
				cudaMalloc(&dKeys, totalValues * sizeof(int));
				cudaMalloc(&dIds, totalValues * sizeof(int));

				std::mt19937 rng(12345);
				std::uniform_int_distribution<int> dist(0, 1'000'000);

				for (int it = 0; it < iters; ++it) {
					for (int i = 0; i < totalValues; ++i) {
						hKeys[i] = dist(rng);
						hIds[i] = i;
					}

					refKeys = hKeys;
					refIds = hIds;

					// CPU reference, per half, per bin
					for (int half = 0; half < 2; ++half) {
						const int base = half * 32;
						for (int b = 0; b < nBins; ++b) {
							const int off = base + b * nValuesPerBin;
							std::vector<std::pair<int, int>> tmp;
							for (int i = 0; i < nValuesPerBin; ++i)
								tmp.emplace_back(refKeys[off + i], refIds[off + i]);

							std::sort(tmp.begin(), tmp.end(),
								[](auto a, auto b) { return a.first < b.first; });

							for (int i = 0; i < nValuesPerBin; ++i) {
								refKeys[off + i] = tmp[i].first;
								refIds[off + i] = tmp[i].second;
							}
						}
					}

					cudaMemcpy(dKeys, hKeys.data(), totalValues * sizeof(int), cudaMemcpyHostToDevice);
					cudaMemcpy(dIds, hIds.data(), totalValues * sizeof(int), cudaMemcpyHostToDevice);

					// first 32
					if (nBins == 1)
						TestSortKernel32<1, 32> << <1, 32 >> > (dKeys, dIds);
					else if (nBins == 4)
						TestSortKernel32<4, 8> << <1, 32 >> > (dKeys, dIds);
					else if (nBins == 16)
						TestSortKernel32<16, 2> << <1, 32 >> > (dKeys, dIds);

					// second 32
					if (nBins == 1)
						TestSortKernel32<1, 32> << <1, 32 >> > (dKeys + 32, dIds + 32);
					else if (nBins == 4)
						TestSortKernel32<4, 8> << <1, 32 >> > (dKeys + 32, dIds + 32);
					else if (nBins == 16)
						TestSortKernel32<16, 2> << <1, 32 >> > (dKeys + 32, dIds + 32);

					cudaDeviceSynchronize();

					cudaMemcpy(hKeys.data(), dKeys, totalValues * sizeof(int), cudaMemcpyDeviceToHost);
					cudaMemcpy(hIds.data(), dIds, totalValues * sizeof(int), cudaMemcpyDeviceToHost);

					for (int i = 0; i < totalValues; ++i) {
						ASSERT(hKeys[i] == refKeys[i], "key mismatch");
						ASSERT(hIds[i] == refIds[i], "id mismatch");
					}
				}

				cudaFree(dKeys);
				cudaFree(dIds);
			};

		// 64 values total, tested via 2x32
		runCase(1, 32);   // effectively 1×64
		runCase(4, 8);    // effectively 4×16
		runCase(16, 2);   // effectively 16×4

		return LimaUnittestResult{ true, "Warp sort correct for all bin configs", true };
	}
}



