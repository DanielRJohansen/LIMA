#pragma once

#include "LimaTypes.cuh"

#include <assert.h>
#include <cmath>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace LIMA_UTILS {

	static int roundUp(int numToRound, int multiple)
	{
		assert(multiple);
		return ((numToRound + multiple - 1) / multiple) * multiple;
	}

	// Throw on CUDA errors. The variants with a text synchronize first
	void genericErrorCheck(const char* text);
	void genericErrorCheck(cudaStream_t stream, const char* text);
	void genericErrorCheckNoSync(const char* text);
	void genericErrorCheck(cudaError_t cuda_status);
}

namespace StringUtils {
    // Formats into "%%.%% [s/min/hr/days/weeks/months/years]". Takes seconds rather than a std::chrono::duration, as <chrono> is slow to compile
    std::string FormatTime(double seconds, int decimalPlacesBeforePoint, int decimalPlacesAfterPoint);

    std::vector<std::string> SplitWords(std::string_view line);
}


// Lima Algorithm Library
namespace LAL {

    // TODO: Make a lowest-level file for these type agnostic algo's so we can use them in limatypes.cuh
    //template <typename T>
    //__device__ __host__ static T max(const T l, const T r) {
    //    return r > l ? r : l;
    //}
    
    bool constexpr isPowerOf2(int n) {
        return (n > 0) && ((n & (n - 1)) == 0);
    }


    constexpr int powi(int base, int exp) {
        int res = 1;
        for (int i = 0; i < exp; i++) {
			res *= base;
		}
        return res;
    }

    constexpr void RotatePoint(Float3& point, const Float3& rotationCenter, const Float3& rotation) {
        point = point - rotationCenter;
        point = Float3::rodriguesRotatation(point, Float3(1, 0, 0), rotation.x);
        point = Float3::rodriguesRotatation(point, Float3(0, 1, 0), rotation.y);
        point = Float3::rodriguesRotatation(point, Float3(0, 0, 1), rotation.z);
        point = point + rotationCenter;
    }


    //float LargestDiff(const Float3 queryPoint, const std::span<Float3>& points);

    template <typename ContainerType>
    float LargestDiff(const Float3 queryPoint, const ContainerType& points) {
        float maxDiff = 0.f;
        for (const Float3& p : points) {
            maxDiff = std::max(maxDiff, (p - queryPoint).len()); // Compute lenSq if optimizing
        }
        return maxDiff;
    }
}
