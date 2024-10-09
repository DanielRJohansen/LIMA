#pragma once

#include "vector"
#include "LimaTypes.cuh"
#include <span>

namespace Statistics {

    // Returns <slope, intersect>
    std::pair<float, float> linearFit(const std::vector<float>& x, const std::vector<float>& y);

    Float3 Mean(const std::span<const Float3>& data);

    float calculateR2(const std::vector<float>& x, const std::vector<float>& y, float slope, float intercept);

    Float3 CalculateMinimaxPoint(const std::span<const Float3>& data);

    float Mean(const std::vector<float>& vec);

	float StdDev(const std::vector<float>& vec);

    float Max(const float* const data, size_t size);

    float MaxLen(const Float3* const data, size_t n);
}