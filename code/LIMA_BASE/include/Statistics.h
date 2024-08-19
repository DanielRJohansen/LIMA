#pragma once

#include "vector"
#include "LimaTypes.cuh"
#include <span>

namespace Statistics {

    // Returns <slope, intersect>
    std::pair<float, float> linearFit(const std::vector<float>& x, const std::vector<float>& y);

    Float3 Mean(const std::span<Float3>& data);

    float calculateR2(const std::vector<float>& x, const std::vector<float>& y, float slope, float intercept);

    Float3 CalculateMinimaxPoint(const std::span<Float3>& data);
}