#pragma once

#include "vector"
#include "LimaTypes.cuh"

namespace Statistics {

    // Returns <slope, intersect>
    std::pair<float, float> linearFit(const std::vector<float>& x, const std::vector<float>& y);

    Float3 Mean(const std::vector<Float3>& data);

    float calculateR2(const std::vector<float>& x, const std::vector<float>& y, float slope, float intercept);


}