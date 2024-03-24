#pragma once

#include "vector"

namespace Statistics {

    // Returns <slope, intersect>
    std::pair<float, float> linearFit(const std::vector<float>& x, const std::vector<float>& y);


    float calculateR2(const std::vector<float>& x, const std::vector<float>& y, float slope, float intercept);
}