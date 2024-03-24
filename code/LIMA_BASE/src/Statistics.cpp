#pragma once

#include "Statistics.h"

#include <assert.h>
#include <cmath>

std::pair<float, float> Statistics::linearFit(const std::vector<float>& x, const std::vector<float>& y) {
    assert(x.size() == y.size());
    float sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
    int n = x.size();

    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumXX += x[i] * x[i];
    }

    const float meanX = sumX / n;
    const float meanY = sumY / n;
    const float slope = (sumXY - n * meanX * meanY) / (sumXX - n * meanX * meanX);
    const float intercept = meanY - slope * meanX;

    return { slope, intercept };
}

float Statistics::calculateR2(const std::vector<float>& x, const std::vector<float>& y, float slope, float intercept) {
    assert(x.size() == y.size());
    float totalVar = 0, explainedVar = 0;

    float meanY = 0;
    for (auto val : y) meanY += val;
    meanY /= y.size();

    for (int i = 0; i < x.size(); i++) {
        float predictedY = slope * x[i] + intercept;
        explainedVar += std::powf(predictedY - meanY, 2);
        totalVar += std::powf(y[i] - meanY, 2);
    }

    return explainedVar / totalVar;
}
