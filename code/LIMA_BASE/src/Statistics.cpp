#pragma once

#include "Statistics.h"

#include <assert.h>
#include <cmath>
#include <math.h>

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
        explainedVar += std::pow(predictedY - meanY, 2.f);
        totalVar += std::pow(y[i] - meanY, 2.f);
    }

    return explainedVar / totalVar;
}

Float3 Statistics::Mean(const std::span<const Float3>& values) {
	Float3 sum = { 0, 0, 0 };
	for (auto val : values) {
		sum += val;
	}
	return sum / values.size();
}

Float3 Statistics::CalculateMinimaxPoint(const std::span<const Float3>& points) {
    float epsilon = 1e-4;

    // Initial guess for the new point as the centroid
    Float3 currentPoint = Mean(points);


    while (true) {
        Float3 downGradient = { 0, 0, 0 };
        float maxDist = 0.0f;
        float secondMaxDist = 0.0f;
        float prevVelocity = std::numeric_limits<float>::max();

        for (const auto& point : points) {
            float dist = (currentPoint-point).len();
            if (dist > maxDist) {
                secondMaxDist = maxDist;
                maxDist = dist;
                downGradient = (point - currentPoint);
            }
        }

        const float velocity = std::min((maxDist - secondMaxDist) / 2.f, prevVelocity);
        prevVelocity = velocity;

        currentPoint = currentPoint + downGradient * velocity;


        if (velocity < epsilon) {
            break;
        }
    }

    return currentPoint;
}