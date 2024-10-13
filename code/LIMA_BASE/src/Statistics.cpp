#pragma once

#include "Statistics.h"
#include <numeric>
#include <assert.h>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <execution>


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

//Float3 Statistics::Mean(const std::span<const Float3>& values) {
//    Float3 sum = std::reduce(std::execution::par, values.begin(), values.end(), Float3{ 0, 0, 0 },
//        [](const Float3& a, const Float3& b) -> Float3 {
//            return { a.x + b.x, a.y + b.y, a.z + b.z };
//        });
//    return sum / values.size();
//}

Float3 Statistics::Mean(const std::span<const Float3>& values) {
    // Use Double3 for more precision in the accumulation
    Double3 sum = std::reduce(std::execution::par, values.begin(), values.end(), Double3{ 0, 0, 0 },
        [](const Double3& a, const Double3& b) -> Double3 {
            return { a.x + b.x, a.y + b.y, a.z + b.z };
        });

    // Convert back to Float3 for the final result
    return (sum / static_cast<double>(values.size()) ).toFloat3();
}

Float3 Statistics::CalculateMinimaxPoint(const std::span<const Float3>& points) {
    float epsilon = 1e-4;

    // Initial guess for the new point as the centroid
    //Float3 currentPoint = Mean(points);
    Float3 currentPoint = points[0];

    float prevVelocity = std::numeric_limits<float>::max();

    const int maxIterations = 5;
    for (int i = 0; i < maxIterations; i++) {
        Float3 downGradient{};
        float maxDist = 0.0f;
        float secondMaxDist = 0.0f;

        for (const auto& point : points) {
            float dist = (currentPoint-point).len();
            if (dist > maxDist) {
                secondMaxDist = maxDist;
                maxDist = dist;
                downGradient = (point - currentPoint);
            }
        }

        const float velocity = std::min((maxDist) / 2.f, prevVelocity * 0.9f);
        prevVelocity = velocity;

        currentPoint = currentPoint + downGradient * velocity;


        if (velocity < epsilon) {
            break;
        }
    }

    return currentPoint;
}

float Statistics::Mean(const std::vector<float>& vec)
{
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    return static_cast<float>(sum / static_cast<double>(vec.size()));
}

float Statistics::StdDev(const std::vector<float>& vec) {
    if (vec.empty()) { return 0.f; }

    const double mean = Mean(vec);

    double variance = std::accumulate(vec.begin(), vec.end(), 0.0,
        [mean](double acc, float elem) {
            double diff = elem - mean;
            return acc + diff * diff;
        });

    const double deviation = variance / static_cast<double>(vec.size());
    return static_cast<float>(std::sqrt(deviation));
}

float Statistics::Max(const float* const data, size_t n) {
    return *std::max_element(std::execution::par, data, data + n);
}

float Statistics::MaxLen(const Float3* const data, size_t n) {
    const Float3 maxElem = *std::max_element(std::execution::par, data, data + n,
        [](const Float3& a, const Float3& b) {
            return a.lenSquared() < b.lenSquared();
        }
    );
    return maxElem.len();
}

double Statistics::Sumd(const float* const data, size_t n) {
	return std::reduce(std::execution::par, data, data + n, 0.0);
}
