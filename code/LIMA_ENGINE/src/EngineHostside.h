#include "KernelConstants.cuh"

constexpr std::array<float, 2 * BSPLINE_LUT_SIZE> PrecomputeBsplineTable()
{
    const int N = BSPLINE_LUT_SIZE;
    std::array<float, 2 * N> result{};

    for (int i = 0; i < N; ++i)
    {
        const double f = static_cast<double>(i) / static_cast<double>(N-1);

        // w0 = (1 - f)^3 / 6
        const double w0 = (1 - f) * (1 - f) * (1. - f) / 6.;
        const double w1 = (4. - 6. * f * f + 3. * f * f * f) / 6.;

        // Store in array: [w0, w1]
        result[i] = static_cast<float>(w0);
        result[N + i] = static_cast<float>(w1);
    }

    return result;
}