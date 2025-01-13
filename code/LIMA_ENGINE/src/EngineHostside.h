#include "KernelConstants.cuh"
#include "PhysicsUtils.cuh"

constexpr std::array<float, 2 * BSPLINE_LUT_SIZE> PrecomputeBsplineTable()
{
    const int N = BSPLINE_LUT_SIZE;
    std::array<float, 2 * N> result{};

    for (int i = 0; i < N; i++)
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

// Precomputes ERFC-related scalars from 0 to cutoffNM
std::array<float, ERFC_LUT_SIZE> PrecomputeErfcForcescalarTable(float cutoffNM) {
	const float ewaldKappa = PhysicsUtils::CalcEwaldkappa(cutoffNM);

	std::array<float, ERFC_LUT_SIZE> result{};

	for (int i = 0; i < ERFC_LUT_SIZE; i++)	{
		const double fraction = static_cast<double>(i) / static_cast<double>(ERFC_LUT_SIZE - 1);
		const double correspondingDistance = fraction * cutoffNM;

        const double erfcTerm = erfc(correspondingDistance * ewaldKappa);
        const float scalar = erfcTerm + 2. * ewaldKappa / sqrt(PI) * correspondingDistance * exp(-ewaldKappa * ewaldKappa * (correspondingDistance* correspondingDistance));

		result[i] = scalar;
	}
	return result;
}

std::array<float, ERFC_LUT_SIZE> PrecomputeErfcPotentialscalarTable(float cutoffNM) {
	const float ewaldKappa = PhysicsUtils::CalcEwaldkappa(cutoffNM);

	std::array<float, ERFC_LUT_SIZE> result{};

	for (int i = 0; i < ERFC_LUT_SIZE; i++)	{
		const double fraction = static_cast<double>(i) / static_cast<double>(ERFC_LUT_SIZE - 1);
		const double correspondingDistance = fraction * cutoffNM;

		const float scalar = erfc(correspondingDistance * ewaldKappa);
		result[i] = scalar;
	}
	return result;
}

