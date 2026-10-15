#pragma once

namespace Virial {
	// Scalar configurational virial [J/mol]. Integer partials make their later reduction independent of block scheduling.
	constexpr double fixedPointScale = 1048576.0; // 2^20
	constexpr double fixedPointScaleInv = 1.0 / fixedPointScale;

	__device__ inline long long ToFixed(float value) {
		return llrint(static_cast<double>(value) * fixedPointScale);
	}
}
