#pragma once

#include "Simulation.cuh"
#include <optional>
#include <string_view>
#include <utility>

namespace EngineBatch {
	// Keep scheduling and Engine validation in agreement as SimParams grows.
	inline std::optional<std::string_view> FindIncompatibility(const Simulation& lhs, const Simulation& rhs) {
		if (!lhs.box || !rhs.box) return "missing box";
		if (lhs.box->boxparams.boxSize != rhs.box->boxparams.boxSize) return "box size";
		const auto reference = lhs.simParams.Params();
		const auto candidate = rhs.simParams.Params();
		std::optional<std::string_view> mismatch;
		[&]<size_t... indices>(std::index_sequence<indices...>) {
			([&] {
				const auto& param = std::get<indices>(candidate);
				if (!mismatch && param.name != "dt" && param.name != "n_steps"
					&& param.value != std::get<indices>(reference).value)
					mismatch = param.name;
			}(), ...);
		}(std::make_index_sequence<std::tuple_size_v<decltype(reference)>>{});
		return mismatch;
	}
}
