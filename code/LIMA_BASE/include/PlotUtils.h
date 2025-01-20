#pragma once

#include <vector>
#include <string>
#include <optional>


namespace PlotUtils {
	void PlotData(const std::vector<std::vector<float>>& vectors, const std::vector<std::string>& labels,
		const std::optional<std::vector<float>>& x_axis = std::nullopt);
};