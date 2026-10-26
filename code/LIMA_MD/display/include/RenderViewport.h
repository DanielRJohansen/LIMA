#pragma once

#include <algorithm>
#include <cmath>
#include <glm.hpp>

// Both rectangles use a top-left origin. Convert Y only at the OpenGL boundary.
struct RenderViewport {
	glm::dvec2 origin{};
	glm::dvec2 size{};
	glm::ivec2 pixelOrigin{};
	glm::ivec2 pixelSize{};

	bool Contains(glm::dvec2 point) const {
		return point.x >= origin.x && point.y >= origin.y
			&& point.x < origin.x + size.x && point.y < origin.y + size.y;
	}

	static RenderViewport Tile(int index, int count, glm::ivec2 windowSize,
		glm::ivec2 framebufferSize, double top) {
		if (count <= 0 || index < 0 || index >= count || windowSize.x <= 0 || windowSize.y <= 0
			|| framebufferSize.x <= 0 || framebufferSize.y <= 0)
			return {};
		top = std::clamp(top, 0., static_cast<double>(windowSize.y));
		const glm::dvec2 scale = glm::dvec2(framebufferSize) / glm::dvec2(windowSize);
		const int pixelTop = static_cast<int>(std::ceil(top * scale.y));
		const int height = framebufferSize.y - pixelTop;

		// All tiles have the same size. The grid is the one that fits the largest square in each tile, preferring a grid
		// the tiles fill exactly (32 tiles: 8x4, not 6x6 with 4 empty) unless that makes the tiles much smaller (7: 3x3, not 7x1)
		const auto SquareSide = [&](int candidateColumns) {
			const int candidateRows = (count + candidateColumns - 1) / candidateColumns;
			return (std::min)(static_cast<double>(framebufferSize.x) / candidateColumns, static_cast<double>(height) / candidateRows);
		};
		int columns = 1, exactColumns = 1;
		for (int candidate = 1; candidate <= count; ++candidate) {
			if (SquareSide(candidate) > SquareSide(columns)) columns = candidate;
			if (count % candidate == 0 && SquareSide(candidate) > SquareSide(exactColumns)) exactColumns = candidate;
		}
		if (SquareSide(exactColumns) >= 0.75 * SquareSide(columns))
			columns = exactColumns;
		const int rows = (count + columns - 1) / columns;
		const int column = index % columns;
		const int row = index / columns;
		const glm::ivec2 pixelStart{ column * framebufferSize.x / columns, pixelTop + row * height / rows };
		const glm::ivec2 pixelEnd{ (column + 1) * framebufferSize.x / columns, pixelTop + (row + 1) * height / rows };
		return { glm::dvec2(pixelStart) / scale, glm::dvec2(pixelEnd - pixelStart) / scale,
			pixelStart, pixelEnd - pixelStart };
	}
};
