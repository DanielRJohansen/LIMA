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
		const int columns = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(count))));
		const int rows = (count + columns - 1) / columns;
		top = std::clamp(top, 0., static_cast<double>(windowSize.y));
		const glm::dvec2 scale = glm::dvec2(framebufferSize) / glm::dvec2(windowSize);
		const int pixelTop = static_cast<int>(std::ceil(top * scale.y));
		const int height = framebufferSize.y - pixelTop;
		const int column = index % columns;
		const int row = index / columns;
		const glm::ivec2 pixelStart{ column * framebufferSize.x / columns, pixelTop + row * height / rows };
		const glm::ivec2 pixelEnd{ (column + 1) * framebufferSize.x / columns, pixelTop + (row + 1) * height / rows };
		return { glm::dvec2(pixelStart) / scale, glm::dvec2(pixelEnd - pixelStart) / scale,
			pixelStart, pixelEnd - pixelStart };
	}
};
