#pragma once

#include <glm.hpp>

struct Vertex {
	glm::vec3 position;
	glm::vec3 normal;
};


enum class UniqueRenderElementIds {
	gizmoArrowX = 1 << 30 + 0,
	gizmoArrowY = 1 << 30 + 1,
	gizmoArrowZ = 1 << 30 + 2
};