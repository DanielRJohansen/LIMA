#pragma once

#include <glm.hpp>

struct Vertex {
	glm::vec3 position;
	glm::vec3 normal;
};


enum class UniqueRenderElementIds {
	gizmoArrowX = 1 << 30 + 0,
	gizmoArrowY,
	gizmoArrowZ
};

static bool ElementIdIsAtomid(int id) {
	return (id & (1 << 30)) == 0; // If the highest bit is not set, it's an atom id
}