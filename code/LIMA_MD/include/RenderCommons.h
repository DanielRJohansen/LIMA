#pragma once

#include <glm.hpp>

struct Vertex {
	glm::vec3 position;
	glm::vec3 normal;
};


enum class UniqueRenderElementIds {
	gizmoArrowX = 1 << 30 + 0,
	gizmoArrowY,
	gizmoArrowZ,
	gizmoRotateX,
	gizmoRotateY,
	gizmoRotateZ
};

static bool ElementIdIsAtomid(int id) {
	return (id & (1 << 30)) == 0; // If the highest bit is not set, it's an atom id
}




namespace ColorScheme {
	static const glm::vec4 backgroundBot(0.07f, 0.065f, 0.06f, 1.f);
	static const glm::vec4 backgroundTop(0.10f, 0.095f, 0.12f, 1.f);
}