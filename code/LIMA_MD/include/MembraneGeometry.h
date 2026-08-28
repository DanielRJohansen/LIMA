#pragma once

#include "LimaTypes.cuh"

#include <variant>

namespace MembraneGeometry {
	// The z coordinate describes the bilayer mid-plane.
	struct Plane {
		float z;
	};

	// The radius describes the bilayer mid-surface. Lipids are placed on both
	// sides of it, with their heads pointing away from the hydrophobic core.
	struct Sphere {
		Float3 center;
		float radius;
	};

	using Figure = std::variant<Plane, Sphere>;
}
