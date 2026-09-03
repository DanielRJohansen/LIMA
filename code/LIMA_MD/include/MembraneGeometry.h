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

	// The three radii describe the bilayer mid-surface along x, y, and z.
	struct Ellipsoid {
		Float3 center;
		Float3 radii;
	};

	using Figure = std::variant<Plane, Sphere, Ellipsoid>;
}
