#include "LimaTypes.cuh"
#include "RenderUtilities.cuh"
#include <limits>

RenderAtom::RenderAtom(Float3 pos, Float3 boxSize, char atomLetter) {
	const float radius = RenderUtilities::getRadius(RenderUtilities::RAS_getTypeFromAtomletter(atomLetter));
	position = float4{ pos.x, pos.y, pos.z, radius };
	color = RenderUtilities::getColor(RenderUtilities::RAS_getTypeFromAtomletter(atomLetter));
}

BoundingBox::BoundingBox(const std::vector<Float3>& points) {
	min = Float3{std::numeric_limits<float>::max()};
	max = Float3{std::numeric_limits<float>::min()};
	for (const Float3& p : points) {
		for (int dim = 0; dim < 3; dim++) {
			min[dim] = std::min(min[dim], p[dim]);
			max[dim] = std::max(max[dim], p[dim]);
		}
	}
}