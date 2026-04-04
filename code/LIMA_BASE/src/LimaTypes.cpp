#include "LimaTypes.cuh"
#include "RenderUtilities.cuh"
#include <limits>

RenderAtom::RenderAtom(Float3 pos, Float3 boxSize, char atomLetter) {
	const float radius = RenderUtilities::getRadius(RenderUtilities::RAS_getTypeFromAtomletter(atomLetter));
	position = float4{ pos.x, pos.y, pos.z, radius };
	color = RenderUtilities::getColor(RenderUtilities::RAS_getTypeFromAtomletter(atomLetter));
}

constexpr BoundingBox::BoundingBox(const std::vector<Float3>& points) {
	min = Float3{std::numeric_limits<float>::max()};
	max = Float3{std::numeric_limits<float>::min()};
	for (const Float3& p : points) {
		min.x = std::min(min.x, p.x);
		min.y = std::min(min.y, p.y);
		min.z = std::min(min.z, p.z);
		max.x = std::max(max.x, p.x);
		max.y = std::max(max.y, p.y);
		max.z = std::max(max.z, p.z);
	}
}
