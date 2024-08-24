#include "LimaTypes.cuh"
#include "RenderUtilities.cuh"
#include <limits>

RenderAtom::RenderAtom(Float3 pos, Float3 boxSize, char atomLetter) {
	const Float3 posNorm = pos / boxSize - 0.5f;
	const float radiusNorm = RenderUtilities::getRadius(RenderUtilities::RAS_getTypeFromAtomletter(atomLetter)) / boxSize.x;
	position = float4{ posNorm.x, posNorm.y, posNorm.z, radiusNorm };
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