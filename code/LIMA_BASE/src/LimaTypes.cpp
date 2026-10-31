#include "LimaTypes.cuh"
#include "RenderUtilities.cuh"
#include <limits>

RenderAtom::RenderAtom(Float3 pos, Float3 boxSize, char atomLetter) {
	const float radius = RenderUtilities::getRadius(RenderUtilities::RAS_getTypeFromAtomletter(atomLetter));
	position = float4{ pos.x, pos.y, pos.z, radius };
	color = RenderUtilities::getColor(RenderUtilities::RAS_getTypeFromAtomletter(atomLetter));
}