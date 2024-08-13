#include "LimaTypes.cuh"
#include "RenderUtilities.cuh"

RenderAtom::RenderAtom(Float3 pos, Float3 boxSize, char atomLetter) {
	const Float3 posNorm = pos / boxSize - 0.5f;
	position = float4{ posNorm.x, posNorm.y, posNorm.z, RenderUtilities::getRadius(RenderUtilities::RAS_getTypeFromAtomletter(atomLetter)) };
	color = RenderUtilities::getColor(RenderUtilities::RAS_getTypeFromAtomletter(atomLetter));
}