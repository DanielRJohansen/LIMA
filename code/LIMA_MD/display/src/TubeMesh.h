#pragma once

#include "Bodies.cuh"
#include "RenderCommons.h"

#include <glm.hpp>
#include <vector>

// Smooth tubes through a list of points, shared by the cartoon and lipid representations
namespace TubeMesh {

struct Sample {
	glm::vec3 position{};
	glm::vec3 tangent{};
	glm::vec3 widthDirection{};
};

glm::vec3 ToVec3(const Float3& value);
glm::vec3 MinimumImage(glm::vec3 delta, const Float3& box);
glm::vec3 SafeNormalize(const glm::vec3& value, const glm::vec3& fallback);
glm::vec3 AnyPerpendicular(const glm::vec3& direction);

// Catmull-Rom spline through positions[firstPoint..lastPoint], with a parallel transported frame
std::vector<Sample> SamplePath(
	const std::vector<glm::vec3>& positions,
	std::size_t firstPoint,
	std::size_t lastPoint,
	int samplesPerSegment);

void AppendTube(std::vector<Vertex>& vertices, const std::vector<Sample>& samples, float radius, int sides);
std::vector<Vertex> BuildTube(const std::vector<Sample>& samples, float radius, int sides);
// Closes a tube end with a hemisphere bulging along the sample's tangent
void AppendRoundCap(std::vector<Vertex>& vertices, const Sample& end, float radius, int sides, int rings);

} // namespace TubeMesh
