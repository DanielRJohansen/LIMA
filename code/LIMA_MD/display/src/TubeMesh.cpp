#include "TubeMesh.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace TubeMesh {
namespace {

constexpr float pi = 3.14159265358979323846f;

glm::vec3 CatmullRom(
	const glm::vec3& p0,
	const glm::vec3& p1,
	const glm::vec3& p2,
	const glm::vec3& p3,
	float t)
{
	const float t2 = t * t;
	const float t3 = t2 * t;
	return 0.5f * ((2.f * p1)
		+ (-p0 + p2) * t
		+ (2.f * p0 - 5.f * p1 + 4.f * p2 - p3) * t2
		+ (-p0 + 3.f * p1 - 3.f * p2 + p3) * t3);
}

} // namespace

glm::vec3 ToVec3(const Float3& value)
{
	return { value.x, value.y, value.z };
}

glm::vec3 MinimumImage(glm::vec3 delta, const Float3& box)
{
	const std::array<float, 3> lengths{ box.x, box.y, box.z };
	for (int axis = 0; axis < 3; ++axis) {
		if (lengths[axis] > 0.f)
			delta[axis] -= std::round(delta[axis] / lengths[axis]) * lengths[axis];
	}
	return delta;
}

glm::vec3 SafeNormalize(const glm::vec3& value, const glm::vec3& fallback)
{
	return glm::dot(value, value) > 1e-10f ? glm::normalize(value) : fallback;
}

glm::vec3 AnyPerpendicular(const glm::vec3& direction)
{
	const glm::vec3 helper = std::abs(direction.z) < 0.9f
		? glm::vec3(0.f, 0.f, 1.f)
		: glm::vec3(0.f, 1.f, 0.f);
	return SafeNormalize(glm::cross(helper, direction), glm::vec3(1.f, 0.f, 0.f));
}

std::vector<Sample> SamplePath(
	const std::vector<glm::vec3>& positions,
	std::size_t firstPoint,
	std::size_t lastPoint,
	int samplesPerSegment)
{
	std::vector<Sample> samples;
	if (lastPoint <= firstPoint || lastPoint >= positions.size())
		return samples;

	samples.reserve((lastPoint - firstPoint) * samplesPerSegment + 1);
	for (std::size_t point = firstPoint; point < lastPoint; ++point) {
		const glm::vec3& p0 = positions[point > 0 ? point - 1 : point];
		const glm::vec3& p1 = positions[point];
		const glm::vec3& p2 = positions[point + 1];
		const glm::vec3& p3 = positions[std::min(point + 2, positions.size() - 1)];
		for (int subdivision = 0; subdivision < samplesPerSegment; ++subdivision) {
			const float t = static_cast<float>(subdivision) / samplesPerSegment;
			samples.push_back({ CatmullRom(p0, p1, p2, p3, t) });
		}
	}
	samples.push_back({ positions[lastPoint] });

	for (std::size_t i = 0; i < samples.size(); ++i) {
		const glm::vec3 before = samples[i > 0 ? i - 1 : i].position;
		const glm::vec3 after = samples[std::min(i + 1, samples.size() - 1)].position;
		samples[i].tangent = SafeNormalize(after - before, glm::vec3(0.f, 0.f, 1.f));

		const glm::vec3 previousWidth = i > 0
			? samples[i - 1].widthDirection
			: AnyPerpendicular(samples[i].tangent);
		const glm::vec3 transportedWidth = previousWidth
			- glm::dot(previousWidth, samples[i].tangent) * samples[i].tangent;
		samples[i].widthDirection = SafeNormalize(transportedWidth, AnyPerpendicular(samples[i].tangent));
	}
	return samples;
}

void AppendTube(std::vector<Vertex>& vertices, const std::vector<Sample>& samples, float radius, int sides)
{
	if (samples.size() < 2)
		return;
	vertices.reserve(vertices.size() + (samples.size() - 1) * sides * 6);

	for (std::size_t i = 0; i + 1 < samples.size(); ++i) {
		const glm::vec3 b0 = SafeNormalize(glm::cross(samples[i].tangent, samples[i].widthDirection), AnyPerpendicular(samples[i].tangent));
		const glm::vec3 b1 = SafeNormalize(glm::cross(samples[i + 1].tangent, samples[i + 1].widthDirection), b0);
		for (int side = 0; side < sides; ++side) {
			const float a0 = 2.f * pi * side / sides;
			const float a1 = 2.f * pi * (side + 1) / sides;
			const glm::vec3 n00 = std::cos(a0) * samples[i].widthDirection + std::sin(a0) * b0;
			const glm::vec3 n01 = std::cos(a1) * samples[i].widthDirection + std::sin(a1) * b0;
			const glm::vec3 n10 = std::cos(a0) * samples[i + 1].widthDirection + std::sin(a0) * b1;
			const glm::vec3 n11 = std::cos(a1) * samples[i + 1].widthDirection + std::sin(a1) * b1;

			vertices.push_back({ samples[i].position + n00 * radius, n00 });
			vertices.push_back({ samples[i].position + n01 * radius, n01 });
			vertices.push_back({ samples[i + 1].position + n11 * radius, n11 });
			vertices.push_back({ samples[i].position + n00 * radius, n00 });
			vertices.push_back({ samples[i + 1].position + n11 * radius, n11 });
			vertices.push_back({ samples[i + 1].position + n10 * radius, n10 });
		}
	}
}

std::vector<Vertex> BuildTube(const std::vector<Sample>& samples, float radius, int sides)
{
	std::vector<Vertex> vertices;
	AppendTube(vertices, samples, radius, sides);
	return vertices;
}

void AppendRoundCap(std::vector<Vertex>& vertices, const Sample& end, float radius, int sides, int rings)
{
	// Same frame as AppendTube's last ring, so the cap meets the tube wall exactly
	const glm::vec3 b = SafeNormalize(glm::cross(end.tangent, end.widthDirection), AnyPerpendicular(end.tangent));
	const auto Normal = [&](int ring, int side) {
		const float polar = 0.5f * pi * ring / rings;	// 0 at the tube wall, pi/2 at the tip
		const float azimuth = 2.f * pi * side / sides;
		const glm::vec3 radial = std::cos(azimuth) * end.widthDirection + std::sin(azimuth) * b;
		return std::cos(polar) * radial + std::sin(polar) * end.tangent;
	};
	for (int ring = 0; ring < rings; ++ring) {
		for (int side = 0; side < sides; ++side) {
			const glm::vec3 n00 = Normal(ring, side);
			const glm::vec3 n01 = Normal(ring, side + 1);
			const glm::vec3 n10 = Normal(ring + 1, side);
			const glm::vec3 n11 = Normal(ring + 1, side + 1);
			vertices.push_back({ end.position + n00 * radius, n00 });
			vertices.push_back({ end.position + n01 * radius, n01 });
			vertices.push_back({ end.position + n11 * radius, n11 });
			vertices.push_back({ end.position + n00 * radius, n00 });
			vertices.push_back({ end.position + n11 * radius, n11 });
			vertices.push_back({ end.position + n10 * radius, n10 });
		}
	}
}

} // namespace TubeMesh
