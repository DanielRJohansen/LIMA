#include "NewCartoonRenderer.h"

#include "Shaders.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace NewCartoon {
namespace {

constexpr int samplesPerResidue = 6;
constexpr int tubeSides = 10;
constexpr float pi = 3.14159265358979323846f;

struct Sample {
	glm::vec3 position{};
	glm::vec3 tangent{};
	glm::vec3 widthDirection{};
};

glm::vec3 ToCartoonVec3(const Float3& value)
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

std::vector<Sample> SamplePath(
	const std::vector<glm::vec3>& positions,
	std::size_t firstPoint,
	std::size_t lastPoint)
{
	std::vector<Sample> samples;
	if (lastPoint <= firstPoint || lastPoint >= positions.size())
		return samples;

	samples.reserve((lastPoint - firstPoint) * samplesPerResidue + 1);
	for (std::size_t point = firstPoint; point < lastPoint; ++point) {
		const glm::vec3& p0 = positions[point > 0 ? point - 1 : point];
		const glm::vec3& p1 = positions[point];
		const glm::vec3& p2 = positions[point + 1];
		const glm::vec3& p3 = positions[std::min(point + 2, positions.size() - 1)];
		for (int subdivision = 0; subdivision < samplesPerResidue; ++subdivision) {
			const float t = static_cast<float>(subdivision) / samplesPerResidue;
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

void AppendTriangle(
	std::vector<Vertex>& vertices,
	const glm::vec3& a,
	const glm::vec3& b,
	const glm::vec3& c,
	const glm::vec3& normal)
{
	vertices.push_back({ a, normal });
	vertices.push_back({ b, normal });
	vertices.push_back({ c, normal });
}

void AppendQuad(
	std::vector<Vertex>& vertices,
	const glm::vec3& a,
	const glm::vec3& b,
	const glm::vec3& c,
	const glm::vec3& d,
	const glm::vec3& normal)
{
	AppendTriangle(vertices, a, b, c, normal);
	AppendTriangle(vertices, a, c, d, normal);
}

std::vector<Vertex> BuildTube(const std::vector<Sample>& samples, float radius)
{
	std::vector<Vertex> vertices;
	if (samples.size() < 2)
		return vertices;
	vertices.reserve((samples.size() - 1) * tubeSides * 6);

	for (std::size_t i = 0; i + 1 < samples.size(); ++i) {
		const glm::vec3 b0 = SafeNormalize(glm::cross(samples[i].tangent, samples[i].widthDirection), AnyPerpendicular(samples[i].tangent));
		const glm::vec3 b1 = SafeNormalize(glm::cross(samples[i + 1].tangent, samples[i + 1].widthDirection), b0);
		for (int side = 0; side < tubeSides; ++side) {
			const float a0 = 2.f * pi * side / tubeSides;
			const float a1 = 2.f * pi * (side + 1) / tubeSides;
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
	return vertices;
}

float SheetArrowScale(float fraction)
{
	if (fraction < 0.70f)
		return 1.f;
	if (fraction < 0.82f)
		return glm::mix(1.f, 1.55f, (fraction - 0.70f) / 0.12f);
	return glm::mix(1.55f, 0.06f, (fraction - 0.82f) / 0.18f);
}

std::vector<Vertex> BuildRibbon(
	const std::vector<Sample>& samples,
	float halfWidth,
	float halfThickness,
	bool arrow)
{
	std::vector<Vertex> vertices;
	if (samples.size() < 2)
		return vertices;
	vertices.reserve((samples.size() - 1) * 24 + 12);

	struct CrossSection {
		glm::vec3 topLeft;
		glm::vec3 topRight;
		glm::vec3 bottomLeft;
		glm::vec3 bottomRight;
		glm::vec3 width;
		glm::vec3 thickness;
	};
	std::vector<CrossSection> sections;
	sections.reserve(samples.size());
	for (std::size_t i = 0; i < samples.size(); ++i) {
		const float fraction = static_cast<float>(i) / static_cast<float>(samples.size() - 1);
		const float width = halfWidth * (arrow ? SheetArrowScale(fraction) : 1.f);
		const glm::vec3 widthVector = samples[i].widthDirection * width;
		const glm::vec3 thicknessDirection = SafeNormalize(
			glm::cross(samples[i].tangent, samples[i].widthDirection), AnyPerpendicular(samples[i].tangent));
		const glm::vec3 thicknessVector = thicknessDirection * halfThickness;
		sections.push_back({
			samples[i].position - widthVector + thicknessVector,
			samples[i].position + widthVector + thicknessVector,
			samples[i].position - widthVector - thicknessVector,
			samples[i].position + widthVector - thicknessVector,
			samples[i].widthDirection,
			thicknessDirection,
		});
	}

	for (std::size_t i = 0; i + 1 < sections.size(); ++i) {
		const CrossSection& a = sections[i];
		const CrossSection& b = sections[i + 1];
		AppendQuad(vertices, a.topLeft, a.topRight, b.topRight, b.topLeft, SafeNormalize(a.thickness + b.thickness, a.thickness));
		AppendQuad(vertices, a.bottomRight, a.bottomLeft, b.bottomLeft, b.bottomRight, -SafeNormalize(a.thickness + b.thickness, a.thickness));
		AppendQuad(vertices, a.bottomLeft, a.topLeft, b.topLeft, b.bottomLeft, -SafeNormalize(a.width + b.width, a.width));
		AppendQuad(vertices, a.topRight, a.bottomRight, b.bottomRight, b.topRight, SafeNormalize(a.width + b.width, a.width));
	}

	const CrossSection& first = sections.front();
	const CrossSection& last = sections.back();
	AppendQuad(vertices, first.bottomLeft, first.bottomRight, first.topRight, first.topLeft, -samples.front().tangent);
	AppendQuad(vertices, last.topLeft, last.topRight, last.bottomRight, last.bottomLeft, samples.back().tangent);
	return vertices;
}

glm::vec4 StructureColor(SecondaryStructure structure)
{
	switch (structure) {
	case SecondaryStructure::Helix: return { 0.88f, 0.20f, 0.30f, 1.f };
	case SecondaryStructure::Sheet: return { 0.95f, 0.72f, 0.16f, 1.f };
	case SecondaryStructure::Coil: return { 0.24f, 0.66f, 0.86f, 1.f };
	}
	return { 0.8f, 0.8f, 0.8f, 1.f };
}

} // namespace

void Renderer::Prepare(
	const BackboneChains& backboneChains,
	const std::vector<Float3>& positions,
	Float3 newBoxSize)
{
	Clear();
	boxSize = newBoxSize;
	positionCount = positions.size();

	for (const BackboneChain& sourceChain : backboneChains) {
		BoundChain chain;
		chain.points.reserve(sourceChain.points.size());
		for (const BackbonePoint& point : sourceChain.points) {
			if (point.particleId < 0 || static_cast<std::size_t>(point.particleId) >= positionCount) {
				if (chain.points.size() >= 2)
					boundChains.push_back(std::move(chain));
				chain = {};
				continue;
			}
			chain.points.push_back({ point.particleId, point.secondaryStructure });
		}
		if (chain.points.size() >= 2)
			boundChains.push_back(std::move(chain));
	}

	BuildDrawableRuns();
	Update(positions);
}

void Renderer::BuildDrawableRuns()
{
	for (std::size_t chainIndex = 0; chainIndex < boundChains.size(); ++chainIndex) {
		const BoundChain& chain = boundChains[chainIndex];
		const auto TypeOfEdge = [&chain](std::size_t edge) {
			const SecondaryStructure left = chain.points[edge].secondaryStructure;
			return left == chain.points[edge + 1].secondaryStructure ? left : SecondaryStructure::Coil;
		};

		std::size_t firstEdge = 0;
		while (firstEdge + 1 < chain.points.size()) {
			const SecondaryStructure type = TypeOfEdge(firstEdge);
			std::size_t lastEdge = firstEdge;
			while (lastEdge + 2 < chain.points.size() && TypeOfEdge(lastEdge + 1) == type)
				++lastEdge;

			drawableRuns.push_back({
				chainIndex, firstEdge, lastEdge + 1, type, {}, chain.points[firstEdge].globalParticleId
			});
			firstEdge = lastEdge + 1;
		}
	}
}

void Renderer::Update(const std::vector<Float3>& positions)
{
	if (positions.empty() || positionCount == 0)
		return;

	for (BoundChain& chain : boundChains) {
		chain.positions.resize(chain.points.size());
		glm::vec3 previousRawPosition{};
		for (std::size_t i = 0; i < chain.points.size(); ++i) {
			const int globalParticleId = chain.points[i].globalParticleId;
			if (globalParticleId < 0 || static_cast<std::size_t>(globalParticleId) >= positions.size())
				continue;

			const glm::vec3 rawPosition = ToCartoonVec3(positions[globalParticleId]);
			if (i == 0)
				chain.positions[i] = rawPosition;
			else
				chain.positions[i] = chain.positions[i - 1] + MinimumImage(rawPosition - previousRawPosition, boxSize);
			previousRawPosition = rawPosition;
		}
	}
	RebuildMeshes();
}

void Renderer::RebuildMeshes()
{
	for (DrawableRun& run : drawableRuns) {
		const std::vector<Sample> samples = SamplePath(
			boundChains[run.chainIndex].positions, run.firstPoint, run.lastPoint);
		switch (run.secondaryStructure) {
		case SecondaryStructure::Helix:
			run.vertices = BuildRibbon(samples, 0.085f, 0.025f, false);
			break;
		case SecondaryStructure::Sheet:
			run.vertices = BuildRibbon(samples, 0.075f, 0.012f, true);
			break;
		case SecondaryStructure::Coil:
			run.vertices = BuildTube(samples, 0.030f);
			break;
		}
	}
}

void Renderer::Clear()
{
	positionCount = 0;
	boundChains.clear();
	drawableRuns.clear();
}

void Renderer::Draw(DrawTrianglesShader& shader, const glm::mat4& viewProjection) const
{
	const glm::mat4 identity(1.f);
	for (const DrawableRun& run : drawableRuns)
		shader.Draw(run.vertices, viewProjection, identity, StructureColor(run.secondaryStructure), run.objectId);
}

} // namespace NewCartoon
