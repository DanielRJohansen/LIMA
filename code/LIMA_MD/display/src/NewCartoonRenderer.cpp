#include "NewCartoonRenderer.h"

#include "Shaders.h"
#include "TubeMesh.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace NewCartoon {

using TubeMesh::AnyPerpendicular;
using TubeMesh::BuildTube;
using TubeMesh::MinimumImage;
using TubeMesh::Sample;
using TubeMesh::SafeNormalize;
using TubeMesh::SamplePath;

namespace {

constexpr int samplesPerResidue = 6;
constexpr int tubeSides = 10;

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

			const glm::vec3 rawPosition = TubeMesh::ToVec3(positions[globalParticleId]);
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
			boundChains[run.chainIndex].positions, run.firstPoint, run.lastPoint, samplesPerResidue);
		switch (run.secondaryStructure) {
		case SecondaryStructure::Helix:
			run.vertices = BuildRibbon(samples, 0.085f, 0.025f, false);
			break;
		case SecondaryStructure::Sheet:
			run.vertices = BuildRibbon(samples, 0.075f, 0.012f, true);
			break;
		case SecondaryStructure::Coil:
			run.vertices = BuildTube(samples, 0.030f, tubeSides);
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
