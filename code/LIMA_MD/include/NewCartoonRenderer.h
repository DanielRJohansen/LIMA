#pragma once

#include "Backbone.h"
#include "Bodies.cuh"
#include "RenderCommons.h"

#include <glm.hpp>
#include <vector>

class DrawTrianglesShader;

namespace NewCartoon {

class Renderer {
public:
	// Bind immutable global backbone ids to renderer-owned global positions.
	void Prepare(
		const BackboneChains& backboneChains,
		const std::vector<Float3>& positions,
		Float3 boxSize);

	// Position-only fast path. The secondary-structure assignment is preserved.
	void Update(const std::vector<Float3>& positions);
	void Clear();
	void Draw(DrawTrianglesShader& shader, const glm::mat4& viewProjection) const;

	[[nodiscard]] bool HasGeometry() const { return !drawableRuns.empty(); }

private:
	struct BoundPoint {
		int globalParticleId = -1;
		SecondaryStructure secondaryStructure = SecondaryStructure::Coil;
	};

	struct BoundChain {
		std::vector<BoundPoint> points;
		std::vector<glm::vec3> positions;
	};

	struct DrawableRun {
		std::size_t chainIndex = 0;
		std::size_t firstPoint = 0;
		std::size_t lastPoint = 0;
		SecondaryStructure secondaryStructure = SecondaryStructure::Coil;
		std::vector<Vertex> vertices;
		int objectId = -1;
	};

	void BuildDrawableRuns();
	void RebuildMeshes();

	Float3 boxSize{};
	std::size_t positionCount = 0;
	std::vector<BoundChain> boundChains;
	std::vector<DrawableRun> drawableRuns;
};

} // namespace NewCartoon
