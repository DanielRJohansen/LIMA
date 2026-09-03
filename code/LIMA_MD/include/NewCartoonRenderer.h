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
	// Bind immutable global backbone ids to the simulation's packed position layout.
	void Prepare(
		const BackboneChains& backboneChains,
		const std::vector<PersistentCluster>& pclusters,
		const std::vector<PersistentClusterMeta>& pcMeta,
		Float3 boxSize);

	// Position-only fast path. The secondary-structure assignment is preserved.
	void Update(const Float3* packedPositions);
	void Clear();
	void Draw(DrawTrianglesShader& shader, const glm::mat4& viewProjection) const;

	[[nodiscard]] bool HasGeometry() const { return !drawableRuns.empty(); }

private:
	struct BoundPoint {
		int globalParticleId = -1;
		int packedPositionIndex = -1;
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
	std::size_t packedPositionCount = 0;
	std::vector<BoundChain> boundChains;
	std::vector<DrawableRun> drawableRuns;
};

} // namespace NewCartoon
