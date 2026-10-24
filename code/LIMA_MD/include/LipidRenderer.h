#pragma once

#include "Bodies.cuh"
#include "LipidModel.h"
#include "RenderCommons.h"

#include <glm.hpp>
#include <vector>

class DrawTrianglesShader;

// The lipid representation: tails are drawn here as tubes, the head is drawn as an enlarged atom by the atom shader
class LipidRenderer {
public:
	static constexpr float headRadius = 0.24f; // [nm]
	static glm::vec4 HeadColor(LipidHeadgroup);

	void Prepare(const LipidModels& lipids, const std::vector<Float3>& positions, Float3 boxSize);

	// Position-only fast path
	void Update(const std::vector<Float3>& positions);
	void Clear();
	void Draw(DrawTrianglesShader& shader, const glm::mat4& viewProjection) const;

	[[nodiscard]] bool HasGeometry() const { return !meshes.empty(); }

private:
	struct Mesh {
		std::size_t firstVertex = 0;	// Into vertices
		std::size_t vertexCount = 0;
		glm::vec4 color{};
		int objectId = -1;
	};

	Float3 boxSize{};
	LipidModels lipids;
	std::vector<Mesh> meshes;	// One per lipid
	std::vector<Vertex> vertices;	// Of all lipids
};
