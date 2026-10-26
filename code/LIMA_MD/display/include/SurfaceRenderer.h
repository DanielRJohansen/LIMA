#pragma once

#include "Bodies.cuh"

#include <glm.hpp>
#include <chrono>
#include <vector>

class DrawSurfaceShader;

// The surface representation: a smooth molecular surface, in the style of VMD's QuickSurf. Each atom adds a Gaussian
// density on a grid, and the mesh is the density's isosurface, which for a lone atom lies at its radius
class SurfaceRenderer {
public:
	struct Atom {
		int id = -1;		// Global atom id, for picking
		float radius = 0;	// [nm]
		glm::vec4 color{};
	};

	struct Vertex {
		glm::vec3 position{};
		glm::vec3 normal{};
		glm::vec4 color{};
		int atomId = -1;	// The atom contributing most density here
	};

	struct Mesh {
		std::vector<Vertex> vertices;
		std::vector<unsigned int> indices;	// Triangles
	};

	static constexpr float gridSpacing = 0.1f;	// [nm]

	// Pure geometry, without OpenGL, so it can be tested
	static Mesh BuildMesh(const std::vector<Atom>& atoms, const std::vector<Float3>& positions);

	SurfaceRenderer();
	~SurfaceRenderer();
	SurfaceRenderer(const SurfaceRenderer&) = delete;
	SurfaceRenderer& operator=(const SurfaceRenderer&) = delete;

	void Prepare(std::vector<Atom> atoms, const std::vector<Float3>& positions);
	// Rebuilds the mesh for new positions, but skips rebuilds of large systems that come too often
	void Update(const std::vector<Float3>& positions);
	void Clear();
	void Draw(DrawSurfaceShader& shader, const glm::mat4& view, const glm::mat4& projection) const;

	[[nodiscard]] bool HasGeometry() const { return indexCount > 0; }

private:
	std::vector<Atom> atoms;
	unsigned int vao = 0;
	unsigned int vertexBuffer = 0;
	unsigned int indexBuffer = 0;
	int indexCount = 0;
	std::chrono::steady_clock::time_point lastBuildEnd{};
	std::chrono::steady_clock::duration lastBuildDuration{};

	void Rebuild(const std::vector<Float3>& positions);
};
