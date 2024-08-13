#pragma once

#include "LimaTypes.cuh"


struct Tri {
	std::array<Float3, 3> vertices;
	Float3 normal;
};
static_assert(sizeof(Tri) % 16 == 0);

struct Facet {
	std::array<Float3, 3> vertices;
	Float3 normal;
	double_t distFromOrigo;
	char _[8];

	Float3 intersectionPoint(Float3 p1, Float3 p2) const {
		//Return the intersection point of a line passing two points and this plane
		return p1 + (p2 - p1) * (-distance(p1) / normal.dot(p2 - p1));
	};
	void invert() { normal *= -1.f; }
	double_t distance(Float3 point) const { return normal.dot(point) + distFromOrigo; }
};
static_assert(sizeof(Facet) % 16 == 0);

template <int maxFacets>
class ConvexHull {
	std::array<Float3, maxFacets * 2> vertices;	// [nm] Only used to determine CoM
	std::array<Facet, maxFacets> facets;


public:
	int numFacets = 0;
	int numVertices = 0;

	ConvexHull() {}
	ConvexHull(const std::vector<Float3>& points);

	Float3 CenterOfMass() const {
		Float3 sum{};
		for (const Float3& v : vertices) {
			sum += v;
		}
		return sum / static_cast<float>(vertices.size());
	}

	const std::span<const Facet>& GetFacets() const { return { facets.begin(), facets.begin() + numFacets }; }
	const std::span<const Float3>& GetVertices() const {
		return { vertices.begin(), vertices.begin() + numVertices };
	}

	void Add(const Float3& vertex) {
		if (numVertices >= maxFacets * 2) {
			throw std::runtime_error("Too many vertices in ConvexHull");
		}

		if (std::count(vertices.begin(), vertices.begin() + numVertices, vertex) == 0) {
			vertices[numVertices++] = vertex;
		}
	}
	void Add(const Facet& plane) {
		if (numFacets >= maxFacets) {
			throw std::runtime_error("Too many facets in ConvexHull");
		}
		facets[numFacets++] = plane;
		for (const auto& v : plane.vertices) {
			Add(v);
		}
	}
};

template <int maxFacets, int maxParticles>
class MoleculeContainer {
	std::array<Float3, maxParticles> particlePositions; // [nm]

public:
	int nParticles = 0;

	ConvexHull<maxFacets> convexHull;

	void AddParticle(const Float3& particle) {
		if (nParticles >= maxParticles) {
			throw std::runtime_error("Too many particles in MoleculeContainer");
		}
		particlePositions[nParticles++] = particle;
	}
	const std::span<const Float3>& GetParticles() const {
		return { particlePositions.begin(), particlePositions.begin() + nParticles };
	}
};

using MoleculeContainerSmall = MoleculeContainer<64, 256>;	// should be 32 after triangle culling


struct MoleculeHull {
	int indexOfFirstFacetInBuffer = -1;
	int nFacets = -1;
	int indexOfFirstParticleInBuffer = -1;
	int nParticles = -1;
};

// A collection that stays on ho st, but the buffers are on device
struct MoleculeHullCollection 
{
	// Device memory
	Facet* facets = nullptr;
	RenderAtom* particles = nullptr; // [normalized -0.5 -> 0.5] Only to be used for rendering
	MoleculeHull* moleculeHulls = nullptr;

	// Hostside
	int nFacets = 0;
	int nParticles = 0;
	int nMoleculeHulls = 0;

	MoleculeHullCollection(const std::vector<MoleculeContainerSmall>& moleculeHullFactories, Float3 boxsize);
};


















