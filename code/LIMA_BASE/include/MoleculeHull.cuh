#pragma once

#include "LimaTypes.cuh"

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

class ConvexHull {
	std::vector<Float3> vertices;	// [nm] Only used to determine CoM
	std::vector<Facet> facets;


public:
	ConvexHull() {}
	ConvexHull(const std::vector<Float3>& points);

	Float3 CenterOfMass() const {
		Float3 sum{};
		for (const Float3& v : vertices) {
			sum += v;
		}
		return sum / static_cast<float>(vertices.size() * 3);
	}

	const std::vector<Facet>& GetFacets() const { return facets; }
	const std::vector<Float3>& GetVertices() const { return vertices; }

	void Add(const Float3& vertex) {
		if (std::count(vertices.begin(), vertices.end(), vertex) == 0) {
			vertices.emplace_back(vertex);
		}
	}
	void Add(const Facet& plane) {

		facets.emplace_back(plane);
		for (const auto& v : plane.vertices) {
			Add(v);
		}
	}
};

class MoleculeHullFactory {
	std::vector<Float3> particlePositions; // [nm]
	std::vector<char> atomLetters;

public:

	ConvexHull convexHull;

	void AddParticle(const Float3& particle, char atomType) {
		particlePositions.emplace_back(particle);
		atomLetters.emplace_back(atomType);
	}
	const std::vector<Float3>& GetParticles() const {
		return particlePositions;
	}
	const std::vector<char>& GetAtomLetters() const {
		return atomLetters;
	}
};


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

	MoleculeHullCollection(const std::vector<MoleculeHullFactory>& moleculeHullFactories, Float3 boxsize);
};


















