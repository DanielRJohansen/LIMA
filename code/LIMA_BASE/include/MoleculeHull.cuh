#pragma once

#include "LimaTypes.cuh"
#include <functional>


enum FacetDrawMode {
	FACES = 0,
	EDGES = 1
};

struct Facet {
	std::array<Float3, 3> vertices;
	Float3 normal;
	float D; // distance from origo along the normal, so it can be negative
	char __[4];

	char _[8];

	Float3 intersectionPoint(Float3 p1, Float3 p2) const {
		//Return the intersection point of a line passing two points and this plane
		return p1 + (p2 - p1) * (-distance(p1) / normal.dot(p2 - p1));
	};
	void invert() { normal *= -1.f; }
	float distance(Float3 point) const { return normal.dot(point) - D; }
};
static_assert(sizeof(Facet) % 16 == 0);


struct Facet2 {
	std::array<int,3> verticesIds;
	Float3 normal;
};

class ConvexHull {
	std::vector<Float3> vertices;	// [nm] Only used to determine CoM
	//std::vector<Facet> facets;

	std::vector<Facet2> facets;

	void Add(const Float3& vertex) {
		if (std::count(vertices.begin(), vertices.end(), vertex) == 0) {
			vertices.emplace_back(vertex);
		}
	}

	int getIndex(const Float3& vertex) const {
		for (int i = 0; i < vertices.size(); i++) {
			if (vertices[i] == vertex) {
				return i;
			}
		}
		throw std::runtime_error("Vertex not found");
	}

	void CullSmallTriangles();

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

	//const std::vector<Facet>& GetFacets() const { return facets; }
	const std::vector<Facet> GetFacets() const { 
		std::vector<Facet> facets2;
		for (const auto& f : facets) {
			Facet newFacet;
			for (int i = 0; i < 3; i++) {
				newFacet.vertices[i] = vertices[f.verticesIds[i]];
			}
			newFacet.normal = f.normal;
			newFacet.D = newFacet.normal.dot(newFacet.vertices[0]);
			facets2.emplace_back(newFacet);
		}
		return facets2;
	}


	const std::vector<Float3>& GetVertices() const { return vertices; }

	void Add(const Facet& plane) {

		Facet2 newFacet;
		for (int i = 0; i < 3; i++) {
			Add(plane.vertices[i]);
			newFacet.verticesIds[i] = getIndex(plane.vertices[i]);
		}
		newFacet.normal = plane.normal;

		facets.emplace_back(newFacet);
		//facets.emplace_back(plane);
	/*	for (const auto& v : plane.vertices) {
			Add(v);
		}*/
	}
};

class MoleculeHullFactory {
	std::vector<Float3> particlePositions; // [nm]
	std::vector<char> atomLetters;

public:

	ConvexHull convexHull;

	void CreateConvexHull();
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





/// <summary>
/// Implementation of SutherlandHodgman algorithm to find the intersection of two convex hulls
/// https://github.com/Alamot/code-snippets/blob/master/graphics/SutherlandHodgman/Linux/SutherlandHodgman.cpp
/// </summary>
/// <param name="ch1"></param>
/// <param name="ch2"></param>
/// <returns></returns>
std::vector<Float3> FindIntersectionConvexhullFrom2Convexhulls(const ConvexHull& ch1, const ConvexHull& ch2,
	std::function<bool(const std::vector<Facet>&, const std::vector<Float3>&, bool)> callback
);












