#include "MoleculeHull.cuh"
#include "RenderUtilities.cuh"
#include "MDFiles.h"
#include "TimeIt.h"

#include <quickhull/quickhull.hpp>
#include <numeric>

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>

std::array<Float3,3> ConvertFacet(const orgQhull::QhullFacet& f) {
    std::array<Float3,3> vertices;

    vertices[0] = Float3{ f.vertices()[0].point().coordinates()[0], f.vertices()[0].point().coordinates()[1], f.vertices()[0].point().coordinates()[2] };
    vertices[1] = Float3{ f.vertices()[1].point().coordinates()[0], f.vertices()[1].point().coordinates()[1], f.vertices()[1].point().coordinates()[2] };
    vertices[2] = Float3{ f.vertices()[2].point().coordinates()[0], f.vertices()[2].point().coordinates()[1], f.vertices()[2].point().coordinates()[2] };

    const Float3 normal{ f.outerplane().coordinates()[0], f.outerplane().coordinates()[1], f.outerplane().coordinates()[2] };

    const Float3 expectedNormal = (vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]).norm();
    if (expectedNormal.dot(normal) < 0.f) {
        std::swap(vertices[1], vertices[2]);
    }

    return vertices;
}



class VertexMapping {
    std::vector<int> map;

public:
    VertexMapping(int nVertices) : map(nVertices) {
        for (int i = 0; i < nVertices; i++) {
            map[i] = i;
        }
    }

    void Add(int from, int to) {
        map[from] = to;
    }

    int Get(int from) {
        int prev = from;
        while (true) {
            if (map[prev] == prev) {
                return prev;
            }
            prev = map[prev];
        }
    }  

    bool MapsToSelf(int from) {
        return map[from] == from;
    }
};


void ConvexHull::CullSmallTriangles() {


    for (auto& facet : facets) {
        Float3 expecetedNormal = (vertices[facet.verticesIds[1]] - vertices[facet.verticesIds[0]]).cross(vertices[facet.verticesIds[2]] - vertices[facet.verticesIds[0]]).norm();
        if (expecetedNormal.dot(facet.normal) < .9f) {
			throw std::runtime_error("A facet's normal changed too much during optimization");
		}
    }

    std::vector<bool> facetMarkedForDeletion(facets.size(), false);
    VertexMapping vertexMapping(vertices.size());


    for (int facetId = 0; facetId < facets.size(); facetId++) {
        Facet2& facet = facets[facetId];

        // If > 0 of the facets vertices are moved, we update the facet, but do no further work
        //if (vertexMarkedForDeletion[facet.verticesIds[0]] || vertexMarkedForDeletion[facet.verticesIds[1]] || vertexMarkedForDeletion[facet.verticesIds[2]]) 
        if (!vertexMapping.MapsToSelf(facet.verticesIds[0]) || !vertexMapping.MapsToSelf(facet.verticesIds[1]) || !vertexMapping.MapsToSelf(facet.verticesIds[2]))            
        {
            facet.verticesIds[0] = vertexMapping.Get(facet.verticesIds[0]);
            facet.verticesIds[1] = vertexMapping.Get(facet.verticesIds[1]);
            facet.verticesIds[2] = vertexMapping.Get(facet.verticesIds[2]);
            continue;
        }


	    const float maxEdgeLen = std::max(std::max(
            (vertices[facet.verticesIds[0]] - vertices[facet.verticesIds[1]]).len(), 
            (vertices[facet.verticesIds[1]] - vertices[facet.verticesIds[2]]).len()), 
            (vertices[facet.verticesIds[2]] - vertices[facet.verticesIds[0]]).len());

        // Face is tiny. Delete V1 and and V2. V0 is moved to the center of the facet. V1 and V2 points to V0
        if (maxEdgeLen < 0.2f) {
            facetMarkedForDeletion[facetId] = true;

            vertices[facet.verticesIds[0]] = (vertices[facet.verticesIds[0]] + vertices[facet.verticesIds[1]] + vertices[facet.verticesIds[2]]) / 3.f;
            
            vertexMapping.Add(facet.verticesIds[1], facet.verticesIds[0]);
            vertexMapping.Add(facet.verticesIds[2], facet.verticesIds[0]);

            continue;
		}

        for (int i = 0; i < 3; i++) {
            int v1 = facet.verticesIds[i];
            int v2 = facet.verticesIds[(i + 1) % 3];

            if ((vertices[v1] - vertices[v2]).len() < 0.1f) {
                facetMarkedForDeletion[facetId] = true;
                vertices[v1] = (vertices[v1] + vertices[v2]) / 2.f;
                vertexMapping.Add(v2, v1);
                break;
            }
        }


	}


    // Finally go through all facets again and update the vertices
    for (auto& facet : facets) {
        facet.verticesIds[0] = vertexMapping.Get(facet.verticesIds[0]);
        facet.verticesIds[1] = vertexMapping.Get(facet.verticesIds[1]);
        facet.verticesIds[2] = vertexMapping.Get(facet.verticesIds[2]);

        Float3 newNormal = (vertices[facet.verticesIds[1]] - vertices[facet.verticesIds[0]]).cross(vertices[facet.verticesIds[2]] - vertices[facet.verticesIds[0]]).norm();
        if (newNormal.dot(facet.normal) < 0.f) {
            //throw std::runtime_error("A facet's normal changed too much during optimization");
        }
        facet.normal = newNormal;
    }






    std::vector<Facet2> compressedFacets;
    for (int i = 0; i < facets.size(); i++) {
        if (!facetMarkedForDeletion[i]) {
            compressedFacets.push_back(facets[i]);          
        }
    }


    facets = compressedFacets;
}




// TODO: If we receive a set of points from the same molecule as we have seen here before, just with a new rotation and translation, 
// we should be able to reuse the old convex hull, simply by applying the old rotation+translation as used before

ConvexHull::ConvexHull(const std::vector<Float3>& points) {
    constexpr size_t dim = 3;
    constexpr float eps = std::numeric_limits<float>::epsilon();




    std::vector<double3> pointsDouble(points.size());
    for (int i = 0; i < points.size(); i++) {
		pointsDouble[i] = { points[i].x, points[i].y, points[i].z };
	}

    orgQhull::Qhull qhull{};
    qhull.runQhull("", 3, points.size(), (double*)pointsDouble.data(), "Qt");


    for (const auto& f : qhull.facetList()) {
        Add(ConvertFacet(f));
    }

    if (vertices.size() != 2 + facets.size() / 2) {
        throw std::runtime_error("ConvexHull algo failed - resulted structure is not convex (generally due to precision errors)");
    }

    //for (int i = 0; i < 2; i++)
}




// TODO: Move to a GeometryUtils file
std::vector<Float3> GenerateSpherePoints(const Float3& origo, float radius, int numPoints) {
    std::vector<Float3> points;
    points.reserve(numPoints);

    const float goldenRatio = (1.0f + std::sqrt(5.0f)) / 2.0f;
    const float angleIncrement = 2.0f * PI * goldenRatio;

    for (int i = 0; i < numPoints; ++i) {
        float t = static_cast<float>(i) / static_cast<float>(numPoints - 1);
        float inclination = std::acos(1.0f - 2.0f * t);
        float azimuth = angleIncrement * i;

        float x = radius * std::sin(inclination) * std::cos(azimuth);
        float y = radius * std::sin(inclination) * std::sin(azimuth);
        float z = radius * std::cos(inclination);

        points.push_back({ origo.x + x, origo.y + y, origo.z + z });
    }

    return points;
}

void MoleculeHullFactory::CreateConvexHull() {
    //TimeIt timer{"CH", true};


    // To ensure the atoms are covered including their VDW dist, we pad them first
    std::vector<Float3> paddedPoints;
    for (const Float3& p : particlePositions) {
        const std::vector<Float3> spherePoints = GenerateSpherePoints(p, .1f, 16);
        paddedPoints.insert(paddedPoints.end(), spherePoints.begin(), spherePoints.end());
    }
    ConvexHull paddedCH{ paddedPoints };
    paddedCH.CullSmallTriangles();


    std::vector<Float3> uniquePrunedHullVertices;
    for (const auto& facet : paddedCH.GetSelfcontainedFacets()) {
        for (int i = 0; i < 3; i++) {
            if (std::count(uniquePrunedHullVertices.begin(), uniquePrunedHullVertices.end(), facet.vertices[i]) == 0)
                uniquePrunedHullVertices.emplace_back(facet.vertices[i]);
		}
	}

    convexHull = ConvexHull(std::vector<Float3>(uniquePrunedHullVertices.begin(), uniquePrunedHullVertices.end()));
}



void Facet::ApplyTransformation(const glm::mat4& transformationMatrix) {
	for (int i = 0; i < 3; i++) {
		const glm::vec4 transformedVertex = transformationMatrix * glm::vec4{ vertices[i].x, vertices[i].y, vertices[i].z, 1.f };
		vertices[i] = Float3{ transformedVertex.x, transformedVertex.y, transformedVertex.z };
	}
	normal = (vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]).norm();
    D = normal.dot(vertices[0]);
}

void ConvexHull::ApplyTransformation(const glm::mat4& transformationMatrix) {
    for (Float3& v : vertices) {
        const glm::vec4 transformedVertex = transformationMatrix * glm::vec4{ v.x, v.y, v.z, 1.f };
		v = Float3{ transformedVertex.x, transformedVertex.y, transformedVertex.z };
	}
    
	for (Facet2& f : facets) {
		f.normal = (vertices[f.verticesIds[1]] - vertices[f.verticesIds[0]]).cross(vertices[f.verticesIds[2]] - vertices[f.verticesIds[0]]).norm();
	}
}
    

void MoleculeHullFactory::ApplyTransformation(const glm::mat4& transformationMatrix) {
	for (Float3& v : particlePositions) {
		const glm::vec4 transformedVertex = transformationMatrix * glm::vec4{ v.x, v.y, v.z, 1.f };
		v = Float3{ transformedVertex.x, transformedVertex.y, transformedVertex.z };
	}
	convexHull.ApplyTransformation(transformationMatrix);
}


void MoleculeHull::ApplyTransformation(const glm::mat4& transformationMatrix, Facet* const facetsInCollection, RenderAtom* const particlesInCollection, Float3 boxSize) {
    for (int facetId = indexOfFirstFacetInBuffer; facetId < indexOfFirstFacetInBuffer + nFacets; facetId++) {
		facetsInCollection[facetId].ApplyTransformation(transformationMatrix);
	}

    const glm::vec4 transformedCenter = transformationMatrix * glm::vec4{ center.x, center.y, center.z, 1.f };
    center = Float3{ transformedCenter.x, transformedCenter.y, transformedCenter.z };

    for (int particleId = indexOfFirstParticleInBuffer; particleId < indexOfFirstParticleInBuffer + nParticles; particleId++) {
        const Float3 screenNormalizedPosition = Float3{ particlesInCollection[particleId].position.x, particlesInCollection[particleId].position.y, particlesInCollection[particleId].position.z };
        const Float3 unNormalizedPosition = (screenNormalizedPosition + Float3{ 0.5f }) * boxSize;

        const glm::vec4 transformedVertex = transformationMatrix * glm::vec4{ unNormalizedPosition.x, unNormalizedPosition.y, unNormalizedPosition.z, 1.f };

        Float3 normalizedTransformedPosition = Float3{ transformedVertex.x, transformedVertex.y, transformedVertex.z} / boxSize - Float3{0.5f};

        particlesInCollection[particleId].position = normalizedTransformedPosition.Tofloat4(particlesInCollection[particleId].position.w);
    }

}






#include "Statistics.h"
#include "Utilities.h"

MoleculeHullCollection::MoleculeHullCollection(const std::vector<MoleculeHullFactory>& molecules, Float3 boxSize) {
    
    std::vector<int> hullFacetcountPrefixsum(molecules.size(), 0);
    std::vector<int> hullParticlecountPrefixsum(molecules.size(), 0);

    for (int i = 0; i < molecules.size(); i++) {
        hullFacetcountPrefixsum[i] = molecules[i].convexHull.GetFacets().size();
        hullParticlecountPrefixsum[i] = molecules[i].GetParticles().size();
	}
    std::exclusive_scan(hullFacetcountPrefixsum.begin(), hullFacetcountPrefixsum.end(), hullFacetcountPrefixsum.begin(), 0);
    std::exclusive_scan(hullParticlecountPrefixsum.begin(), hullParticlecountPrefixsum.end(), hullParticlecountPrefixsum.begin(), 0);


    std::vector<MoleculeHull> moleculehullsHost(molecules.size());
    for (int i = 0; i < molecules.size(); i++) {
        const Float3 center = Statistics::CalculateMinimaxPoint(molecules[i].convexHull.GetVertices());
        const float radius = LAL::LargestDiff(center, molecules[i].convexHull.GetVertices());
        moleculehullsHost[i] = MoleculeHull{ // Todo, just make a contructor taking the factory as input
            hullFacetcountPrefixsum[i], 
            static_cast<int>(molecules[i].convexHull.GetFacets().size()), 
            hullParticlecountPrefixsum[i], 
            static_cast<int>(molecules[i].GetParticles().size()),
            center, 
            radius
        };
	}

    nFacets = hullFacetcountPrefixsum.back() + molecules.back().convexHull.GetFacets().size();
    nParticles = hullParticlecountPrefixsum.back() + molecules.back().GetParticles().size();
    nMoleculeHulls = molecules.size();

    std::vector<Facet> facetsHost(nFacets);
    std::vector<RenderAtom> particlesHost(nParticles);

    for (int i = 0; i < molecules.size(); i++) {
        for (int j = 0; j < molecules[i].convexHull.GetFacets().size(); j++) {
            facetsHost[hullFacetcountPrefixsum[i] + j] = molecules[i].convexHull.GetSelfcontainedFacets()[j];
        }
        for (int j = 0; j < molecules[i].GetParticles().size(); j++) {
            particlesHost[hullParticlecountPrefixsum[i] + j] = RenderAtom{ molecules[i].GetParticles()[j] , boxSize, molecules[i].GetAtomLetters()[j] };
		}
    }

    cudaMalloc(&facets, sizeof(Facet) * facetsHost.size());
    cudaMemcpy(facets, facetsHost.data(), sizeof(Facet) * facetsHost.size(), cudaMemcpyHostToDevice);
    cudaMalloc(&particles, sizeof(RenderAtom) * particlesHost.size());
    cudaMemcpy(particles, particlesHost.data(), sizeof(RenderAtom) * particlesHost.size(), cudaMemcpyHostToDevice);
    cudaMalloc(&moleculeHulls, sizeof(MoleculeHull) * moleculehullsHost.size());
    cudaMemcpy(moleculeHulls, moleculehullsHost.data(), sizeof(MoleculeHull) * moleculehullsHost.size(), cudaMemcpyHostToDevice);
}

MoleculeHullCollection::~MoleculeHullCollection() {
    cudaFree(facets);
    cudaFree(particles);
    cudaFree(moleculeHulls);
}

void MoleculeHullCollection::WritePositionsToGrofile(GroFile& grofile) {
    if (nParticles != grofile.atoms.size())
        throw ("GroFile does not seem to belong to this MoleculeHullCollection");

    std::vector<RenderAtom> renderatoms(nParticles);
    cudaMemcpy(renderatoms.data(), particles, sizeof(RenderAtom) * nParticles, cudaMemcpyDeviceToHost);

    for (int i = 0; i < nParticles; i++) {
        grofile.atoms[i].position = Float3{ renderatoms[i].position.x, renderatoms[i].position.y, renderatoms[i].position.z };
	}
}
