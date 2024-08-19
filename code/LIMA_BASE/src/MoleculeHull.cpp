#include "MoleculeHull.cuh"
#include "RenderUtilities.cuh"

#include <quickhull/quickhull.hpp>
#include <numeric>

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>

Facet ConvertFacet(const orgQhull::QhullFacet& f) {
	Facet facet;
	facet.vertices[0].x = f.vertices()[0].point().coordinates()[0];
	facet.vertices[0].y = f.vertices()[0].point().coordinates()[1];
	facet.vertices[0].z = f.vertices()[0].point().coordinates()[2];

	facet.vertices[1].x = f.vertices()[1].point().coordinates()[0];
	facet.vertices[1].y = f.vertices()[1].point().coordinates()[1];
	facet.vertices[1].z = f.vertices()[1].point().coordinates()[2];

	facet.vertices[2].x = f.vertices()[2].point().coordinates()[0];
	facet.vertices[2].y = f.vertices()[2].point().coordinates()[1];
	facet.vertices[2].z = f.vertices()[2].point().coordinates()[2];

    facet.normal = Float3{ f.outerplane().coordinates()[0], f.outerplane().coordinates()[1], f.outerplane().coordinates()[2] };
    

    const Float3 expectedNormal = (facet.vertices[1] - facet.vertices[0]).cross(facet.vertices[2] - facet.vertices[0]).norm();
    if (expectedNormal.dot(facet.normal) < 0.f) {
        std::swap(facet.vertices[1], facet.vertices[2]);
    }


    facet.D = facet.normal.dot(facet.vertices[0]);

	return facet;
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

#include "TimeIt.h"
void MoleculeHullFactory::CreateConvexHull() {
    TimeIt timer{"CH", true};

 //   convexHull = ConvexHull(particlePositions);

 //   ConvexHull closelyFittingHull{ particlePositions };

 //   std::vector<Float3> hullVerticesPadded;
 //   for (const Float3& v : closelyFittingHull.GetVertices()) {
 //       const std::vector<Float3> spherePoints = GenerateSpherePoints(v, .1f, 16);

 //       hullVerticesPadded.insert(hullVerticesPadded.end(), spherePoints.begin(), spherePoints.end());
 //   }


 //   float min = 99999.f;
 //   for(int i = 0; i < hullVerticesPadded.size(); i++) {
	//	for (int j = i+1; j < hullVerticesPadded.size(); j++) {
	//		float dist = (hullVerticesPadded[i] - hullVerticesPadded[j]).len();
	//		if (dist < min) {
	//			min = dist;
	//		}
	//	}
	//}
 //   convexHull = ConvexHull(hullVerticesPadded);
 //   

	//convexHull.CullSmallTriangles();

    // To ensure the atoms are covered including their VDW dist, we pad them first
    std::vector<Float3> paddedPoints;
    for (const Float3& p : particlePositions) {
        const std::vector<Float3> spherePoints = GenerateSpherePoints(p, .1f, 16);
        paddedPoints.insert(paddedPoints.end(), spherePoints.begin(), spherePoints.end());
    }
    ConvexHull paddedCH{ paddedPoints };
    paddedCH.CullSmallTriangles();


    std::vector<Float3> uniquePrunedHullVertices;
    for (const auto& facet : paddedCH.GetFacets()) {
        for (int i = 0; i < 3; i++) {

            if (std::count(uniquePrunedHullVertices.begin(), uniquePrunedHullVertices.end(), facet.vertices[i]) == 0)
                uniquePrunedHullVertices.emplace_back(facet.vertices[i]);
		}
	}

    convexHull = ConvexHull(std::vector<Float3>(uniquePrunedHullVertices.begin(), uniquePrunedHullVertices.end()));

    int a = 0;
}



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
        moleculehullsHost[i] = MoleculeHull{ hullFacetcountPrefixsum[i], (int)molecules[i].convexHull.GetFacets().size(), hullParticlecountPrefixsum[i], (int)molecules[i].GetParticles().size() };
	}

    nFacets = hullFacetcountPrefixsum.back() + molecules.back().convexHull.GetFacets().size();
    nParticles = hullParticlecountPrefixsum.back() + molecules.back().GetParticles().size();
    nMoleculeHulls = molecules.size();

    std::vector<Facet> facetsHost(nFacets);
    std::vector<RenderAtom> particlesHost(nParticles);

    for (int i = 0; i < molecules.size(); i++) {
        for (int j = 0; j < molecules[i].convexHull.GetFacets().size(); j++) {
            facetsHost[hullFacetcountPrefixsum[i] + j] = molecules[i].convexHull.GetFacets()[j];
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


