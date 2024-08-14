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

    facet.normal = -Float3{ f.outerplane().coordinates()[0], f.outerplane().coordinates()[1], f.outerplane().coordinates()[2] };
    
    facet.D = facet.normal.dot(facet.vertices[0]);

	return facet;
}


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
    convexHull = ConvexHull(particlePositions);



    ConvexHull closelyFittingHull{ particlePositions };

    std::vector<Float3> hullVerticesPadded;
    for (const Float3& v : closelyFittingHull.GetVertices()) {
        const std::vector<Float3> spherePoints = GenerateSpherePoints(v, .1f, 16);

        hullVerticesPadded.insert(hullVerticesPadded.end(), spherePoints.begin(), spherePoints.end());
    }


    float min = 99999.f;
    for(int i = 0; i < hullVerticesPadded.size(); i++) {
		for (int j = i+1; j < hullVerticesPadded.size(); j++) {
			float dist = (hullVerticesPadded[i] - hullVerticesPadded[j]).len();
			if (dist < min) {
				min = dist;
			}
		}
	}


    convexHull = ConvexHull(hullVerticesPadded);
    
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





ConvexHull FindIntersectionConvexhullFrom2Convexhulls(const ConvexHull& ch1, const ConvexHull& ch2) {
    float D1, D2 = 0;
    ConvexHull polygon = ConvexHull(ch1);
    const float EPSILON = 0.00001f;

    const auto& clippingPlanes = ch2.GetFacets();

    for (unsigned int c = 0; c < clippingPlanes.size(); c++) 
    {
        ConvexHull clippedPolygon{};

        std::vector<Float3> points = polygon.GetVertices();
        for (unsigned int i = 0; i < points.size() - 1; i++) {
            D1 = clippingPlanes[c].distance(points[i]);
            D2 = clippingPlanes[c].distance(points[i + 1]);
            if ((D1 <= 0) && (D2 <= 0))
            {
                //clippedPolygon.add(points[i + 1]);
                clippedPolygon.Add(points[i + 1]);
            }
            else if ((D1 > 0) && ((D2 > -EPSILON) && (D2 < EPSILON)))
            {
                //clippedPolygon.add(points[i + 1]);
                clippedPolygon.Add(points[i + 1]);

            }
            else if (((D1 > -EPSILON) && (D1 < EPSILON)) && (D2 > 0))
            {
                continue;
            }
            else if ((D1 <= 0) && (D2 > 0))
            {
                clippedPolygon.Add(clippingPlanes[c].intersectionPoint(points[i], points[i + 1]));
            }
            else if ((D1 > 0) && (D2 <= 0))
            {                
                clippedPolygon.Add(clippingPlanes[c].intersectionPoint(points[i], points[i + 1]));
                clippedPolygon.Add(points[i + 1]);
            }
        }
        polygon = clippedPolygon; // keep on working with the new polygon
    }
    return polygon;
}