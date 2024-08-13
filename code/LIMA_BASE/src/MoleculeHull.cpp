#include "MoleculeHull.cuh"
#include "RenderUtilities.cuh"

#include <quickhull/quickhull.hpp>
#include <numeric>

ConvexHull::ConvexHull(const std::vector<Float3>& points) {
    constexpr size_t dim = 3;
    constexpr float eps = std::numeric_limits<float>::epsilon();

    quick_hull<std::vector<Float3>::const_iterator> qh{ dim, eps };

    qh.add_points(std::cbegin(points), std::cend(points));


    auto initial_simplex = qh.get_affine_basis();
    if (initial_simplex.size() < dim + 1) {
        throw std::runtime_error("Degenerated input set in ConvexHull algo");
    }

    qh.create_initial_simplex(std::cbegin(initial_simplex), std::prev(std::cend(initial_simplex)));
    qh.create_convex_hull();
    //if (!qh.check()) {
    //    throw std::runtime_error("ConvexHull algo failed - resulted structure is not convex (generally due to precision errors)");
    //}


    for (const auto& f : qh.facets_) {
        Facet facet;
        facet.vertices[0] = *f.vertices_[0];
        facet.vertices[1] = *f.vertices_[1];
        facet.vertices[2] = *f.vertices_[2];

        facet.distFromOrigo = f.D;

        facet.normal = Float3{ f.normal_[0], f.normal_[1], f.normal_[2] };

        Add(facet);
    }

    if (vertices.size() != 2 + facets.size() / 2) {
        throw std::runtime_error("ConvexHull algo failed - resulted structure is not convex (generally due to precision errors)");
    }
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