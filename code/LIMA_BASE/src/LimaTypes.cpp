#include "LimaTypes.cuh"
#include "Statistics.h"
#include "quickhull/quickhull.hpp"

template <int maxFacets>
ConvexHull<maxFacets>::ConvexHull(const std::vector<Float3>& points) {
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
        Plane facet;
        facet.vertices[0] = *f.vertices_[0];
        facet.vertices[1] = *f.vertices_[1];
        facet.vertices[2] = *f.vertices_[2];

        facet.distFromOrigo = f.D;

        facet.normal = Float3{ f.normal_[0], f.normal_[1], f.normal_[2] };

        Add(facet);
    }

    if (numVertices != 2+numFacets / 2) {
		throw std::runtime_error("ConvexHull algo failed - resulted structure is not convex (generally due to precision errors)");
	}
}
template class ConvexHull<64>;