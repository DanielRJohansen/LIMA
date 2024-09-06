#include "Geometry.cuh"
#include <optional>

namespace Geometry {
	// Function to get the intersection of the plane with an edge between two points
	std::optional<Float3> IntersectPlaneWithEdge(const Plane& plane, const Float3& pointA, const Float3& pointB) {
		Float3 AB = pointB - pointA;
		float denom = plane.normal.dot(AB);
		if (std::abs(denom) > 1e-6f) {
			float t = -(plane.normal.dot(pointA) + plane.signedDistance) / denom;
			if (t >= 0 && t <= 1) {
				Float3 intersection = pointA + AB * t;
				return intersection;
			}
		}
		return std::nullopt;
	}


	float PolygonArea(const std::vector<Float3>& vertices) {
		if (vertices.size() < 3) return 0.0f;

		Float3 refPoint = vertices[0];
		float area = 0.0f;

		for (size_t i = 1; i < vertices.size() - 1; ++i) {
			Float3 v0 = vertices[i] - refPoint;
			Float3 v1 = vertices[i + 1] - refPoint;
			area += v0.cross(v1).len() * 0.5f;
		}

		return area;
	}

	float CalcAreaOfPlaneInBoundingBox(const Plane& plane, const Float3& boxSize) {
		std::vector<Float3> intersections;

		Float3 vertices[8] = {
			Float3{0},
			{0.f, 0.f, boxSize.z},
			{0.f, boxSize.y, 0.f},
			{0.f, boxSize.y, boxSize.z},
			{boxSize.x, 0.f, 0.f},
			{boxSize.x, 0.f, boxSize.z},
			{boxSize.x, boxSize.y, 0.f},
			boxSize
		};

		int edges[12][2] = {
			{0, 1}, {0, 2}, {0, 4}, {1, 3}, {1, 5}, {2, 3},
			{2, 6}, {3, 7}, {4, 5}, {4, 6}, {5, 7}, {6, 7}
		};

		for (int i = 0; i < 12; i++) {
			auto intersection = IntersectPlaneWithEdge(plane, vertices[edges[i][0]], vertices[edges[i][1]]);
			if (intersection) {
				intersections.push_back(intersection.value());
			}
		}


		float area = 0.0f;

		if (intersections.size() > 2) {
			Float3 refPoint = intersections[0];

			for (size_t i = 1; i < intersections.size() - 1; ++i) {
				Float3 v0 = intersections[i] - refPoint;
				Float3 v1 = intersections[i + 1] - refPoint;
				area += v0.cross(v1).len() * 0.5f;
			}
		}

		return area;	
	}
}
