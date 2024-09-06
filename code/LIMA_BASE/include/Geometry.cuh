#include "LimaTypes.cuh"





namespace Geometry {
	struct Plane {
		Float3 normal;
		float signedDistance;

		Plane(Float3 normal, float signedDistance) : normal(normal), signedDistance(signedDistance) {};

	};

	// Function to get the intersection of the plane with an edge between two points
	bool IntersectPlaneWithEdge(const Plane& plane, const Float3& pointA, const Float3& pointB, Float3& intersection);



	float CalcAreaOfPlaneInBoundingBox(const Plane& plane, const Float3& boxSize);
}

