#include "ConvexHullEngine.cuh"
#include "Statistics.h"
#include "Utilities.h"

#include <glm.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <gtx/rotate_vector.hpp>
#undef GLM_ENABLE_EXPERIMENTAL

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/tuple.h>

#include "DeviceAlgorithms.cuh"
#include "Utilities.h"

#include "TimeIt.h"
#include<span>

//#define WITHOUT_NUMPY
//#include "matplotlib-cpp/matplotlibcpp.h"
//#undef WITHOUT_NUMPY

const int maxFacetsInCH = 256 + 128;
const int maxCollisionsPerMH = 32;

struct Overlap {
	bool isOverlapping = false;
	Float3 intersectionCenter{};
	float depth{};
};


struct TriList {
	Float3 triVertices[maxFacetsInCH * 3];
	int nTris = 0;

	TriList() {}
	TriList(const TriList& other) {
		for (int i = 0; i < other.nTris; i++) {
			triVertices[i * 3 + 0] = other.triVertices[i * 3 + 0];
			triVertices[i * 3 + 1] = other.triVertices[i * 3 + 1];
			triVertices[i * 3 + 2] = other.triVertices[i * 3 + 2];
		}
		nTris = other.nTris;
	}	
};

__device__ void FindMaxDepth(const Float3& queryPoint, const Float3* const points, int nPoints, float& result) {
	float localMaxDepth = 0.f;
	for (int vertexIndex = threadIdx.x; vertexIndex < nPoints; vertexIndex += blockDim.x) {
		const auto& point = points[vertexIndex];
		const float dist = (queryPoint - point).len();
		if (dist > localMaxDepth) {
			localMaxDepth = dist;
		}
	}
	for (int i = 0; i < blockDim.x; i++) {
		if (i == threadIdx.x) {
			if (localMaxDepth > result)
				result = localMaxDepth;
		}
		__syncthreads();
	}
}


__device__ void FindMinimaxPoint(Float3& sharedDownGradient, float& sharedMaxDist, float& sharedSecondMaxDist, const Float3* const points, const int nPoints, Float3& sharedMinimaxPoint, float& sharedPrevVelocity) {
	if (threadIdx.x == 0) {
		sharedMinimaxPoint = points[0];
		sharedPrevVelocity = FLT_MAX;
	}

	const int maxIterations = 5;
	for (int i = 0; i < maxIterations; i++) {


		if (threadIdx.x == 0) {
			sharedDownGradient = Float3{};
			sharedMaxDist = 0.0f;
			sharedSecondMaxDist = 0.0f;
		}
		__syncthreads();

		Float3 localDownGradient{};
		float localMaxDist = 0.0f;
		float localSecondMaxDist = 0.0f;

		// Each thread processes a subset of the triangleVertices
		for (int vertexIndex = threadIdx.x; vertexIndex < nPoints; vertexIndex += blockDim.x) {
			const auto& point = points[vertexIndex];
			const float dist = (sharedMinimaxPoint - point).len();

			if (dist > localMaxDist) {
				localSecondMaxDist = localMaxDist;
				localMaxDist = dist;
				localDownGradient = (point - sharedMinimaxPoint);
			}
			else if (dist > localSecondMaxDist) {
				localSecondMaxDist = dist;
			}
		}


		// Reduction across threads within a block to find the global max and second max distances
		for (int i = 0; i < blockDim.x; i++) {
			if (i == threadIdx.x) {
				if (localMaxDist > sharedMaxDist) {
					sharedSecondMaxDist = sharedMaxDist;
					sharedMaxDist = localMaxDist;
					sharedDownGradient = localDownGradient;

					if (localSecondMaxDist > sharedSecondMaxDist)
						sharedSecondMaxDist = localSecondMaxDist;
				}
				else if (localMaxDist > sharedSecondMaxDist) {
					sharedSecondMaxDist = localMaxDist;
				}
			}
			__syncthreads();
		}


		__syncthreads();
		if (localMaxDist == sharedMaxDist) {
			sharedDownGradient = localDownGradient;
		}
		__syncthreads();

		if (threadIdx.x == 0) {
			// Calculate velocity
			float velocity = std::min((sharedMaxDist) / 2.0f, sharedPrevVelocity * 0.9f);
			sharedPrevVelocity = velocity;

			// Normalize the gradient
			sharedDownGradient = sharedDownGradient;
			sharedMinimaxPoint += sharedDownGradient.norm() * velocity;
		}
		__syncthreads();
	}
}


// Input is in global memory
__global__ void CalculateIntersectionCenterAndDepth(const Float3* const triVerticesBuffer, const int* const nVerticesBuffer, 
	MoleculeHull* hullsBuffer, Float3* pivotPoints, Float3* forceApplicationPoints, Float3* forceDirections, float* intersectionRadii)
{
	__shared__ float prevVelocity;
	__shared__ int nVertices;

	__shared__ Float3 sharedDownGradient; // Shared memory to store gradient direction
	__shared__ float sharedMaxDist;        // Shared memory to store maximum distance
	__shared__ float sharedSecondMaxDist;  // Shared memory to store second maximum distance

	__shared__ Overlap overlap;


	// First figure out if there is any work to do
	if (threadIdx.x == 0) {
		nVertices = nVerticesBuffer[blockIdx.x];
	}
	__syncthreads();
	if (nVertices < 4)
		return;


	const float epsilon = 1e-4;
	const Float3* const triVertices = &(triVerticesBuffer[blockIdx.x * maxFacetsInCH * 3]);

	if (threadIdx.x == 0) {
		overlap.intersectionCenter = triVertices[0];
		overlap.depth = 0.0f;
		prevVelocity = FLT_MAX;
	}
	__syncthreads();
	
	FindMinimaxPoint(sharedDownGradient, sharedMaxDist, sharedSecondMaxDist, triVertices, nVertices, overlap.intersectionCenter, prevVelocity);

	FindMaxDepth(overlap.intersectionCenter, triVertices, nVertices, overlap.depth);

	if (threadIdx.x == 0) {

		const int collisionId = hullsBuffer[blockIdx.x].nCollisions-1; // -1 Because the previous kernel did the "++" 
		if (collisionId < 0 || collisionId >= maxCollisionsPerMH) {
			printf("Invalid collision id %d nVert %d\n", collisionId, nVertices);
		}
		forceApplicationPoints[blockIdx.x * maxCollisionsPerMH + collisionId] = overlap.intersectionCenter;
		intersectionRadii[blockIdx.x * maxCollisionsPerMH + collisionId] = overlap.depth;
	}
}



const int nThreadsInFindIntersectionKernel = 32;
__global__ void FindIntersectionConvexhullFrom2Convexhulls(MoleculeHull* hullsBuffer, const Facet* const facets, Float3 boxSize, int targetHull, 
	Float3* verticesOutBuffer, int* nVerticesOutBuffer, Float3* pivotPoints
	//,bool* collisionsBuffer
) {
	__shared__ TriList clippedFacets;
	__shared__ Facet currentClippingFacet;
	__shared__ int prefixSumBuffer[nThreadsInFindIntersectionKernel];

	if (threadIdx.x == 0)
		nVerticesOutBuffer[blockIdx.x] = 0;

	if (targetHull == blockIdx.x) {
		if (threadIdx.x == 0) {
			//collisionsBuffer[blockIdx.x] = false;
		}
		return;
	}
		

	MoleculeHull mh1 = hullsBuffer[blockIdx.x];
	MoleculeHull mh2 = hullsBuffer[targetHull];

	if ((mh1.center - mh2.center).len() > mh1.radius + mh2.radius)
		return;

	// First initialize the facets which vertices we are clipping
	for (int facetIdOffset = threadIdx.x; facetIdOffset < mh2.nFacets; facetIdOffset += blockDim.x) {
		const int facetId = mh2.indexOfFirstFacetInBuffer + facetIdOffset;

		if (facetIdOffset >= maxFacetsInCH)
			printf("Too many facets: %d\n", facetIdOffset);

		for (int vertexId = 0; vertexId < 3; vertexId++) 
			clippedFacets.triVertices[facetIdOffset * 3 + vertexId] = facets[facetId].vertices[vertexId];
	}
	if (threadIdx.x == 0) {
		clippedFacets.nTris = mh2.nFacets;
	}
	prefixSumBuffer[threadIdx.x] = 0;
	__syncthreads();



	for (int clippingFacetId = mh1.indexOfFirstFacetInBuffer; clippingFacetId < mh1.indexOfFirstFacetInBuffer + mh1.nFacets; clippingFacetId++) {
		if (threadIdx.x == 0) {
			currentClippingFacet = facets[clippingFacetId];
		}
		__syncthreads();

		Float3 clippedFacetsToBeAdded[8][3];
		int nClippedFacetsToBeAdded = 0;

		for (int facetId = threadIdx.x; facetId < clippedFacets.nTris; facetId += blockDim.x) {
			const Float3 queryFacet[3] = { clippedFacets.triVertices[facetId * 3 + 0], clippedFacets.triVertices[facetId * 3 + 1], clippedFacets.triVertices[facetId * 3 + 2] };

			bool vertexIsInsideFacet[] = {
				(queryFacet[0] - currentClippingFacet.vertices[0]).dot(currentClippingFacet.normal) <= 0,
				(queryFacet[1] - currentClippingFacet.vertices[0]).dot(currentClippingFacet.normal) <= 0,
				(queryFacet[2] - currentClippingFacet.vertices[0]).dot(currentClippingFacet.normal) <= 0
			};

			if (!vertexIsInsideFacet[0] && !vertexIsInsideFacet[1] && !vertexIsInsideFacet[2]) {
				// Entire facet is outside the clipping plane
				continue;
			}

			Float3 clippedFacet[3];

			for (int vertexIndex = 0; vertexIndex < 3; vertexIndex++) {
				if (vertexIsInsideFacet[vertexIndex]) {
					clippedFacet[vertexIndex] = queryFacet[vertexIndex];
				}
				else if (vertexIsInsideFacet[(vertexIndex + 1) % 3]) {
					clippedFacet[vertexIndex] = currentClippingFacet.intersectionPoint(queryFacet[vertexIndex], queryFacet[(vertexIndex + 1) % 3]);
				}
				else {
					clippedFacet[vertexIndex] = currentClippingFacet.intersectionPoint(queryFacet[vertexIndex], queryFacet[(vertexIndex + 2) % 3]);
				}
			}

			if (nClippedFacetsToBeAdded == 8)
				printf("Too many facets to be added\n");
			clippedFacetsToBeAdded[nClippedFacetsToBeAdded][0] = clippedFacet[0];
			clippedFacetsToBeAdded[nClippedFacetsToBeAdded][1] = clippedFacet[1];
			clippedFacetsToBeAdded[nClippedFacetsToBeAdded][2] = clippedFacet[2];
			nClippedFacetsToBeAdded++;

		}

		prefixSumBuffer[threadIdx.x] = nClippedFacetsToBeAdded;
		__syncthreads();
		LAL::ExclusiveScan(prefixSumBuffer, blockDim.x);
		__syncthreads();

		for (int i = 0; i < nClippedFacetsToBeAdded; i++) {
			clippedFacets.triVertices[(prefixSumBuffer[threadIdx.x] + i) * 3 + 0] = clippedFacetsToBeAdded[i][0];
			clippedFacets.triVertices[(prefixSumBuffer[threadIdx.x] + i) * 3 + 1] = clippedFacetsToBeAdded[i][1];
			clippedFacets.triVertices[(prefixSumBuffer[threadIdx.x] + i) * 3 + 2] = clippedFacetsToBeAdded[i][2];
		}
		if (threadIdx.x == blockDim.x - 1) {
			clippedFacets.nTris = prefixSumBuffer[threadIdx.x] + nClippedFacetsToBeAdded;
		}
		__syncthreads();
	}

	if (threadIdx.x == 0)
		nVerticesOutBuffer[blockIdx.x] = clippedFacets.nTris * 3;

	// If we dont have a proper overlap, we can exit now
	if (clippedFacets.nTris < 2)
		return;

	for (int vertexIndex = threadIdx.x; vertexIndex < clippedFacets.nTris * 3; vertexIndex += blockDim.x) {
		const int moleculeOffset = blockIdx.x * maxFacetsInCH * 3;
		const int outputIndex = moleculeOffset + vertexIndex;
		verticesOutBuffer[outputIndex] = clippedFacets.triVertices[vertexIndex];
	}
	if (threadIdx.x == 0) {
		pivotPoints[blockIdx.x * maxCollisionsPerMH + mh1.nCollisions] = mh1.center;
		hullsBuffer[blockIdx.x].nCollisions++;
	}
}

__global__ void FindForceDirection(const MoleculeHull* const hullsBuffer, const Facet* const facets, const Float3* const intersectionCenters, 
	Float3* forceDirections, const int* const nVerticesInCollision, int targetMoleculeId) {
	__shared__ float minDistToFacet;
	__shared__ Float3 intersectionCenter;
	__shared__ Float3 forceDirection;
	__shared__ bool intersectionFound;
	__shared__ MoleculeHull otherMolecule;

	// First figure out if there is any work to do
	{
		if (threadIdx.x == 0) {
			intersectionFound = nVerticesInCollision[blockIdx.x] > 3;
		}
		__syncthreads();
		if (!intersectionFound)
			return;
	}

	// Setup shared data
	if (threadIdx.x == 0) {
		otherMolecule = hullsBuffer[targetMoleculeId];
		intersectionCenter = intersectionCenters[blockIdx.x];
		minDistToFacet = FLT_MAX;
	}
	__syncthreads();

	// Run algo, local first then sequential reduciton
	float localMinDist = FLT_MAX;
	Float3 localBestForceDirection{};
	for (int clippingFacetId = otherMolecule.indexOfFirstFacetInBuffer; clippingFacetId < otherMolecule.indexOfFirstFacetInBuffer + otherMolecule.nFacets; clippingFacetId++) {
		const float dist = std::abs(facets[clippingFacetId].signedDistance(intersectionCenter));

		if (dist < localMinDist) {
			localMinDist = dist;
			localBestForceDirection = facets[clippingFacetId].normal;
		}
	}
	for (int i = 0; i < blockDim.x; i++) {
		if (i == threadIdx.x) {
			if (localMinDist < minDistToFacet) {
				minDistToFacet = localMinDist;
				forceDirection = localBestForceDirection;
			}
		}
		__syncthreads();
	}

	// Finally push result 
	if (threadIdx.x == 0) {
		const int localCollisionId = hullsBuffer[blockIdx.x].nCollisions - 1;

		if (forceDirection.dot(hullsBuffer[blockIdx.x].center - otherMolecule.center) < 0) {
			forceDirection = hullsBuffer[blockIdx.x].center - otherMolecule.center;
		}

		forceDirections[blockIdx.x * maxCollisionsPerMH + localCollisionId] = forceDirection;
		//forceDirections[blockIdx.x * maxCollisionsPerMH + localCollisionId] = hullsBuffer[blockIdx.x].center - otherMolecule.center;
	}
}




__global__ void ApplyTransformations(MoleculeHull* moleculeHulls, Facet* facetsBuffer, RenderAtom* renderatomsBuffer, const glm::mat4* const transformMatrices, Float3 boxSize) {
	__shared__ glm::mat4 transformMatrix;
	__shared__ MoleculeHull moleculeHull;



	if (threadIdx.x == 0)
		moleculeHull = moleculeHulls[blockIdx.x];
	__syncthreads();


	for (int collisionId = 0; collisionId < moleculeHull.nCollisions; collisionId++) {

		const int matrixIndex = blockIdx.x * maxCollisionsPerMH + collisionId;

		if (threadIdx.x < 16)
			reinterpret_cast<float*>(&transformMatrix)[threadIdx.x] = reinterpret_cast<const float*>(&transformMatrices[matrixIndex])[threadIdx.x];




		for (int facetIdOffset = threadIdx.x; facetIdOffset < moleculeHull.nFacets; facetIdOffset += blockDim.x) {
			const int facetId = moleculeHull.indexOfFirstFacetInBuffer + facetIdOffset;

			Facet facet = facetsBuffer[facetId];

			for (int i = 0; i < 3; i++) {
				const glm::vec4 transformedVertex = transformMatrix * glm::vec4{ facet.vertices[i].x, facet.vertices[i].y, facet.vertices[i].z, 1.f };
				facet.vertices[i] = Float3{ transformedVertex.x, transformedVertex.y, transformedVertex.z };
			}

			facet.normal = ((facet.vertices[1] - facet.vertices[0]).norm()).cross(facet.vertices[2] - facet.vertices[0]).norm();
			facet.D = facet.normal.dot(facet.vertices[0]);

			facetsBuffer[facetId] = facet;
		}




		if (threadIdx.x == 0) {
			const glm::vec4 transformedCenter = transformMatrix * ToVec4(moleculeHull.center,1.f);
			moleculeHulls[blockIdx.x].center = Float3{ transformedCenter.x, transformedCenter.y, transformedCenter.z };
		}

		for (int particleId = moleculeHull.indexOfFirstParticleInBuffer; particleId < moleculeHull.indexOfFirstParticleInBuffer + moleculeHull.nParticles; particleId++) {
			const glm::vec4 pos{ renderatomsBuffer[particleId].position.x, renderatomsBuffer[particleId].position.y, renderatomsBuffer[particleId].position.z , 1.f };
			glm::vec4 transformedVertex = transformMatrix * pos;
			transformedVertex.w = renderatomsBuffer[particleId].position.w;

			renderatomsBuffer[particleId].position = float4{ transformedVertex.x, transformedVertex.y, transformedVertex.z, transformedVertex.w };// normalizedTransformedPosition.Tofloat4(renderatomsBuffer[particleId].position.w);
		}	
	}


	if (threadIdx.x == 0)
		moleculeHulls[blockIdx.x].nCollisions = 0;
}

// Define a functor for the transformation matrix computation
struct ComputeTransformationMatrix
{
	__host__ __device__
		glm::mat4 operator()(const thrust::tuple<glm::vec3, glm::vec3, glm::vec3, float>& t) const
	{
		const float magnitudeScalar = 0.003f;

		glm::vec3 rotationPivotPoint = thrust::get<0>(t);
		glm::vec3 forceApplicationPoint = thrust::get<1>(t);
		glm::vec3 forceDirection = thrust::get<2>(t);
		float forceMagnitude = thrust::get<3>(t) * magnitudeScalar;

		// Step 2: Compute the torque axis
		glm::vec3 torqueAxis = glm::cross(forceApplicationPoint - rotationPivotPoint, forceDirection);
		torqueAxis = glm::normalize(torqueAxis);

		// Step 3: Determine the rotation angle
		const float translationRatio = std::abs(glm::dot(glm::normalize(forceApplicationPoint - rotationPivotPoint), forceDirection));
		const float translationMagnitude = forceMagnitude * translationRatio;
		const float rotationMagnitude = forceMagnitude * (1.f - translationMagnitude);
		//float rotationAngle = forceMagnitude;

		// Step 4: Create the rotation matrix
		glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), rotationMagnitude, torqueAxis);

		// Step 5: Create translation matrices
		glm::mat4 translationMatrixToPivot = glm::translate(glm::mat4(1.0f), -rotationPivotPoint);
		glm::mat4 translationMatrixBack = glm::translate(glm::mat4(1.0f), rotationPivotPoint);

		// Step 6: Compute translation caused by the force
		glm::vec3 translationVector = forceMagnitude * forceDirection;
		glm::mat4 translationMatrix = glm::translate(glm::mat4(1.0f), translationVector);

		// Step 7: Combine the matrices
		glm::mat4 transformationMatrix = translationMatrix * translationMatrixBack * rotationMatrix * translationMatrixToPivot;

		return transformationMatrix;
	}
};



__global__ void MeasureOverallOverlapRadius(const MoleculeHull* const moleculeHulls, const float* const intersectionRadii, int nMolsTotal, float* results) {
	__shared__ double sharedSum;
	__shared__ float sharedMaxRatio;

	if (threadIdx.x == 0) {
		sharedSum = 0.f;
		sharedMaxRatio = 0.f;
	}
	__syncthreads();

	double localSum = 0.f;
	float localMaxRatio = 0.f;
	for (int index = threadIdx.x; index < nMolsTotal * maxCollisionsPerMH; index+=blockDim.x) {
		const float moleculeRadius = moleculeHulls[index / maxCollisionsPerMH].radius;
		localSum += intersectionRadii[index] / moleculeRadius;
		localMaxRatio = std::max(localMaxRatio, intersectionRadii[index] / moleculeRadius);
	}

	for (int i = 0; i < blockDim.x; i++) {
		if (i == threadIdx.x) {
			sharedSum += localSum;
			sharedMaxRatio = std::max(sharedMaxRatio, localMaxRatio);
		}
		__syncthreads();
	}

	if (threadIdx.x == 0) {
		results[0] = sharedSum;
		results[1] = sharedMaxRatio;
	}
}




void ConvexHullEngine::MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize, std::function<void()> callEachIteration)
{
	TimeIt timer{ "FindIntersect", false };

	Float3* verticesOutDev = nullptr;
	int* nVerticesOutDev = nullptr;

	cudaMalloc(&verticesOutDev, mhCol.nMoleculeHulls * maxFacetsInCH * 3 * sizeof(Float3));
	cudaMalloc(&nVerticesOutDev, mhCol.nMoleculeHulls * sizeof(int));




	// Setup buffers for computing transforms. These are filled when a MH finds a collision
	Float3* pivotPoints, *forceApplicationPoints, *forceDirections;
	float* forceMagnitudes;
	glm::mat4* transformMatrices;

	cudaMalloc(&pivotPoints, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(Float3));
	cudaMalloc(&forceApplicationPoints, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(Float3));
	cudaMalloc(&forceDirections, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(Float3));
	cudaMalloc(&forceMagnitudes, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(float));
	cudaMalloc(&transformMatrices, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(glm::mat4));


	thrust::device_ptr<glm::vec3> rotationPivotPointsPtr(reinterpret_cast<glm::vec3*>(pivotPoints));
	thrust::device_ptr<glm::vec3> forceApplicationPointsPtr(reinterpret_cast<glm::vec3*>(forceApplicationPoints));
	thrust::device_ptr<glm::vec3> forceDirectionsPtr(reinterpret_cast<glm::vec3*>(forceDirections));
	thrust::device_ptr<float> forceMagnitudesPtr(forceMagnitudes);
	thrust::device_ptr<glm::mat4> transformationMatricesPtr(transformMatrices);

	float* totalOverlapDev;
	cudaMalloc(&totalOverlapDev, sizeof(float)*2);

	

	std::vector<float> avgOverlapRatioHistory, maxOverlapHistory;
	const int maxIterations = 1000;
	for (int iteration = 0; iteration < maxIterations; iteration++) {
		TimeIt timer{ "FindIntersectIteration", false };

		bool anyIntersecting = false;

		// We set this to 0, so we can count the total sum in an iteration without accounting for some molecules not having all collisions filled
		cudaMemset(forceMagnitudes, 0, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(float)); 

		for (int queryMoleculeId = 0; queryMoleculeId < mhCol.nMoleculeHulls; queryMoleculeId++) {

			//ResetCollisionCounters<< <(mhCol.nMoleculeHulls + 64 -1)/ 64, 64 >> > (mhCol.moleculeHulls, mhCol.nMoleculeHulls);

			FindIntersectionConvexhullFrom2Convexhulls << <mhCol.nMoleculeHulls, nThreadsInFindIntersectionKernel >> > (mhCol.moleculeHulls, mhCol.facets, boxSize, queryMoleculeId,
				verticesOutDev, nVerticesOutDev, pivotPoints);
			LIMA_UTILS::genericErrorCheck("FindIntersectionConvexhullFrom2Convexhulls");
			cudaDeviceSynchronize();

			CalculateIntersectionCenterAndDepth<<<mhCol.nMoleculeHulls, 32>>>(
				verticesOutDev, nVerticesOutDev, mhCol.moleculeHulls, pivotPoints, forceApplicationPoints, forceDirections, forceMagnitudes
				);
			LIMA_UTILS::genericErrorCheck("CalculateIntersectionCenterAndDepth");
			
			FindForceDirection << <mhCol.nMoleculeHulls, 32 >> > (mhCol.moleculeHulls, mhCol.facets, forceApplicationPoints, forceDirections, nVerticesOutDev, queryMoleculeId);
		}

		// Apply the transformation using thrust::transform
		const int totalElements = mhCol.nMoleculeHulls * maxCollisionsPerMH;
		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(rotationPivotPointsPtr, forceApplicationPointsPtr, forceDirectionsPtr, forceMagnitudesPtr)),
			thrust::make_zip_iterator(thrust::make_tuple(rotationPivotPointsPtr + totalElements, forceApplicationPointsPtr + totalElements, forceDirectionsPtr + totalElements, forceMagnitudesPtr + totalElements)),
			transformationMatricesPtr,
			ComputeTransformationMatrix()
		);
		LIMA_UTILS::genericErrorCheck("Thrust error: ");

		cudaDeviceSynchronize();

		ApplyTransformations << <mhCol.nMoleculeHulls, 32 >> > (mhCol.moleculeHulls, mhCol.facets, mhCol.particles, transformMatrices, boxSize);
		LIMA_UTILS::genericErrorCheck("Apply transform: ");




		callEachIteration();

		MeasureOverallOverlapRadius<<<1, 128>>>( mhCol.moleculeHulls, forceMagnitudes, mhCol.nMoleculeHulls, totalOverlapDev);

		std::array<float,2> totalOverlapHost; // totalRatio, maxRatio;
		cudaMemcpy(totalOverlapHost.data(), totalOverlapDev, sizeof(float) * 2, cudaMemcpyDeviceToHost);

		const float avgOverlapRatio = totalOverlapHost[0] / static_cast<float>(mhCol.nMoleculeHulls);
		avgOverlapRatioHistory.push_back(avgOverlapRatio);
		maxOverlapHistory.push_back(totalOverlapHost[1]);
		if (avgOverlapRatio < 0.01f && totalOverlapHost[1] < 0.05f)
			break;
	}

	cudaFree(verticesOutDev);
	cudaFree(nVerticesOutDev);


	cudaFree(pivotPoints);
	cudaFree(forceApplicationPoints);
	cudaFree(forceDirections);
	cudaFree(forceMagnitudes);

	LIMA_UTILS::genericErrorCheck("Algo Exit");

	//matplotlibcpp::named_plot("Avg overlap", avgOverlapRatioHistory);
	//matplotlibcpp::named_plot("Max overlap", maxOverlapHistory);
	//matplotlibcpp::show();
}