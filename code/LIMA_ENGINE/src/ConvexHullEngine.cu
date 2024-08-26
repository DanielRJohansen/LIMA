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


#include "DeviceAlgorithms.cuh"
#include "Utilities.h"

const int maxFacetsInCH = 128;
const int maxCollisionsPerMH = 8;

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
			sharedDownGradient = { 0, 0, 0 };
			sharedMaxDist = 0.0f;
			sharedSecondMaxDist = 0.0f;
		}
		__syncthreads();

		Float3 localDownGradient = { 0, 0, 0 };
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


__global__ void ResetCollisionCounters(MoleculeHull* hullsBuffer, int nMolecules) {
	const int mhIndex = blockIdx.x * blockDim.x + threadIdx.x;
	if (mhIndex < nMolecules)
		hullsBuffer[mhIndex].nCollisions = 0;
}


// Input is in global memory
__global__ void CalculateIntersectionCenterAndDepth(const Float3* const triVerticesBuffer, const int* const nVerticesBuffer, Overlap* overlapsBuffer, 
	MoleculeHull* hullsBuffer, Float3* pivotPoints, Float3* forceApplicationPoints, Float3* forceDirections, float* forceMagnitudes)
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
		if (nVertices < 4)
			overlapsBuffer[blockIdx.x].isOverlapping = false;
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
		overlap.isOverlapping = true;
		overlapsBuffer[blockIdx.x] = overlap;

		const float magnitudeScalar = 0.003f;

		const int collisionId = hullsBuffer[blockIdx.x].nCollisions-1; // -1 Because the previous kernel did the "++" 

		forceApplicationPoints[blockIdx.x * maxCollisionsPerMH + collisionId] = overlap.intersectionCenter;
		forceMagnitudes[blockIdx.x * maxCollisionsPerMH + collisionId] = overlap.depth * magnitudeScalar;
	}
}

float LargestDiff(Float3 queryPoint, TriList triList) {
	float largestDiff = 0.f;

	for (int triIndex = 0; triIndex < triList.nTris; triIndex++) {
		for (int vertexIndex = 0; vertexIndex < 3; vertexIndex++) {
			const auto& point = triList.triVertices[triIndex*3+vertexIndex];
			const float diff = (queryPoint - point).len();
			if (diff > largestDiff) {
				largestDiff = diff;
			}
		}
	}

	return largestDiff;
}



const int nThreadsInFindIntersectionKernel = 32;
__global__ void FindIntersectionConvexhullFrom2Convexhulls(MoleculeHull* hullsBuffer, const Facet* const facets, Float3 boxSize, int targetHull, 
	Float3* verticesOutBuffer, int* nVerticesOutBuffer, Float3* pivotPoints, Float3* forceDirections) {
	__shared__ TriList clippedFacets;
	__shared__ Facet currentClippingFacet;
	__shared__ int prefixSumBuffer[nThreadsInFindIntersectionKernel];

	if (targetHull == blockIdx.x) {
		if (threadIdx.x == 0)
			nVerticesOutBuffer[blockIdx.x] = 0;
		return;
	}
		

	MoleculeHull mh1 = hullsBuffer[blockIdx.x];
	MoleculeHull mh2 = hullsBuffer[targetHull];

	if ((mh1.center - mh2.center).len() > mh1.radius + mh2.radius)
		return;

	// First initialize the facets which vertices we are clipping
	for (int facetIdOffset = threadIdx.x; facetIdOffset < mh2.nFacets; facetIdOffset += blockDim.x) {
		const int facetId = mh2.indexOfFirstFacetInBuffer + facetIdOffset;

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
		forceDirections[blockIdx.x * maxCollisionsPerMH + mh1.nCollisions] = (mh1.center - mh2.center).norm();
		hullsBuffer[blockIdx.x].nCollisions++;
	}
}


glm::mat4 computeTransformationMatrix(const glm::vec3& rotationPivotPoint,
	const glm::vec3& forceApplicationPoint,
	const glm::vec3& forceDirection,
	float forceMagnitude) {


	// Step 2: Compute the torque axis (cross product between the force direction and the vector from pivot to application point)
	glm::vec3 torqueAxis = glm::cross(forceApplicationPoint - rotationPivotPoint, forceDirection);
	torqueAxis = glm::normalize(torqueAxis);

	// Step 3: Determine the rotation angle based on the force magnitude
	float rotationAngle = forceMagnitude;  // Scale this based on your needs

	// Step 4: Create the rotation matrix
	glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), rotationAngle, torqueAxis);

	// Step 5: Create the translation matrix to move the object to and from the pivot point
	glm::mat4 translationMatrixToPivot = glm::translate(glm::mat4(1.0f), -rotationPivotPoint);
	glm::mat4 translationMatrixBack = glm::translate(glm::mat4(1.0f), rotationPivotPoint);

	// Step 6: Compute translation caused by the force
	glm::vec3 translationVector = forceMagnitude * forceDirection;
	glm::mat4 translationMatrix = glm::translate(glm::mat4(1.0f), translationVector);
	// Step 7: Combine the matrices to get the final transformation matrix
	glm::mat4 transformationMatrix = translationMatrix * translationMatrixBack * rotationMatrix * translationMatrixToPivot;

	return transformationMatrix;
}


// Define a functor for the transformation matrix computation
struct ComputeTransformationMatrix
{
	__host__ __device__
		glm::mat4 operator()(const glm::vec3& rotationPivotPoint, const glm::vec3& forceApplicationPoint, const glm::vec3& forceDirection, float forceMagnitude)
	{
		// Step 2: Compute the torque axis
		glm::vec3 torqueAxis = glm::cross(forceApplicationPoint - rotationPivotPoint, forceDirection);
		torqueAxis = glm::normalize(torqueAxis);

		// Step 3: Determine the rotation angle
		float rotationAngle = forceMagnitude;

		// Step 4: Create the rotation matrix
		glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), rotationAngle, torqueAxis);

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



#include "TimeIt.h"
#include<span>

void ConvexHullEngine::MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize) 
{
	TimeIt timer{ "FindIntersect", false };

	Facet* facets = new Facet[mhCol.nFacets];
	RenderAtom* atoms = new RenderAtom[mhCol.nParticles];
	MoleculeHull* moleculeHulls = new MoleculeHull[mhCol.nMoleculeHulls];

	cudaMemcpy(facets, mhCol.facets, mhCol.nFacets * sizeof(Facet), cudaMemcpyDeviceToHost);
	cudaMemcpy(atoms, mhCol.particles, mhCol.nParticles * sizeof(RenderAtom), cudaMemcpyDeviceToHost);
	cudaMemcpy(moleculeHulls, mhCol.moleculeHulls, mhCol.nMoleculeHulls * sizeof(MoleculeHull), cudaMemcpyDeviceToHost);

	Float3* verticesOutDev = nullptr;
	int* nVerticesOutDev = nullptr;
	Overlap* overlapDev = nullptr;

	cudaMalloc(&verticesOutDev, mhCol.nMoleculeHulls * maxFacetsInCH * 3 * sizeof(Float3));
	cudaMalloc(&nVerticesOutDev, mhCol.nMoleculeHulls * sizeof(int));
	cudaMalloc(&overlapDev, mhCol.nMoleculeHulls * sizeof(Overlap));

	Overlap* overlaps = new Overlap[mhCol.nMoleculeHulls];





	// Setup buffers for computing transforms. These are filled when a MH finds a collision
	Float3* pivotPoints, *forceApplicationPoints, *forceDirections;
	float* forceMagnitudes;
	cudaMalloc(&pivotPoints, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(Float3));
	cudaMalloc(&forceApplicationPoints, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(Float3));
	cudaMalloc(&forceDirections, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(Float3));
	cudaMalloc(&forceMagnitudes, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(float));


	Float3* pivotPointsHost = new Float3[mhCol.nMoleculeHulls * maxCollisionsPerMH];
	Float3* forceApplicationPointsHost = new Float3[mhCol.nMoleculeHulls * maxCollisionsPerMH];
	Float3* forceDirectionsHost = new Float3[mhCol.nMoleculeHulls * maxCollisionsPerMH];
	float* forceMagnitudesHost = new float[mhCol.nMoleculeHulls * maxCollisionsPerMH];

	while (true) {

		bool anyIntersecting = false;

		for (int j = 0; j < mhCol.nMoleculeHulls; j++) {

			ResetCollisionCounters<< <(mhCol.nMoleculeHulls + 64 -1)/ 64, 64 >> > (mhCol.moleculeHulls, mhCol.nMoleculeHulls);

			FindIntersectionConvexhullFrom2Convexhulls << <mhCol.nMoleculeHulls, nThreadsInFindIntersectionKernel >> > (mhCol.moleculeHulls, mhCol.facets, boxSize, j, 
				verticesOutDev, nVerticesOutDev, pivotPoints, forceDirections);
			LIMA_UTILS::genericErrorCheck("FindIntersectionConvexhullFrom2Convexhulls");
			cudaDeviceSynchronize();

			CalculateIntersectionCenterAndDepth<<<mhCol.nMoleculeHulls, 32>>>(
				verticesOutDev, nVerticesOutDev, overlapDev, mhCol.moleculeHulls, pivotPoints, forceApplicationPoints, forceDirections, forceMagnitudes
				);
			LIMA_UTILS::genericErrorCheck("CalculateIntersectionCenterAndDepth");

			// Now move the outputBuffers to host and compute Overlap

			cudaMemcpy(overlaps, overlapDev, mhCol.nMoleculeHulls * sizeof(Overlap), cudaMemcpyDeviceToHost);

			cudaMemcpy(pivotPointsHost, pivotPoints, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(Float3), cudaMemcpyDeviceToHost);
			cudaMemcpy(forceApplicationPointsHost, forceApplicationPoints, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(Float3), cudaMemcpyDeviceToHost);
			cudaMemcpy(forceDirectionsHost, forceDirections, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(Float3), cudaMemcpyDeviceToHost);
			cudaMemcpy(forceMagnitudesHost, forceMagnitudes, mhCol.nMoleculeHulls * maxCollisionsPerMH * sizeof(float), cudaMemcpyDeviceToHost);









			for (int i = 0; i < mhCol.nMoleculeHulls; i++) {
				Overlap devOverlap = overlaps[i];
				if (!devOverlap.isOverlapping)
					continue;

				if (i == j)
					continue;

				const Float3 pivotPoint = pivotPointsHost[i * maxCollisionsPerMH + 0]; // Same as moleculeCenter
				const Float3 forceApplicationPoint = forceApplicationPointsHost[i * maxCollisionsPerMH + 0]; // center of intersect
				const float forceMagnitude = forceMagnitudesHost[i * maxCollisionsPerMH + 0] * 2.f;	// 
				const Float3 forceDirection = forceDirectionsHost[i * maxCollisionsPerMH + 0];

				anyIntersecting = true;

				const glm::mat4 leftTransformMatrix = computeTransformationMatrix(pivotPoint.ToVec3(), forceApplicationPoint.ToVec3(), forceDirection.ToVec3(), forceMagnitude);

				moleculeHulls[i].ApplyTransformation(leftTransformMatrix, facets, atoms, boxSize);

				cudaMemcpy(mhCol.facets, facets, mhCol.nFacets * sizeof(Facet), cudaMemcpyHostToDevice);
				cudaMemcpy(mhCol.particles, atoms, mhCol.nParticles * sizeof(RenderAtom), cudaMemcpyHostToDevice);
				cudaMemcpy(mhCol.moleculeHulls, moleculeHulls, mhCol.nMoleculeHulls * sizeof(MoleculeHull), cudaMemcpyHostToDevice);	// Only needed as long as we need to move center from host to device

			}

		}



		//}
		//cudaMemcpy(mhCol.facets, facets, mhCol.nFacets * sizeof(Facet), cudaMemcpyHostToDevice);
		//cudaMemcpy(mhCol.particles, atoms, mhCol.nParticles * sizeof(RenderAtom), cudaMemcpyHostToDevice);
		break;
	}

	cudaFree(verticesOutDev);
	cudaFree(nVerticesOutDev);
	cudaFree(overlapDev);

	delete[] overlaps;

	cudaFree(pivotPoints);
	cudaFree(forceApplicationPoints);
	cudaFree(forceDirections);
	cudaFree(forceMagnitudes);


	delete[] pivotPointsHost;
	delete[] forceApplicationPointsHost;
	delete[] forceDirectionsHost;
	delete[] forceMagnitudesHost;

	LIMA_UTILS::genericErrorCheck("Algo Exit");

}