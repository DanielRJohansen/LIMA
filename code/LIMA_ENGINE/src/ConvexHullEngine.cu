#include "ConvexHullEngine.cuh"
#include "Statistics.h"
#include "Utilities.h"

#include <glm.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <gtx/rotate_vector.hpp>
#undef GLM_ENABLE_EXPERIMENTAL

const int maxFacetsInCH = 128;


struct Overlap {
	bool isOverlapping = false;
	Float3 intersectionCenter{};
	float depth{};
};

struct Tri {
	Float3 vertices[3];
};

struct TriList {
	//Tri tris[maxFacetsInCH];
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

Float3 MeanPoint(const TriList& triList) {
	Float3 sum{};
	for (int i = 0; i < triList.nTris; i++) {
		for (int j = 0; j < 3; j++) {
			sum += triList.triVertices[i*3+j];
		}
	}
	return sum / (triList.nTris * 3);
}


Float3 CalculateMinimaxPoint(const Float3* const triVertices, int nVertices) {
	float epsilon = 1e-4;

	// Initial guess for the new point as the centroid
	//Float3 currentPoint = MeanPoint(triList);
	Float3 currentPoint = triVertices[0];

	float prevVelocity = std::numeric_limits<float>::max();

	const int maxIterations = 5;
	for (int i = 0; i < maxIterations; i++) {
		Float3 downGradient = { 0, 0, 0 };
		float maxDist = 0.0f;
		float secondMaxDist = 0.0f;

		for (int vertexIndex = 0; vertexIndex < nVertices; vertexIndex++) {

			const auto& point = triVertices[vertexIndex];

			float dist = (currentPoint - point).len();
			if (dist > maxDist) {
				secondMaxDist = maxDist;
				maxDist = dist;
				downGradient = (point - currentPoint);
			}
			else if (dist > secondMaxDist) {
				secondMaxDist = dist;
			}

		}

		const float velocity = std::min((maxDist - secondMaxDist) / 2.f, prevVelocity * 0.9f);
		prevVelocity = velocity;

		currentPoint += downGradient.norm() * velocity;


		//if (velocity < epsilon) {
		//	break;
		//}
	}

	return currentPoint;
}

// Input is in global memory
__global__ void CalculateIntersectionCenterAndDepth(const Float3* const triVerticesBuffer, const int* const nVerticesBuffer, Overlap* overlapsBuffer) {
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
		for (int vertexIndex = threadIdx.x; vertexIndex < nVertices; vertexIndex += blockDim.x) {
			const auto& point = triVertices[vertexIndex];
			const float dist = (overlap.intersectionCenter - point).len();

			if (dist > localMaxDist) {
				localSecondMaxDist = localMaxDist;
				localMaxDist = dist;
				localDownGradient = (point - overlap.intersectionCenter);
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
			float velocity = std::min((sharedMaxDist - sharedSecondMaxDist) / 2.0f, prevVelocity * 0.9f);
			prevVelocity = velocity;

			// Normalize the gradient
			sharedDownGradient = sharedDownGradient;
			overlap.intersectionCenter += sharedDownGradient.norm() * velocity;
		}
		__syncthreads();
	}


	float localMaxDepth = 0.f;
	for (int vertexIndex = threadIdx.x; vertexIndex < nVertices; vertexIndex += blockDim.x) {
		const auto& point = triVertices[vertexIndex];
		const float dist = (overlap.intersectionCenter - point).len();
		if (dist > localMaxDepth) {
			localMaxDepth = dist;
		}
	}
	for (int i = 0; i < blockDim.x; i++) {
		if (i == threadIdx.x) {
			if (localMaxDepth > overlap.depth)
				overlap.depth = localMaxDepth;
		}
		__syncthreads();
	}
	__syncthreads();
	
	if (threadIdx.x == 0) {
		overlap.isOverlapping = true;
		overlapsBuffer[blockIdx.x] = overlap;
	}
	__syncthreads();

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

__device__ void LargestDiffDev(Float3 queryPoint, TriList triList, float* result) {
	float largestDiffLocal = 0.f;

	for (int vertexId = threadIdx.x; vertexId < triList.nTris * 3; vertexId += blockDim.x) {
		const float diff = (queryPoint - triList.triVertices[vertexId]).len();
		if (diff > largestDiffLocal) {
			largestDiffLocal = diff;
		}
	}

	for (int i = 0; i < blockDim.x; i++) {
		if (i == threadIdx.x)
			if (largestDiffLocal > *result)
				*result = largestDiffLocal;
		__syncthreads();
	}
}

#include "DeviceAlgorithms.cuh"
#include "Utilities.h"
const int nThreadsInFindIntersectionKernel = 32;
__global__ void FindIntersectionConvexhullFrom2Convexhulls(MoleculeHull* hullsBuffer, Facet* facets, Float3 boxSize, int targetHull, Float3* verticesOutBuffer, int* nVerticesOutBuffer) {
	__shared__ TriList clippedFacets;
	//__shared__ TriList newClippedFacets;
	__shared__ Facet currentClippingFacet;
	__shared__ int prefixSumBuffer[nThreadsInFindIntersectionKernel];

	if (targetHull == blockIdx.x)
		return;

	MoleculeHull mh1 = hullsBuffer[blockIdx.x];
	MoleculeHull mh2 = hullsBuffer[targetHull];

	// First initialize the facets which vertices we are clipping
	for (int facetIdOffset = threadIdx.x; facetIdOffset < mh2.nFacets; facetIdOffset += blockDim.x) {
		const int facetId = mh2.indexOfFirstFacetInBuffer + facetIdOffset;
		Tri tri{ facets[facetId].vertices[0], facets[facetId].vertices[1], facets[facetId].vertices[2] };

		for (int vertexId = 0; vertexId < 3; vertexId++) 
			clippedFacets.triVertices[facetIdOffset * 3 + vertexId] = tri.vertices[vertexId];
	}
	if (threadIdx.x == 0) {
		clippedFacets.nTris = mh2.nFacets;
	}
	prefixSumBuffer[threadIdx.x] = 0;
	__syncthreads();



	for (int clippingFacetId = mh1.indexOfFirstFacetInBuffer; clippingFacetId < mh1.indexOfFirstFacetInBuffer + mh1.nFacets; clippingFacetId++) {
	//for (int clippingFacetIdOffset = 0; clippingFacetIdOffset < mh1.nFacets; clippingFacetIdOffset += blockDim.x) {

		if (threadIdx.x == 0) {
			//newClippedFacets.nTris = 0;
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
			//newClippedFacets.Add({ clippedFacet[0], clippedFacet[1], clippedFacet[2] });
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
			//clippedFacets.tris[prefixSumBuffer[threadIdx.x] + i] = { clippedFacetsToBeAdded[i][0], clippedFacetsToBeAdded[i][1], clippedFacetsToBeAdded[i][2] };
			clippedFacets.triVertices[(prefixSumBuffer[threadIdx.x] + i) * 3 + 0] = clippedFacetsToBeAdded[i][0];
			clippedFacets.triVertices[(prefixSumBuffer[threadIdx.x] + i) * 3 + 1] = clippedFacetsToBeAdded[i][1];
			clippedFacets.triVertices[(prefixSumBuffer[threadIdx.x] + i) * 3 + 2] = clippedFacetsToBeAdded[i][2];
		}
		if (threadIdx.x == blockDim.x - 1) {
			clippedFacets.nTris = prefixSumBuffer[threadIdx.x] + nClippedFacetsToBeAdded;

		}
		__syncthreads();
	}

	for (int vertexIndex = threadIdx.x; vertexIndex < clippedFacets.nTris * 3; vertexIndex += blockDim.x) {
		const int moleculeOffset = blockIdx.x * maxFacetsInCH * 3;
		const int outputIndex = moleculeOffset + vertexIndex;
		verticesOutBuffer[outputIndex] = clippedFacets.triVertices[vertexIndex];
	}
	if (threadIdx.x == 0) {
		nVerticesOutBuffer[blockIdx.x] = clippedFacets.nTris * 3;
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

#include "TimeIt.h"
#include<span>

void ConvexHullEngine::MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize) 
{
	TimeIt timer{ "FindIntersect", true };

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

	while (true) {

		bool anyIntersecting = false;

		for (int j = 0; j < mhCol.nMoleculeHulls; j++) {
			const Float3 rightCenter = Statistics::CalculateMinimaxPoint(moleculeHulls[j].GetFacetVertices(facets));
			const float rightRadius = LAL::LargestDiff(rightCenter, moleculeHulls[j].GetFacetVertices(facets) );



			FindIntersectionConvexhullFrom2Convexhulls << <mhCol.nMoleculeHulls, nThreadsInFindIntersectionKernel >> > (mhCol.moleculeHulls, mhCol.facets, boxSize, j, verticesOutDev, nVerticesOutDev);
			LIMA_UTILS::genericErrorCheck("FindIntersectionConvexhullFrom2Convexhulls");
			cudaDeviceSynchronize();

			CalculateIntersectionCenterAndDepth<<<mhCol.nMoleculeHulls, 32>>>(verticesOutDev, nVerticesOutDev, overlapDev);
			LIMA_UTILS::genericErrorCheck("CalculateIntersectionCenterAndDepth");

			// Now move the outputBuffers to host and compute Overlap

			cudaMemcpy(overlaps, overlapDev, mhCol.nMoleculeHulls * sizeof(Overlap), cudaMemcpyDeviceToHost);


			for (int i = 0; i < mhCol.nMoleculeHulls; i++) {

				Overlap devOverlap = overlaps[i];
				if (!devOverlap.isOverlapping)
					continue;

				const Float3 leftCenter = Statistics::CalculateMinimaxPoint(moleculeHulls[i].GetFacetVertices(facets));
				const float leftRadius = LAL::LargestDiff(leftCenter, moleculeHulls[i].GetFacetVertices(facets));

				if (i == j)
					continue;	


				anyIntersecting = true;

				const glm::mat4 leftTransformMatrix = computeTransformationMatrix(leftCenter.ToVec3(), devOverlap.intersectionCenter.ToVec3(), (leftCenter - rightCenter).norm().ToVec3(), devOverlap.depth * 0.05f);
				const glm::mat4 rightTransformMatrix = computeTransformationMatrix(rightCenter.ToVec3(), devOverlap.intersectionCenter.ToVec3(), (rightCenter - leftCenter).norm().ToVec3(), devOverlap.depth * 0.05f);

				moleculeHulls[i].ApplyTransformation(leftTransformMatrix, facets, atoms, boxSize);
				moleculeHulls[j].ApplyTransformation(rightTransformMatrix, facets, atoms, boxSize);

				cudaMemcpy(mhCol.facets, facets, mhCol.nFacets * sizeof(Facet), cudaMemcpyHostToDevice);
				cudaMemcpy(mhCol.particles, atoms, mhCol.nParticles * sizeof(RenderAtom), cudaMemcpyHostToDevice);

			}

		}

		cudaFree(verticesOutDev);
		cudaFree(nVerticesOutDev);
		cudaFree(overlapDev);

		delete[] overlaps;

		//}
		//cudaMemcpy(mhCol.facets, facets, mhCol.nFacets * sizeof(Facet), cudaMemcpyHostToDevice);
		//cudaMemcpy(mhCol.particles, atoms, mhCol.nParticles * sizeof(RenderAtom), cudaMemcpyHostToDevice);
		break;
	}


}