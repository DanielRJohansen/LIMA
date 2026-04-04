#include "CompoundBuilder.h"






//std::array<int, CompoundInteractionBoundary::k> kMeansClusterCenters(const Float3* const positions, int n_elems, const Float3& boxlen_nm, BoundaryConditionSelect bc) {
//	const int k = CompoundInteractionBoundary::k;
//	std::array<int, k> center_indices{};
//
//
//	// Initialize k centers randomly
//	// Randomly pick k particles and set them as initial centers
//	for (int i = 0; i < k; ++i) {
//		//center_indices[i] = std::rand() % n_elems;
//		center_indices[i] = std::min(i * n_elems, n_elems);
//	}
//
//	// Holds the index of the center that each point is closest to
//	std::vector<int> labels(n_elems);
//
//	// Max number of iterations or convergence criteria
//	const int max_iter = 50;
//
//	// Main k-means loop
//	for (int iter = 0; iter < max_iter; ++iter) {
//		// Assignment step
//		// Go through each particle and assign it to the closest center
//		for (int i = 0; i < n_elems; ++i) {
//			float min_dist = std::numeric_limits<float>::infinity();
//			for (int j = 0; j < k; ++j) {
//				//float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM<PeriodicBoundaryCondition>(&positions[i], &positions[center_indices[j]]);
//				float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(positions[i], positions[center_indices[j]], boxlen_nm, bc);
//				if (dist < min_dist) {
//					min_dist = dist;
//					labels[i] = j; // Assign this particle to cluster j
//				}
//			}
//		}
//
//		// Update step
//		// Calculate new centers as the mean of all particles assigned to each center
//		std::vector<Float3> new_centers(k, Float3{});
//		std::vector<int> counts(k, 0);
//
//		for (int i = 0; i < n_elems; ++i) {
//			int label = labels[i]; // Cluster label of this particle
//			new_centers[label] += positions[i]; // Summing up for mean calculation
//			counts[label] += 1; // Counting particles in each cluster for mean calculation
//		}
//
//		// Divide by the number of particles in each cluster to get the new center
//		for (int j = 0; j < k; ++j) {
//			if (counts[j] > 0) {
//				new_centers[j] *= 1.f / static_cast<float>(counts[j]);
//			}
//		}
//
//		// Find the index in the original positions array that is closest to the new centers
//		for (int j = 0; j < k; ++j) {
//			float min_dist = std::numeric_limits<float>::infinity();
//			for (int i = 0; i < n_elems; ++i) {
//				float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(positions[i], new_centers[j], boxlen_nm, bc);
//				if (dist < min_dist) {
//					min_dist = dist;
//					center_indices[j] = i; // Update the center index to this particle
//				}
//			}
//		}
//	}
//
//	return center_indices; // Return the indices of particles that are final centers
//}



Float3 calcCOM(const Float3* positions, int n_elems, const Float3& boxlen_nm, BoundaryConditionSelect bc) {
	Float3 com{};
	const Float3& designatedCenterPosition = positions[0];
	for (int i = 0; i < n_elems; i++) {
		Float3 pos = positions[i];
		BoundaryConditionPublic::applyHyperposNM(designatedCenterPosition, pos, boxlen_nm, bc);	// Hyperpos around particle 0, since we dont know key position yet 
		com += pos;
	}
	return com / static_cast<float>(n_elems);
}


int indexOfParticleClosestToCom(const Float3* positions, int n_elems, const Float3& com, const Float3& boxlen_nm, BoundaryConditionSelect bc) {
	int closest_particle_index = 0;
	float closest_particle_distance = std::numeric_limits<float>::infinity();
	for (int i = 0; i < n_elems; i++) {
		const float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(positions[i], com, boxlen_nm, bc);
		if (dist < closest_particle_distance) {
			closest_particle_distance = dist;
			closest_particle_index = i;
		}
	}
	return closest_particle_index;
}



