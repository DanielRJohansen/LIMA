#include "CompoundBuilder.h"


void CompoundFactory::addParticle(const ParticleFactory& particle, int global_id, float boxlen_nm, BoundaryConditionSelect bc) {
	if (n_particles >= MAX_COMPOUND_PARTICLES) {
		throw std::runtime_error("Failed to add particle to compound");
	}

	// Hyperposition each compound relative to particle 0, so we can find key_particle and radius later
	Float3 hyperpos = particle.position;
	if (n_particles > 0) {
		BoundaryConditionPublic::applyHyperposNM(positions[0], hyperpos, boxlen_nm, bc); // i dont LOVE this.. :(
	}

	// Variables only present in factory
	positions[n_particles] = hyperpos;
	global_ids[n_particles] = global_id;

	// Variables present in Compound
	atom_types[n_particles] = particle.activeLJParamIndex;
	atomLetters[n_particles] = particle.topologyAtom.type[0];

	atom_charges[n_particles] = particle.topologyAtom.charge * elementaryChargeToKiloCoulombPerMole;
	atomMasses[n_particles] = particle.topologyAtom.mass / KILO;

	indicesInGrofile[n_particles] = particle.indexInGrofile;

	n_particles++;
}




std::array<int, CompoundInteractionBoundary::k> kMeansClusterCenters(const Float3* const positions, int n_elems, float boxlen_nm, BoundaryConditionSelect bc) {
	const int k = CompoundInteractionBoundary::k;
	std::array<int, k> center_indices{};


	// Initialize k centers randomly
	// Randomly pick k particles and set them as initial centers
	for (int i = 0; i < k; ++i) {
		//center_indices[i] = std::rand() % n_elems;
		center_indices[i] = std::min(i * n_elems, n_elems);
	}

	// Holds the index of the center that each point is closest to
	std::vector<int> labels(n_elems);

	// Max number of iterations or convergence criteria
	const int max_iter = 50;

	// Main k-means loop
	for (int iter = 0; iter < max_iter; ++iter) {
		// Assignment step
		// Go through each particle and assign it to the closest center
		for (int i = 0; i < n_elems; ++i) {
			float min_dist = std::numeric_limits<float>::infinity();
			for (int j = 0; j < k; ++j) {
				//float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM<PeriodicBoundaryCondition>(&positions[i], &positions[center_indices[j]]);
				float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(positions[i], positions[center_indices[j]], boxlen_nm, bc);
				if (dist < min_dist) {
					min_dist = dist;
					labels[i] = j; // Assign this particle to cluster j
				}
			}
		}

		// Update step
		// Calculate new centers as the mean of all particles assigned to each center
		std::vector<Float3> new_centers(k, Float3{});
		std::vector<int> counts(k, 0);

		for (int i = 0; i < n_elems; ++i) {
			int label = labels[i]; // Cluster label of this particle
			new_centers[label] += positions[i]; // Summing up for mean calculation
			counts[label] += 1; // Counting particles in each cluster for mean calculation
		}

		// Divide by the number of particles in each cluster to get the new center
		for (int j = 0; j < k; ++j) {
			if (counts[j] > 0) {
				new_centers[j] *= 1.f / static_cast<float>(counts[j]);
			}
		}

		// Find the index in the original positions array that is closest to the new centers
		for (int j = 0; j < k; ++j) {
			float min_dist = std::numeric_limits<float>::infinity();
			for (int i = 0; i < n_elems; ++i) {
				float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(positions[i], new_centers[j], boxlen_nm, bc);
				if (dist < min_dist) {
					min_dist = dist;
					center_indices[j] = i; // Update the center index to this particle
				}
			}
		}
	}

	return center_indices; // Return the indices of particles that are final centers
}


std::array<float, CompoundInteractionBoundary::k> clusterRadii(const Float3* const positions, int n_particles, const std::array<Float3, CompoundInteractionBoundary::k>& key_positions,
	float boxlen_nm, BoundaryConditionSelect bc) {
	std::array<float, CompoundInteractionBoundary::k> radii;
	std::vector<int> labels(n_particles);  // Holds the index of the closest key particle for each particle

	// Assign each particle to the closest key particle
	for (int i = 0; i < n_particles; ++i) {
		float min_dist = std::numeric_limits<float>::infinity();
		for (size_t j = 0; j < CompoundInteractionBoundary::k; ++j) {
			float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(positions[i], key_positions[j], boxlen_nm, bc);
			if (dist < min_dist) {
				min_dist = dist;
				labels[i] = j;  // Assign this particle to cluster j
			}
		}
	}

	// Calculate the furthest distance for each cluster
	for (size_t j = 0; j < CompoundInteractionBoundary::k; ++j) {
		float max_radius = 0.0f;  // Initialize maximum radius for this cluster

		for (int i = 0; i < n_particles; ++i) {
			if (labels[i] == j) {  // If this particle belongs to the current cluster
				float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(positions[i], key_positions[j], boxlen_nm, bc);
				if (dist > max_radius) {
					max_radius = dist; // Update maximum radius
				}
			}
		}

		// Save the maximum radius for this cluster
		radii[j] = max_radius;
	}

	return radii;
}

Float3 calcCOM(const Float3* positions, int n_elems, float boxlen_nm, BoundaryConditionSelect bc) {
	Float3 com{};
	const Float3& designatedCenterPosition = positions[0];
	for (int i = 0; i < n_elems; i++) {
		Float3 pos = positions[i];
		BoundaryConditionPublic::applyHyperposNM(designatedCenterPosition, pos, boxlen_nm, bc);	// Hyperpos around particle 0, since we dont know key position yet 
		com += pos;
	}
	return com / static_cast<float>(n_elems);
}


int indexOfParticleClosestToCom(const Float3* positions, int n_elems, const Float3& com, float boxlen_nm, BoundaryConditionSelect bc) {
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


void CompoundFactory::CalcCompoundMetaInfo(float boxlen_nm, std::vector<CompoundFactory>& compounds, BoundaryConditionSelect bc_select) {
	//#pragma omp parallel for // TODO Add OMP here, 
		//for (CompoundFactory& compound : compounds) {
	for (int cid = 0; cid < compounds.size(); cid++) {
		CompoundFactory& compound = compounds[cid];

		//for (auto& p : compound.atomLetters) { p = 'C'; }

		const int k = CompoundInteractionBoundary::k;
		std::array<int, k> key_indices = kMeansClusterCenters(compound.positions, compound.n_particles, boxlen_nm, bc_select);
		std::array<Float3, k> key_positions;
		for (int i = 0; i < k; i++) { key_positions[i] = compound.positions[key_indices[i]]; }

		std::array<float, k> radii = clusterRadii(compound.positions, compound.n_particles, key_positions, boxlen_nm, bc_select);

		for (int i = 0; i < k; i++) {
			compound.interaction_boundary.key_particle_indices[i] = key_indices[i];
			compound.interaction_boundary.radii[i] = radii[i] * 1.1f;	// Add 10% leeway
		}

		const Float3 com = calcCOM(compound.positions, compound.n_particles, boxlen_nm, bc_select);
		compound.centerparticle_index = indexOfParticleClosestToCom(compound.positions, compound.n_particles, com, boxlen_nm, bc_select);

		//char letter[3] = { 'O', 'N', 'P' };

		//for (auto index : key_indices)
		//	compound.atomLetters[index] = letter[cid];
	}

	// Calc absolute ids of particles in compounds 
	if (compounds.empty()) { return; }
	compounds[0].absoluteIndexOfFirstParticle = 0;
	for (int i = 1; i < compounds.size(); i++) {
		compounds[i].absoluteIndexOfFirstParticle = compounds[i - 1].absoluteIndexOfFirstParticle + compounds[i - 1].n_particles;
	}
}


