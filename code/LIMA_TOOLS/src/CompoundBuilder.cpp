#include "CompoundBuilder.h"
#include "Forcefield.h"

#include <unordered_set>
#include <format>
#include <array>

using namespace LIMA_MOLECULEBUILD;

bool compoundsAreAlreadyConnected(const std::vector<int> compound_ids, const std::vector<BridgeFactory>& compoundbriges) 
{
	for (auto& bridge : compoundbriges) {

		bool allCompoundIdsPresent = true;
		for (auto& compound_id : compound_ids) {
			if (!bridge.containsCompound(compound_id)) {
				allCompoundIdsPresent = false;
			}
		}

		if (allCompoundIdsPresent) {
			return true;
		}
	}

	return false;
}


bool bridgeContainsTheseTwoCompounds(const BridgeFactory& bridge, const std::array<int, 2> compound_ids) {
	for (const int id : compound_ids) {
		bool found = false;

		for (auto compoundId : bridge.compound_ids) {
			if (compoundId == id) {
				found = true;
				break;
			}
		}

		if (!found)
			return false;
	}

	return true;
}
















std::array<int, CompoundInteractionBoundary::k> kMeansClusterCenters(const Float3* const positions, int n_elems, float boxlen_nm, BoundaryConditionSelect bc) {
	const int k = CompoundInteractionBoundary::k;
	std::array<int, k> center_indices{};


	// Initialize k centers randomly
	// Randomly pick k particles and set them as initial centers
	for (int i = 0; i < k; ++i) {
		//center_indices[i] = std::rand() % n_elems;
		center_indices[i] = std::min(i*n_elems, n_elems);
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
				new_centers[j] *= 1.f/static_cast<float>(counts[j]);
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




// --------------------------------------------------------------- Factory Functions --------------------------------------------------------------- //
void CompoundFactory::addParticle(const Float3& position, int atomtype_id, char atomLetter,
	int global_id, float boxlen_nm, BoundaryConditionSelect bc, float charge) {
	if (n_particles >= MAX_COMPOUND_PARTICLES) {
		throw std::runtime_error("Failed to add particle to compound");
	}

	// Hyperposition each compound relative to particle 0, so we can find key_particle and radius later
	Float3 hyperpos = position;
	if (n_particles > 0) {
		BoundaryConditionPublic::applyHyperposNM(positions[0], hyperpos, boxlen_nm, bc);
	}

	// Variables only present in factory
	positions[n_particles] = hyperpos;
	global_ids[n_particles] = global_id;

	// Variables present in Compound
	atom_types[n_particles] = atomtype_id;
	atomLetters[n_particles] = atomLetter;

	atom_charges[n_particles] = charge;

	n_particles++;
}

void CompoundFactory::AddBond(const std::vector<ParticleInfo>& atoms, const SingleBondFactory& bond) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add singlebond to compound");
	}
	singlebonds[n_singlebonds++] = SingleBond(
		{
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[0]].localIdInCompound),
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[1]].localIdInCompound)
		},
		bond.params.b0,
		bond.params.kb
		);
}
void CompoundFactory::AddBond(const std::vector<ParticleInfo>& atoms, const AngleBondFactory& bond) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add anglebond to compound");
	}
	anglebonds[n_anglebonds++] = AngleUreyBradleyBond(
		{
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[0]].localIdInCompound),
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[1]].localIdInCompound),
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[2]].localIdInCompound)
		},
		bond.params.theta0,
		bond.params.kTheta,
		bond.params.ub0,
		bond.params.kUB
		);
}
void CompoundFactory::AddBond(const std::vector<ParticleInfo>& atoms, const DihedralBondFactory& bond) {
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add dihedralbond to compound");
	}
	dihedrals[n_dihedrals++] = DihedralBond(
		{
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[0]].localIdInCompound),
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[1]].localIdInCompound),
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[2]].localIdInCompound),
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[3]].localIdInCompound)
		},
		bond.params.phi_0,
		bond.params.k_phi,
		bond.params.n
		);
}
void CompoundFactory::AddBond(const std::vector<ParticleInfo>& atoms, const ImproperDihedralBondFactory& bond) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add improperdihedralbond to compound");
	}
	impropers[n_improperdihedrals++] = ImproperDihedralBond(
		{
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[0]].localIdInCompound),
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[1]].localIdInCompound),
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[2]].localIdInCompound),
			static_cast<uint8_t>(atoms[bond.global_atom_indexes[3]].localIdInCompound)
		},
		bond.params.psi_0,
		bond.params.k_psi
		);
}

void CompoundFactory::addIdOfBondedCompound(int id) {
	if (n_bonded_compounds == max_bonded_compounds) { throw std::runtime_error("Failed to add bonded compound id to compound"); }

	for (int i = 0; i < n_bonded_compounds; i++) {
		// If the other compound is already saved, move do nothing
		if (bonded_compound_ids[i] == id)
			return;
	}
	bonded_compound_ids[n_bonded_compounds++] = id;
}

template <typename BondFactory_type>
std::array<uint8_t, BondFactory_type::nAtoms> BridgeFactory::ConvertGlobalIdsToCompoundlocalIds (std::vector<ParticleInfo>& particle_info, const BondFactory_type& bond) {
	std::array<uint8_t, BondFactory_type::nAtoms> localIds;
	for (int i = 0; i < BondFactory_type::nAtoms; i++) {
		localIds[i] = getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[i]]);
	}
	return localIds;
}

void BridgeFactory::AddBond(std::vector<ParticleInfo>& particle_info, const SingleBondFactory& bond) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add singlebond to bridge"); }
	singlebonds[n_singlebonds++] = SingleBond{
		ConvertGlobalIdsToCompoundlocalIds(particle_info, bond),
		bond.params.b0,
		bond.params.kb
	};
}
void BridgeFactory::AddBond(std::vector<ParticleInfo>& particle_info, const AngleBondFactory& bond) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add anglebond to bridge"); }
	anglebonds[n_anglebonds++] = AngleUreyBradleyBond{
		ConvertGlobalIdsToCompoundlocalIds(particle_info, bond),
		bond.params.theta0,
		bond.params.kTheta,
		bond.params.ub0,
		bond.params.kUB
	};
}
void BridgeFactory::AddBond(std::vector<ParticleInfo>& particle_info, const DihedralBondFactory& bond) {
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_BRIDGE) {
		throw std::runtime_error("Failed to add dihedralbond to bridge");
	}
	dihedrals[n_dihedrals++] = DihedralBond{
		ConvertGlobalIdsToCompoundlocalIds(particle_info, bond),
		bond.params.phi_0,
		bond.params.k_phi,
		bond.params.n
	};
}
void BridgeFactory::AddBond(std::vector<ParticleInfo>& particle_info, const ImproperDihedralBondFactory& bond) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add improperdihedralbond to bridge"); }
	impropers[n_improperdihedrals++] = ImproperDihedralBond{
		ConvertGlobalIdsToCompoundlocalIds(particle_info, bond),
		bond.params.psi_0,
		bond.params.k_psi,
	};
}


uint8_t BridgeFactory::getBridgelocalIdOfParticle(ParticleInfo& particle_info) {

	// Assign bridge id to particle if it doesnt have one
	if (particle_info.bridgeId == -1) {
		particle_info.bridgeId = this->bridge_id;
	}
	// If it has another id, fail as we dont support that for now
	else if (particle_info.bridgeId != this->bridge_id) {
		//printf("Particle global id %d\n", particle_info.globalId);
		throw std::runtime_error(std::format("Cannot add particle to this bridge ({}) as it already has another bridge ({})", bridge_id, particle_info.bridgeId).c_str());
	}

	if (particle_info.localIdInBridge == -1) {	// This limits the amount of particles in a bridge
		addParticle(particle_info);
	}

	return particle_info.localIdInBridge;
}

void BridgeFactory::addParticle(ParticleInfo& particle_info) {	// This can be made const, and return it's iD in bridge
	if (n_particles == MAX_PARTICLES_IN_BRIDGE) { throw std::runtime_error("Failed to add particle to bridge"); }

	particle_info.localIdInBridge = n_particles;

	int compoundlocalid_in_bridge = -1;
	for (int i = 0; i < this->n_compounds; i++) {
		if (particle_info.compoundId == compound_ids[i])
			compoundlocalid_in_bridge = i;
	}
	if (compoundlocalid_in_bridge == -1) {
		throw std::runtime_error("Could not find compoundlocalid_in_bridge");
	}

	particle_refs[n_particles] = ParticleReference{
		particle_info.compoundId,
		particle_info.localIdInCompound,
#ifdef LIMAKERNELDEBUGMODE
		particle_info.global_id,
#endif
		static_cast<uint8_t>(compoundlocalid_in_bridge)
	};

	n_particles++;
}





struct TopologyFileRef {
	TopologyFileRef(int offset, const TopologyFile& topol) : atomsOffset(offset), topology(topol) {}
	int atomsOffset = 0;	// TODO this logic is being moved into the topology file
	const TopologyFile& topology;
};

std::vector<TopologyFileRef> PrepareTopologies(const TopologyFile& topol) {
	std::vector<TopologyFileRef> topologies;
	int offset = 0;

	topologies.push_back({ offset, topol });
	offset += topol.GetLocalAtoms().size();
	for (const auto& molecule : topol.GetAllSubMolecules()) {
		topologies.emplace_back(TopologyFileRef{offset, *molecule.includeTopologyFile });
		offset += molecule.includeTopologyFile->GetLocalAtoms().size();
	}
	return topologies;
}



std::vector<ParticleInfo> PrepareAtoms(const std::vector<TopologyFileRef>& topologyFiles, const GroFile& grofile, ForcefieldManager& forcefield, int nNonsolventAtoms) {
	std::vector<ParticleInfo> atomRefs(nNonsolventAtoms);
	
	int globalIndex = 0;
	int uniqueResId = 0;
	for (const auto& topology : topologyFiles) {
		bool atomsAddedThisTop = false;
		for (const auto& atom : topology.topology.GetLocalAtoms()) {
			atomsAddedThisTop = true;
			// TODO: Handle ignore atoms (H) somehow?			
			// TODO: When accessing gro file, first check that it is not a solvent

			if (globalIndex > 0 && atomRefs[globalIndex - 1].topAtom->resnr != atom.resnr) {
				uniqueResId++;
			}

			if (grofile.atoms[globalIndex].atomName != atom.atomname || grofile.atoms[globalIndex].residueName.substr(0, 3) != atom.residue.substr(0,3))
				throw std::runtime_error("Atom names do not match between gro and topology file");			

			const int activeLJParamIndex = forcefield.GetActiveLjParameterIndex(topology.topology.forcefieldIncludes, atom.type);
			atomRefs[globalIndex] = ParticleInfo{ &grofile.atoms[globalIndex], &atom, activeLJParamIndex, uniqueResId };
			atomRefs[globalIndex].sourceLine = grofile.atoms[globalIndex].sourceLine;
			
			globalIndex++;		
		}

		if (atomsAddedThisTop) { uniqueResId++; }
	}

	// TODO: At this point we still need to load the solvents, they should come now
	assert(globalIndex == atomRefs.size());

	return atomRefs;
}
std::vector<Float3> LoadSolventPositions(const GroFile& grofile) {
	std::vector<Float3> solventPositions;// (grofile.atoms.size() - nNonsolventAtoms);
	solventPositions.reserve(grofile.atoms.size());

	for (const auto& atom : grofile.atoms)
		if (atom.residueName == "SOL" && (atom.atomName[0] == 'O'))
			solventPositions.emplace_back(atom.position);
	
	return solventPositions;
}


template <int n, typename GenericBond, typename BondtypeFactory>
void LoadBondIntoTopology(const int bondIdsRelativeToTopolFile[n], int atomIdOffset, ForcefieldManager& forcefield,
	const std::vector<ParticleInfo>& atomRefs, std::vector<BondtypeFactory>& topology, const std::vector<fs::path>& forcefieldNames, std::string sourceLine = "")
{
	std::array<uint32_t, n> globalIds{};
	for (int i = 0; i < n; i++) {
		globalIds[i] = bondIdsRelativeToTopolFile[i] + atomIdOffset;
		if (globalIds[i] >= atomRefs.size())
			return;
			//throw std::runtime_error("Bond refers to atom that has not been loaded");
	}

	std::array<std::string, n> atomTypenames{};
	for (int i = 0; i < n; i++) {
		atomTypenames[i] = atomRefs[globalIds[i]].topAtom->type;
	}

	// TODO: This is the correct place to make this check, but we cant implement it here untill we do a graph-based search for LJ ignores,
	// since the bonds here are currently used for that
	//if (forcefield.BondHasZeroParameter<typename GenericBond::Parameters>(atomTypenames))
	//	return;

	auto bondParams = forcefield.GetBondParameters<GenericBond>(forcefieldNames, atomTypenames);

	for (const auto& param : bondParams) {
		topology.emplace_back(BondtypeFactory{ globalIds, param });
		topology.back().sourceLine = sourceLine;
	}
}

void ReserveSpaceForAllBonds(Topology& topology, const std::vector<TopologyFileRef>& topologyFiles) {
	int expectedSinglebonds = 0;
	int expectedAnglebonds = 0;
	int expectedDihedralbonds = 0;
	int expectedImproperDihedralbonds = 0;
	for (const auto& topologyFile : topologyFiles) {
		expectedSinglebonds += topologyFile.topology.GetLocalSinglebonds().size();
		expectedAnglebonds += topologyFile.topology.GetLocalAnglebonds().size();
		expectedDihedralbonds += topologyFile.topology.GetLocalDihedralbonds().size();
		expectedImproperDihedralbonds += topologyFile.topology.GetLocalImproperDihedralbonds().size();
	}
	topology.singlebonds.reserve(expectedSinglebonds);
	topology.anglebonds.reserve(expectedAnglebonds);
	topology.dihedralbonds.reserve(expectedDihedralbonds);
	topology.improperdihedralbonds.reserve(expectedImproperDihedralbonds);
}

Topology LoadTopology(const std::vector<TopologyFileRef>& topologyFiles, ForcefieldManager& forcefield, const std::vector<ParticleInfo>& atomRefs) {
	Topology topology;	

	ReserveSpaceForAllBonds(topology, topologyFiles);

	for (const auto& topologyFile : topologyFiles) {

		for (const auto& bondTopol : topologyFile.topology.GetLocalSinglebonds()) {			
			LoadBondIntoTopology<2, SingleBond, SingleBondFactory>(
				bondTopol.ids, topologyFile.atomsOffset, forcefield, atomRefs, topology.singlebonds, topologyFile.topology.forcefieldIncludes, bondTopol.sourceLine);
		}

		for (const auto& bondTopol : topologyFile.topology.GetLocalAnglebonds()) {			
			LoadBondIntoTopology<3, AngleUreyBradleyBond, AngleBondFactory>(
				bondTopol.ids, topologyFile.atomsOffset, forcefield, atomRefs, topology.anglebonds, topologyFile.topology.forcefieldIncludes);
		}

		for (const auto& bondTopol : topologyFile.topology.GetLocalDihedralbonds()) {			
			LoadBondIntoTopology<4, DihedralBond, DihedralBondFactory>(
				bondTopol.ids, topologyFile.atomsOffset, forcefield, atomRefs, topology.dihedralbonds, topologyFile.topology.forcefieldIncludes);
		}

		for (const auto& bondTopol : topologyFile.topology.GetLocalImproperDihedralbonds()) {
			LoadBondIntoTopology<4, ImproperDihedralBond, ImproperDihedralBondFactory>(
				bondTopol.ids, topologyFile.atomsOffset, forcefield, atomRefs, topology.improperdihedralbonds, topologyFile.topology.forcefieldIncludes);
		}
	}

	return topology;
}



std::vector<std::vector<int>> MapAtomsToSinglebonds(const std::vector<ParticleInfo>& preparedAtoms, const Topology& topology) {
	std::vector<std::vector<int>> atomIdToSinglebondsMap(preparedAtoms.size());


	for (int i = 0; i < topology.singlebonds.size(); i++) {
		const auto& bond = topology.singlebonds[i];
		atomIdToSinglebondsMap[bond.global_atom_indexes[0]].push_back(i);
		atomIdToSinglebondsMap[bond.global_atom_indexes[1]].push_back(i);
	}

	return atomIdToSinglebondsMap;
}


// Usually just a residue, but we sometimes also need to split lipids into smaller groups 
struct AtomGroup {
	std::vector<int> atomIds;
};

std::vector<AtomGroup> GroupAtoms(const std::vector<ParticleInfo>& preparedAtoms, const Topology& topology) {
	const int nResidues = preparedAtoms.empty() ? 0 : preparedAtoms.back().uniqueResId + 1;
	std::vector<AtomGroup> atomGroups;
	atomGroups.reserve(nResidues);

	int currentResidueId = -1;

	for (int particleId = 0; particleId < preparedAtoms.size(); particleId++) {
		const ParticleInfo& atom = preparedAtoms[particleId];
		if (atomGroups.empty() || atom.topAtom->section_name.has_value() || atom.uniqueResId != preparedAtoms[particleId-1].uniqueResId)
			atomGroups.push_back(AtomGroup{});
		atomGroups.back().atomIds.push_back(particleId);
		//atomGroups[atom.uniqueResId].atomIds.push_back(particleId);
	}

	return atomGroups;
}	


bool areBonded(const AtomGroup& left, const AtomGroup& right, const std::vector<ParticleInfo>& preparedAtoms, const std::vector<std::vector<int>>& atomIdToSinglebondsMap) {
	//assert(left != right);

	for (auto& atomleft_gid : left.atomIds) {
		for (auto& atomright_gid : right.atomIds) {
			const std::vector<int>& atomleft_singlesbonds = atomIdToSinglebondsMap[atomleft_gid];
			const std::vector<int>& atomright_singlesbonds = atomIdToSinglebondsMap[atomright_gid];
			
			// Check if any singlebond id from left matches any in right
			for (int i = 0; i < atomleft_singlesbonds.size(); i++) {
				for (int j = 0; j < atomright_singlesbonds.size(); j++) {
					if (atomleft_singlesbonds[i] == atomright_singlesbonds[j]) {
						return true;
					}
				}
			}
		}
	}
	return false;
}

std::vector<CompoundFactory> CreateCompounds(const Topology& topology, float boxlen_nm, const std::vector<AtomGroup>& residues,
	std::vector<ParticleInfo>& preparedAtoms, const std::vector<std::vector<int>>& atomIdToSinglebondsMap, BoundaryConditionSelect bc_select)
{
	std::vector<CompoundFactory> compounds;
	int current_residue_id = -1;

	for (int residue_index = 0; residue_index < residues.size(); residue_index++) {
		const AtomGroup& residue = residues[residue_index];

		const bool is_bonded_with_previous_residue = residue_index > 0 && areBonded(residues[residue_index - 1], residue, preparedAtoms, atomIdToSinglebondsMap);
		const bool compound_has_room_for_residue = residue_index > 0 && compounds.back().hasRoomForRes(residue.atomIds.size());

		// If we are either a new molecule, or same molecule but the current compound has no more room, make new compound
		if (residue_index == 0 || !is_bonded_with_previous_residue || !compound_has_room_for_residue) {
			if (compounds.size() >= MAX_COMPOUNDS) {
				throw std::runtime_error(std::format("Cannot handle more than {} compounds", MAX_COMPOUNDS).c_str());
			}

			compounds.push_back(CompoundFactory{ static_cast<int>(compounds.size()) });

			if (is_bonded_with_previous_residue) {
				compounds.back().indicesOfBondconnectedCompounds.push_back(compounds.size() - 2);
				compounds[compounds.size() - 2].indicesOfBondconnectedCompounds.push_back(compounds.size() - 1);
			}
		}
		

		// Add all atoms of residue to current compound
		for (int i = 0; i < residue.atomIds.size(); i++) {
			const int atom_gid = residue.atomIds[i];
			compounds.back().addParticle(
				preparedAtoms[atom_gid].groAtom->position,
				preparedAtoms[atom_gid].activeLjtypeParameterIndex,
				preparedAtoms[atom_gid].topAtom->type[0],
				atom_gid,
				boxlen_nm, 
				bc_select,
				preparedAtoms[atom_gid].topAtom->charge * elementaryChargeToKiloCoulombPerMole
			);

			preparedAtoms[atom_gid].compoundId = compounds.size() - 1;
			preparedAtoms[atom_gid].localIdInCompound = compounds.back().n_particles - 1;
		}
	}

	return compounds;
}


std::vector<BridgeFactory> CreateBridges(const std::vector<SingleBondFactory>& singlebonds, const std::vector<CompoundFactory>& compounds, const std::vector<ParticleInfo>& atoms)
{
	std::vector<BridgeFactory> compoundBridges;


	// First find which ompounds are bonded to other compounds. I assume singlebonds should find all these, 
	// since angles and dihedrals make no sense if the singlebond is not present, as the distances between the atoms can then be extremely large
	std::vector<std::unordered_set<int>> compound_to_bondedcompound_table(compounds.size());
	for (const auto& bond : singlebonds) {
		const int particle_gid_left = bond.global_atom_indexes[0];
		const int particle_gid_right = bond.global_atom_indexes[1];

		const int compound_gid_left = atoms[particle_gid_left].compoundId;
		const int compound_gid_right = atoms[particle_gid_right].compoundId;

		if (compound_gid_left != compound_gid_right) {
			compound_to_bondedcompound_table[compound_gid_left].insert(compound_gid_right);
			compound_to_bondedcompound_table[compound_gid_right].insert(compound_gid_left);
		}
	}

	// Knowing which compounds are bonded we can create the bridges.
	for (int compound_id_self = 0; compound_id_self < compound_to_bondedcompound_table.size(); compound_id_self++) {
		const std::unordered_set<int>& bonded_compounds_set = compound_to_bondedcompound_table[compound_id_self];

		// if a compound is the entire molecule, this set will be empty, and that is alright
		if (bonded_compounds_set.empty()) {
			continue;
		}

		// A compound can not be bonded to itself
		if (compound_id_self == *std::min_element(bonded_compounds_set.begin(), bonded_compounds_set.end())) {
			throw std::runtime_error("Something went wrong in the createBridges algorithm");
		}


		// Each compound is only allowed to create a bridge to compounds with a larger id
		std::vector<int> bridge_compound_ids{ compound_id_self };
		for (auto compound_id_other : bonded_compounds_set) {
			if (compound_id_other > compound_id_self) {
				bridge_compound_ids.push_back(compound_id_other);
			}
		}

		// Note that this bridge might contain more compounds than this compound proposes to bridge, as it was made by a compound with a lower id, that is compound does not include!
		if (!compoundsAreAlreadyConnected(bridge_compound_ids, compoundBridges)) {
			compoundBridges.push_back(BridgeFactory{ static_cast<int>(compoundBridges.size()), bridge_compound_ids });

		}
	}

	return compoundBridges;
}









template <int n_ids>
bool SpansTwoCompounds(const std::array<uint32_t, n_ids> bond_globalids, const std::vector<ParticleInfo>& atoms) {
	const int compoundid_left = atoms[bond_globalids[0]].compoundId;

	for (auto globalId : bond_globalids) {
		if (compoundid_left != atoms[globalId].compoundId) {
			return true;
		}
	}
	return false;
}

template <int n_ids>
void DistributeLJIgnores(BondedParticlesLUTManager* bplut_man, const std::vector<ParticleInfo>& atoms, const std::array<uint32_t, n_ids>& global_ids) {
	for (auto gid_self : global_ids) {
		for (auto gid_other : global_ids) {

			if (gid_self == gid_other) { continue; }
		
			const int compoundIndexSelf = atoms[gid_self].compoundId;
			const int compoundIndexOther = atoms[gid_other].compoundId;


			bplut_man->addNewConnectedCompoundIfNotAlreadyConnected(compoundIndexSelf, compoundIndexOther);
			BondedParticlesLUT* lut = bplut_man->get(compoundIndexSelf, compoundIndexOther, true);
			lut->set(atoms[gid_self].localIdInCompound, atoms[gid_other].localIdInCompound, true);
		}
	}
}

// Returns two compound id's of a bond in a bridge. The id's are sorted with lowest first
template <int N>
std::array<int, 2> getTheTwoDifferentIds(const std::array<uint32_t, N>& particle_global_ids, const std::vector<ParticleInfo>& atoms) {
	std::array<int, 2> compound_ids = { atoms[particle_global_ids[0]].compoundId, -1 };

	for (auto globalId : particle_global_ids) {
		int compoundid_of_particle = atoms[globalId].compoundId;
		if (compoundid_of_particle != compound_ids[0]) {
			compound_ids[1] = compoundid_of_particle;
			break;
		}
	}

	if (compound_ids[1] == -1) {
		throw std::runtime_error("Failed to find the second compound of bridge");
	}

	if (compound_ids[0] > compound_ids[1]) {
		std::swap(compound_ids[0], compound_ids[1]);
	}

	return compound_ids;
}

template <int N>
BridgeFactory& getBridge(std::vector<BridgeFactory>& bridges, const std::array<uint32_t, N>& particle_global_ids, const std::vector<ParticleInfo>& atoms,
	std::unordered_map<int, std::vector<BridgeFactory*>>& compoundToBridges
)
{
	std::array<int, 2> compound_ids = getTheTwoDifferentIds<N>(particle_global_ids, atoms);

	// Check if the compound_id is already associated with bridges
	auto it = compoundToBridges.find(compound_ids[0]);
	if (it != compoundToBridges.end()) {
		for (BridgeFactory* bridge : it->second) {
			if (bridgeContainsTheseTwoCompounds(*bridge, compound_ids)) {
				return *bridge;
			}
		}
	}


	// If not found in the lookup or no matching bridge, search through all bridges
	for (BridgeFactory& bridge : bridges) {
		if (bridgeContainsTheseTwoCompounds(bridge, compound_ids)) {
			// Cache the bridge reference for the first compound_id
			compoundToBridges[compound_ids[0]].push_back(&bridge);
			return bridge;
		}
	}

	throw std::runtime_error(std::format("Failed to find the bridge ({})", N).c_str());
}

template <typename BondTypeFactory>
void DistributeBondsToCompoundsAndBridges(const std::vector<BondTypeFactory>& bonds, std::vector<ParticleInfo>& particleinfotable, std::vector<CompoundFactory>& compounds, std::vector<BridgeFactory>& bridges,
	BondedParticlesLUTManager& bpLutManager, std::unordered_map<int, std::vector<BridgeFactory*>>& compoundToBridges) {
	
	constexpr int atoms_in_bond = BondTypeFactory::nAtoms;

	for (const BondTypeFactory& bond : bonds) {
		
		DistributeLJIgnores<atoms_in_bond>(&bpLutManager, particleinfotable, bond.global_atom_indexes);

		// TODO Remove this once we switch to a graph based LJignores search
		if (bond.params.HasZeroParam())
			continue;
		//


		if (SpansTwoCompounds<atoms_in_bond>(bond.global_atom_indexes, particleinfotable)) {
			BridgeFactory& bridge = getBridge<atoms_in_bond>(bridges, bond.global_atom_indexes, particleinfotable, compoundToBridges);
			bridge.AddBond(particleinfotable, bond);
		}
		else {
			const int compound_id = particleinfotable[bond.global_atom_indexes[0]].compoundId;	// Pick compound using first particle in bond
			CompoundFactory& compound = compounds.at(compound_id);
			compound.AddBond(particleinfotable, bond);
		}

	}
}

void DistributeBondsToCompoundsAndBridges(const Topology& topology, float boxlen_nm, std::vector<ParticleInfo>& atoms, BoundaryConditionSelect bc_select, 
	std::vector<CompoundFactory>& compounds, std::vector<BridgeFactory>& bridges, BondedParticlesLUTManager& bpLutManager)
{
	std::unordered_map<int, std::vector<BridgeFactory*>> compoundToBridges{};

	DistributeBondsToCompoundsAndBridges(topology.singlebonds, atoms, compounds, bridges, bpLutManager, compoundToBridges);
	//m_logger->print(std::format("Added {} singlebonds to molecule\n", topology.singlebonds.size()));
	DistributeBondsToCompoundsAndBridges(topology.anglebonds, atoms, compounds, bridges, bpLutManager, compoundToBridges);
	//m_logger->print(std::format("Added {} anglebonds to molecule\n", topology.anglebonds.size()));
	DistributeBondsToCompoundsAndBridges(topology.dihedralbonds, atoms, compounds, bridges, bpLutManager, compoundToBridges);
	//m_logger->print(std::format("Added {} dihedralbonds to molecule\n", topology.dihedralbonds.size()));
	DistributeBondsToCompoundsAndBridges(topology.improperdihedralbonds, atoms, compounds, bridges, bpLutManager, compoundToBridges);
	//m_logger->print(std::format("Added {} improper dihedralbonds to molecule\n", topology.improperdihedralbonds.size()));


	// While are particle is never bonded to itself, we are not allowed to calc LJ with itself
	// so we can doubly use this lut to avoid doing that

	for (int com_id = 0; com_id < compounds.size(); com_id++) {
		auto* lut = bpLutManager.get(com_id, com_id);
		for (int pid = 0; pid < MAX_COMPOUND_PARTICLES; pid++) {
			lut->set(pid, pid, true);
		}
	}

	// Finally tell each compound which other compounds they are bonded to for faster lookup of LUTs
	for (const auto& bridge : bridges) {
		for (int i = 0; i < bridge.n_compounds; i++) {
			for (int ii = i + 1; ii < bridge.n_compounds; ii++) {
				if (bridge.compound_ids[i] == bridge.compound_ids[ii]) {
					throw std::runtime_error("A bridge contains the same compound twice");
				}
				compounds[bridge.compound_ids[i]].addIdOfBondedCompound(bridge.compound_ids[ii]);
				compounds[bridge.compound_ids[ii]].addIdOfBondedCompound(bridge.compound_ids[i]);
			}
		}
	}
}

void CalcCompoundMetaInfo(float boxlen_nm, std::vector<CompoundFactory>& compounds, BoundaryConditionSelect bc_select) {
	for (CompoundFactory& compound : compounds) {

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
	}

	// Calc absolute ids of particles in compounds 
	if (compounds.empty()) { return; }
	compounds[0].absoluteIndexOfFirstParticle = 0;
	for (int i = 1; i < compounds.size(); i++) {
		compounds[i].absoluteIndexOfFirstParticle = compounds[i - 1].absoluteIndexOfFirstParticle + compounds[i - 1].n_particles;
	}
}


void VerifyBondsAreStable(const std::vector<SingleBondFactory>& singlebonds, const std::vector<ParticleInfo> atoms, float boxlen_nm, BoundaryConditionSelect bc_select, bool energyMinimizationMode) {
	// First check that we dont have any unrealistic bonds, and warn immediately.
	for (const auto& bond : singlebonds) {
		auto atom1 = atoms[bond.global_atom_indexes[0]];
		auto atom2 = atoms[bond.global_atom_indexes[1]];

		const Float3 pos1 = atom1.groAtom->position;
		const Float3 pos2 = atom2.groAtom->position;
		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos1, pos2, boxlen_nm, bc_select);
		const float bondRelaxedDist = bond.params.b0 * LIMA_TO_NANO;
		const float allowedScalar = energyMinimizationMode ? 7.f : 1.8f;

		if (hyper_dist > bondRelaxedDist * allowedScalar) {
			throw std::runtime_error(std::format("Loading singlebond with illegally large dist ({}). b0: {}", hyper_dist, bond.params.b0 * LIMA_TO_NANO).c_str());
		}
	}
}


std::unique_ptr<BoxImage> LIMA_MOLECULEBUILD::buildMolecules(
	const GroFile& gro_file,
	const TopologyFile& topol_file,
	VerbosityLevel vl,
	std::unique_ptr<LimaLogger> logger,
	bool ignore_hydrogens,
	const SimParams& simparams
) 
{
	// Solvents are not present in top file, so we can't require these to match
	const int nNonsolventAtoms = topol_file.GetElementCount<TopologyFile::AtomsEntry>();
	assert(gro_file.atoms.size() >= nNonsolventAtoms);
	const std::vector<TopologyFileRef> preparedTopologyFiles = PrepareTopologies(topol_file);

	ForcefieldManager forcefieldManager{};

	std::vector<ParticleInfo> preparedAtoms = PrepareAtoms(preparedTopologyFiles, gro_file, forcefieldManager, nNonsolventAtoms);

	Topology topology = LoadTopology(preparedTopologyFiles, forcefieldManager, preparedAtoms);

	VerifyBondsAreStable(topology.singlebonds, preparedAtoms, gro_file.box_size.x, simparams.bc_select, simparams.em_variant);

	std::vector<std::vector<int>> atomIdToSinglebondsMap = MapAtomsToSinglebonds(preparedAtoms, topology);

	std::vector<AtomGroup> atomGroups = GroupAtoms(preparedAtoms, topology);


	std::vector<CompoundFactory> compounds = CreateCompounds(topology, gro_file.box_size.x, atomGroups, preparedAtoms, atomIdToSinglebondsMap, simparams.bc_select);
	std::vector<BridgeFactory> bridges = CreateBridges(topology.singlebonds, compounds, preparedAtoms);
	auto bpLutManager = std::make_unique<BondedParticlesLUTManager>();
	DistributeBondsToCompoundsAndBridges(topology, gro_file.box_size.x, preparedAtoms, simparams.bc_select, compounds, bridges, *bpLutManager);

	//bpLutManager->get(0, 0)->printMatrix(compounds.begin()->n_particles);

	CalcCompoundMetaInfo(gro_file.box_size.x, compounds, simparams.bc_select);


	auto bridges_compact = std::make_unique<CompoundBridgeBundleCompact>(
		std::vector<CompoundBridge>(bridges.begin(), bridges.end())
	);


	const std::vector<Float3> solventPositions = LoadSolventPositions(gro_file);

	auto boxImage = std::make_unique<BoxImage>(
		std::move(compounds),
		static_cast<int>(preparedAtoms.size()),
		std::move(bpLutManager),
		std::move(bridges_compact),
		std::move(solventPositions),
		gro_file.box_size.x,	// TODO: Find a better way..
		std::move(preparedAtoms),
		GroFile{ gro_file },	// TODO: wierd ass copy here
		forcefieldManager.GetActiveLjParameters(), 
		std::move(topology)
	);

	return boxImage;
}
