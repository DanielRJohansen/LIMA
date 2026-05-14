#include "CompoundBuilder.h"
#include "Forcefield.h"
#include "MoleculeGraph.h"
#include "TimeIt.h"

#include <unordered_set>
#include <format>
#include <array>
#include <numeric>
#include <set>

//#include "Display.h"
using namespace LIMA_MOLECULEBUILD;
using namespace LimaMoleculeGraph;




template <typename BondType, typename BondtypeFactory, typename BondTypeTopologyfile>
void SuperTopology::LoadBondsIntoTopology(const std::vector<BondTypeTopologyfile>& bondsInTopfile, int atomIdOffset, LIMAForcefield& forcefield,
	std::vector<BondtypeFactory>& topology)
{
	for (const auto& bondTopol : bondsInTopfile) {
		std::array<int, BondType::nAtoms> globalIds;
		std::array<std::string, BondType::nAtoms> atomTypenames;

		bool bondExists = true;
		for (int i = 0; i < BondType::nAtoms; ++i) {
			if (bondTopol.ids[i] + atomIdOffset >= particles.size()) {
				bondExists = false;
				break;
			}
		}
		if (!bondExists)
			continue;


		for (int i = 0; i < BondType::nAtoms; ++i) {
			globalIds[i] = bondTopol.ids[i] + atomIdOffset;
			atomTypenames[i] = particles[globalIds[i]].topologyAtom.type;
		}


		// Solvent's bonds are defined in the forcefield, rather the params are directly in the topology... Not sure how to deal with that rn
		//const bool getParamsFromForcefield = particles[globalIds[0]].topologyAtom.residue != "SOL" && particles[globalIds[0]].topologyAtom.residue != "TIP3";
		/*const bool getParamsFromForcefield = true;
		const bool getParamsFromForcefield1 = bondTopol.parameters.has_value();*/

		// In rare cases, the bond parameters are directly in the topology file
		if (bondTopol.parameters.has_value()) {
			topology.emplace_back(BondtypeFactory{ globalIds, bondTopol.parameters.value() });
		}
		else {
			// A bond may be described as multiple bonds, so this is a vector
			const std::vector<typename BondType::Parameters>& bondParams = forcefield.GetBondParameters<BondType>(atomTypenames);

			for (const auto& param : bondParams) {
				topology.emplace_back(BondtypeFactory{ globalIds, param });
			}
		}
	}
}

SuperTopology::SuperTopology(const TopologyFile::System& system, const GroFile& grofile, LIMAForcefield& forcefield) {

	int nextUniqueParticleId = 0;
	int indexInGrofile = 0;



	for (int topologyMoleculeIndex = 0; topologyMoleculeIndex < system.molecules.size(); topologyMoleculeIndex++) {
		const TopologyFile::MoleculeEntry& molecule = system.molecules[topologyMoleculeIndex];

#if ENABLE_SOLVENTS != 1
		if (molecule.name == "SOL" || molecule.name == "TIP3") {// TODO: Add the other Solvent labels
			continue;
		}
#endif 

		const int particleIdOffset = nextUniqueParticleId;
		const TopologyFile::Moleculetype& molType = *molecule.moleculetype;

		if (molType.atoms.empty())
			throw std::runtime_error("Molecule has no atoms");

		for (int localId = 0; localId < molType.atoms.size(); localId++) {

			// Here's we fetch an LJ param, but we dont yet know if this particle is in a tinyMol. So this is a waste...
			const int activeLJParamIndex = forcefield.GetActiveLjParameterIndex(molType.atoms[localId].type);

			particles.push_back(ParticleFactory{ molType.atoms[localId], grofile.atoms[indexInGrofile].position, indexInGrofile, activeLJParamIndex });
			nextUniqueParticleId++;
			indexInGrofile++;
		}

		LoadBondsIntoTopology<SingleBond, SingleBondFactory, TopologyFile::SingleBond>(molType.singlebonds, particleIdOffset, forcefield, singlebonds);
		LoadBondsIntoTopology<PairBond, PairBondFactory, TopologyFile::PairBond>(molType.pairbonds, particleIdOffset, forcefield, pairbonds);
		LoadBondsIntoTopology<AngleUreyBradleyBond, AngleBondFactory, TopologyFile::AngleBond>(molType.anglebonds, particleIdOffset, forcefield, anglebonds);
		LoadBondsIntoTopology<DihedralBond, DihedralBondFactory, TopologyFile::DihedralBond>(molType.dihedralbonds, particleIdOffset, forcefield, dihedralbonds);
		LoadBondsIntoTopology<ImproperDihedralBond, ImproperDihedralBondFactory, TopologyFile::ImproperDihedralBond>(molType.improperdihedralbonds, particleIdOffset, forcefield, improperdihedralbonds);
	}
}

void SuperTopology::VerifyBondsAreStable(const Float3& boxlen_nm, BoundaryConditionSelect bc_select, bool energyMinimizationMode) const {
	const float allowedScalar = 1.9f;

	for (const auto& bond : singlebonds)
	{
		const Float3 pos1 = particles[bond.global_atom_indexes[0]].position;
		const Float3 pos2 = particles[bond.global_atom_indexes[1]].position;
		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos1, pos2, boxlen_nm, bc_select);
		const float bondRelaxedDist = bond.params.b0;

		if (hyper_dist > bondRelaxedDist * allowedScalar) {
			/*throw std::runtime_error(std::format("Loading singlebond with illegally large dist ({}). b0: {}. AtomIndices: {} {}",
				hyper_dist, bond.params.b0, bond.global_atom_indexes[0], bond.global_atom_indexes[1]));*/
		}
		if (hyper_dist < bondRelaxedDist * 0.001) {
			/*throw std::runtime_error(std::format("Loading singlebond with illegally small dist ({}). b0: {}. AtomIndices: {} {}",
				hyper_dist, bond.params.b0, bond.global_atom_indexes[0], bond.global_atom_indexes[1]));*/
		}
	}
	for (const auto& bond : anglebonds)
	{
		const Float3 pos1 = particles[bond.global_atom_indexes[0]].position;
		const Float3 pos2 = particles[bond.global_atom_indexes[1]].position;
		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos1, pos2, boxlen_nm, bc_select);
		if (hyper_dist < 0.001)
			throw std::runtime_error(std::format("Loading singlebond with illegally small dist ({}). b0: {}", hyper_dist, bond.params.ub0));
	}
}



// --------------------------------------------------------------- Factory Functions --------------------------------------------------------------- //


struct ParticleBondedToParticlesLookup {
	ParticleBondedToParticlesLookup(const SuperTopology& system) {
		particleBondedToParticle.resize(system.particles.size());



		// Then go through all bonds to make the particle nointeraction matrix
		for (const auto& singlebond : system.singlebonds)
			AddBond(singlebond.global_atom_indexes);
		for (const auto& anglebond : system.anglebonds)
			AddBond(anglebond.global_atom_indexes);
		for (const auto& dihedralbond : system.dihedralbonds)
			AddBond(dihedralbond.global_atom_indexes);
		for (const auto& improperdihedralbond : system.improperdihedralbonds)
			AddBond(improperdihedralbond.global_atom_indexes);
	}

	bool AreAllBonded(std::span<int>ids) const {
		std::vector<bool> bonded(ids.size(), false);

		for (int i = 0; i < ids.size(); i++) {
			const int pid_self = ids[i];
			for (int j = i + 1; j < ids.size(); j++) {
				const int pid_other = ids[j];

				if (particleBondedToParticle[pid_self].find(pid_other) != particleBondedToParticle[pid_self].end()) {
					bonded[i] = true;
					bonded[j] = true;
					break;
				}
			}
		}
		return std::any_of(bonded.begin(), bonded.end(), [](bool v) { return v; });
	}

	std::vector<bool> BondedToFirst(std::span<int> ids) const {
		std::vector<bool> bonded(ids.size(), false);
		bonded[0] = true;
		const int pid_self = ids[0];
		for (int j = 1; j < ids.size(); j++) {
			const int pid_other = ids[j];
			if (particleBondedToParticle[pid_self].find(pid_other) != particleBondedToParticle[pid_self].end()) {
				bonded[j] = true;
			}
		}
		return bonded;
	}

private:

	void AddBond(std::span<const int> particleIdsInBond) {
		for (int i = 0; i < particleIdsInBond.size(); i++) {
			const int pid_self = particleIdsInBond[i];

			for (int j = i + 1; j < particleIdsInBond.size(); j++) {
				const int pid_other = particleIdsInBond[j];

				particleBondedToParticle[pid_self].insert(pid_other);
				particleBondedToParticle[pid_other].insert(pid_self);
			}
		}
	}
	std::vector<std::set<int>> particleBondedToParticle;
};

void SplitClusters(std::span<int> ids, const ParticleBondedToParticlesLookup& particleBondedToParticlesLookup, std::vector<std::array<int, 4>>& outClusters) 
{
	std::vector<bool> bondedToFirst = particleBondedToParticlesLookup.BondedToFirst(ids);

	std::array<int, 4> thisCluster{ -1, -1, -1, -1 };
	std::vector<int> remainingIds;

	for (int i = 0; i < ids.size(); i++) {
		if (bondedToFirst[i])
			thisCluster[i] = ids[i];
		else
			remainingIds.push_back(ids[i]);
	}

	outClusters.push_back(thisCluster);

	if (!remainingIds.empty())
		SplitClusters(std::span<int>(remainingIds), particleBondedToParticlesLookup, outClusters);
}

std::shared_ptr<MoleculeGraph> MakeMoleculeGraph(const SuperTopology& system) {
	std::vector<std::pair<int, std::string>> atoms;
	atoms.reserve(system.particles.size());
	for (int pid = 0; pid < system.particles.size(); pid++) {
		atoms.push_back({ pid, system.particles[pid].topologyAtom.type });
	}
	std::vector<std::array<int, 2>> edges;
	edges.reserve(system.singlebonds.size());
	for (const auto& bond : system.singlebonds) {
		edges.push_back(bond.global_atom_indexes);
	}

	// TODO: This is under the assumption that we get a ideally sorted graph back, ill need to verify that
	auto systemGraph = std::make_shared<MoleculeGraph>(atoms, edges);
	return systemGraph;
}

std::vector<std::array<int, 4>> SplitIntoPersistentClusters(const SuperTopology& system, const MoleculeGraph& systemGraph, const ParticleBondedToParticlesLookup& particleBondedToParticlesLookup, Float3 box_size) {

	const std::vector<std::vector<int>> particleidCollectionsOfMolecules = systemGraph.GetListOfListsofConnectedNodeids();

	std::vector<std::array<int, 4>> persistentClusters;
	persistentClusters.reserve(system.particles.size()); // A bit too big..


	
	auto StoreCurrentCluster = [&](std::array<int, 4>& cluster, int& nextIndex) {
		persistentClusters.push_back(cluster);
		cluster = { -1,-1,-1,-1 };
		nextIndex = 0;
		};

	auto CanAppendToCluster = [&systemGraph](int particleId, std::array<int, 4>& cluster, int nextIndex) {
		if (nextIndex == 0)
			return true;

		std::optional<int> distanceToPreviousNode = systemGraph.DistanceBetweenNodes(cluster[nextIndex - 1], particleId, 5);
		std::optional<int> distanceToFirstNode = systemGraph.DistanceBetweenNodes(cluster[0], particleId, 5);

		const bool canAppend = distanceToFirstNode.value_or(INT_MAX) < 3 ||
			distanceToFirstNode.value_or(INT_MAX) <= 4 && distanceToPreviousNode.value_or(INT_MAX) <= 2;

		return canAppend;
		};


	


	for (const std::vector<int>& collection : particleidCollectionsOfMolecules) {

		// Quick handling of small molecules (water/ions)
		if (collection.size() <= 4) {
			auto& cluster = persistentClusters.emplace_back(std::array{ -1, -1, -1, -1 });
			std::ranges::copy(collection, cluster.begin());
			continue;
		}
		


		// Complex compression of larger molecules
		std::array<int, 4> cluster{ -1,-1,-1,-1 };
		std::unordered_set<int> addedByLookahead;
		addedByLookahead.clear();
		int nextIndex = 0;
		for (int i = 0; i < collection.size(); i++) {		
			const int particleId = collection[i];
			if (addedByLookahead.contains(particleId))
				continue;

			if (nextIndex == 0) {}
			else {
				const bool canAppend = CanAppendToCluster(particleId, cluster, nextIndex);

				if (!canAppend) {
					// Look ahead and add other particles if possible
					const int lookaheadCnt = 6;
					for (int lookaheadIndex = i + 1; (lookaheadIndex <= std::min(i + lookaheadCnt, (int)collection.size() - 2)) && nextIndex < 4; lookaheadIndex++) {
						const int lookaheadId = collection[lookaheadIndex];
						if (CanAppendToCluster(lookaheadId, cluster, nextIndex)) {
							addedByLookahead.insert(lookaheadId);
							cluster[nextIndex++] = lookaheadId;
						}

					}

					StoreCurrentCluster(cluster, nextIndex);
				}
			}

			cluster[nextIndex++] = particleId;

			// If full or finished, push
			if (i == collection.size() - 1 || nextIndex == 4) {
				StoreCurrentCluster(cluster, nextIndex);
			}
		}
	}


	// Compute cluster vacancy
	int vacantCount = 0;
	for (auto& cluster : persistentClusters) {
		for (int i = 0; i < PersistentCluster::maxParticles; i++) {
			if (cluster[i] == -1)
				vacantCount++;
		}
	}
	double vacancyFraction = static_cast<double>(vacantCount) / (double)(persistentClusters.size() * PersistentCluster::maxParticles);

	float largestDistInsidePcluster = 0.f;
	for (auto& cluster : persistentClusters) {
		std::vector<Float3> positions;
		for (int i = 0; i < PersistentCluster::maxParticles; i++) {
			if (cluster[i] != -1) {
				positions.push_back(system.particles[cluster[i]].position);
			}
		}
		for (int i = 0; i < positions.size(); i++) {
			for (int j = i + 1; j < positions.size(); j++) {
				const float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(positions[i], positions[j],  box_size, BoundaryConditionSelect::PBC);
				if (dist > largestDistInsidePcluster)
					largestDistInsidePcluster = dist;
			}
		}
	}
	// Next compute the largest distances inside clusters
	
	//either the pclusters are made with particles far from eachother??


	return persistentClusters;
}
std::tuple<std::vector<PersistentCluster>, std::vector<PersistentClusterMeta>, ParticleToPclusterMap> MakePersistentClusters(const std::vector<std::array<int, 4>>& clustersParticleIds, const SuperTopology& system, LIMAForcefield& forcefield) {

	std::vector<PersistentCluster> pClusters(clustersParticleIds.size());
	std::vector<PersistentClusterMeta> pClusterMetas(clustersParticleIds.size());
	ParticleToPclusterMap particleToPclusterMap(system.particles.size());

	std::vector<std::string> solventResNames{ "SOL", "SPC", "SPCE", "TIP3", "TIP3P" };

	for (int pcId = 0; pcId < clustersParticleIds.size(); pcId++) {
		for (int pidRel = 0; pidRel < PersistentCluster::maxParticles; pidRel++) {
			const int pId = clustersParticleIds[pcId][pidRel];			

			if (pId == -1) {
				pClusters[pcId].pqd[pidRel] = PData{};
				pClusterMetas[pcId].particleIdsGlobal[pidRel] = -1;
				continue;
			}
			else {
				const auto topAtom = system.particles[pId].topologyAtom;
				const std::string& atomType = topAtom.type;
				const Float3 pos = system.particles[pId].position;
				NBParams nbParams = forcefield.GetLjParameters(atomType);
				if (topAtom.charge.has_value())
					nbParams.charge = topAtom.charge.value() * elementaryChargeToKiloCoulombPerMole;

				pClusters[pcId].pqd[pidRel] = PData{ pos,  nbParams };
				pClusterMetas[pcId].particleIdsGlobal[pidRel] = pId;
				if (topAtom.mass.has_value())
					pClusterMetas[pcId].mass[pidRel] = topAtom.mass.value() / KILO;	// TODO: I dont like this conversion here. Actually we should get the mass from the forcefield, which already does the conversion??
				else if (forcefield.GetAtomtype(atomType).has_value())
					pClusterMetas[pcId].mass[pidRel] = forcefield.GetAtomtype(atomType)->mass;
				pClusterMetas[pcId].atomLetter[pidRel] = !topAtom.atomname.empty() ? topAtom.atomname[0] : ' ';
				assert(pClusterMetas[pcId].mass[pidRel] > 0.f );

				pClusterMetas[pcId].isSolvent = std::find(solventResNames.begin(), solventResNames.end(), topAtom.residue) != solventResNames.end();
				pClusterMetas[pcId].nParticles++;

				// Also set mapping
				particleToPclusterMap[pId] = ParticleToPclusterMapping{ pcId, pidRel };
			}
		}
	}

	return { pClusters, pClusterMetas, particleToPclusterMap };
}

// returns bondedParticles, bondedPclusters
std::pair<std::vector<std::set<int>>, std::vector<std::set<int>>> GetBondedPersistentClusters(const std::vector<std::array<int, 4>>& clustersParticleIds, const SuperTopology& system) {
	// First make a particle-2-pcluster map	
	std::vector<int> particleIdToPclusterIdMap(system.particles.size(), -1);
	for (int pcId = 0; pcId < clustersParticleIds.size(); pcId++) {
		for (int pidRel = 0; pidRel < PersistentCluster::maxParticles; pidRel++) {
			const int pId = clustersParticleIds[pcId][pidRel];
			if (pId == -1) { continue; }
			particleIdToPclusterIdMap[pId] = pcId;
		}
	}

	std::vector<std::set<int>> particleBondedToParticle(system.particles.size());
	std::vector<std::set<int>> pclusterBondedToPcluster(clustersParticleIds.size());



	auto AddBond = [&](const auto& particleIdsInBond) {
		for (int i = 0; i < particleIdsInBond.size(); i++) {
			const int pid_self = particleIdsInBond[i];
			const int pcid_self = particleIdToPclusterIdMap[pid_self];
			//if (pcid_self == -1) { continue; }
			assert(pcid_self != -1);

			for (int j = i + 1; j < particleIdsInBond.size(); j++) {
				const int pid_other = particleIdsInBond[j];
				const int pcid_other = particleIdToPclusterIdMap[pid_other];
				assert(pcid_other != -1);

				particleBondedToParticle[pid_self].insert(pid_other);
				particleBondedToParticle[pid_other].insert(pid_self);
				pclusterBondedToPcluster[pcid_self].insert(pcid_other);
				pclusterBondedToPcluster[pcid_other].insert(pcid_self);
			}
		}
		};

	// Then go through all bonds to make the particle nointeraction matrix
	for (const auto& singlebond : system.singlebonds)
		AddBond(singlebond.global_atom_indexes);
	for (const auto& anglebond : system.anglebonds)
		AddBond(anglebond.global_atom_indexes);
	for (const auto& dihedralbond : system.dihedralbonds)
		AddBond(dihedralbond.global_atom_indexes);
	for (const auto& improperdihedralbond : system.improperdihedralbonds)
		AddBond(improperdihedralbond.global_atom_indexes);

	return { std::move(particleBondedToParticle), std::move(pclusterBondedToPcluster )};
}



template <typename BondType, typename BondtypeFactory, typename BondTypeTopologyfile>
void LoadBondIntoTopology(const std::vector<BondTypeTopologyfile>& bondsInTopfile,	int atomIdOffset, LIMAForcefield& forcefield,
	const std::vector<ParticleFactory>& particles, std::vector<BondtypeFactory>& topology)
{
	for (const auto& bondTopol : bondsInTopfile) {
		std::array<int, BondType::nAtoms> globalIds;
		std::array<std::string, BondType::nAtoms> atomTypenames;

		for (int i = 0; i < BondType::nAtoms; ++i) {
			globalIds[i] = bondTopol.ids[i] + atomIdOffset;
			atomTypenames[i] = particles[globalIds[i]].topAtom->type;
		}

		auto bondParams = forcefield.GetBondParameters<BondType>(atomTypenames);

		for (const auto& param : bondParams) {
			topology.emplace_back(BondtypeFactory{ globalIds, param });
		}
	}
}



// Usually just a residue, but we sometimes also need to split lipids into smaller groups 
struct AtomGroup {
	std::vector<int> atomIds;
	std::set<int> idsOfBondedAtomgroups;

	bool operator != (const AtomGroup& other) const {
		return atomIds != other.atomIds || idsOfBondedAtomgroups != other.idsOfBondedAtomgroups;
	}
};


bool AreBonded(const AtomGroup& left, const AtomGroup& right, const std::vector<std::unordered_set<int>>& atomIdToSinglebondsMap) {
	for (auto& atomleft_gid : left.atomIds) {
		for (auto& atomright_gid : right.atomIds) {
			const std::unordered_set<int>& atomleft_singlesbonds = atomIdToSinglebondsMap[atomleft_gid];
			const std::unordered_set<int>& atomright_singlesbonds = atomIdToSinglebondsMap[atomright_gid];

			for (int singlebondId : atomleft_singlesbonds) {
				if (atomright_singlesbonds.contains(singlebondId)) {
					return true;
				}
			}
		}
	}
	return false;
}



std::vector<int> ReorderSubchains(const std::vector<int>& ids, const std::unordered_map<int,int>& nodeIdToNumDownstream, int spaceLeft) {
	std::vector<int> bestOrder = ids;
	int maxElements = 0;

	// Start with the initial permutation of ids
	std::vector<int> currentOrder = ids;

	do {
		int currentElements = 0;

		// Try to fit as many elements as possible in the current permutation
		for (int id : currentOrder) {
			if (currentElements + nodeIdToNumDownstream.at(id) <= spaceLeft)
				currentElements += nodeIdToNumDownstream.at(id);
		}

		// Update the best order if this permutation fits more elements
		if (currentElements > maxElements) {
			maxElements = currentElements;
			bestOrder = currentOrder;
		}
	} while (std::next_permutation(currentOrder.begin(), currentOrder.end()));

	// Update the original ids to reflect the best order found
	return bestOrder;
}



//const std::vector<AtomGroup> GroupAtoms(const std::vector<std::vector<int>>& particleidsInMolecules, const SuperTopology& topology) {
//	std::vector<AtomGroup> atomGroups;
//
//
//	std::vector<std::unordered_set<int>> pidToSinglebondidMap(topology.particles.size());
//	for (int bid = 0; bid < topology.singlebonds.size(); bid++) {
//		pidToSinglebondidMap[topology.singlebonds[bid].global_atom_indexes[0]].insert(bid);
//		pidToSinglebondidMap[topology.singlebonds[bid].global_atom_indexes[1]].insert(bid);
//	}
//
//
//	for (const auto& particleIdsInMolecule : particleidsInMolecules) {
//
//		std::vector<std::pair<int, std::string>> atoms;
//		atoms.reserve(particleIdsInMolecule.size());
//		for (int pid : particleIdsInMolecule) {
//			atoms.emplace_back( pid, topology.particles[pid].topologyAtom.type );
//		}
//
//		std::unordered_set<int> bondIdsInMolecule;
//		for (int pid : particleIdsInMolecule) {
//			for (int bid : pidToSinglebondidMap[pid]) {
//				bondIdsInMolecule.insert(bid);
//			}
//		}
//
//		std::vector<std::array<int, 2>> edges;
//		edges.reserve(bondIdsInMolecule.size());
//		for (int bid : bondIdsInMolecule) {
//			edges.emplace_back(topology.singlebonds[bid].global_atom_indexes);
//		}
//
//
//
//
//		const MoleculeGraph molGraph(atoms, edges);
//		const MoleculeTree moleculeTree = molGraph.ConstructMoleculeTree();
//		const std::unordered_map<int, int> nodeIdNumDownstreamNodes = molGraph.ComputeNumDownstreamNodes(moleculeTree);
//
//		std::stack<const MoleculeGraph::Node*> nodeStack;
//		nodeStack.push(molGraph.root);
//
//		atomGroups.emplace_back();
//
//		while (!nodeStack.empty()) {
//			const MoleculeGraph::Node* node = nodeStack.top();
//			nodeStack.pop();
//
//			if (MAX_COMPOUND_PARTICLES - atomGroups.back().atomIds.size() == 0)
//				atomGroups.emplace_back();
//			atomGroups.back().atomIds.emplace_back(node->atomid);
//
//			std::vector<int> nodeChildren = moleculeTree.GetChildIds(node->atomid);
//
//			if (nodeChildren.empty()) {
//				// finished
//			}
//			else {
//				// Add the longest childchain to our stack, and remove it from the current children
//				const int indexOfLongestChain = std::max_element(nodeChildren.begin(), nodeChildren.end(),
//					[&nodeIdNumDownstreamNodes](const int& a, const int& b) { return nodeIdNumDownstreamNodes.at(a) < nodeIdNumDownstreamNodes.at(b); }
//				) - nodeChildren.begin();
//				nodeStack.push(&molGraph.nodes.at(nodeChildren[indexOfLongestChain]));
//				nodeChildren[indexOfLongestChain] = nodeChildren.back();
//				nodeChildren.pop_back();
//
//				const std::vector<int> nodeChildrenIdsIdealOrder = ReorderSubchains(nodeChildren, nodeIdNumDownstreamNodes, MAX_COMPOUND_PARTICLES - atomGroups.back().atomIds.size());
//
//				for (int id : nodeChildrenIdsIdealOrder) {
//					if (nodeIdNumDownstreamNodes.at(id) > MAX_COMPOUND_PARTICLES - atomGroups.back().atomIds.size())
//						atomGroups.emplace_back();
//
//					moleculeTree.ForSelfAndAllChildrenIds(id,
//						[&atomGroups](int _id) { atomGroups.back().atomIds.emplace_back(_id); }
//					);
//				}
//			}
//		}
//	}
//
//	// Now figure out which groups are bonded to others
//	for (int i = 1; i < atomGroups.size(); i++) {
//		if (AreBonded(atomGroups[i - 1], atomGroups[i], pidToSinglebondidMap)) {
//			atomGroups[i].idsOfBondedAtomgroups.insert(i - 1);
//			atomGroups[i - 1].idsOfBondedAtomgroups.insert(i);
//		}
//	}
//
//	return atomGroups;
//}






//void VerifyBondsAreLegal(const std::vector<AngleBondFactory>& anglebonds, const std::vector<DihedralBondFactory>& dihedralbonds) {
//	for (const auto& bond : anglebonds) {
//		if (bond.global_atom_indexes[0] == bond.global_atom_indexes[1] || bond.global_atom_indexes[1] == bond.global_atom_indexes[2] || bond.global_atom_indexes[0] == bond.global_atom_indexes[2])
//			throw ("Bond contains the same index multiple times: %d %d %d");
//	}
//	for (const auto& bond : dihedralbonds) {
//		std::unordered_set<int> ids;
//		for (int i = 0; i < 4; i++) {
//			if (ids.contains(bond.global_atom_indexes[i]))
//				throw ("Dihedral bond contains the same index multiple times");
//			ids.insert(bond.global_atom_indexes[i]);
//		}
//	}
//}


std::unique_ptr<BoxImage> LIMA_MOLECULEBUILD::buildMolecules(
	const GroFile& grofile,
	const TopologyFile& topol_file,
	VerbosityLevel vl,
	std::unique_ptr<LimaLogger> logger,
	bool ignore_hydrogens,
	const SimParams& simparams
)
{
	TimeIt timer("Build Molecules", true);
	LIMAForcefield forcefield{ topol_file.forcefieldInclude ? topol_file.forcefieldInclude->contents : GenericItpFile{} };

	SuperTopology superTopology(topol_file.GetSystem(), grofile, forcefield);
	superTopology.VerifyBondsAreStable(grofile.box_size, simparams.bc_select, simparams.em_variant);




	const ParticleBondedToParticlesLookup particleBondedToParticlesLookup(superTopology);

	// Make PersistenClusters
	std::shared_ptr<MoleculeGraph> systemGraph = MakeMoleculeGraph(superTopology);
	std::vector<std::array<int,4>> pClustersParticleids = SplitIntoPersistentClusters(superTopology, *systemGraph, particleBondedToParticlesLookup, grofile.box_size);

	auto [pClusters, pClusterMetas, particleToPclusterMap] = MakePersistentClusters(pClustersParticleids, superTopology, forcefield);

	auto [particleBondedToParticle, pclusterBondedToPcluster] = GetBondedPersistentClusters(pClustersParticleids, superTopology);



	std::vector<BondGroupFactory> bondGroups = BondGroupFactory::MakeBondgroups(superTopology, particleToPclusterMap, pClusters.data());
	const auto particleToBondgroupMap = BondGroupFactory::MakeParticleToBondgroupsMap(bondGroups, superTopology.particles.size());

	//std::swap(bondGroups.front().singlebonds[0], bondGroups.front().singlebonds[2]);
	//bondGroups.front().nSinglebonds = 2;
	/*bondGroups.front().nAnglebonds = 0;
	bondGroups.front().nDihedralbonds = 0;*/

	//{
	//	std::vector<Float3> bgPositions;
	//	for (int pid = 0; pid < bondGroups.front().nParticles; pid++) {
	//		auto pref = bondGroups.front().particles[pid];
	//		bgPositions.push_back(pClusters[pref.pcid].pqd[pref.pid].position);
	//	}

	//	for (int i = 0; i < bondGroups.front().nSinglebonds; i++) {
	//		int id0 = bondGroups.front().singlebonds[i].idInBondgroup[0];
	//		int id1 = bondGroups.front().singlebonds[i].idInBondgroup[1];
	//		Float3 p0 = bgPositions[id0];
	//		Float3 p1 = bgPositions[id1];
	//		float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(p0, p1, grofile.box_size, simparams.bc_select);
	//		//printf("p0 %d %f %f %f p1 %d %f %f %f dist %f\n", id0, p0.x, p0.y, p0.z, id1, p1.x, p1.y, p1.z, dist);
	//		int a = 0;

	//		int id0Global = pClusterMetas[bondGroups.front().particles[id0].pcid].particleIdsGlobal[bondGroups.front().particles[id0].pid];
	//		int id1Global = pClusterMetas[bondGroups.front().particles[id1].pcid].particleIdsGlobal[bondGroups.front().particles[id1].pid];
	//		int id0Groid = grofile.atoms[superTopology.particles[id0Global].indexInGrofile].gro_id;
	//		int id1Groid = grofile.atoms[superTopology.particles[id1Global].indexInGrofile].gro_id;

	//		printf("id0Global %d id1Global %d id0Groid %d id1Groid %d dist %f\n", id0Global, id1Global, id0Groid, id1Groid, dist);
	//	}
	//}


	for (int i = 0; i < particleToPclusterMap.size(); i++) {
		const auto pcRef = particleToPclusterMap[i];
		const std::set<BondgroupRef>& bgRefs = particleToBondgroupMap[i];

		for (const BondgroupRef& bgRef : bgRefs) {
			pClusterMetas[pcRef.pcid].bondgroupReferences[pcRef.pid].Add(bgRef);
			//compounds[pcRef.compoundId].AddBondgroupReference(pcRef.localIdInCompound, bgRef);
		}
	}



	int nParticles = 0;
	//int nSolvents = 0;
	for (const auto& pc : pClusterMetas) {
		for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
			if (pc.particleIdsGlobal[pid] == -1)
				continue;
			nParticles++;
		}
	}

	std::vector<std::tuple<int, int>> gpidToPcidAndPid(nParticles);
	for (int pcid = 0; pcid < pClusterMetas.size(); pcid++) {
		for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
			int gpid = pClusterMetas[pcid].particleIdsGlobal[pid];
			if (gpid != -1)
				gpidToPcidAndPid[gpid] = { pcid, pid };
		}
	}

	return std::make_unique<BoxImage>(
		grofile,	// TODO: wierd ass copy here. Probably make the input a sharedPtr?
		forcefield.GetActiveLjParameters(),
		superTopology,
		systemGraph,
		forcefield.GetNonbondedInteractionParams(),
		BondGroupFactory::FinishBondgroups(bondGroups),
		pClusters,
		pClusterMetas,
		particleBondedToParticle,
		pclusterBondedToPcluster,
		gpidToPcidAndPid,
		nParticles
	);

}
