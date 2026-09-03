#include "BoxImageBuilder.h"
#include "Forcefield.h"
#include "MoleculeGraph.h"
#include "TimeIt.h"

#include <unordered_set>
#include <unordered_map>
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
	moleculeInstances.reserve(system.molecules.size());



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

		moleculeInstances.push_back(MoleculeInstance{ &molType, particleIdOffset });

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

namespace {
	struct PersistentParticleTemplate {
		NBParams nbParams{};
		float mass = 0.f;
		char atomLetter = ' ';
		bool isSolvent = false;
	};

	struct PersistentClusterTemplate {
		std::vector<std::array<int, PersistentCluster::maxParticles>> clusters;
		std::vector<int> particleToPcluster;
		std::vector<std::set<int>> particleBondedToParticle;
		std::vector<std::set<int>> pclusterBondedToPcluster;
		std::vector<PersistentParticleTemplate> particles;
	};

	bool IsSolventResidue(const std::string& residue) {
		return residue == "SOL" || residue == "SPC" || residue == "SPCE" || residue == "TIP3" || residue == "TIP3P";
	}

	std::vector<std::array<int, PersistentCluster::maxParticles>> MakeLocalPersistentClusters(
		const TopologyFile::Moleculetype& molecule
	) {
		const MoleculeGraph moleculeGraph(molecule);
		const std::vector<std::vector<int>> connectedComponents = moleculeGraph.GetListOfListsofConnectedNodeids();

		std::vector<std::array<int, PersistentCluster::maxParticles>> clusters;
		clusters.reserve((molecule.atoms.size() + PersistentCluster::maxParticles - 1) / PersistentCluster::maxParticles);

		auto CanAppendToCluster = [&moleculeGraph](int particleId, const auto& cluster, int nextIndex) {
			if (nextIndex == 0)
				return true;

			const std::optional<int> distanceToPreviousNode = moleculeGraph.DistanceBetweenNodes(cluster[nextIndex - 1], particleId, 5);
			const std::optional<int> distanceToFirstNode = moleculeGraph.DistanceBetweenNodes(cluster[0], particleId, 5);

			return distanceToFirstNode.value_or(INT_MAX) < 3 ||
				distanceToFirstNode.value_or(INT_MAX) <= 4 && distanceToPreviousNode.value_or(INT_MAX) <= 2;
		};

		std::unordered_set<int> addedByLookahead;
		for (const std::vector<int>& component : connectedComponents) {
			if (component.size() <= PersistentCluster::maxParticles) {
				auto& cluster = clusters.emplace_back(std::array{ -1, -1, -1, -1 });
				std::ranges::copy(component, cluster.begin());
				continue;
			}

			addedByLookahead.clear();
			clusters.emplace_back(std::array{ -1, -1, -1, -1 });
			int nextIndex = 0;

			for (int i = 0; i < component.size(); i++) {
				const int particleId = component[i];
				if (addedByLookahead.contains(particleId))
					continue;

				if (nextIndex != 0 && !CanAppendToCluster(particleId, clusters.back(), nextIndex)) {
					constexpr int lookaheadCount = 6;
					const int lastLookaheadIndex = std::min(i + lookaheadCount, static_cast<int>(component.size()) - 2);
					for (int lookaheadIndex = i + 1; lookaheadIndex <= lastLookaheadIndex && nextIndex < PersistentCluster::maxParticles; lookaheadIndex++) {
						const int lookaheadId = component[lookaheadIndex];
						if (CanAppendToCluster(lookaheadId, clusters.back(), nextIndex)) {
							addedByLookahead.insert(lookaheadId);
							clusters.back()[nextIndex++] = lookaheadId;
						}
					}

					clusters.emplace_back(std::array{ -1, -1, -1, -1 });
					nextIndex = 0;
				}

				clusters.back()[nextIndex++] = particleId;
				if (nextIndex == PersistentCluster::maxParticles && i != component.size() - 1) {
					clusters.emplace_back(std::array{ -1, -1, -1, -1 });
					nextIndex = 0;
				}
			}
		}

		return clusters;
	}

	PersistentClusterTemplate BuildPersistentClusterTemplate(
		const TopologyFile::Moleculetype& molecule,
		LIMAForcefield& forcefield
	) {
		PersistentClusterTemplate result;
		result.clusters = MakeLocalPersistentClusters(molecule);
		result.particleToPcluster.resize(molecule.atoms.size(), -1);
		result.particleBondedToParticle.resize(molecule.atoms.size());
		result.pclusterBondedToPcluster.resize(result.clusters.size());
		result.particles.resize(molecule.atoms.size());

		for (int pcid = 0; pcid < result.clusters.size(); pcid++) {
			for (const int particleId : result.clusters[pcid]) {
				if (particleId != -1)
					result.particleToPcluster[particleId] = pcid;
			}
		}

		auto AddBond = [&](const auto& bond) {
			const auto& ids = bond.ids;
			for (const int id : ids) {
				if (id < 0 || id >= result.particleToPcluster.size())
					return;
			}

			for (int i = 0; i < ids.size(); i++) {
				const int pidSelf = ids[i];
				const int pcidSelf = result.particleToPcluster[pidSelf];
				assert(pcidSelf != -1);

				for (int j = i + 1; j < ids.size(); j++) {
					const int pidOther = ids[j];
					const int pcidOther = result.particleToPcluster[pidOther];
					assert(pcidOther != -1);

					result.particleBondedToParticle[pidSelf].insert(pidOther);
					result.particleBondedToParticle[pidOther].insert(pidSelf);
					result.pclusterBondedToPcluster[pcidSelf].insert(pcidOther);
					result.pclusterBondedToPcluster[pcidOther].insert(pcidSelf);
				}
			}
		};

		for (const auto& bond : molecule.singlebonds)
			AddBond(bond);
		for (const auto& bond : molecule.anglebonds)
			AddBond(bond);
		for (const auto& bond : molecule.dihedralbonds)
			AddBond(bond);
		for (const auto& bond : molecule.improperdihedralbonds)
			AddBond(bond);

		for (int particleId = 0; particleId < molecule.atoms.size(); particleId++) {
			const TopologyFile::AtomsEntry& atom = molecule.atoms[particleId];
			PersistentParticleTemplate& particle = result.particles[particleId];

			particle.nbParams = forcefield.GetLjParameters(atom.type);
			if (atom.charge.has_value())
				particle.nbParams.charge = atom.charge.value() * elementaryChargeToKiloCoulombPerMole;

			if (atom.mass.has_value()) {
				particle.mass = atom.mass.value() / KILO;
			}
			else {
				const std::optional<AtomType> atomType = forcefield.GetAtomtype(atom.type);
				if (atomType.has_value())
					particle.mass = atomType->mass;
			}

			particle.atomLetter = !atom.atomname.empty() ? atom.atomname[0] : ' ';
			particle.isSolvent = IsSolventResidue(atom.residue);
			assert(particle.mass > 0.f);
		}

		return result;
	}
}

PersistentClusterFactory MakePersistentClusters(const SuperTopology& system, LIMAForcefield& forcefield) {
	std::unordered_map<const TopologyFile::Moleculetype*, PersistentClusterTemplate> templates;
	std::vector<const PersistentClusterTemplate*> instanceTemplates;
	instanceTemplates.reserve(system.moleculeInstances.size());

	size_t totalClusterCount = 0;
	for (const SuperTopology::MoleculeInstance& instance : system.moleculeInstances) {
		auto templateIt = templates.find(instance.type);
		if (templateIt == templates.end())
			templateIt = templates.emplace(instance.type, BuildPersistentClusterTemplate(*instance.type, forcefield)).first;

		instanceTemplates.push_back(&templateIt->second);
		totalClusterCount += templateIt->second.clusters.size();
	}

	PersistentClusterFactory pcFactory{};
	pcFactory.pClusters.resize(totalClusterCount);
	pcFactory.pClusterMetas.resize(totalClusterCount);
	pcFactory.particleToPclusterMap.resize(system.particles.size());
	pcFactory.particleBondedToParticle.resize(system.particles.size());
	pcFactory.pclusterBondedToPcluster.resize(totalClusterCount);

	int pclusterOffset = 0;
	for (int instanceId = 0; instanceId < system.moleculeInstances.size(); instanceId++) {
		const SuperTopology::MoleculeInstance& instance = system.moleculeInstances[instanceId];
		const PersistentClusterTemplate& persistentTemplate = *instanceTemplates[instanceId];

		for (int localPcid = 0; localPcid < persistentTemplate.clusters.size(); localPcid++) {
			const int globalPcid = pclusterOffset + localPcid;
			const auto& localCluster = persistentTemplate.clusters[localPcid];

			for (int pidRel = 0; pidRel < PersistentCluster::maxParticles; pidRel++) {
				const int localParticleId = localCluster[pidRel];
				if (localParticleId == -1) {
					pcFactory.pClusters[globalPcid].pqd[pidRel] = PData{};
					continue;
				}

				const int globalParticleId = instance.particleOffset + localParticleId;
				const PersistentParticleTemplate& particle = persistentTemplate.particles[localParticleId];
				PersistentClusterMeta& meta = pcFactory.pClusterMetas[globalPcid];

				pcFactory.pClusters[globalPcid].pqd[pidRel] = PData{ system.particles[globalParticleId].position, particle.nbParams };
				meta.particleIdsGlobal[pidRel] = globalParticleId;
				meta.mass[pidRel] = particle.mass;
				meta.atomLetter[pidRel] = particle.atomLetter;
				meta.isSolvent = particle.isSolvent;
				meta.nParticles++;

				pcFactory.particleToPclusterMap[globalParticleId] = ParticleToPclusterMapping{ globalPcid, pidRel };
			}
		}

		for (int localParticleId = 0; localParticleId < persistentTemplate.particleBondedToParticle.size(); localParticleId++) {
			pcFactory.particleBondedToParticle[instance.particleOffset + localParticleId] =
				ParticlesBondedToParticle::Create(persistentTemplate.particleBondedToParticle[localParticleId], instance.particleOffset);
		}

		for (int localPcid = 0; localPcid < persistentTemplate.pclusterBondedToPcluster.size(); localPcid++) {
			pcFactory.pclusterBondedToPcluster[pclusterOffset + localPcid] =
				PclustersBondedToPcluster::Create(persistentTemplate.pclusterBondedToPcluster[localPcid], pclusterOffset);
		}

		pclusterOffset += static_cast<int>(persistentTemplate.clusters.size());
	}

	assert(pclusterOffset == totalClusterCount);
	return pcFactory;
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
	//TimeIt timer("buildMolecules", true);
	LIMAForcefield forcefield{ topol_file.forcefieldInclude ? topol_file.forcefieldInclude->contents : GenericItpFile{} };

	SuperTopology superTopology(topol_file.GetSystem(), grofile, forcefield);
	superTopology.VerifyBondsAreStable(grofile.box_size, simparams.bc_select, simparams.em_variant);


	std::future<BondGroupFactory> bgfFuture = std::async(std::launch::async, [&]{ return BondGroupFactory(superTopology); });

	// Make PersistenClusters
	std::shared_ptr<MoleculeGraph> systemGraph = MakeMoleculeGraph(superTopology);
	std::future<PersistentClusterFactory> pcFactoryFuture = std::async(std::launch::async, [&]{ return MakePersistentClusters(superTopology, forcefield); });
	
	BondGroupFactory bgFactory = bgfFuture.get();
	const auto particleToBondgroupMap = bgFactory.MakeParticleToBondgroupsMap(superTopology.particles.size());


	PersistentClusterFactory pcFactory = pcFactoryFuture.get();
	bgFactory.AddPclusterRefs(pcFactory.particleToPclusterMap);

	std::vector<BondGroup> bondGroups = bgFactory.GetBondgroups();
	

	for (int i = 0; i < pcFactory.particleToPclusterMap.size(); i++) {
		const auto pcRef = pcFactory.particleToPclusterMap[i];
		const std::set<BondgroupRef>& bgRefs = particleToBondgroupMap[i];

		for (const BondgroupRef& bgRef : bgRefs) {
			pcFactory.pClusterMetas[pcRef.pcid].bondgroupReferences[pcRef.pid].Add(bgRef);
		}
	}



	int nParticles = 0;
	//int nSolvents = 0;
	for (const auto& pc : pcFactory.pClusterMetas) {
		for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
			if (pc.particleIdsGlobal[pid] == -1)
				continue;
			nParticles++;
		}
	}

	std::vector<std::tuple<int, int>> gpidToPcidAndPid(nParticles);
	for (int pcid = 0; pcid < pcFactory.pClusterMetas.size(); pcid++) {
		for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
			int gpid = pcFactory.pClusterMetas[pcid].particleIdsGlobal[pid];
			if (gpid != -1)
				gpidToPcidAndPid[gpid] = { pcid, pid };
		}
	}

	return std::make_unique<BoxImage>(
		grofile,	// TODO: wierd ass copy here. Probably make the input a sharedPtr?
		std::move(superTopology),
		systemGraph,
		std::move(bondGroups),
		std::move(pcFactory.pClusters),
		std::move(pcFactory.pClusterMetas),
		std::move(pcFactory.particleBondedToParticle),
		std::move(pcFactory.pclusterBondedToPcluster),
		std::move(gpidToPcidAndPid),
		nParticles
	);

}
