#include "CompoundBuilder.h"
#include "Forcefield.h"
#include "MoleculeGraph.h"

#include <unordered_set>
#include <format>
#include <array>
#include <numeric>
#include <set>

#include "Display.h"
using namespace LIMA_MOLECULEBUILD;
using namespace LimaMoleculeGraph;



class BondedParticlesLUTManagerFactory {
	std::vector<BondedParticlesLUT> luts;
	const int nCompounds;

	void DistributeLJIgnores(BondedParticlesLUTManagerFactory* bplut_man, const std::vector<ParticleToCompoundMapping>& particleToCompoundMap, std::span<const int> global_ids) {
		for (auto gid_self : global_ids) {
			for (auto gid_other : global_ids) {
				if (gid_self == gid_other) { continue; }

				const ParticleToCompoundMapping& mappingSelf = particleToCompoundMap[gid_self];
				const ParticleToCompoundMapping& mappingOther = particleToCompoundMap[gid_other];

				if (mappingOther.compoundId == -1 || mappingSelf.compoundId == -1)
					throw std::runtime_error("compoundId is -1");


				BondedParticlesLUT& lut = bplut_man->get(mappingSelf.compoundId, mappingOther.compoundId);
				lut.set(mappingSelf.localIdInCompound, mappingOther.localIdInCompound, true);
			}
		}
	}

public:
	BondedParticlesLUTManagerFactory(int nCompounds, const SuperTopology& topology, const std::vector<ParticleToCompoundMapping>& p2cMap) :
		nCompounds(nCompounds) ,
		luts(nCompounds * BondedParticlesLUTHelpers::max_bonded_compounds, BondedParticlesLUT(false)) 
	{
		for (const auto& bond : topology.singlebonds)
			DistributeLJIgnores(this, p2cMap, bond.global_atom_indexes);
		for (const auto& bond : topology.anglebonds)
			DistributeLJIgnores(this, p2cMap, bond.global_atom_indexes);
		for (const auto& bond : topology.dihedralbonds)
			DistributeLJIgnores(this, p2cMap, bond.global_atom_indexes);
		for (const auto& bond : topology.improperdihedralbonds)
			DistributeLJIgnores(this, p2cMap, bond.global_atom_indexes);

		// Set LUT so no particle ever interacts with itself
		for (int com_id = 0; com_id < nCompounds; com_id++) {
			auto& lut = get(com_id, com_id);
			for (int pid = 0; pid < MAX_COMPOUND_PARTICLES; pid++) {
				lut.set(pid, pid, true);
			}
		}
	}


	BondedParticlesLUT& get(int id_self, int id_other) {
		if (std::abs(id_self - id_other > BondedParticlesLUTHelpers::maxDiff)) {
			throw std::runtime_error(std::format("Cannot get BPLUT for compounds with distances > 2 in id-space {} {}", id_self, id_other));
		}
		if (id_self >= nCompounds || id_other >= nCompounds) {
			throw std::runtime_error("Cannot get BPLUT for compounds with id >= nCompounds");
		}

		const int local_index = BondedParticlesLUTHelpers::getLocalIndex(id_self, id_other);
		return luts[BondedParticlesLUTHelpers::getGlobalIndex(local_index, id_self)];
	}

	// Returns device ptr to the bpLUTS
	std::vector<BondedParticlesLUT> Finish() {
		return std::move(luts);
	}

	// So actually returns ID's of compounds with which queryCompound has bonded interactions
	std::vector<int> GetIdsOfCompoundsWithEntriesInLut(int queryCompoundId) {
		std::vector<int> ids;
		for (int relIndex = -BondedParticlesLUTHelpers::maxDiff; relIndex <= BondedParticlesLUTHelpers::maxDiff; relIndex++) {
			const int otherIndex = queryCompoundId + relIndex;
			if (otherIndex < 0 || otherIndex >= nCompounds || otherIndex == queryCompoundId) { continue; }

			if (get(queryCompoundId, otherIndex).HasEntries()) {
				ids.push_back(otherIndex);
			}
		}
		return ids;
	}
};



template <int n>
std::array<int, n> TransformBondIds(const std::array<int, n>& ids, int offset) {
	std::array<int, n> out;
	for (int i = 0; i < n; i++) {
		out[i] = ids[i] + offset;
	}
	return out;
}

template <typename BondType, typename BondtypeFactory, typename BondTypeTopologyfile>
void SuperTopology::LoadBondsIntoTopology(const std::vector<BondTypeTopologyfile>& bondsInTopfile, int atomIdOffset, LIMAForcefield& forcefield, std::vector<BondtypeFactory>& topology)
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
		const bool getParamsFromForcefield = particles[globalIds[0]].topologyAtom.residue != "SOL";
			

		// A bond may be describe as multiple bonds, so this is a vector
		const std::vector<typename BondType::Parameters>& bondParams = getParamsFromForcefield
			? forcefield.GetBondParameters<BondType>(atomTypenames)
			: std::vector<typename BondType::Parameters>{ {typename BondType::Parameters{}} };

		for (const auto& param : bondParams) {
			topology.emplace_back(BondtypeFactory{ globalIds, param });
		}
	}
}

SuperTopology::SuperTopology(const TopologyFile::System& system, const GroFile& grofile, LIMAForcefield& forcefield) {

	int nextUniqueParticleId = 0;
	int indexInGrofile = 0;
	for (int topologyMoleculeIndex = 0; topologyMoleculeIndex < system.molecules.size(); topologyMoleculeIndex++) {
		const TopologyFile::MoleculeEntry& molecule = system.molecules[topologyMoleculeIndex];
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

			if (molType.atoms[0].residue == "SOL" || molType.atoms[0].residue == "TIP3") {
				indexInGrofile += molType.atoms.size() - 1;
				break;
			}
		}

		LoadBondsIntoTopology<SingleBond, SingleBondFactory, TopologyFile::SingleBond>(molType.singlebonds, particleIdOffset, forcefield, singlebonds);
		LoadBondsIntoTopology<AngleUreyBradleyBond, AngleBondFactory, TopologyFile::AngleBond>(molType.anglebonds, particleIdOffset, forcefield, anglebonds);
		LoadBondsIntoTopology<DihedralBond, DihedralBondFactory, TopologyFile::DihedralBond>(molType.dihedralbonds, particleIdOffset, forcefield, dihedralbonds);
		LoadBondsIntoTopology<ImproperDihedralBond, ImproperDihedralBondFactory, TopologyFile::ImproperDihedralBond>(molType.improperdihedralbonds, particleIdOffset, forcefield, improperdihedralbonds);
	}
}

void SuperTopology::VerifyBondsAreStable(float boxlen_nm, BoundaryConditionSelect bc_select, bool energyMinimizationMode) const {
	const float allowedScalar = energyMinimizationMode ? 7.f : 3.f;//1.9999f;

	for (const auto& bond : singlebonds)
	{
		const Float3 pos1 = particles[bond.global_atom_indexes[0]].position;
		const Float3 pos2 = particles[bond.global_atom_indexes[1]].position;
		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos1, pos2, boxlen_nm, bc_select);
		const float bondRelaxedDist = bond.params.b0 * LIMA_TO_NANO;

		if (hyper_dist > bondRelaxedDist * allowedScalar) {
			throw std::runtime_error(std::format("Loading singlebond with illegally large dist ({}). b0: {}", hyper_dist, bond.params.b0 * LIMA_TO_NANO));
		}
		if (hyper_dist < bondRelaxedDist * 0.001)
			throw std::runtime_error(std::format("Loading singlebond with illegally small dist ({}). b0: {}", hyper_dist, bond.params.b0 * LIMA_TO_NANO));
	}
	for (const auto& bond : anglebonds)
	{
		const Float3 pos1 = particles[bond.global_atom_indexes[0]].position;
		const Float3 pos2 = particles[bond.global_atom_indexes[1]].position;
		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos1, pos2, boxlen_nm, bc_select);
		if (hyper_dist < 0.001)
			throw std::runtime_error(std::format("Loading singlebond with illegally small dist ({}). b0: {}", hyper_dist, bond.params.ub0 * LIMA_TO_NANO));
	}
}











// --------------------------------------------------------------- Factory Functions --------------------------------------------------------------- //


template <int n>
std::array<uint8_t, n> ConvertGlobalAtomidsToCompoundlocalIds(const std::vector<ParticleToCompoundMapping>& p2cMap, const std::array<uint32_t, n>& global_atom_ids) {
	std::array<uint8_t, n> localIds;
	for (int i = 0; i < n; i++) {
		localIds[i] = static_cast<uint8_t>(p2cMap[global_atom_ids[i]].localIdInCompound);
	}
	return localIds;
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

void CompoundFactory::AddBondgroupReference(int particleId, const BondgroupRef& bgRef) {
	if (bondgroupReferences[particleId].nBondgroupApperances >= bondgroupReferences[particleId].maxBondgroupApperances)
		throw std::runtime_error("Failed to add bondgroup reference to compound");

	bondgroupReferences[particleId].bondgroupApperances[bondgroupReferences[particleId].nBondgroupApperances++] = bgRef;
}



std::pair<const std::vector<std::vector<int>>, const std::vector<std::vector<int>>> SeparateMolecules(const SuperTopology& system) {
	std::vector<std::vector<int>> molecules;
	std::vector<std::vector<int>> tinyMolecules;

	std::vector<std::pair<int, std::string>> atoms;
	atoms.reserve(system.particles.size());
	for (int pid = 0; pid < system.particles.size(); pid++) {
		atoms.push_back({ pid, system.particles[pid].topologyAtom.type});
	}
	std::vector<std::array<int, 2>> edges;
	edges.reserve(system.singlebonds.size());
	for (const auto& bond : system.singlebonds) {
		edges.push_back(bond.global_atom_indexes);
	}


	const auto systemGraph = std::make_shared<MoleculeGraph>(atoms, edges);

	const std::vector<std::vector<int>> particleidCollectionsOfMolecules = systemGraph->GetListOfListsofConnectedNodeids();

	for (const std::vector<int>& collection : particleidCollectionsOfMolecules) {

		const bool collectionIsCustomLimaMolecule = system.particles[collection[0]].topologyAtom.residue == "lxx";

		if (collection.size() > 3 || collectionIsCustomLimaMolecule)
			molecules.emplace_back(collection);
		else
			tinyMolecules.emplace_back(collection);
	}

	return { molecules, tinyMolecules };
}



std::vector<TinyMolFactory> LoadTinyMols(const std::vector<std::vector<int>>& particleidsInTinymols, const SuperTopology& topology, LIMAForcefield& forcefield) {
	std::vector<TinyMolFactory> tinyMols;
	tinyMols.reserve(topology.particles.size()); // not accurate

	for (const std::vector<int>& particleIds : particleidsInTinymols) {

		// FOr now ignore that there are multiple atoms in a tinymol. LOOONG TODO
		const int onlyParticleToTake = particleIds[0];

		tinyMols.emplace_back(TinyMolFactory{ 
			topology.particles[onlyParticleToTake].position,
			forcefield.GetActiveTinymoltypeIndex(topology.particles[onlyParticleToTake].topologyAtom.type),
			topology.particles[onlyParticleToTake].topologyAtom.type
			});
	}

	return tinyMols;
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

	

std::vector<int> ReorderSubchains(const std::vector<int>& ids, const std::unordered_map<int,int> nodeIdToNumDownstream, int spaceLeft) {
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



const std::vector<AtomGroup> GroupAtoms2(const std::vector<std::vector<int>>& particleidsInMolecules, const SuperTopology& topology) {
	std::vector<AtomGroup> atomGroups;


	std::vector<std::unordered_set<int>> pidToSinglebondidMap(topology.particles.size());
	for (int bid = 0; bid < topology.singlebonds.size(); bid++) {
		pidToSinglebondidMap[topology.singlebonds[bid].global_atom_indexes[0]].insert(bid);
		pidToSinglebondidMap[topology.singlebonds[bid].global_atom_indexes[1]].insert(bid);
	}


	for (const auto& particleIdsInMolecule : particleidsInMolecules) {

		std::vector<std::pair<int, std::string>> atoms;
		atoms.reserve(particleIdsInMolecule.size());
		for (int pid : particleIdsInMolecule) {
			atoms.push_back({ pid, topology.particles[pid].topologyAtom.type });
		}

		std::unordered_set<int> bondIdsInMolecule;
		for (int pid : particleIdsInMolecule) {
			for (int bid : pidToSinglebondidMap[pid]) {
				bondIdsInMolecule.insert(bid);
			}
		}

		std::vector<std::array<int, 2>> edges;
		edges.reserve(bondIdsInMolecule.size());
		for (int bid : bondIdsInMolecule) {
			edges.push_back(topology.singlebonds[bid].global_atom_indexes);
		}




		const MoleculeGraph molGraph(atoms, edges);
		const std::unordered_map<int, int> nodeIdNumDownstreamNodes = molGraph.ComputeNumDownstreamNodes();
		const MoleculeTree moleculeTree = molGraph.ConstructMoleculeTree();

		std::stack<const MoleculeGraph::Node*> nodeStack;
		nodeStack.push(molGraph.root);

		atomGroups.push_back({});

		while (!nodeStack.empty()) {
			const MoleculeGraph::Node* node = nodeStack.top();
			nodeStack.pop();

			if (MAX_COMPOUND_PARTICLES - atomGroups.back().atomIds.size() == 0)
				atomGroups.push_back({});
			atomGroups.back().atomIds.emplace_back(node->atomid);

			std::vector<int> nodeChildren = moleculeTree.GetChildIds(node->atomid);

			if (nodeChildren.empty()) {
				// finished
			}
			else {
				// Add the longest childchain to our stack, and remove it from the current children
				const int indexOfLongestChain = std::max_element(nodeChildren.begin(), nodeChildren.end(),
					[&nodeIdNumDownstreamNodes](const int& a, const int& b) { return nodeIdNumDownstreamNodes.at(a) < nodeIdNumDownstreamNodes.at(b); }
				) - nodeChildren.begin();
				nodeStack.push(&molGraph.nodes.at(nodeChildren[indexOfLongestChain]));
				nodeChildren[indexOfLongestChain] = nodeChildren.back();
				nodeChildren.pop_back();

				const std::vector<int> nodeChildrenIdsIdealOrder = ReorderSubchains(nodeChildren, nodeIdNumDownstreamNodes, MAX_COMPOUND_PARTICLES - atomGroups.back().atomIds.size());

				for (int id : nodeChildrenIdsIdealOrder) {
					if (nodeIdNumDownstreamNodes.at(id) > MAX_COMPOUND_PARTICLES - atomGroups.back().atomIds.size())
						atomGroups.push_back({});

					std::vector<int> allChildIds = moleculeTree.GetAllChildIdsAndSelf(id, nodeIdNumDownstreamNodes);
					for (int childId : allChildIds) {
						atomGroups.back().atomIds.emplace_back(childId);
					}
				}
			}
		}
	}

	// Now figure out which groups are bonded to others
	for (int i = 1; i < atomGroups.size(); i++) {
		if (AreBonded(atomGroups[i - 1], atomGroups[i], pidToSinglebondidMap)) {
			atomGroups[i].idsOfBondedAtomgroups.insert(i - 1);
			atomGroups[i - 1].idsOfBondedAtomgroups.insert(i);
		}
	}

	return atomGroups;
}


std::vector<CompoundFactory> CreateCompounds(const SuperTopology& topology, float boxlen_nm, const std::vector<AtomGroup>& atomGroups, BoundaryConditionSelect bc_select)
{
	std::vector<CompoundFactory> compounds;
	std::vector<int> atomGroupToCompoundIdMap(atomGroups.size());

	for (int atomgroupIndex = 0; atomgroupIndex < atomGroups.size(); atomgroupIndex++) {
		const AtomGroup& atomGroup = atomGroups[atomgroupIndex];

		const bool is_bonded_with_previous_residue = atomgroupIndex > 0 && atomGroup.idsOfBondedAtomgroups.contains(atomgroupIndex - 1);
		const bool compound_has_room_for_residue = atomgroupIndex > 0 && compounds.back().hasRoomForRes(atomGroup.atomIds.size());

		// If we are either a new molecule, or same molecule but the current compound has no more room, make new compound
		if (atomgroupIndex == 0 || !is_bonded_with_previous_residue || !compound_has_room_for_residue) {
			if (compounds.size() >= MAX_COMPOUNDS) {
				throw std::runtime_error(std::format("Cannot handle more than {} compounds", MAX_COMPOUNDS).c_str());
			}

			compounds.push_back(CompoundFactory{ static_cast<int>(compounds.size()) });
		}
		
		atomGroupToCompoundIdMap[atomgroupIndex] = compounds.size() - 1;

		// Add all atoms of residue to current compound
		for (int i = 0; i < atomGroup.atomIds.size(); i++) {
			const int atom_gid = atomGroup.atomIds[i];
			compounds.back().addParticle(topology.particles[atom_gid], atom_gid, boxlen_nm, bc_select);
		}
	}

	// Now find all compounds that are bonded to each other
	for (int atomGroupId = 0; atomGroupId < atomGroups.size(); atomGroupId++) {
		for (const int& bondedAtomGroupId : atomGroups[atomGroupId].idsOfBondedAtomgroups) {

			const int cidLeft = atomGroupToCompoundIdMap[atomGroupId];
			const int cidRight = atomGroupToCompoundIdMap[bondedAtomGroupId];

			// If two bonded atomGroups map to 2 different compounds, those compounds must also be bonded
			if (cidLeft != cidRight) {
				compounds[cidLeft].addIdOfBondedCompound(cidRight);
				compounds[cidRight].addIdOfBondedCompound(cidLeft);
			}
		}
	}
	
	return compounds;
}



const std::vector<ParticleToCompoundMapping> MakeParticleToCompoundidMap(const std::vector<CompoundFactory>& compounds, int nParticlesTotal) {
	std::vector<ParticleToCompoundMapping> particleToCompoundidMap(nParticlesTotal);
	for (int cid = 0; cid < compounds.size(); cid++) {
		for (int pid = 0; pid < compounds[cid].n_particles; pid++) {
			particleToCompoundidMap[compounds[cid].global_ids[pid]] = ParticleToCompoundMapping{ cid, pid };  // cid;
		}
	}
	return particleToCompoundidMap;
}






//// Check that we dont have any unrealistic bonds, and warn immediately.
//void VerifyBondsAreStable(const Topology& topology, float boxlen_nm, BoundaryConditionSelect bc_select, bool energyMinimizationMode) {	
//	const float allowedScalar = energyMinimizationMode ? 7.f : 3.f;//1.9999f;
//
//	for (const auto& bond : topology.singlebonds) 
//	{		
//		const Float3 pos1 = topology.particles[bond.global_atom_indexes[0]].position;
//		const Float3 pos2 = topology.particles[bond.global_atom_indexes[1]].position;
//		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos1, pos2, boxlen_nm, bc_select);
//		const float bondRelaxedDist = bond.params.b0 * LIMA_TO_NANO;
//
//		if (hyper_dist > bondRelaxedDist * allowedScalar) {
//			throw std::runtime_error(std::format("Loading singlebond with illegally large dist ({}). b0: {}", hyper_dist, bond.params.b0 * LIMA_TO_NANO));
//		}
//		if (hyper_dist < bondRelaxedDist * 0.001)
//			throw std::runtime_error(std::format("Loading singlebond with illegally small dist ({}). b0: {}", hyper_dist, bond.params.b0 * LIMA_TO_NANO));
//	}
//	for (const auto& bond : topology.anglebonds)
//	{
//		const Float3 pos1 = topology.particles[bond.global_atom_indexes[0]].position;
//		const Float3 pos2 = topology.particles[bond.global_atom_indexes[1]].position;
//		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos1, pos2, boxlen_nm, bc_select);
//		if (hyper_dist < 0.001)
//			throw std::runtime_error(std::format("Loading singlebond with illegally small dist ({}). b0: {}", hyper_dist, bond.params.ub0 * LIMA_TO_NANO));
//	}
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
	LIMAForcefield forcefield{ topol_file.forcefieldInclude->contents };

	const SuperTopology superTopology(topol_file.GetSystem(), grofile, forcefield);
	superTopology.VerifyBondsAreStable(grofile.box_size.x, simparams.bc_select, simparams.em_variant);

	auto [molecules, tinyMolecules] = SeparateMolecules(superTopology);

	const std::vector<AtomGroup> atomGroups = GroupAtoms2(molecules, superTopology);


	std::vector<CompoundFactory> compounds = CreateCompounds(superTopology, grofile.box_size.x, atomGroups, simparams.bc_select);

	//Display d;
	//std::vector<std::array<Float3, MAX_COMPOUND_PARTICLES>> positions;
	//for (const auto& compound : compounds) {
	//	positions.push_back({});
	//	for (int i = 0; i < MAX_COMPOUND_PARTICLES; i++) {
	//		positions.back()[i] = compound.positions[i];
	//	}
	//}
	//std::vector<Compound> compounds2;
	//for (const auto& compound : compounds) {
	//	compounds2.push_back(compound);
	//}	
	//d.Render(std::make_unique<Rendering::CompoundsTask>(compounds2, positions, grofile.box_size), true);

	//printf("%d compounds\n", compounds.size());

	const std::vector<ParticleToCompoundMapping> particleToCompoundidMap = MakeParticleToCompoundidMap(compounds, superTopology.particles.size());


	auto bpLutManager = std::make_unique<BondedParticlesLUTManagerFactory>(compounds.size(), superTopology, particleToCompoundidMap);

	for (int cid = 0; cid < compounds.size(); cid++) {
		std::vector<int> compoundIdsWithBondedInteractions = bpLutManager->GetIdsOfCompoundsWithEntriesInLut(cid);
		for (int id : compoundIdsWithBondedInteractions) {

			bool found = false;
			for (int i = 0; i < compounds[cid].n_bonded_compounds; i++)
				if (compounds[cid].bonded_compound_ids[i] == id)
					found = true;
			if (!found)
				int a = 0;
			compounds[cid].addIdOfBondedCompound(id);
		}

	}

	std::vector<BondGroupFactory> bondGroups = BondGroupFactory::MakeBondgroups(superTopology, particleToCompoundidMap);
	const auto particleToBondgroupMap = BondGroupFactory::MakeParticleToBondgroupsMap(bondGroups, superTopology.particles.size());


	for (int i = 0; i < particleToCompoundidMap.size(); i++) {
		const auto cRef = particleToCompoundidMap[i];
		const std::set<BondgroupRef>& bgRefs = particleToBondgroupMap[i];

		for (const BondgroupRef& bgRef : bgRefs) {
			compounds[cRef.compoundId].AddBondgroupReference(cRef.localIdInCompound, bgRef);
		}
	}



	//bpLutManager->get(0, 0)->printMatrix(compounds.begin()->n_particles);

	CompoundFactory::CalcCompoundMetaInfo(grofile.box_size.x, compounds, simparams.bc_select);

	const std::vector<TinyMolFactory> tinyMols = LoadTinyMols(tinyMolecules, superTopology, forcefield);
	const int totalCompoundParticles = std::accumulate(compounds.begin(), compounds.end(), 0, [](int sum, const auto& compound) { return sum + compound.n_particles; });

	return std::make_unique<BoxImage>(
		std::move(compounds),
		static_cast<int>(totalCompoundParticles),
		bpLutManager->Finish(),
		std::move(tinyMols),
		grofile,	// TODO: wierd ass copy here. Probably make the input a sharedPtr?
		forcefield.GetActiveLjParameters(),
		forcefield.GetTinymolTypes(),
		superTopology,
		forcefield.GetNonbondedInteractionParams(),
		BondGroupFactory::FinishBondgroups(bondGroups)
	);

}
