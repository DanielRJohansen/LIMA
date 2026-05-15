#include "CompoundBuilder.h"

#include "map"
#include "queue"
#include "set"
#include <algorithm>
#include "TimeIt.h"


class RandomAccessDeleteSet {
private:
	std::vector<bool> deleted;    // Tracks deleted elements.
	size_t first_valid_index;     // Caches the index of the first valid element.
	const size_t nElements;       // Number of elements in the set.

	// Helper function to find the next valid index
	void updateFirstValidIndex() {
		while (first_valid_index < deleted.size() && deleted[first_valid_index]) {
			++first_valid_index;
		}
	}

public:
	// Constructor that initializes with a list of elements.
	RandomAccessDeleteSet(int nElements)
		: nElements(nElements), deleted(nElements, false), first_valid_index(0) {}

	// Returns the first not deleted element.
	int front() const {
		if (first_valid_index >= nElements) {
			throw std::out_of_range("No valid elements in the queue");
		}
		return first_valid_index;
	}

	// Marks an element at the given index as deleted.
	void erase(size_t id) {
		if (id >= nElements) {
			throw std::out_of_range("Index out of range");
		}
		deleted[id] = true;
		if (id == first_valid_index) {
			updateFirstValidIndex();
		}
	}

	void pop() {
		deleted[first_valid_index] = true;
		updateFirstValidIndex();
	}

	// Checks if the element at the given index is not deleted.
	bool contains(size_t id) const {
		if (id >= nElements) {
			return false;
		}
		return !deleted[id];
	}

	bool empty() const {
		return first_valid_index >= nElements;
	}
};



















bool MoreWorkToBeDone(
	const RandomAccessDeleteSet& availableSinglebondIds,
	const RandomAccessDeleteSet& availablePairbondIds,
	const RandomAccessDeleteSet& availableAnglebondIds,
	const RandomAccessDeleteSet& availableDihedralbondIds,
	const RandomAccessDeleteSet& availableImproperDihedralbondIds)
{
	return !availableSinglebondIds.empty() || !availablePairbondIds.empty() || !availableAnglebondIds.empty() || !availableDihedralbondIds.empty() || !availableImproperDihedralbondIds.empty();
}


enum Bondtype
{
	single, pair, angle, dihedral, improper
};

bool ReplaceIfSmaller(int& min, int newVal) {
	if (newVal < min) {
		min = newVal;
		return true;
	}
	return false;
}

Bondtype GetBondtypeWithLowestAvailableParticleId(
	const LIMA_MOLECULEBUILD::SuperTopology& topology,
	const RandomAccessDeleteSet& availableSinglebondIds,
	const RandomAccessDeleteSet& availablePairbondIds,
	const RandomAccessDeleteSet& availableAnglebondIds,
	const RandomAccessDeleteSet& availableDihedralbondIds,
	const RandomAccessDeleteSet& availableImproperDihedralbondIds)
{
	int minParticleId = std::numeric_limits<int>::max();

	Bondtype type = single;

	if (!availableSinglebondIds.empty()) {
		const auto& ids = topology.singlebonds[availableSinglebondIds.front()].global_atom_indexes;
		if (ReplaceIfSmaller(minParticleId, *std::min_element(ids.begin(), ids.end())))
			type = single;		
	}
	if (!availablePairbondIds.empty()) {
		const auto& ids = topology.pairbonds[availablePairbondIds.front()].global_atom_indexes;
		if (ReplaceIfSmaller(minParticleId, *std::min_element(ids.begin(), ids.end())))
			type = pair;
	}
	if (!availableAnglebondIds.empty()) {
		const auto& ids = topology.anglebonds[availableAnglebondIds.front()].global_atom_indexes;
		if (ReplaceIfSmaller(minParticleId, *std::min_element(ids.begin(), ids.end())))
			type = angle;
	}
	if (!availableDihedralbondIds.empty()) {
		const auto& ids = topology.dihedralbonds[availableDihedralbondIds.front()].global_atom_indexes;
		if (ReplaceIfSmaller(minParticleId, *std::min_element(ids.begin(), ids.end())))
			type = dihedral;
	}
	if (!availableImproperDihedralbondIds.empty()) {
		const auto& ids = topology.improperdihedralbonds[availableImproperDihedralbondIds.front()].global_atom_indexes;
		if (ReplaceIfSmaller(minParticleId, *std::min_element(ids.begin(), ids.end())))
			type = improper;
	}

	return type;
}

template <typename BondType>
std::vector<std::vector<int>> mapParticleToBondIds(const std::vector<BondType>& bonds, size_t particleCount) {
	std::vector<std::vector<int>> particleToBondMap(particleCount);
	for (int bondId = 0; bondId < bonds.size(); ++bondId) {
		for (uint32_t particleId : bonds[bondId].global_atom_indexes) {
			particleToBondMap[particleId].push_back(bondId);
		}
	}
	return particleToBondMap;
}






BondGroupFactory::BondGroupFactory(const LIMA_MOLECULEBUILD::SuperTopology& topology) {
	if (topology.singlebonds.empty()) return;
	

	const std::vector<std::vector<int>> pid2SinglebondIdMap = mapParticleToBondIds(topology.singlebonds, topology.particles.size());
	const std::vector<std::vector<int>> pid2PairbondIdMap = mapParticleToBondIds(topology.pairbonds, topology.particles.size());
	const std::vector<std::vector<int>> pid2AnglebondIdMap = mapParticleToBondIds(topology.anglebonds, topology.particles.size());
	const std::vector<std::vector<int>> pid2DihedralbondIdMap = mapParticleToBondIds(topology.dihedralbonds, topology.particles.size());
	const std::vector<std::vector<int>> pid2ImproperDihedralbondIdMap = mapParticleToBondIds(topology.improperdihedralbonds, topology.particles.size());


	RandomAccessDeleteSet availableSinglebondIds(topology.singlebonds.size());
	RandomAccessDeleteSet availablePairbondIds(topology.pairbonds.size());
	RandomAccessDeleteSet availableAnglebondIds(topology.anglebonds.size());
	RandomAccessDeleteSet availableDihedralbondIds(topology.dihedralbonds.size());
	RandomAccessDeleteSet availableImproperDihedralbondIds(topology.improperdihedralbonds.size());




	const int expectedNumGroups = static_cast<int>((static_cast<float>(topology.singlebonds.size()) / static_cast<float>(BondGroup::maxSinglebonds)) * 2.f);
	bondgroups.reserve(expectedNumGroups);
	particleGlobalIds.reserve(expectedNumGroups);

	int currentParticleIndexInGroup = 0;

	while (MoreWorkToBeDone(availableSinglebondIds, availablePairbondIds, availableAnglebondIds, availableDihedralbondIds, availableImproperDihedralbondIds)) {
		
		// To start or continue a group, simple add the bond containing the next lowest particleId
		const Bondtype typeOfBondWithLowestId = GetBondtypeWithLowestAvailableParticleId(topology, availableSinglebondIds, availablePairbondIds, availableAnglebondIds, availableDihedralbondIds, availableImproperDihedralbondIds);
		

		// Because if we start a new group with a zero-param bond, we skip that bond so this is to avoid empty groups.. 
		//Annoying to deal with here, maybe discard the bonds as we make the topology instead?
		if (bondgroups.empty() || bondgroups.back().nParticles != 0) {
			bondgroups.push_back({});
			particleGlobalIds.push_back({});
		}
		currentParticleIndexInGroup = 0;
		

		const int currentBondgroupIndex = static_cast<int>(bondgroups.size()) - 1;

		switch (typeOfBondWithLowestId) {
		case single:
			AddBond(currentBondgroupIndex, topology.singlebonds[availableSinglebondIds.front()]);
			availableSinglebondIds.erase(availableSinglebondIds.front());
			break;
		case pair:
			AddBond(currentBondgroupIndex, topology.pairbonds[availablePairbondIds.front()]);
			availablePairbondIds.erase(availablePairbondIds.front());
			break;
		case angle:
			AddBond(currentBondgroupIndex, topology.anglebonds[availableAnglebondIds.front()]);
			availableAnglebondIds.erase(availableAnglebondIds.front());
			break;
		case dihedral:
			AddBond(currentBondgroupIndex, topology.dihedralbonds[availableDihedralbondIds.front()]);
			availableDihedralbondIds.erase(availableDihedralbondIds.front());
			break;
		case improper:
			AddBond(currentBondgroupIndex, topology.improperdihedralbonds[availableImproperDihedralbondIds.front()]);
			availableImproperDihedralbondIds.erase(availableImproperDihedralbondIds.front());
			break;
		}


		// Now go from the current index of particles in the group to the last. 
		// For each particle find all bonds that contain said particle, and add those bonds to this group IF we have room
		// Then move to the next particle and repeat
		// We exit when, either we have no more bonds in the chain, or the group has no more room
		for (; currentParticleIndexInGroup < bondgroups.back().nParticles; currentParticleIndexInGroup++) {
			const int currentParticleId = particleGlobalIds.back()[currentParticleIndexInGroup];
			

			AddBondsFromMap(currentBondgroupIndex, pid2SinglebondIdMap[currentParticleId], availableSinglebondIds, topology.singlebonds);
			AddBondsFromMap(currentBondgroupIndex, pid2PairbondIdMap[currentParticleId], availablePairbondIds, topology.pairbonds);
			AddBondsFromMap(currentBondgroupIndex, pid2AnglebondIdMap[currentParticleId], availableAnglebondIds, topology.anglebonds);
			AddBondsFromMap(currentBondgroupIndex, pid2DihedralbondIdMap[currentParticleId], availableDihedralbondIds, topology.dihedralbonds);
			AddBondsFromMap(currentBondgroupIndex, pid2ImproperDihedralbondIdMap[currentParticleId], availableImproperDihedralbondIds, topology.improperdihedralbonds);
		}
	}


	//return bondgroups;
}
std::vector<BondGroup> BondGroupFactory::GetBondgroups() {
	// TODO Add safety checks here that addpclusterinfo has been called, and that this hasnt been called before.
	return std::move(bondgroups);
}

void BondGroupFactory::AddPclusterRefs(const ParticleToPclusterMap& particleToPclusterMap) {
	for (int bgIndex = 0; bgIndex < bondgroups.size(); bgIndex++) {
		BondGroup& group = bondgroups[bgIndex];
		for (int i = 0; i < group.nParticles; i++) {
			const int globalId = particleGlobalIds[bgIndex][i];
			group.particles[i] = BondGroup::ParticleRef{ particleToPclusterMap[globalId].pcid, particleToPclusterMap[globalId].pid };
		}
	}
}



std::vector<std::set<BondgroupRef>> BondGroupFactory::MakeParticleToBondgroupsMap(int nParticlesTotal) const {
	std::vector<std::set<BondgroupRef>> particleToBondgroupsMap(nParticlesTotal);

	for (int bgIndex = 0; bgIndex < bondgroups.size(); bgIndex++) {
		const BondGroup& group = bondgroups[bgIndex];

		for (int particleLocalId = 0; particleLocalId < group.nParticles; particleLocalId++) {
			const uint32_t particleGlobalId = particleGlobalIds[bgIndex][particleLocalId];
			particleToBondgroupsMap[particleGlobalId].insert({ bgIndex, particleLocalId });
		}
	}
	
	int maxGroupsForAParticle = 0;
	for (const auto& set : particleToBondgroupsMap) {
		maxGroupsForAParticle = std::max(maxGroupsForAParticle, static_cast<int>(set.size()));
	}
//	printf("%d\n", maxGroupsForAParticle);

	return particleToBondgroupsMap;
}

//std::vector<BondGroup> BondGroupFactory::FinishBondgroups(const std::vector<BondGroupFactory>& in) {
//	std::vector<BondGroup> out(in.size());
//	for (int i = 0; i < in.size(); i++) {
//		out[i] = in[i];
//	}
//	return out;
//}










template <int n>
std::tuple<int, std::array<uint8_t, n>> BondGroupFactory::TryAssignLocalIds(int bgIndex, const std::array<int, n>& particleIds) const {
	int nNewParticles = 0;
	std::array<uint8_t, n> localIds;

	for (int i = 0; i < n; i++) {
		int localId = FindLocalParticleId(bgIndex, particleIds[i]);		
		if (localId == -1) {
			localIds[i] = bondgroups[bgIndex].nParticles + nNewParticles;
			nNewParticles++;			
		}
		else {
			localIds[i] = localId;
		}

	}

	return { nNewParticles, localIds };
}
template std::tuple<int, std::array<uint8_t, 2>> BondGroupFactory::TryAssignLocalIds(int bgIndex, const std::array<int, 2>& particleIds) const;
template std::tuple<int, std::array<uint8_t, 3>> BondGroupFactory::TryAssignLocalIds(int bgIndex, const std::array<int, 3>& particleIds) const;
template std::tuple<int, std::array<uint8_t, 4>> BondGroupFactory::TryAssignLocalIds(int bgIndex, const std::array<int, 4>& particleIds) const;


//template <int n>
//std::array<uint8_t, n> BondGroupFactory::GetLocalIds(const std::array<int, n>& globalIds) const {
//	std::array<uint8_t, n> localIds;
//
//	for (int i = 0; i < n; i++) {
//		const int localId = FindLocalParticleId(globalIds[i]);
//		if (localId == -1)
//			throw std::runtime_error("Global id not found in bondgroup");
//
//		localIds[i] = static_cast<uint8_t>(localId);
//	}
//
//	return localIds;
//}
//template std::array<uint8_t, 2> BondGroupFactory::GetLocalIds(const std::array<int, 2>& globalIds) const;
//template std::array<uint8_t, 3> BondGroupFactory::GetLocalIds(const std::array<int, 3>& globalIds) const;
//template std::array<uint8_t, 4> BondGroupFactory::GetLocalIds(const std::array<int, 4>& globalIds) const;

void BondGroupFactory::AddBondParticles(int bgIndex, std::span<const int> bondGlobalIds, std::span<const uint8_t> localIds) {
	if (bondGlobalIds.front() == bondGlobalIds.back()) {
		assert(false);
	}

	for (int i = 0; i < bondGlobalIds.size(); i++) {		
		if (localIds[i] >= bondgroups[bgIndex].nParticles){
			const int id = bondGlobalIds[i];

			particleGlobalIds[bgIndex][bondgroups[bgIndex].nParticles] = id;
			//particles[nParticles] = ParticleRef{ particleToPclustermap.at(id).pcid, particleToPclustermap.at(id).pid };
			bondgroups[bgIndex].nParticles++;
		}
	}
}

int BondGroupFactory::FindLocalParticleId(int bgIndex, int globalId) const {
	for (int i = 0; i < bondgroups[bgIndex].nParticles; i++) {
		if (particleGlobalIds[bgIndex][i] == globalId)
			return i;
	}

	return -1;
}


bool BondGroupFactory::AddBond(int bgIndex, const SingleBondFactory& bond) {
	if (bond.params.HasZeroParam())
		return true;

	auto [nNewParticles, localIds] = TryAssignLocalIds(bgIndex, bond.global_atom_indexes);
	if (bondgroups[bgIndex].nParticles + nNewParticles > BondGroup::maxParticles) {
		return false;
	}

	AddBondParticles(bgIndex, bond.global_atom_indexes, localIds);
	if (bondgroups[bgIndex].nSinglebonds >= bondgroups[bgIndex].maxSinglebonds) {
		throw std::runtime_error("Too many bonds in bondgroup");
	}
	bondgroups[bgIndex].singlebonds[bondgroups[bgIndex].nSinglebonds++] = SingleBond{ localIds, bond.params };
	return true;
}

bool BondGroupFactory::AddBond(int bgIndex, const PairBondFactory& bond) {
	if (bond.params.HasZeroParam())
		return true;

	auto [nNewParticles, localIds] = TryAssignLocalIds(bgIndex, bond.global_atom_indexes);
	if (bondgroups[bgIndex].nParticles + nNewParticles > BondGroup::maxParticles) {
		return false;
	}

	AddBondParticles(bgIndex, bond.global_atom_indexes, localIds);

	if (bondgroups[bgIndex].nPairbonds >= bondgroups[bgIndex].maxPairbonds) {
		throw std::runtime_error("Too many bonds in bondgroup");
	}
	bondgroups[bgIndex].pairbonds[bondgroups[bgIndex].nPairbonds++] = PairBond{ localIds, bond.params };
	return true;
}

bool BondGroupFactory::AddBond(int bgIndex, const AngleBondFactory& bond) {
	if (bond.params.HasZeroParam())
		return true;

	auto [nNewParticles, localIds] = TryAssignLocalIds(bgIndex, bond.global_atom_indexes);
	if (bondgroups[bgIndex].nParticles + nNewParticles > BondGroup::maxParticles) {
		return false;
	}
	AddBondParticles(bgIndex, bond.global_atom_indexes, localIds);

	if (bondgroups[bgIndex].nAnglebonds >= bondgroups[bgIndex].maxAnglebonds) {
		throw std::runtime_error("Too many bonds in bondgroup");
	}
	bondgroups[bgIndex].anglebonds[bondgroups[bgIndex].nAnglebonds++] = AngleUreyBradleyBond{ localIds, bond.params };
	return true;
}

bool BondGroupFactory::AddBond(int bgIndex, const DihedralBondFactory& bond) {
	if (bond.params.HasZeroParam())
		return true;

	auto [nNewParticles, localIds] = TryAssignLocalIds(bgIndex, bond.global_atom_indexes);
	if (bondgroups[bgIndex].nParticles + nNewParticles > BondGroup::maxParticles) {
		return false;
	}
	AddBondParticles(bgIndex, bond.global_atom_indexes, localIds);

	if (bondgroups[bgIndex].nDihedralbonds >= bondgroups[bgIndex].maxDihedralbonds) {
		throw std::runtime_error("Too many bonds in bondgroup");
	}
	bondgroups[bgIndex].dihedralbonds[bondgroups[bgIndex].nDihedralbonds++] = DihedralBond{ localIds, bond.params };
	return true;
}

bool BondGroupFactory::AddBond(int bgIndex, const ImproperDihedralBondFactory& bond) {
	if (bond.params.HasZeroParam())
		return true;

	auto [nNewParticles, localIds] = TryAssignLocalIds(bgIndex, bond.global_atom_indexes);
	if (bondgroups[bgIndex].nParticles + nNewParticles > BondGroup::maxParticles) {
		return false;
	}
	AddBondParticles(bgIndex, bond.global_atom_indexes, localIds);

	if (bondgroups[bgIndex].nImproperdihedralbonds >= bondgroups[bgIndex].maxImproperdihedralbonds) {
		throw std::runtime_error("Too many bonds in bondgroup");
	}
	bondgroups[bgIndex].improperdihedralbonds[bondgroups[bgIndex].nImproperdihedralbonds++] = ImproperDihedralBond{ localIds, bond.params };
	return true;
}
