#include "CompoundBuilder.h"

#include "map"
#include "queue"
#include "set"
#include <algorithm>

class WorkQueue {
	struct ParticleEntry {
		int particleId;
		std::queue<int> singlebondIds;
		std::queue<int> anglebondIds;
		std::queue<int> dihedralbondIds;
		std::queue<int> improperDihedralbondIds;

		bool empty() const {
			return singlebondIds.empty() && anglebondIds.empty() && dihedralbondIds.empty() && improperDihedralbondIds.empty();
		}
	};


	std::queue<ParticleEntry> queue;
};


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
	const RandomAccessDeleteSet& availableAnglebondIds,
	const RandomAccessDeleteSet& availableDihedralbondIds,
	const RandomAccessDeleteSet& availableImproperDihedralbondIds)
{
	return !availableSinglebondIds.empty() || !availableAnglebondIds.empty() || !availableDihedralbondIds.empty() || !availableImproperDihedralbondIds.empty();
}


enum Bondtype
{
	single, angle, dihedral, improper
};

Bondtype GetBondtypeWithLowestAvailableParticleId(
	const LIMA_MOLECULEBUILD::SuperTopology& topology,
	const RandomAccessDeleteSet& availableSinglebondIds,
	const RandomAccessDeleteSet& availableAnglebondIds,
	const RandomAccessDeleteSet& availableDihedralbondIds,
	const RandomAccessDeleteSet& availableImproperDihedralbondIds)
{
	int minParticleId = std::numeric_limits<int>::max();

	Bondtype type = single;

	if (!availableSinglebondIds.empty()) {
		const auto& ids = topology.singlebonds[availableSinglebondIds.front()].global_atom_indexes;
		minParticleId = std::min(minParticleId, *std::min_element(ids.begin(), ids.end()));
		type = single;
	}
	if (!availableAnglebondIds.empty()) {
		const auto& ids = topology.anglebonds[availableAnglebondIds.front()].global_atom_indexes;
		minParticleId = std::min(minParticleId, *std::min_element(ids.begin(), ids.end()));
		type = angle;
	}
	if (!availableDihedralbondIds.empty()) {
		const auto& ids = topology.dihedralbonds[availableDihedralbondIds.front()].global_atom_indexes;
		minParticleId = std::min(minParticleId, *std::min_element(ids.begin(), ids.end()));
		type = dihedral;
	}
	if (!availableImproperDihedralbondIds.empty()) {
		const auto& ids = topology.improperdihedralbonds[availableImproperDihedralbondIds.front()].global_atom_indexes;
		minParticleId = std::min(minParticleId, *std::min_element(ids.begin(), ids.end()));
		type = improper;
	}

	return type;
}


std::vector<BondGroupFactory> BondGroupFactory::MakeBondgroups(const LIMA_MOLECULEBUILD::SuperTopology& topology, const std::vector<ParticleToCompoundMapping>& particlesToCompoundIdMap) {

	if (topology.singlebonds.empty()) return {};

	std::vector<std::vector<int>> pid2SinglebondIdMap(topology.particles.size());
	for (int bid = 0; bid < topology.singlebonds.size(); bid++)
		for (uint32_t pid : topology.singlebonds[bid].global_atom_indexes)
			pid2SinglebondIdMap[pid].push_back(bid);
	std::vector<std::vector<int>> pid2AnglebondIdMap(topology.particles.size());
	for (int bid = 0; bid < topology.anglebonds.size(); bid++)
		for (uint32_t pid : topology.anglebonds[bid].global_atom_indexes)
			pid2AnglebondIdMap[pid].push_back(bid);
	std::vector<std::vector<int>> pid2DihedralbondIdMap(topology.particles.size());
	for (int bid = 0; bid < topology.dihedralbonds.size(); bid++)
		for (uint32_t pid : topology.dihedralbonds[bid].global_atom_indexes)
			pid2DihedralbondIdMap[pid].push_back(bid);
	std::vector<std::vector<int>> pid2ImproperDihedralbondIdMap(topology.particles.size());
	for (int bid = 0; bid < topology.improperdihedralbonds.size(); bid++)
		for (uint32_t pid : topology.improperdihedralbonds[bid].global_atom_indexes)
			pid2ImproperDihedralbondIdMap[pid].push_back(bid);


	RandomAccessDeleteSet availableSinglebondIds(topology.singlebonds.size());
	RandomAccessDeleteSet availableAnglebondIds(topology.anglebonds.size());
	RandomAccessDeleteSet availableDihedralbondIds(topology.dihedralbonds.size());
	RandomAccessDeleteSet availableImproperDihedralbondIds(topology.improperdihedralbonds.size());


	std::vector<BondGroupFactory> bondgroups;
	//std::vector<BondGroupFactory> bondgroups(1, {});
	const int expectedNumGroups = static_cast<int>((static_cast<float>(topology.singlebonds.size()) / static_cast<float>(BondGroup::maxSinglebonds)) * 2.f);
	bondgroups.reserve(expectedNumGroups);


	int currentParticleIndexInGroup = 0;

	while (MoreWorkToBeDone(availableSinglebondIds, availableAnglebondIds, availableDihedralbondIds, availableImproperDihedralbondIds)) {
		
		// We can have disconnected molecules in the same bondgroup. But if we reach that point, and we only have room for < 4 particles we start a new group.
		// OR NOT! Since we need relative positions within a bondgroup
		/*if (BondGroup::maxParticles - bondgroups.back().nParticles < 4) {
			bondgroups.push_back({});			
			currentParticleIndexInGroup = 0;
		}*/


		// Because if we start a new group with a zero-param bond, we skip that bond so this is to avoid empty groups.. 
			//Annoying to deal with here, maybe discard the bonds as we make the topology instead?
		if (bondgroups.empty() || bondgroups.back().nParticles != 0) 
			bondgroups.push_back({});
		currentParticleIndexInGroup = 0;

		// To start or continue a group, simple add the bond containing the next lowest particleId
		const Bondtype typeOfBondWithLowestId = GetBondtypeWithLowestAvailableParticleId(topology, availableSinglebondIds, availableAnglebondIds, availableDihedralbondIds, availableImproperDihedralbondIds);
		switch (typeOfBondWithLowestId) {
		case single:
			bondgroups.back().AddBond(particlesToCompoundIdMap, topology.singlebonds[availableSinglebondIds.front()]);
			availableSinglebondIds.erase(availableSinglebondIds.front());
			break;
		case angle:
			bondgroups.back().AddBond(particlesToCompoundIdMap, topology.anglebonds[availableAnglebondIds.front()]);
			availableAnglebondIds.erase(availableAnglebondIds.front());
			break;
		case dihedral:
			bondgroups.back().AddBond(particlesToCompoundIdMap, topology.dihedralbonds[availableDihedralbondIds.front()]);
			availableDihedralbondIds.erase(availableDihedralbondIds.front());
			break;
		case improper:
			bondgroups.back().AddBond(particlesToCompoundIdMap, topology.improperdihedralbonds[availableImproperDihedralbondIds.front()]);
			availableImproperDihedralbondIds.erase(availableImproperDihedralbondIds.front());
			break;
		}


		// Now go from the current index of particles in the group to the last. 
		// For each particle find all bonds that contain said particle, and add those bonds to this group IF we have room
		// Then move to the next particle and repeat
		// We exit when, either we have no more bonds in the chain, or the group has no more room
		for (; currentParticleIndexInGroup < bondgroups.back().nParticles; currentParticleIndexInGroup++) {
			const int currentParticleId = bondgroups.back().particleGlobalIds[currentParticleIndexInGroup];

			for (const int& singlebondId : pid2SinglebondIdMap[currentParticleId]) {
				const bool bondExists = availableSinglebondIds.contains(singlebondId);
				const bool atomGroupHasSpace = bondgroups.back().HasSpaceForParticlesInBond(topology.singlebonds[singlebondId].global_atom_indexes);
				if (bondExists && atomGroupHasSpace) {
					bondgroups.back().AddBond(particlesToCompoundIdMap, topology.singlebonds[singlebondId]);
					availableSinglebondIds.erase(singlebondId);
				}
			}

			for (const int& anglebondId : pid2AnglebondIdMap[currentParticleId]) {
				const bool bondExists = availableAnglebondIds.contains(anglebondId);
				const bool atomGroupHasSpace = bondgroups.back().HasSpaceForParticlesInBond(topology.anglebonds[anglebondId].global_atom_indexes);
				if (bondExists && atomGroupHasSpace) {
					bondgroups.back().AddBond(particlesToCompoundIdMap, topology.anglebonds[anglebondId]);
					availableAnglebondIds.erase(anglebondId);
				}
			}


			for (const int& dihedralbondId : pid2DihedralbondIdMap[currentParticleId]) {
				const bool bondExists = availableDihedralbondIds.contains(dihedralbondId);
				const bool atomGroupHasSpace = bondgroups.back().HasSpaceForParticlesInBond(topology.dihedralbonds[dihedralbondId].global_atom_indexes);
				if (bondExists && atomGroupHasSpace) {
					bondgroups.back().AddBond(particlesToCompoundIdMap, topology.dihedralbonds[dihedralbondId]);
					availableDihedralbondIds.erase(dihedralbondId);
				}
			}

			for (const int& improperDihedralbondId : pid2ImproperDihedralbondIdMap[currentParticleId]) {
				const bool bondExists = availableImproperDihedralbondIds.contains(improperDihedralbondId);
				const bool atomGroupHasSpace = bondgroups.back().HasSpaceForParticlesInBond(topology.improperdihedralbonds[improperDihedralbondId].global_atom_indexes);
				if (bondExists && atomGroupHasSpace) {
					bondgroups.back().AddBond(particlesToCompoundIdMap, topology.improperdihedralbonds[improperDihedralbondId]);
					availableImproperDihedralbondIds.erase(improperDihedralbondId);
				}
			}
		}
	}

	return bondgroups;
}





std::vector<std::set<BondgroupRef>> BondGroupFactory::MakeParticleToBondgroupsMap(const std::vector<BondGroupFactory>& bondGroups, int nParticlesTotal) {
	std::vector<std::set<BondgroupRef>> particleToBondgroupsMap(nParticlesTotal);

	for (int groupId = 0; groupId < bondGroups.size(); groupId++) {
		const BondGroupFactory& group = bondGroups[groupId];

		for (int particleLocalId = 0; particleLocalId < group.nParticles; particleLocalId++) {
			const uint32_t particleGlobalId = group.particleGlobalIds[particleLocalId];
			particleToBondgroupsMap[particleGlobalId].insert({ groupId, particleLocalId });
		}
	}
	
	int maxGroupsForAParticle = 0;
	for (const auto& set : particleToBondgroupsMap) {
		maxGroupsForAParticle = std::max(maxGroupsForAParticle, static_cast<int>(set.size()));
	}
//	printf("%d\n", maxGroupsForAParticle);

	return particleToBondgroupsMap;
}

std::vector<BondGroup> BondGroupFactory::FinishBondgroups(const std::vector<BondGroupFactory>& in) {
	std::vector<BondGroup> out(in.size());
	for (int i = 0; i < in.size(); i++) {
		out[i] = in[i];
	}
	return out;
}












bool BondGroupFactory::HasSpaceForParticlesInBond(const std::span<const int>& particleIds) const {
	int nNewParticles = 0;
	for (const int id : particleIds) {
		if (!particleGlobalToLocalId.contains(id)) {
			nNewParticles++;
		}
	}
	return nNewParticles < maxParticles - nParticles;
}

template <int n>
std::array<uint8_t, n> GetLocalIds(const std::unordered_map<int, uint8_t>& particleGlobalToLocalId, const std::array<int, n>& globalIds) {
	std::array<uint8_t, n> localIds;
	for (int i = 0; i < n; i++) {
		if (!particleGlobalToLocalId.contains(globalIds[i])) {
			throw std::runtime_error("Global id not found in particleGlobalToLocalId");
		}
		localIds[i] = particleGlobalToLocalId.at(globalIds[i]);
	}
	return localIds;
}

void BondGroupFactory::AddBondParticles(const ParticleToCompoundMap& particleToCompoundMap, std::span<const int> bondGlobalIds) {
	for (int id : bondGlobalIds) {
		if (!particleGlobalToLocalId.contains(id)) {

			if (nParticles >= maxParticles) {
				throw std::runtime_error("Too many particles in bondgroup");
			};

			particleGlobalToLocalId.insert({ id, nParticles });
			particleGlobalIds[nParticles] = id;
			particles[nParticles] = ParticleRef{ particleToCompoundMap.at(id).compoundId, particleToCompoundMap.at(id).localIdInCompound };
			nParticles++;
		}
	}
}



void BondGroupFactory::AddBond(const ParticleToCompoundMap& particleToCompoundMap, const SingleBondFactory& bond) {
	if (bond.params.HasZeroParam())
		return;

	AddBondParticles(particleToCompoundMap, bond.global_atom_indexes);

	if (nSinglebonds >= maxSinglebonds) {
		throw std::runtime_error("Too many bonds in bondgroup");
	}
    singlebonds[nSinglebonds++] = SingleBond{ GetLocalIds<SingleBond::nAtoms>(particleGlobalToLocalId, bond.global_atom_indexes), bond.params };
}

void BondGroupFactory::AddBond(const ParticleToCompoundMap& particleToCompoundMap, const AngleBondFactory& bond) {
	if (bond.params.HasZeroParam())
		return;

	AddBondParticles(particleToCompoundMap, bond.global_atom_indexes);

	if (nAnglebonds >= maxAnglebonds) {
		throw std::runtime_error("Too many bonds in bondgroup");
	}
    anglebonds[nAnglebonds++] = AngleUreyBradleyBond{ GetLocalIds<AngleUreyBradleyBond::nAtoms>(particleGlobalToLocalId, bond.global_atom_indexes), bond.params };
}

void BondGroupFactory::AddBond(const ParticleToCompoundMap& particleToCompoundMap, const DihedralBondFactory& bond) {
	if (bond.params.HasZeroParam())
		return;

	AddBondParticles(particleToCompoundMap, bond.global_atom_indexes);

	if (nDihedralbonds >= maxDihedralbonds) {
		throw std::runtime_error("Too many bonds in bondgroup");
	}
    dihedralbonds[nDihedralbonds++] = DihedralBond{ GetLocalIds<DihedralBond::nAtoms>(particleGlobalToLocalId, bond.global_atom_indexes), bond.params };
}

void BondGroupFactory::AddBond(const ParticleToCompoundMap& particleToCompoundMap, const ImproperDihedralBondFactory& bond) {
	if (bond.params.HasZeroParam())
		return;

	AddBondParticles(particleToCompoundMap, bond.global_atom_indexes);

	if (nImproperdihedralbonds >= maxImproperdihedralbonds) {
		throw std::runtime_error("Too many bonds in bondgroup");
	}
    improperdihedralbonds[nImproperdihedralbonds++] = ImproperDihedralBond{ GetLocalIds<ImproperDihedralBond::nAtoms>(particleGlobalToLocalId, bond.global_atom_indexes), bond.params };
}
