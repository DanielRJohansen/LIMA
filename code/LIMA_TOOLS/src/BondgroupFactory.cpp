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
		minParticleId = std::min(minParticleId, *std::min_element(ids.begin(), ids.end()));
		type = single; // TODO: these are all wrong, only set type if we set new min
	}
	if (!availablePairbondIds.empty()) {
		const auto& ids = topology.pairbonds[availablePairbondIds.front()].global_atom_indexes;
		minParticleId = std::min(minParticleId, *std::min_element(ids.begin(), ids.end()));
		type = pair;
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

// Add bonds from a specific type to the bond group
void AddBondsFromMap(auto& bondGroup, const auto& bondMap, auto& availableBondIds, const auto& bonds, const auto& particlesToCompoundIdMap) {
	for (const int bondId : bondMap) {
		if (availableBondIds.contains(bondId) && bondGroup.HasSpaceForParticlesInBond(bonds[bondId].global_atom_indexes)) {
			bondGroup.AddBond(particlesToCompoundIdMap, bonds[bondId]);
			availableBondIds.erase(bondId);
		}
	}
};

std::vector<BondgroupTinymol> TinyMolFactory::MakeBondgroups(const LIMA_MOLECULEBUILD::SuperTopology& topology, const std::vector<std::vector<int>>& tinymolParticlesIds) {
	std::vector<BondgroupTinymol> bondgroups(tinymolParticlesIds.size());

	// TODO: OPTIM: Waste that we make these maps both for compounds and tinymols, they could reuse the same
	const std::vector<std::vector<int>> pid2SinglebondIdMap = mapParticleToBondIds(topology.singlebonds, topology.particles.size());
	const std::vector<std::vector<int>> pid2AnglebondIdMap = mapParticleToBondIds(topology.anglebonds, topology.particles.size());

	for (int tid = 0; tid < tinymolParticlesIds.size(); tid++) {
		bondgroups[tid].nParticles = tinymolParticlesIds[tid].size();

		// First find the bonds that are in this tinymol
		std::set<int> availableSinglebondIds;
		std::set<int> availableAnglebondIds;
		for (int pid = 0; pid < tinymolParticlesIds[tid].size(); pid++) {
			const int particleId = tinymolParticlesIds[tid][pid];
			for (int singlebondId : pid2SinglebondIdMap[particleId]) {
				availableSinglebondIds.insert(singlebondId);
			}
			for (int anglebondId : pid2AnglebondIdMap[particleId]) {
				availableAnglebondIds.insert(anglebondId);
			}
		}

		assert(availableSinglebondIds.size() <= BondGroup::maxSinglebonds);
		assert(availableAnglebondIds.size() <= BondGroup::maxAnglebonds);

		for (int singlebondId : availableSinglebondIds) {
			std::array<uint8_t, 2> idsLocalToTinymol = {
				static_cast<uint8_t>(topology.singlebonds[singlebondId].global_atom_indexes[0] - tinymolParticlesIds[tid][0]),
				static_cast<uint8_t>(topology.singlebonds[singlebondId].global_atom_indexes[1] - tinymolParticlesIds[tid][0])
			};
			bondgroups[tid].singlebonds[bondgroups[tid].nSinglebonds++] = SingleBond(idsLocalToTinymol, topology.singlebonds[singlebondId].params );
		}
		for (int anglebondId : availableAnglebondIds) {
			std::array<uint8_t, 3> idsLocalToTinymol = {
				static_cast<uint8_t>(topology.anglebonds[anglebondId].global_atom_indexes[0] - tinymolParticlesIds[tid][0]),
				static_cast<uint8_t>(topology.anglebonds[anglebondId].global_atom_indexes[1] - tinymolParticlesIds[tid][0]),
				static_cast<uint8_t>(topology.anglebonds[anglebondId].global_atom_indexes[2] - tinymolParticlesIds[tid][0])
			};
			bondgroups[tid].anglebonds[bondgroups[tid].nAnglebonds++] = AngleUreyBradleyBond(idsLocalToTinymol, topology.anglebonds[anglebondId].params);
		}
	}

	return bondgroups;
}


std::vector<BondGroupFactory> BondGroupFactory::MakeBondgroups(const LIMA_MOLECULEBUILD::SuperTopology& topology, const std::vector<ParticleToCompoundMapping>& particlesToCompoundIdMap) {

	if (topology.singlebonds.empty()) return {};

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


	std::vector<BondGroupFactory> bondgroups;
	const int expectedNumGroups = static_cast<int>((static_cast<float>(topology.singlebonds.size()) / static_cast<float>(BondGroup::maxSinglebonds)) * 2.f);
	bondgroups.reserve(expectedNumGroups);


	int currentParticleIndexInGroup = 0;

	while (MoreWorkToBeDone(availableSinglebondIds, availablePairbondIds, availableAnglebondIds, availableDihedralbondIds, availableImproperDihedralbondIds)) {

		// Because if we start a new group with a zero-param bond, we skip that bond so this is to avoid empty groups.. 
			//Annoying to deal with here, maybe discard the bonds as we make the topology instead?
		if (bondgroups.empty() || bondgroups.back().nParticles != 0) 
			bondgroups.push_back({});
		currentParticleIndexInGroup = 0;

		// To start or continue a group, simple add the bond containing the next lowest particleId
		const Bondtype typeOfBondWithLowestId = GetBondtypeWithLowestAvailableParticleId(topology, availableSinglebondIds, availablePairbondIds, availableAnglebondIds, availableDihedralbondIds, availableImproperDihedralbondIds);
		switch (typeOfBondWithLowestId) {
		case single:
			bondgroups.back().AddBond(particlesToCompoundIdMap, topology.singlebonds[availableSinglebondIds.front()]);
			availableSinglebondIds.erase(availableSinglebondIds.front());
			break;
		case pair:
			bondgroups.back().AddBond(particlesToCompoundIdMap, topology.pairbonds[availablePairbondIds.front()]);
			availablePairbondIds.erase(availablePairbondIds.front());
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
			
			AddBondsFromMap(bondgroups.back(), pid2SinglebondIdMap[currentParticleId], availableSinglebondIds, topology.singlebonds, particlesToCompoundIdMap);
			AddBondsFromMap(bondgroups.back(), pid2PairbondIdMap[currentParticleId], availablePairbondIds, topology.pairbonds, particlesToCompoundIdMap);
			AddBondsFromMap(bondgroups.back(), pid2AnglebondIdMap[currentParticleId], availableAnglebondIds, topology.anglebonds, particlesToCompoundIdMap);
			AddBondsFromMap(bondgroups.back(), pid2DihedralbondIdMap[currentParticleId], availableDihedralbondIds, topology.dihedralbonds, particlesToCompoundIdMap);
			AddBondsFromMap(bondgroups.back(), pid2ImproperDihedralbondIdMap[currentParticleId], availableImproperDihedralbondIds, topology.improperdihedralbonds, particlesToCompoundIdMap);
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

void BondGroupFactory::AddBond(const ParticleToCompoundMap& particleToCompoundMap, const PairBondFactory& bond) {
	if (bond.params.HasZeroParam())
		return;

	AddBondParticles(particleToCompoundMap, bond.global_atom_indexes);

	if (nPairbonds >= maxPairbonds) {
		throw std::runtime_error("Too many bonds in bondgroup");
	}
	pairbonds[nPairbonds++] = PairBond{ GetLocalIds<PairBond::nAtoms>(particleGlobalToLocalId, bond.global_atom_indexes), bond.params };
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
