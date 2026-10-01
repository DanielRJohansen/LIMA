#include "BoxImageBuilder.h"
#include "ParallelFor.h"

#include <algorithm>


namespace {
	using LIMA_MOLECULEBUILD::SuperTopology;

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
		RandomAccessDeleteSet(size_t nElements)
			: nElements(nElements), deleted(nElements, false), first_valid_index(0) {}

		// Returns the first not deleted element.
		int front() const {
			if (first_valid_index >= nElements) {
				throw std::out_of_range("No valid elements in the queue");
			}
			return static_cast<int>(first_valid_index);
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


	// Maps each particle in a range to the ids of the bonds it is in, in ascending order
	class ParticleToBondsMap {
		int firstParticle = 0;
		std::vector<int> offsets;
		std::vector<int> bondIds;
	public:
		template <typename BondFactoryType>
		ParticleToBondsMap(std::span<const BondFactoryType> bonds, int firstParticle, int nParticles) : firstParticle(firstParticle) {
			offsets.assign(nParticles + 1, 0);
			for (const auto& bond : bonds) {
				for (const int particleId : bond.global_atom_indexes) {
					if (particleId < firstParticle || particleId >= firstParticle + nParticles)
						throw std::runtime_error("Bond references a particle outside its molecule");
					offsets[particleId - firstParticle + 1]++;
				}
			}
			for (int i = 0; i < nParticles; i++)
				offsets[i + 1] += offsets[i];

			bondIds.resize(offsets.back());
			std::vector<int> cursors(offsets.begin(), offsets.end() - 1);
			for (int bondId = 0; bondId < bonds.size(); bondId++)
				for (const int particleId : bonds[bondId].global_atom_indexes)
					bondIds[cursors[particleId - firstParticle]++] = bondId;
		}

		std::span<const int> operator[](int particleId) const {
			const int local = particleId - firstParticle;
			return { bondIds.data() + offsets[local], static_cast<size_t>(offsets[local + 1] - offsets[local]) };
		}
	};


	enum Bondtype { single, pair, angle, dihedral, improper };

	// The bonds of one or more consecutive molecules
	struct BondsChunk {
		int firstParticle = 0;
		int nParticles = 0;
		std::span<const SingleBondFactory> singlebonds;
		std::span<const PairBondFactory> pairbonds;
		std::span<const AngleBondFactory> anglebonds;
		std::span<const DihedralBondFactory> dihedralbonds;
		std::span<const ImproperDihedralBondFactory> improperdihedralbonds;
	};

	// Greedily groups bonds: Start a group with the bond containing the lowest available particle id, then add all bonds
	// containing the particles of the group, as long as there is room. A molecule is always finished before the next is started,
	// so chunks of whole molecules can be grouped independently, and the results concatenated.
	class BondgroupChunkBuilder {
		static constexpr int maxParticlesPerBondgroup = 64;
		const BondsChunk& chunk;

	public:
		BondGroups bondgroups;
		std::vector<int> particleGlobalIds;

		explicit BondgroupChunkBuilder(const BondsChunk& chunk) : chunk(chunk) {
			const ParticleToBondsMap pid2SinglebondIdMap(chunk.singlebonds, chunk.firstParticle, chunk.nParticles);
			const ParticleToBondsMap pid2PairbondIdMap(chunk.pairbonds, chunk.firstParticle, chunk.nParticles);
			const ParticleToBondsMap pid2AnglebondIdMap(chunk.anglebonds, chunk.firstParticle, chunk.nParticles);
			const ParticleToBondsMap pid2DihedralbondIdMap(chunk.dihedralbonds, chunk.firstParticle, chunk.nParticles);
			const ParticleToBondsMap pid2ImproperDihedralbondIdMap(chunk.improperdihedralbonds, chunk.firstParticle, chunk.nParticles);

			RandomAccessDeleteSet availableSinglebondIds(chunk.singlebonds.size());
			RandomAccessDeleteSet availablePairbondIds(chunk.pairbonds.size());
			RandomAccessDeleteSet availableAnglebondIds(chunk.anglebonds.size());
			RandomAccessDeleteSet availableDihedralbondIds(chunk.dihedralbonds.size());
			RandomAccessDeleteSet availableImproperDihedralbondIds(chunk.improperdihedralbonds.size());

			auto MoreWorkToBeDone = [&]() {
				return !availableSinglebondIds.empty() || !availablePairbondIds.empty() || !availableAnglebondIds.empty() || !availableDihedralbondIds.empty() || !availableImproperDihedralbondIds.empty();
				};

			auto GetBondtypeWithLowestAvailableParticleId = [&]() {
				int minParticleId = std::numeric_limits<int>::max();
				Bondtype type = single;
				auto Check = [&](const RandomAccessDeleteSet& available, const auto& bonds, Bondtype bondtype) {
					if (available.empty())
						return;
					const auto& ids = bonds[available.front()].global_atom_indexes;
					const int minId = *std::min_element(ids.begin(), ids.end());
					if (minId < minParticleId) {
						minParticleId = minId;
						type = bondtype;
					}
					};
				Check(availableSinglebondIds, chunk.singlebonds, single);
				Check(availablePairbondIds, chunk.pairbonds, pair);
				Check(availableAnglebondIds, chunk.anglebonds, angle);
				Check(availableDihedralbondIds, chunk.dihedralbonds, dihedral);
				Check(availableImproperDihedralbondIds, chunk.improperdihedralbonds, improper);
				return type;
				};

			auto AddBondsFromMap = [this](std::span<const int> bondIds, RandomAccessDeleteSet& availableBondIds, const auto& bonds) {
				for (const int bondId : bondIds) {
					if (availableBondIds.contains(bondId)) {
						if (AddBond(bonds[bondId]))
							availableBondIds.erase(bondId);
					}
				}
				};

			while (MoreWorkToBeDone()) {
				// To start or continue a group, simple add the bond containing the next lowest particleId
				const Bondtype typeOfBondWithLowestId = GetBondtypeWithLowestAvailableParticleId();

				// Because if we start a new group with a zero-param bond, we skip that bond so this is to avoid empty groups..
				if (bondgroups.empty() || bondgroups.groups.back().nParticles != 0) {
					bondgroups.groups.push_back({
						static_cast<int>(bondgroups.particles.size()), 0,
						static_cast<int>(bondgroups.singlebonds.size()), 0,
						static_cast<int>(bondgroups.pairbonds.size()), 0,
						static_cast<int>(bondgroups.anglebonds.size()), 0,
						static_cast<int>(bondgroups.dihedralbonds.size()), 0,
						static_cast<int>(bondgroups.improperdihedralbonds.size()), 0
						});
				}

				switch (typeOfBondWithLowestId) {
				case single:
					AddBond(chunk.singlebonds[availableSinglebondIds.front()]);
					availableSinglebondIds.erase(availableSinglebondIds.front());
					break;
				case pair:
					AddBond(chunk.pairbonds[availablePairbondIds.front()]);
					availablePairbondIds.erase(availablePairbondIds.front());
					break;
				case angle:
					AddBond(chunk.anglebonds[availableAnglebondIds.front()]);
					availableAnglebondIds.erase(availableAnglebondIds.front());
					break;
				case dihedral:
					AddBond(chunk.dihedralbonds[availableDihedralbondIds.front()]);
					availableDihedralbondIds.erase(availableDihedralbondIds.front());
					break;
				case improper:
					AddBond(chunk.improperdihedralbonds[availableImproperDihedralbondIds.front()]);
					availableImproperDihedralbondIds.erase(availableImproperDihedralbondIds.front());
					break;
				}

				// Now go from the current index of particles in the group to the last.
				// For each particle find all bonds that contain said particle, and add those bonds to this group IF we have room
				// Then move to the next particle and repeat
				// We exit when, either we have no more bonds in the chain, or the group has no more room
				for (int currentParticleIndexInGroup = 0; currentParticleIndexInGroup < bondgroups.groups.back().nParticles; currentParticleIndexInGroup++) {
					const int currentParticleId = particleGlobalIds[bondgroups.groups.back().indexOfFirstParticle + currentParticleIndexInGroup];

					AddBondsFromMap(pid2SinglebondIdMap[currentParticleId], availableSinglebondIds, chunk.singlebonds);
					AddBondsFromMap(pid2PairbondIdMap[currentParticleId], availablePairbondIds, chunk.pairbonds);
					AddBondsFromMap(pid2AnglebondIdMap[currentParticleId], availableAnglebondIds, chunk.anglebonds);
					AddBondsFromMap(pid2DihedralbondIdMap[currentParticleId], availableDihedralbondIds, chunk.dihedralbonds);
					AddBondsFromMap(pid2ImproperDihedralbondIdMap[currentParticleId], availableImproperDihedralbondIds, chunk.improperdihedralbonds);
				}
			}
		}

	private:
		int FindLocalParticleId(const BondGroup& group, int globalId) const {
			for (int i = 0; i < group.nParticles; i++) {
				if (particleGlobalIds[group.indexOfFirstParticle + i] == globalId)
					return i;
			}
			return -1;
		}

		std::vector<SingleBond>& Output(const SingleBondFactory&) { return bondgroups.singlebonds; }
		std::vector<PairBond>& Output(const PairBondFactory&) { return bondgroups.pairbonds; }
		std::vector<AngleUreyBradleyBond>& Output(const AngleBondFactory&) { return bondgroups.anglebonds; }
		std::vector<DihedralBond>& Output(const DihedralBondFactory&) { return bondgroups.dihedralbonds; }
		std::vector<ImproperDihedralBond>& Output(const ImproperDihedralBondFactory&) { return bondgroups.improperdihedralbonds; }

		static int& Count(BondGroup& group, const SingleBondFactory&) { return group.nSinglebonds; }
		static int& Count(BondGroup& group, const PairBondFactory&) { return group.nPairbonds; }
		static int& Count(BondGroup& group, const AngleBondFactory&) { return group.nAnglebonds; }
		static int& Count(BondGroup& group, const DihedralBondFactory&) { return group.nDihedralbonds; }
		static int& Count(BondGroup& group, const ImproperDihedralBondFactory&) { return group.nImproperdihedralbonds; }

		// Adds the bond to the last group. Returns false if there is no room for its particles
		template <typename BondFactoryType>
		bool AddBond(const BondFactoryType& bond) {
			if (bond.params.HasZeroParam())
				return true;

			BondGroup& group = bondgroups.groups.back();
			constexpr int n = BondFactoryType::nAtoms;
			const auto& globalIds = bond.global_atom_indexes;

			std::array<uint8_t, n> localIds;
			int nNewParticles = 0;
			for (int i = 0; i < n; i++) {
				const int localId = FindLocalParticleId(group, globalIds[i]);
				localIds[i] = localId == -1 ? group.nParticles + nNewParticles++ : localId;
			}
			if (group.nParticles + nNewParticles > maxParticlesPerBondgroup)
				return false;

			assert(globalIds.front() != globalIds.back());
			for (int i = 0; i < n; i++) {
				if (localIds[i] >= group.nParticles) {
					particleGlobalIds.push_back(globalIds[i]);
					bondgroups.particles.push_back({});
					group.nParticles++;
				}
			}

			Output(bond).emplace_back(localIds, bond.params);
			Count(group, bond)++;
			return true;
		}
	};


	// Splits the molecules into chunks with roughly the same number of bonds
	std::vector<BondsChunk> MakeChunks(const SuperTopology& topology) {
		const auto& instances = topology.moleculeInstances;
		const size_t totalBonds = topology.singlebonds.size() + topology.pairbonds.size() + topology.anglebonds.size()
			+ topology.dihedralbonds.size() + topology.improperdihedralbonds.size();
		const size_t targetBondsPerChunk = std::max<size_t>(totalBonds / 256, 1);

		auto FirstBonds = [&](size_t instanceId) {
			if (instanceId == instances.size())
				return std::array<size_t, 5>{ topology.singlebonds.size(), topology.pairbonds.size(), topology.anglebonds.size(), topology.dihedralbonds.size(), topology.improperdihedralbonds.size() };
			const auto& instance = instances[instanceId];
			return std::array<size_t, 5>{ size_t(instance.firstSinglebond), size_t(instance.firstPairbond), size_t(instance.firstAnglebond), size_t(instance.firstDihedralbond), size_t(instance.firstImproperdihedralbond) };
			};
		auto FirstParticle = [&](size_t instanceId) {
			return instanceId == instances.size() ? static_cast<int>(topology.particles.size()) : instances[instanceId].particleOffset;
			};

		std::vector<BondsChunk> chunks;
		size_t chunkStart = 0;
		for (size_t instanceId = 1; instanceId <= instances.size(); instanceId++) {
			const auto begin = FirstBonds(chunkStart);
			const auto end = FirstBonds(instanceId);
			size_t nBonds = 0;
			for (int i = 0; i < 5; i++) nBonds += end[i] - begin[i];
			if (nBonds < targetBondsPerChunk && instanceId != instances.size())
				continue;

			BondsChunk chunk;
			chunk.firstParticle = FirstParticle(chunkStart);
			chunk.nParticles = FirstParticle(instanceId) - chunk.firstParticle;
			chunk.singlebonds = std::span(topology.singlebonds).subspan(begin[0], end[0] - begin[0]);
			chunk.pairbonds = std::span(topology.pairbonds).subspan(begin[1], end[1] - begin[1]);
			chunk.anglebonds = std::span(topology.anglebonds).subspan(begin[2], end[2] - begin[2]);
			chunk.dihedralbonds = std::span(topology.dihedralbonds).subspan(begin[3], end[3] - begin[3]);
			chunk.improperdihedralbonds = std::span(topology.improperdihedralbonds).subspan(begin[4], end[4] - begin[4]);
			chunks.push_back(chunk);
			chunkStart = instanceId;
		}
		return chunks;
	}

	template <typename T>
	void Append(std::vector<T>& dst, const std::vector<T>& src) {
		dst.insert(dst.end(), src.begin(), src.end());
	}
}


BondGroupFactory::BondGroupFactory(const LIMA_MOLECULEBUILD::SuperTopology& topology) {
	if (topology.singlebonds.empty()) return;

	const std::vector<BondsChunk> chunks = MakeChunks(topology);
	std::vector<std::unique_ptr<BondgroupChunkBuilder>> results(chunks.size());
	ParallelUtils::ParallelFor(chunks.size(), [&](size_t i) {
		results[i] = std::make_unique<BondgroupChunkBuilder>(chunks[i]);
		});

	// Concatenate the chunks. The sequential algorithm reuses an empty last group for the next molecule,
	// so an empty group only remains if it is the very last
	size_t nGroups = 0;
	for (const auto& result : results) nGroups += result->bondgroups.size();
	bondgroups.groups.reserve(nGroups);

	for (const auto& result : results) {
		const BondGroups& chunkGroups = result->bondgroups;
		if (!bondgroups.groups.empty() && bondgroups.groups.back().nParticles == 0 && !chunkGroups.empty())
			bondgroups.groups.pop_back();

		for (BondGroup group : chunkGroups.groups) {
			group.indexOfFirstParticle += static_cast<int>(bondgroups.particles.size());
			group.indexOfFirstSinglebond += static_cast<int>(bondgroups.singlebonds.size());
			group.indexOfFirstPairbond += static_cast<int>(bondgroups.pairbonds.size());
			group.indexOfFirstAnglebond += static_cast<int>(bondgroups.anglebonds.size());
			group.indexOfFirstDihedralbond += static_cast<int>(bondgroups.dihedralbonds.size());
			group.indexOfFirstImproperdihedralbond += static_cast<int>(bondgroups.improperdihedralbonds.size());
			bondgroups.groups.push_back(group);
		}
		Append(bondgroups.particles, chunkGroups.particles);
		Append(bondgroups.singlebonds, chunkGroups.singlebonds);
		Append(bondgroups.pairbonds, chunkGroups.pairbonds);
		Append(bondgroups.anglebonds, chunkGroups.anglebonds);
		Append(bondgroups.dihedralbonds, chunkGroups.dihedralbonds);
		Append(bondgroups.improperdihedralbonds, chunkGroups.improperdihedralbonds);
		Append(particleGlobalIds, result->particleGlobalIds);
	}
}

BondGroups BondGroupFactory::GetBondgroups() {
	// TODO Add safety checks here that addpclusterinfo has been called, and that this hasnt been called before.
	return std::move(bondgroups);
}

void BondGroupFactory::AddPclusterRefs(const ParticleToPclusterMap& particleToPclusterMap) {
	ParallelUtils::ParallelForBlocked(bondgroups.particles.size(), [&](size_t i) {
		const ParticleToPclusterMapping& mapping = particleToPclusterMap[particleGlobalIds[i]];
		bondgroups.particles[i] = BondGroup::ParticleRef{ mapping.pcid, mapping.pid };
		});
}

void BondGroupFactory::AddBondgroupRefsToPclusters(const ParticleToPclusterMap& particleToPclusterMap, std::vector<PersistentClusterMeta>& pclusterMetas) const {
	for (int bgIndex = 0; bgIndex < bondgroups.size(); bgIndex++) {
		const BondGroup& group = bondgroups.groups[bgIndex];
		for (int particleLocalId = 0; particleLocalId < group.nParticles; particleLocalId++) {
			const int particleIndex = group.indexOfFirstParticle + particleLocalId;
			const ParticleToPclusterMapping& pcRef = particleToPclusterMap[particleGlobalIds[particleIndex]];
			pclusterMetas[pcRef.pcid].bondgroupReferences[pcRef.pid].Add({ bgIndex, particleIndex });
		}
	}
}
