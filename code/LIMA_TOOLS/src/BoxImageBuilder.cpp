#include "BoxImageBuilder.h"
#include "Forcefield.h"
#include "MoleculeGraph.h"
#include "TimeIt.h"
#include "ParallelFor.h"

#include <unordered_set>
#include <unordered_map>
#include <format>
#include <array>
#include <numeric>
#include <set>
#include <mutex>

//#include "Display.h"
using namespace LIMA_MOLECULEBUILD;
using namespace LimaMoleculeGraph;




namespace {
	// The bonds of a moleculetype with their parameters, with particle ids relative to the molecule
	struct MoleculetypeBonds {
		std::vector<SingleBondFactory> singlebonds;
		std::vector<PairBondFactory> pairbonds;
		std::vector<AngleBondFactory> anglebonds;
		std::vector<DihedralBondFactory> dihedralbonds;
		std::vector<ImproperDihedralBondFactory> improperdihedralbonds;
	};

	// The atomtypes of a moleculetype, as small ids so we can make cheap lookup keys for bonds
	struct MoleculetypeAtomtypes {
		std::vector<int> atomToType;
		std::vector<std::string> typeNames;	// In order of first appearance

		explicit MoleculetypeAtomtypes(const TopologyFile::Moleculetype& molType) {
			std::unordered_map<std::string, int> typeIds;
			atomToType.reserve(molType.atoms.size());
			for (const auto& atom : molType.atoms) {
				const auto [it, inserted] = typeIds.try_emplace(atom.type, static_cast<int>(typeNames.size()));
				if (inserted)
					typeNames.push_back(atom.type);
				atomToType.push_back(it->second);
			}
		}
	};

	// Parameters come from the topology if present, otherwise from the forcefield. The forcefield is not thread safe,
	// so it is only accessed under the mutex, but most lookups are served by the local cache
	template <typename BondType, typename BondtypeFactory, typename BondTypeTopologyfile>
	void ResolveBonds(const std::vector<BondTypeTopologyfile>& bondsInTopfile, const MoleculetypeAtomtypes& atomtypes,
		LIMAForcefield& forcefield, std::mutex& forcefieldMutex, std::vector<BondtypeFactory>& bonds)
	{
		std::unordered_map<uint64_t, const std::vector<typename BondType::Parameters>*> cache;
		const int nAtoms = static_cast<int>(atomtypes.atomToType.size());

		for (const auto& bondTopol : bondsInTopfile) {
			if (std::ranges::any_of(bondTopol.ids, [nAtoms](int id) { return id >= nAtoms; }))
				continue;

			// In rare cases, the bond parameters are directly in the topology file
			if (bondTopol.parameters.has_value()) {
				bonds.emplace_back(BondtypeFactory{ bondTopol.ids, bondTopol.parameters.value() });
				continue;
			}

			uint64_t key = 0;
			for (const int id : bondTopol.ids)
				key = (key << 16) | static_cast<uint64_t>(atomtypes.atomToType[id]);

			auto cached = cache.find(key);
			if (cached == cache.end()) {
				std::array<std::string, BondType::nAtoms> atomTypenames;
				for (int i = 0; i < BondType::nAtoms; ++i)
					atomTypenames[i] = atomtypes.typeNames[atomtypes.atomToType[bondTopol.ids[i]]];

				std::scoped_lock lock(forcefieldMutex);
				// A bond may be described as multiple bonds, so this is a vector
				cached = cache.emplace(key, &forcefield.GetBondParameters<BondType>(atomTypenames)).first;
			}

			for (const auto& param : *cached->second)
				bonds.emplace_back(BondtypeFactory{ bondTopol.ids, param });
		}
	}

	template <typename BondtypeFactory>
	void InstantiateBonds(const std::vector<BondtypeFactory>& moleculetypeBonds, int particleOffset, BondtypeFactory* out) {
		for (const auto& bond : moleculetypeBonds) {
			*out = bond;
			for (int& id : out->global_atom_indexes)
				id += particleOffset;
			out++;
		}
	}
}

SuperTopology::SuperTopology(const TopologyFile::System& system, const GroFile& grofile, LIMAForcefield& forcefield) {
	moleculeInstances.reserve(system.MoleculeCount());

	// Find the molecules, and the unique moleculetypes in order of appearance
	std::unordered_map<const TopologyFile::Moleculetype*, int> moleculetypeIndices;
	std::vector<const TopologyFile::Moleculetype*> moleculetypes;
	std::vector<int> instanceMoleculetype;
	int nParticles = 0;
	for (const TopologyFile::MoleculeEntry& molecule : system.Instances()) {

#if ENABLE_SOLVENTS != 1
		if (molecule.name == "SOL" || molecule.name == "TIP3") {// TODO: Add the other Solvent labels
			continue;
		}
#endif

		const TopologyFile::Moleculetype& molType = *molecule.moleculetype;
		if (molType.atoms.empty())
			throw std::runtime_error("Molecule has no atoms");

		const auto [it, inserted] = moleculetypeIndices.try_emplace(&molType, static_cast<int>(moleculetypes.size()));
		if (inserted)
			moleculetypes.push_back(&molType);
		instanceMoleculetype.push_back(it->second);

		moleculeInstances.push_back(MoleculeInstance{ molecule.moleculetype, nParticles, static_cast<int>(molType.atoms.size()) });
		nParticles += static_cast<int>(molType.atoms.size());
	}

	std::vector<std::unique_ptr<MoleculetypeAtomtypes>> atomtypes(moleculetypes.size());
	std::vector<MoleculetypeBonds> moleculetypeBonds(moleculetypes.size());
	std::mutex forcefieldMutex;
	ParallelUtils::ParallelFor(moleculetypes.size(), [&](size_t i) {
		const TopologyFile::Moleculetype& molType = *moleculetypes[i];
		atomtypes[i] = std::make_unique<MoleculetypeAtomtypes>(molType);
		const MoleculetypeAtomtypes& types = *atomtypes[i];
		MoleculetypeBonds& bonds = moleculetypeBonds[i];
		ResolveBonds<SingleBond>(molType.singlebonds, types, forcefield, forcefieldMutex, bonds.singlebonds);
		ResolveBonds<PairBond>(molType.pairbonds, types, forcefield, forcefieldMutex, bonds.pairbonds);
		ResolveBonds<AngleUreyBradleyBond>(molType.anglebonds, types, forcefield, forcefieldMutex, bonds.anglebonds);
		ResolveBonds<DihedralBond>(molType.dihedralbonds, types, forcefield, forcefieldMutex, bonds.dihedralbonds);
		ResolveBonds<ImproperDihedralBond>(molType.improperdihedralbonds, types, forcefield, forcefieldMutex, bonds.improperdihedralbonds);
		});

	// The active LJ indices are assigned in the order the types are first requested, so this must be sequential and in order of appearance
	std::vector<std::vector<int>> activeLJParamIndices(moleculetypes.size());
	for (size_t i = 0; i < moleculetypes.size(); i++) {
		std::vector<int> typeIndices;
		for (const std::string& typeName : atomtypes[i]->typeNames)
			typeIndices.push_back(forcefield.GetActiveLjParameterIndex(typeName));
		for (const int type : atomtypes[i]->atomToType)
			activeLJParamIndices[i].push_back(typeIndices[type]);
	}

#if ENABLE_SOLVENTS == 1	// Otherwise solvents are skipped above, and the counts legitimately differ
	// Every particle below indexes the coordinates, so a mismatch would read past their end
	if (nParticles != static_cast<int>(grofile.atoms.size()))
		throw std::runtime_error(std::format("The topology describes {} atoms, but the coordinates contain {}",
			nParticles, grofile.atoms.size()));
#endif

	particles.reserve(nParticles);
	for (size_t instanceId = 0; instanceId < moleculeInstances.size(); instanceId++) {
		const MoleculeInstance& instance = moleculeInstances[instanceId];
		const int moleculetypeIndex = instanceMoleculetype[instanceId];
		const TopologyFile::Moleculetype& molType = *moleculetypes[moleculetypeIndex];
		for (int localId = 0; localId < molType.atoms.size(); localId++) {
			const int indexInGrofile = instance.particleOffset + localId;
			particles.push_back(ParticleFactory{ molType.atoms[localId], grofile.atoms[indexInGrofile].position, indexInGrofile, activeLJParamIndices[moleculetypeIndex][localId] });
		}
	}

	// Each molecule gets a contiguous range of each bondtype
	size_t nSinglebonds = 0, nPairbonds = 0, nAnglebonds = 0, nDihedralbonds = 0, nImproperdihedralbonds = 0;
	for (size_t instanceId = 0; instanceId < moleculeInstances.size(); instanceId++) {
		MoleculeInstance& instance = moleculeInstances[instanceId];
		const MoleculetypeBonds& bonds = moleculetypeBonds[instanceMoleculetype[instanceId]];
		instance.firstSinglebond = static_cast<int>(nSinglebonds);
		instance.firstPairbond = static_cast<int>(nPairbonds);
		instance.firstAnglebond = static_cast<int>(nAnglebonds);
		instance.firstDihedralbond = static_cast<int>(nDihedralbonds);
		instance.firstImproperdihedralbond = static_cast<int>(nImproperdihedralbonds);
		nSinglebonds += bonds.singlebonds.size();
		nPairbonds += bonds.pairbonds.size();
		nAnglebonds += bonds.anglebonds.size();
		nDihedralbonds += bonds.dihedralbonds.size();
		nImproperdihedralbonds += bonds.improperdihedralbonds.size();
	}
	if (nSinglebonds > INT_MAX || nPairbonds > INT_MAX || nAnglebonds > INT_MAX || nDihedralbonds > INT_MAX || nImproperdihedralbonds > INT_MAX)
		throw std::runtime_error("Too many bonds in system");

	singlebonds.resize(nSinglebonds);
	pairbonds.resize(nPairbonds);
	anglebonds.resize(nAnglebonds);
	dihedralbonds.resize(nDihedralbonds);
	improperdihedralbonds.resize(nImproperdihedralbonds);
	ParallelUtils::ParallelForBlocked(moleculeInstances.size(), [&](size_t instanceId) {
		const MoleculeInstance& instance = moleculeInstances[instanceId];
		const MoleculetypeBonds& bonds = moleculetypeBonds[instanceMoleculetype[instanceId]];
		InstantiateBonds(bonds.singlebonds, instance.particleOffset, singlebonds.data() + instance.firstSinglebond);
		InstantiateBonds(bonds.pairbonds, instance.particleOffset, pairbonds.data() + instance.firstPairbond);
		InstantiateBonds(bonds.anglebonds, instance.particleOffset, anglebonds.data() + instance.firstAnglebond);
		InstantiateBonds(bonds.dihedralbonds, instance.particleOffset, dihedralbonds.data() + instance.firstDihedralbond);
		InstantiateBonds(bonds.improperdihedralbonds, instance.particleOffset, improperdihedralbonds.data() + instance.firstImproperdihedralbond);
		}, 64, 64);
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
		std::vector<ParticlesBondedToParticle> particleBondedToParticle;	// Relative to the molecule
		std::vector<PclustersBondedToPcluster> pclusterBondedToPcluster;	// Relative to the molecule
		std::vector<PersistentParticleTemplate> particles;
	};

	bool IsSolventResidue(const std::string& residue) {
		return residue == "SOL" || residue == "SPC" || residue == "SPCE" || residue == "TIP3" || residue == "TIP3P";
	}

	// The singlebond graph of a moleculetype, as flat arrays. Gives the same results as MoleculeGraph for the queries
	// we need here, but without allocating for each query, as we do millions of them for large systems
	class LocalMoleculeGraph {
		static constexpr int maxNeighbors = 8;
		std::vector<std::array<int, maxNeighbors>> neighbors;
		std::vector<int> nNeighbors;

		// Reused between searches
		mutable std::vector<int> visitedStamp;
		mutable int currentStamp = 0;
		mutable std::vector<std::pair<int, int>> queue;	// { nodeId, depth }

		std::span<const int> Neighbors(int id) const { return { neighbors[id].data(), static_cast<size_t>(nNeighbors[id]) }; }

		void Connect(int a, int b) {
			if (nNeighbors[a] >= maxNeighbors || nNeighbors[b] >= maxNeighbors)
				throw std::runtime_error("Exceeded maximum number of neighbors for node " + std::to_string(nNeighbors[a] >= maxNeighbors ? a : b));
			neighbors[a][nNeighbors[a]++] = b;
			neighbors[b][nNeighbors[b]++] = a;
		}

		// Visits nodes in BFS order from start, until f(nodeId, depth) returns false
		template <typename F>
		void BFS(int start, F&& f) const {
			if (++currentStamp == 0) {	// Overflow, reset
				std::ranges::fill(visitedStamp, 0);
				currentStamp = 1;
			}
			queue.clear();
			queue.emplace_back(start, 0);
			visitedStamp[start] = currentStamp;
			for (size_t head = 0; head < queue.size(); head++) {
				const auto [id, depth] = queue[head];
				if (!f(id, depth))
					return;
				for (const int neighbor : Neighbors(id)) {
					if (visitedStamp[neighbor] != currentStamp) {
						visitedStamp[neighbor] = currentStamp;
						queue.emplace_back(neighbor, depth + 1);
					}
				}
			}
		}

	public:
		explicit LocalMoleculeGraph(const TopologyFile::Moleculetype& molecule) {
			const int n = static_cast<int>(molecule.atoms.size());
			for (int i = 0; i < n; i++)
				if (molecule.atoms[i].id != i)
					throw std::runtime_error(std::format("Moleculetype {}: atom at index {} has id {}", molecule.name, i, molecule.atoms[i].id));

			neighbors.resize(n);
			nNeighbors.resize(n, 0);
			visitedStamp.resize(n, 0);
			for (const auto& bond : molecule.singlebonds) {
				if (bond.ids[0] < 0 || bond.ids[0] >= n || bond.ids[1] < 0 || bond.ids[1] >= n)
					continue;
				Connect(bond.ids[0], bond.ids[1]);
			}
		}

		// Connected components in BFS order, starting from the lowest id of each
		std::vector<std::vector<int>> ConnectedComponents() const {
			std::vector<std::vector<int>> components;
			std::vector<bool> visited(neighbors.size(), false);
			for (int id = 0; id < neighbors.size(); id++) {
				if (visited[id])
					continue;
				auto& component = components.emplace_back();
				BFS(id, [&](int node, int) { visited[node] = true; component.push_back(node); return true; });
			}
			return components;
		}

		// Shortest path length, or nullopt if it is longer than maxDistance
		std::optional<int> Distance(int from, int to, int maxDistance) const {
			std::optional<int> result;
			BFS(from, [&](int node, int depth) {
				if (depth > maxDistance)
					return false;
				if (node == to) {
					result = depth;
					return false;
				}
				return true;
				});
			return result;
		}
	};

	std::vector<std::array<int, PersistentCluster::maxParticles>> MakeLocalPersistentClusters(
		const TopologyFile::Moleculetype& molecule
	) {
		const LocalMoleculeGraph moleculeGraph(molecule);
		const std::vector<std::vector<int>> connectedComponents = moleculeGraph.ConnectedComponents();

		std::vector<std::array<int, PersistentCluster::maxParticles>> clusters;
		clusters.reserve((molecule.atoms.size() + PersistentCluster::maxParticles - 1) / PersistentCluster::maxParticles);

		// Only distances up to 4 affect the result
		auto CanAppendToCluster = [&moleculeGraph](int particleId, const auto& cluster, int nextIndex) {
			if (nextIndex == 0)
				return true;

			const std::optional<int> distanceToPreviousNode = moleculeGraph.Distance(cluster[nextIndex - 1], particleId, 4);
			const std::optional<int> distanceToFirstNode = moleculeGraph.Distance(cluster[0], particleId, 4);

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

	void SortAndRemoveDuplicates(std::vector<int>& values) {
		std::ranges::sort(values);
		values.erase(std::unique(values.begin(), values.end()), values.end());
	}

	PersistentClusterTemplate BuildPersistentClusterTemplate(
		const TopologyFile::Moleculetype& molecule,
		const LIMAForcefield& forcefield
	) {
		PersistentClusterTemplate result;
		result.clusters = MakeLocalPersistentClusters(molecule);
		result.particleToPcluster.resize(molecule.atoms.size(), -1);
		result.particles.resize(molecule.atoms.size());

		for (int pcid = 0; pcid < result.clusters.size(); pcid++) {
			for (const int particleId : result.clusters[pcid]) {
				if (particleId != -1)
					result.particleToPcluster[particleId] = pcid;
			}
		}

		// Collect as vectors and sort afterwards, which is much faster than inserting into sets
		std::vector<std::vector<int>> particleBondedToParticle(molecule.atoms.size());
		std::vector<std::vector<int>> pclusterBondedToPcluster(result.clusters.size());
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

					particleBondedToParticle[pidSelf].push_back(pidOther);
					particleBondedToParticle[pidOther].push_back(pidSelf);
					pclusterBondedToPcluster[pcidSelf].push_back(pcidOther);
					pclusterBondedToPcluster[pcidOther].push_back(pcidSelf);
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

		result.particleBondedToParticle.reserve(particleBondedToParticle.size());
		for (auto& values : particleBondedToParticle) {
			SortAndRemoveDuplicates(values);
			result.particleBondedToParticle.push_back(ParticlesBondedToParticle::CreateFromSorted(values));
		}
		result.pclusterBondedToPcluster.reserve(pclusterBondedToPcluster.size());
		for (auto& values : pclusterBondedToPcluster) {
			SortAndRemoveDuplicates(values);
			result.pclusterBondedToPcluster.push_back(PclustersBondedToPcluster::CreateFromSorted(values));
		}

		// Many atoms share a type, so we only look each type up in the forcefield once
		std::unordered_map<std::string, std::pair<NBParams, std::optional<AtomType>>> forcefieldTypes;
		for (int particleId = 0; particleId < molecule.atoms.size(); particleId++) {
			const TopologyFile::AtomsEntry& atom = molecule.atoms[particleId];
			PersistentParticleTemplate& particle = result.particles[particleId];

			auto type = forcefieldTypes.find(atom.type);
			if (type == forcefieldTypes.end())
				type = forcefieldTypes.emplace(atom.type, std::pair{ forcefield.GetLjParameters(atom.type), forcefield.GetAtomtype(atom.type) }).first;

			particle.nbParams = type->second.first;
			if (atom.charge.has_value())
				particle.nbParams.charge = atom.charge.value() * elementaryChargeToKiloCoulombPerMole;

			if (atom.mass.has_value()) {
				particle.mass = atom.mass.value() / KILO;
			}
			else {
				const std::optional<AtomType>& atomType = type->second.second;
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

PersistentClusterFactory MakePersistentClusters(const SuperTopology& system, const LIMAForcefield& forcefield) {
	// Build a template per moleculetype in parallel, then instantiate the templates for each molecule
	std::unordered_map<const TopologyFile::Moleculetype*, int> templateIndices;
	std::vector<const TopologyFile::Moleculetype*> uniqueMoleculetypes;
	std::vector<int> instanceTemplates;
	instanceTemplates.reserve(system.moleculeInstances.size());
	for (const SuperTopology::MoleculeInstance& instance : system.moleculeInstances) {
		const auto [it, inserted] = templateIndices.try_emplace(instance.type.get(), static_cast<int>(uniqueMoleculetypes.size()));
		if (inserted)
			uniqueMoleculetypes.push_back(instance.type.get());
		instanceTemplates.push_back(it->second);
	}

	std::vector<PersistentClusterTemplate> templates(uniqueMoleculetypes.size());
	ParallelUtils::ParallelFor(uniqueMoleculetypes.size(), [&](size_t i) {
		templates[i] = BuildPersistentClusterTemplate(*uniqueMoleculetypes[i], forcefield);
		});

	// Each instance gets a contiguous range of pclusters
	std::vector<int> instancePclusterOffsets(system.moleculeInstances.size() + 1, 0);
	for (size_t instanceId = 0; instanceId < system.moleculeInstances.size(); instanceId++)
		instancePclusterOffsets[instanceId + 1] = instancePclusterOffsets[instanceId] + static_cast<int>(templates[instanceTemplates[instanceId]].clusters.size());
	const size_t totalClusterCount = instancePclusterOffsets.back();

	PersistentClusterFactory pcFactory{};
	pcFactory.pClusters.resize(totalClusterCount);
	pcFactory.pClusterMetas.resize(totalClusterCount);
	pcFactory.particleToPclusterMap.resize(system.particles.size());
	pcFactory.particleBondedToParticle.resize(system.particles.size());
	pcFactory.pclusterBondedToPcluster.resize(totalClusterCount);

	// Instances write to disjoint ranges, so they can be instantiated in parallel
	ParallelUtils::ParallelForBlocked(system.moleculeInstances.size(), [&](size_t instanceId) {
		const SuperTopology::MoleculeInstance& instance = system.moleculeInstances[instanceId];
		const PersistentClusterTemplate& persistentTemplate = templates[instanceTemplates[instanceId]];
		const int pclusterOffset = instancePclusterOffsets[instanceId];

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
			ParticlesBondedToParticle bonded = persistentTemplate.particleBondedToParticle[localParticleId];
			bonded.AddOffset(instance.particleOffset);
			pcFactory.particleBondedToParticle[instance.particleOffset + localParticleId] = bonded;
		}

		for (int localPcid = 0; localPcid < persistentTemplate.pclusterBondedToPcluster.size(); localPcid++) {
			PclustersBondedToPcluster bonded = persistentTemplate.pclusterBondedToPcluster[localPcid];
			bonded.AddOffset(pclusterOffset);
			pcFactory.pclusterBondedToPcluster[pclusterOffset + localPcid] = bonded;
		}
		}, 64, 16);

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
	std::future<PersistentClusterFactory> pcFactoryFuture = std::async(std::launch::async, [&]{ return MakePersistentClusters(superTopology, forcefield); });

	BondGroupFactory bgFactory = bgfFuture.get();


	PersistentClusterFactory pcFactory = pcFactoryFuture.get();
	bgFactory.AddPclusterRefs(pcFactory.particleToPclusterMap);

	bgFactory.AddBondgroupRefsToPclusters(pcFactory.particleToPclusterMap, pcFactory.pClusterMetas);
	BondGroups bondGroups = bgFactory.GetBondgroups();



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
		std::move(bondGroups),
		std::move(pcFactory.pClusters),
		std::move(pcFactory.pClusterMetas),
		std::move(pcFactory.particleBondedToParticle),
		std::move(pcFactory.pclusterBondedToPcluster),
		std::move(gpidToPcidAndPid),
		nParticles
	);
}
