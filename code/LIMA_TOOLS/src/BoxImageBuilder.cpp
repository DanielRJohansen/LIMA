#include "BoxImageBuilder.h"
#include "Format.h"
#include "Forcefield.h"
#include "MoleculeGraph.h"
#include "TimeIt.h"
#include "ParallelFor.h"
#include "RenderUtilities.cuh"

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
		std::vector<CmapBondFactory> cmapbonds;
	};

	bool IsHydrogen(const TopologyFile::AtomsEntry& atom) {
		return !atom.atomname.empty() && std::toupper(static_cast<unsigned char>(atom.atomname[0])) == 'H';	// As GROMACS
	}

	// The pairs of a moleculetype the nonbonded kernel excludes, as GROMACS: those at most nrexcl bonds apart. Not those sharing an
	// angle or dihedral, as a topology may list only some of them, eg. united-atom lipids one dihedral per torsion.
	// Ids relative to the molecule, each pair once with first < second, sorted
	std::vector<std::pair<int, int>> ExcludedPairs(const TopologyFile::Moleculetype& molType) {
		const int nAtoms = static_cast<int>(molType.atoms.size());
		std::vector<std::vector<int>> neighbors(nAtoms);
		for (const auto& bond : molType.singlebonds) {
			const auto [a, b] = bond.ids;
			if (a < 0 || b < 0 || a >= nAtoms || b >= nAtoms || a == b)
				continue;
			neighbors[a].push_back(b);
			neighbors[b].push_back(a);
		}

		std::vector<std::pair<int, int>> pairs;
		std::vector<int> distance(nAtoms, -1);
		std::vector<int> visited, frontier, next;
		for (int start = 0; start < nAtoms; start++) {
			// Breadth first, nrexcl bonds deep. Only the atoms reached are reset, so each search costs its neighborhood
			distance[start] = 0;
			visited = { start };
			frontier = { start };
			for (int depth = 1; depth <= molType.nrexcl && !frontier.empty(); depth++) {
				next.clear();
				for (const int atom : frontier)
					for (const int neighbor : neighbors[atom])
						if (distance[neighbor] == -1) {
							distance[neighbor] = depth;
							visited.push_back(neighbor);
							next.push_back(neighbor);
							if (neighbor > start)
								pairs.emplace_back(start, neighbor);
						}
				std::swap(frontier, next);
			}
			for (const int atom : visited)
				distance[atom] = -1;
		}
		std::ranges::sort(pairs);
		return pairs;
	}

	// The constraints of a moleculetype: its settled waters ([ settles ]) always, and with constrainHBonds also its bonds to
	// hydrogens at their b0, and in a 3-atom flexible water the distance between the hydrogens at the angle's theta0.
	// The constrained bonds and angles stay: at the constrained lengths their energy is 0 and the constraints absorb their forces,
	// so MD is unaffected, while energy minimization, which ignores the constraints, can use the same box
	MoleculetypeConstraints MakeConstraints(const TopologyFile::Moleculetype& molType, const MoleculetypeBonds& bonds, bool constrainHBonds) {
		MoleculetypeConstraints constraints;
		std::vector<bool> settled(molType.atoms.size(), false);
		for (const auto& settle : molType.settles) {
			const int o = settle.oxygen;
			constraints.pairs.push_back({ o, o + 1 });
			constraints.pairs.push_back({ o, o + 2 });
			constraints.pairs.push_back({ o + 1, o + 2 });
			constraints.distances.insert(constraints.distances.end(), { settle.distanceOH, settle.distanceOH, settle.distanceHH });
			settled[o] = settled[o + 1] = settled[o + 2] = true;
		}
		if (constrainHBonds) {
			for (const auto& bond : bonds.singlebonds) {
				const auto [a, b] = bond.global_atom_indexes;
				if ((IsHydrogen(molType.atoms[a]) || IsHydrogen(molType.atoms[b])) && !settled[a] && !settled[b]) {
					constraints.pairs.push_back({ a, b });
					constraints.distances.push_back(bond.params.b0);
				}
			}
		}
		if (constraints.pairs.empty())
			return constraints;

		if (molType.settles.empty() && molType.atoms.size() == 3 && constraints.pairs.size() == 2) {
			for (size_t i = 0; i < bonds.anglebonds.size(); i++) {
				const auto [h1, o, h2] = bonds.anglebonds[i].global_atom_indexes;
				if (!IsHydrogen(molType.atoms[h1]) || !IsHydrogen(molType.atoms[h2]))
					continue;
				const auto Distance = [&](int h) {
					for (size_t c = 0; c < constraints.pairs.size(); c++)
						if (constraints.pairs[c] == std::array{ o, h } || constraints.pairs[c] == std::array{ h, o })
							return constraints.distances[c];
					throw std::runtime_error(Lima::Format("Moleculetype {}: water angle without matching O-H bonds", molType.name));
				};
				const float d1 = Distance(h1), d2 = Distance(h2);
				const float theta = bonds.anglebonds[i].params.theta0;
				constraints.pairs.push_back({ h1, h2 });
				constraints.distances.push_back(std::sqrt(d1 * d1 + d2 * d2 - 2.f * d1 * d2 * std::cos(theta)));
				break;
			}
		}

		// Connected groups, which must each fit in a pcluster
		const int nAtoms = static_cast<int>(molType.atoms.size());
		std::vector<int> root(nAtoms);
		for (int i = 0; i < nAtoms; i++) root[i] = i;
		const auto Find = [&](int i) { while (root[i] != i) i = root[i] = root[root[i]]; return i; };
		for (const auto& [a, b] : constraints.pairs)
			root[Find(a)] = Find(b);
		constraints.groupOf.assign(nAtoms, -1);
		std::vector<int> groupOfRoot(nAtoms, -1);
		std::vector<int> nGroupConstraints;
		for (const auto& [a, b] : constraints.pairs) {
			for (const int atom : { a, b }) {
				int& group = groupOfRoot[Find(atom)];
				if (group == -1) {
					group = static_cast<int>(constraints.groups.size());
					constraints.groups.emplace_back();
					nGroupConstraints.push_back(0);
				}
				if (constraints.groupOf[atom] == -1) {
					constraints.groupOf[atom] = group;
					constraints.groups[group].push_back(atom);
				}
			}
			nGroupConstraints[groupOfRoot[Find(a)]]++;
		}
		for (size_t g = 0; g < constraints.groups.size(); g++) {
			std::ranges::sort(constraints.groups[g]);
			if (constraints.groups[g].size() > PersistentCluster::maxParticles || nGroupConstraints[g] > PclusterConstraints::maxConstraints)
				throw std::runtime_error(Lima::Format("Moleculetype {}: a group of {} constrained atoms around atom {} ({}), at most {} atoms "
					"and {} constraints are supported (eg. not CH4)", molType.name, constraints.groups[g].size(), constraints.groups[g][0] + 1,
					molType.atoms[constraints.groups[g][0]].atomname, PersistentCluster::maxParticles, PclusterConstraints::maxConstraints));
		}
		return constraints;
	}

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
		std::unordered_map<uint64_t, const std::vector<typename BondType::Parameters>*> periodicImproperCache;	// Dihedrals of funct 4
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

			// GROMACS matches periodic impropers (funct 4) only against their own [ dihedraltypes ]
			bool periodicImproper = false;
			if constexpr (std::is_same_v<BondType, DihedralBond>)
				periodicImproper = bondTopol.funct == 4;
			auto& bondCache = periodicImproper ? periodicImproperCache : cache;
			auto cached = bondCache.find(key);
			if (cached == bondCache.end()) {
				std::array<std::string, BondType::nAtoms> atomTypenames;
				for (int i = 0; i < BondType::nAtoms; ++i)
					atomTypenames[i] = atomtypes.typeNames[atomtypes.atomToType[bondTopol.ids[i]]];

				std::scoped_lock lock(forcefieldMutex);
				// A bond may be described as multiple bonds, so this is a vector
				const std::vector<typename BondType::Parameters>* parameters = nullptr;
				if constexpr (std::is_same_v<BondType, DihedralBond>)
					if (periodicImproper)
						parameters = &forcefield.GetPeriodicImproperParameters(atomTypenames);
				if (!parameters)
					parameters = &forcefield.GetBondParameters<BondType>(atomTypenames);
				cached = bondCache.emplace(key, parameters).first;
			}

			for (const auto& param : *cached->second)
				bonds.emplace_back(BondtypeFactory{ bondTopol.ids, param });
		}
	}

	// CMAP parameters are grids in the forcefield, never in the topology
	void ResolveCmaps(const std::vector<TopologyFile::CmapBond>& cmapsInTopfile, const MoleculetypeAtomtypes& atomtypes,
		const LIMAForcefield& forcefield, std::mutex& forcefieldMutex, std::vector<CmapBondFactory>& cmaps)
	{
		const int nAtoms = static_cast<int>(atomtypes.atomToType.size());
		for (const auto& cmapTopol : cmapsInTopfile) {
			if (std::ranges::any_of(cmapTopol.ids, [nAtoms](int id) { return id >= nAtoms; }))
				continue;
			std::array<std::string, CmapBond::nAtoms> atomTypenames;
			std::array<int, CmapBond::nAtoms> ids;
			for (int i = 0; i < CmapBond::nAtoms; ++i) {
				ids[i] = cmapTopol.ids[i];
				atomTypenames[i] = atomtypes.typeNames[atomtypes.atomToType[ids[i]]];
			}
			std::optional<int> gridIndex;
			{
				std::scoped_lock lock(forcefieldMutex);
				gridIndex = forcefield.GetCmapGridIndex(atomTypenames);
			}
			if (!gridIndex)
				throw std::runtime_error(Lima::Format("No cmaptypes entry for {} {} {} {} {}", atomTypenames[0], atomTypenames[1],
					atomTypenames[2], atomTypenames[3], atomTypenames[4]));
			cmaps.emplace_back(CmapBondFactory{ ids, CmapBond::Parameters{ *gridIndex } });
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

SuperTopology::SuperTopology(const TopologyFile::System& system, const GroFile& grofile, LIMAForcefield& forcefield, ConstraintsSelect constraintsSelect,
	bool electrostatics) {
	moleculeInstances.reserve(system.MoleculeCount());

	// Find the molecules, and the unique moleculetypes in order of appearance
	std::unordered_map<const TopologyFile::Moleculetype*, int> moleculetypeIndices;
	std::vector<const TopologyFile::Moleculetype*> moleculetypes;
	std::vector<int> instanceMoleculetype;
	int nParticles = 0;
	for (const TopologyFile::MoleculeEntry& molecule : system.Instances()) {
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
	std::vector<MoleculetypeConstraints> moleculetypeConstraints(moleculetypes.size());
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
		ResolveCmaps(molType.cmapbonds, types, forcefield, forcefieldMutex, bonds.cmapbonds);
		if (electrostatics) {
			// After all bonded terms are resolved, as the exclusions come from them. The 1-4 Coulomb, scaled by fudgeQQ. A charge
			// missing in the topology comes from the atomtype, as for the nonbonded charges
			std::vector<float> charges(molType.atoms.size());	// [kC/mol]
			{
				std::scoped_lock lock(forcefieldMutex);
				for (size_t a = 0; a < molType.atoms.size(); a++) {
					const auto& atom = molType.atoms[a];
					const std::optional<AtomType> atomtype = atom.charge.has_value() ? std::nullopt : forcefield.GetAtomtype(atom.type);
					charges[a] = atom.charge.has_value() ? atom.charge.value() * elementaryChargeToKiloCoulombPerMole
						: atomtype.has_value() ? atomtype->charge : 0.f;
				}
			}
			for (auto& pair : bonds.pairbonds)
				pair.params.chargeProduct = forcefield.FudgeQQ() * charges[pair.global_atom_indexes[0]] * charges[pair.global_atom_indexes[1]];

			// The Ewald exclusion correction, for every pair the nonbonded kernel excludes (ExcludedPairs, as in
			// MakePersistentClusters). 1-4 pairs carry it in their entry, the others get an entry of their own
			std::vector<uint64_t> excluded;
			for (const auto [i, j] : ExcludedPairs(molType))
				excluded.push_back(uint64_t(i) << 32 | uint32_t(j));
			std::unordered_map<uint64_t, size_t> pairIndices;
			for (size_t p = 0; p < bonds.pairbonds.size(); p++) {
				const auto [i, j] = bonds.pairbonds[p].global_atom_indexes;
				pairIndices.emplace(uint64_t(std::min(i, j)) << 32 | uint32_t(std::max(i, j)), p);
			}
			for (const uint64_t key : excluded) {
				const int i = static_cast<int>(key >> 32), j = static_cast<int>(key & 0xFFFFFFFFu);
				const float exclusionChargeProduct = charges[i] * charges[j];
				if (exclusionChargeProduct == 0.f)
					continue;
				if (const auto found = pairIndices.find(key); found != pairIndices.end())
					bonds.pairbonds[found->second].params.exclusionChargeProduct = exclusionChargeProduct;
				else
					bonds.pairbonds.emplace_back(PairBondFactory{ { i, j }, PairBond::Parameters{ 0.f, 0.f, 0.f, exclusionChargeProduct } });
			}
		}
		if (constraintsSelect == ConstraintsSelect::HBonds || !molType.settles.empty())
			moleculetypeConstraints[i] = MakeConstraints(molType, bonds, constraintsSelect == ConstraintsSelect::HBonds);
		});
	for (size_t i = 0; i < moleculetypes.size(); i++)
		if (!moleculetypeConstraints[i].pairs.empty())
			constraints.emplace(moleculetypes[i], std::move(moleculetypeConstraints[i]));
	for (const int moleculetype : instanceMoleculetype) {
		const auto found = constraints.find(moleculetypes[moleculetype]);
		if (found != constraints.end())
			nConstraints += static_cast<int64_t>(found->second.pairs.size());
	}

	// The active LJ indices are assigned in the order the types are first requested, so this must be sequential and in order of appearance
	std::vector<std::vector<int>> activeLJParamIndices(moleculetypes.size());
	for (size_t i = 0; i < moleculetypes.size(); i++) {
		std::vector<int> typeIndices;
		for (const std::string& typeName : atomtypes[i]->typeNames)
			typeIndices.push_back(forcefield.GetActiveLjParameterIndex(typeName));
		for (const int type : atomtypes[i]->atomToType)
			activeLJParamIndices[i].push_back(typeIndices[type]);
	}

	// Every particle below indexes the coordinates, so a mismatch would read past their end
	if (nParticles != static_cast<int>(grofile.atoms.size()))
		throw std::runtime_error(Lima::Format("The topology describes {} atoms, but the coordinates contain {}",
			nParticles, grofile.atoms.size()));

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
	size_t nSinglebonds = 0, nPairbonds = 0, nAnglebonds = 0, nDihedralbonds = 0, nImproperdihedralbonds = 0, nCmapbonds = 0;
	for (size_t instanceId = 0; instanceId < moleculeInstances.size(); instanceId++) {
		MoleculeInstance& instance = moleculeInstances[instanceId];
		const MoleculetypeBonds& bonds = moleculetypeBonds[instanceMoleculetype[instanceId]];
		instance.firstSinglebond = static_cast<int>(nSinglebonds);
		instance.firstPairbond = static_cast<int>(nPairbonds);
		instance.firstAnglebond = static_cast<int>(nAnglebonds);
		instance.firstDihedralbond = static_cast<int>(nDihedralbonds);
		instance.firstImproperdihedralbond = static_cast<int>(nImproperdihedralbonds);
		instance.firstCmapbond = static_cast<int>(nCmapbonds);
		nSinglebonds += bonds.singlebonds.size();
		nPairbonds += bonds.pairbonds.size();
		nAnglebonds += bonds.anglebonds.size();
		nDihedralbonds += bonds.dihedralbonds.size();
		nImproperdihedralbonds += bonds.improperdihedralbonds.size();
		nCmapbonds += bonds.cmapbonds.size();
	}
	if (nSinglebonds > INT_MAX || nPairbonds > INT_MAX || nAnglebonds > INT_MAX || nDihedralbonds > INT_MAX || nImproperdihedralbonds > INT_MAX || nCmapbonds > INT_MAX)
		throw std::runtime_error("Too many bonds in system");

	singlebonds.resize(nSinglebonds);
	pairbonds.resize(nPairbonds);
	anglebonds.resize(nAnglebonds);
	dihedralbonds.resize(nDihedralbonds);
	improperdihedralbonds.resize(nImproperdihedralbonds);
	cmapbonds.resize(nCmapbonds);
	ParallelUtils::ParallelForBlocked(moleculeInstances.size(), [&](size_t instanceId) {
		const MoleculeInstance& instance = moleculeInstances[instanceId];
		const MoleculetypeBonds& bonds = moleculetypeBonds[instanceMoleculetype[instanceId]];
		InstantiateBonds(bonds.singlebonds, instance.particleOffset, singlebonds.data() + instance.firstSinglebond);
		InstantiateBonds(bonds.pairbonds, instance.particleOffset, pairbonds.data() + instance.firstPairbond);
		InstantiateBonds(bonds.anglebonds, instance.particleOffset, anglebonds.data() + instance.firstAnglebond);
		InstantiateBonds(bonds.dihedralbonds, instance.particleOffset, dihedralbonds.data() + instance.firstDihedralbond);
		InstantiateBonds(bonds.improperdihedralbonds, instance.particleOffset, improperdihedralbonds.data() + instance.firstImproperdihedralbond);
		InstantiateBonds(bonds.cmapbonds, instance.particleOffset, cmapbonds.data() + instance.firstCmapbond);
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
			/*throw std::runtime_error(Lima::Format("Loading singlebond with illegally large dist ({}). b0: {}. AtomIndices: {} {}",
				hyper_dist, bond.params.b0, bond.global_atom_indexes[0], bond.global_atom_indexes[1]));*/
		}
		if (hyper_dist < bondRelaxedDist * 0.001) {
			/*throw std::runtime_error(Lima::Format("Loading singlebond with illegally small dist ({}). b0: {}. AtomIndices: {} {}",
				hyper_dist, bond.params.b0, bond.global_atom_indexes[0], bond.global_atom_indexes[1]));*/
		}
	}
	for (const auto& bond : anglebonds)
	{
		const Float3 pos1 = particles[bond.global_atom_indexes[0]].position;
		const Float3 pos2 = particles[bond.global_atom_indexes[1]].position;
		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos1, pos2, boxlen_nm, bc_select);
		if (hyper_dist < 0.001)
			throw std::runtime_error(Lima::Format("Loading singlebond with illegally small dist ({}). b0: {}", hyper_dist, bond.params.ub0));
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
		std::vector<PclusterConstraints> constraints;	// Per cluster, empty without constraints
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
					throw std::runtime_error(Lima::Format("Moleculetype {}: atom at index {} has id {}", molecule.name, i, molecule.atoms[i].id));

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
		const TopologyFile::Moleculetype& molecule, const MoleculetypeConstraints* constraints
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

		// With constraints, each constrained group is added to a cluster as a whole, so it never spans clusters
		if (constraints) {
			std::vector<bool> added(molecule.atoms.size(), false);
			for (const std::vector<int>& component : connectedComponents) {
				clusters.emplace_back(std::array{ -1, -1, -1, -1 });
				int nextIndex = 0;
				for (const int particleId : component) {
					if (added[particleId])
						continue;
					const int group = constraints->groupOf[particleId];
					const std::vector<int> single{ particleId };
					const std::vector<int>& members = group == -1 ? single : constraints->groups[group];
					if (nextIndex != 0 && (nextIndex + members.size() > PersistentCluster::maxParticles || !CanAppendToCluster(particleId, clusters.back(), nextIndex))) {
						clusters.emplace_back(std::array{ -1, -1, -1, -1 });
						nextIndex = 0;
					}
					for (const int member : members) {
						clusters.back()[nextIndex++] = member;
						added[member] = true;
					}
				}
			}
			return clusters;
		}

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
		const LIMAForcefield& forcefield,
		const MoleculetypeConstraints* constraints
	) {
		PersistentClusterTemplate result;
		result.clusters = MakeLocalPersistentClusters(molecule, constraints);
		result.particleToPcluster.resize(molecule.atoms.size(), -1);
		result.particles.resize(molecule.atoms.size());

		for (int pcid = 0; pcid < result.clusters.size(); pcid++) {
			for (const int particleId : result.clusters[pcid]) {
				if (particleId != -1)
					result.particleToPcluster[particleId] = pcid;
			}
		}

		if (constraints) {
			result.constraints.resize(result.clusters.size());
			const auto IndexInCluster = [&](int particleId) {
				const auto& cluster = result.clusters[result.particleToPcluster[particleId]];
				return static_cast<int>(std::ranges::find(cluster, particleId) - cluster.begin());
			};
			for (size_t c = 0; c < constraints->pairs.size(); c++) {
				const auto [a, b] = constraints->pairs[c];
				if (result.particleToPcluster[a] != result.particleToPcluster[b])
					throw std::runtime_error("A constraint spans two pclusters");
				PclusterConstraints& clusterConstraints = result.constraints[result.particleToPcluster[a]];
				if (clusterConstraints.nConstraints == PclusterConstraints::maxConstraints)
					throw std::runtime_error("Too many constraints in a pcluster");
				clusterConstraints.particles[clusterConstraints.nConstraints] = static_cast<uint8_t>(IndexInCluster(a) | IndexInCluster(b) << 2);
				clusterConstraints.distance[clusterConstraints.nConstraints] = constraints->distances[c];
				clusterConstraints.nConstraints++;
			}
		}

		// Collect as vectors and sort afterwards, which is much faster than inserting into sets
		std::vector<std::vector<int>> particleBondedToParticle(molecule.atoms.size());
		std::vector<std::vector<int>> pclusterBondedToPcluster(result.clusters.size());
		const auto AddPclusterPair = [&](int pidA, int pidB) {
			const int pcidA = result.particleToPcluster[pidA];
			const int pcidB = result.particleToPcluster[pidB];
			assert(pcidA != -1 && pcidB != -1);
			pclusterBondedToPcluster[pcidA].push_back(pcidB);
			pclusterBondedToPcluster[pcidB].push_back(pcidA);
		};
		// The pclusters are bonded where they share a bonded term or an excluded pair
		auto AddBond = [&](const auto& bond) {
			const auto& ids = bond.ids;
			for (const int id : ids) {
				if (id < 0 || id >= result.particleToPcluster.size())
					return;
			}
			for (int i = 0; i < ids.size(); i++)
				for (int j = i + 1; j < ids.size(); j++)
					AddPclusterPair(ids[i], ids[j]);
		};

		for (const auto& bond : molecule.singlebonds)
			AddBond(bond);
		for (const auto& bond : molecule.anglebonds)
			AddBond(bond);
		for (const auto& bond : molecule.dihedralbonds)
			AddBond(bond);
		for (const auto& bond : molecule.improperdihedralbonds)
			AddBond(bond);
		for (const auto [a, b] : ExcludedPairs(molecule)) {
			particleBondedToParticle[a].push_back(b);
			particleBondedToParticle[b].push_back(a);
			AddPclusterPair(a, b);
		}

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

			particle.atomLetter = RenderUtilities::AtomLetter(atom.atomname, atom.residue); // TODO: This is hacky. Boxbuilder must be agnostic to rendering
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
		const auto constraints = system.constraints.find(uniqueMoleculetypes[i]);
		templates[i] = BuildPersistentClusterTemplate(*uniqueMoleculetypes[i], forcefield,
			constraints == system.constraints.end() ? nullptr : &constraints->second);
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
	if (!system.constraints.empty())
		pcFactory.pclusterConstraints.resize(totalClusterCount);

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

		for (int localPcid = 0; localPcid < persistentTemplate.constraints.size(); localPcid++)
			pcFactory.pclusterConstraints[pclusterOffset + localPcid] = persistentTemplate.constraints[localPcid];

		for (int localPcid = 0; localPcid < persistentTemplate.pclusterBondedToPcluster.size(); localPcid++) {
			PclustersBondedToPcluster bonded = persistentTemplate.pclusterBondedToPcluster[localPcid];
			bonded.AddOffset(pclusterOffset);
			pcFactory.pclusterBondedToPcluster[pclusterOffset + localPcid] = bonded;
		}
		}, 64, 16);

	return pcFactory;
}





std::unique_ptr<BoxImage> LIMA_MOLECULEBUILD::buildMolecules(
	const GroFile& grofile,
	const TopologyFile& topol_file,
	const SimParams& simparams
)
{
	//TimeIt timer("buildMolecules", true);
	LIMAForcefield forcefield{ topol_file.forcefieldInclude ? topol_file.forcefieldInclude->contents : GenericItpFile{} };

	SuperTopology superTopology(topol_file.GetSystem(), grofile, forcefield, simparams.constraints.value, simparams.enable_electrostatics);
	superTopology.VerifyBondsAreStable(grofile.box_size, simparams.bc_select, simparams.em_variant);


	std::future<BondGroupFactory> bgfFuture = std::async(std::launch::async, [&]{ return BondGroupFactory(superTopology); });

	// Make PersistenClusters
	std::future<PersistentClusterFactory> pcFactoryFuture = std::async(std::launch::async, [&]{ return MakePersistentClusters(superTopology, forcefield); });

	BondGroupFactory bgFactory = bgfFuture.get();


	PersistentClusterFactory pcFactory = pcFactoryFuture.get();
	bgFactory.AddPclusterRefs(pcFactory.particleToPclusterMap);

	bgFactory.AddBondgroupRefsToPclusters(pcFactory.particleToPclusterMap, pcFactory.pClusterMetas);
	BondGroups bondGroups = bgFactory.GetBondgroups();
	if (!bondGroups.cmapbonds.empty())
		bondGroups.cmapGrids = forcefield.GetCmapGrids();



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

	auto pclusterConstraints = std::move(pcFactory.pclusterConstraints);
	auto boxImage = std::make_unique<BoxImage>(
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
	boxImage->pclusterConstraints = std::move(pclusterConstraints);
	boxImage->builtConstraints = simparams.constraints.value;
	return boxImage;
}
