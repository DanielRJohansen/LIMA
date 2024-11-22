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
public:
	BondedParticlesLUTManagerFactory(int nCompounds) : 
		nCompounds(nCompounds) ,
		luts(nCompounds * BondedParticlesLUTHelpers::max_bonded_compounds, BondedParticlesLUT(false)) 
	{}


	BondedParticlesLUT& get(int id_self, int id_other) {
		if (std::abs(id_self - id_other > 2)) {
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
};


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
void CompoundFactory::addParticle(const ParticleFactory& particle, int global_id, float boxlen_nm, BoundaryConditionSelect bc) {
	if (n_particles >= MAX_COMPOUND_PARTICLES) {
		throw std::runtime_error("Failed to add particle to compound");
	}

	// Hyperposition each compound relative to particle 0, so we can find key_particle and radius later
	Float3 hyperpos = particle.position;
	if (n_particles > 0) {
		BoundaryConditionPublic::applyHyperposNM(positions[0], hyperpos, boxlen_nm, bc); // i dont LOVE this.. :(
	}

	// Variables only present in factory
	positions[n_particles] = hyperpos;
	global_ids[n_particles] = global_id;

	// Variables present in Compound
	atom_types[n_particles] = particle.activeLjtypeParameterIndex;
	atomLetters[n_particles] = particle.topAtom->type[0];

	atom_charges[n_particles] = particle.topAtom->charge * elementaryChargeToKiloCoulombPerMole;
	atomMasses[n_particles] = particle.topAtom->mass / KILO;

	indicesInGrofile[n_particles] = particle.indexInGrofile;

	n_particles++;
}

template <int n> 
std::array<uint8_t, n> ConvertGlobalAtomidsToCompoundlocalIds(const std::vector<ParticleToCompoundMapping>& p2cMap, const std::array<uint32_t, n>& global_atom_ids) {
	std::array<uint8_t, n> localIds;
	for (int i = 0; i < n; i++) {
		localIds[i] = static_cast<uint8_t>(p2cMap[global_atom_ids[i]].localIdInCompound);
	}
	return localIds;
}

void CompoundFactory::AddBond(const std::vector<ParticleToCompoundMapping>& p2cMap, const SingleBondFactory& bond) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add singlebond to compound");
	}
	singlebonds[n_singlebonds++] = SingleBond(ConvertGlobalAtomidsToCompoundlocalIds<2>(p2cMap, bond.global_atom_indexes), bond.params);
}
void CompoundFactory::AddBond(const std::vector<ParticleToCompoundMapping>& p2cMap, const AngleBondFactory& bond) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add anglebond to compound");
	}
	anglebonds[n_anglebonds++] = AngleUreyBradleyBond(ConvertGlobalAtomidsToCompoundlocalIds<3>(p2cMap, bond.global_atom_indexes), bond.params);
}
void CompoundFactory::AddBond(const std::vector<ParticleToCompoundMapping>& p2cMap, const DihedralBondFactory& bond) {
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add dihedralbond to compound");
	}
	dihedrals[n_dihedrals++] = DihedralBond(ConvertGlobalAtomidsToCompoundlocalIds<4>(p2cMap, bond.global_atom_indexes), bond.params);
}
void CompoundFactory::AddBond(const std::vector<ParticleToCompoundMapping>& p2cMap, const ImproperDihedralBondFactory& bond) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add improperdihedralbond to compound");
	}
	impropers[n_improperdihedrals++] = ImproperDihedralBond(ConvertGlobalAtomidsToCompoundlocalIds<4>(p2cMap, bond.global_atom_indexes), bond.params);
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
std::array<uint8_t, BondFactory_type::nAtoms> BridgeFactory::ConvertGlobalIdsToBridgelocalIds (const std::vector<ParticleToCompoundMapping>& p2cMap, 
	ParticleToBridgeMap& p2bMap, const BondFactory_type& bond) 
{
	for (const auto& global_atom_id : bond.global_atom_indexes) {
		if (!p2bMap[global_atom_id].has_value()) {
			addParticle(p2cMap[global_atom_id], p2bMap[global_atom_id]);
		}
		const auto bridgeId = p2bMap[global_atom_id]->bridgeId;
		if (bridgeId != this->bridge_id) {
			throw std::runtime_error("Cannot add bond to bridge, as it contains atoms from another bridge");
		}
	}
	
	std::array<uint8_t, BondFactory_type::nAtoms> localIds;
	std::unordered_set<int> usedIds;
	for (int i = 0; i < BondFactory_type::nAtoms; i++) {
		localIds[i] = p2bMap[bond.global_atom_indexes[i]]->localIdInBridge;
		if (usedIds.contains(localIds[i])) {
			throw std::runtime_error("Failed to add bond to bridge, as it contains multiple atoms from the same compound");
		}
		usedIds.insert(localIds[i]);
	}
	return localIds;
}

void BridgeFactory::AddBond(const ParticleToCompoundMap& p2cMap, ParticleToBridgeMap& p2bMap, const SingleBondFactory& bond) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add singlebond to bridge"); }
	singlebonds[n_singlebonds++] = SingleBond{
		ConvertGlobalIdsToBridgelocalIds(p2cMap, p2bMap, bond),
		bond.params
	};
}
void BridgeFactory::AddBond(const std::vector<ParticleToCompoundMapping>& p2cMap, ParticleToBridgeMap& p2bMap, const AngleBondFactory& bond) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add anglebond to bridge"); }
	anglebonds[n_anglebonds++] = AngleUreyBradleyBond{
		ConvertGlobalIdsToBridgelocalIds(p2cMap, p2bMap, bond),
		bond.params
	};
}
void BridgeFactory::AddBond(const std::vector<ParticleToCompoundMapping>& p2cMap, ParticleToBridgeMap& p2bMap, const DihedralBondFactory& bond) {
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_BRIDGE) {
		throw std::runtime_error("Failed to add dihedralbond to bridge");
	}
	dihedrals[n_dihedrals++] = DihedralBond{
		ConvertGlobalIdsToBridgelocalIds(p2cMap, p2bMap, bond),
		bond.params
	};
}
void BridgeFactory::AddBond(const std::vector<ParticleToCompoundMapping>& p2cMap, ParticleToBridgeMap& p2bMap, const ImproperDihedralBondFactory& bond) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add improperdihedralbond to bridge"); }
	impropers[n_improperdihedrals++] = ImproperDihedralBond{
		ConvertGlobalIdsToBridgelocalIds(p2cMap, p2bMap, bond),
		bond.params
	};
}

void BridgeFactory::addParticle(const ParticleToCompoundMapping& p2cMapping, std::optional<ParticleToBridgeMapping>& p2bMapping) {	// This can be made const, and return it's iD in bridge
	if (n_particles == MAX_PARTICLES_IN_BRIDGE) { throw std::runtime_error("Too many particles in bridge"); }

	p2bMapping = ParticleToBridgeMapping{ bridge_id, n_particles };

	int compoundlocalid_in_bridge = -1;
	for (int i = 0; i < this->n_compounds; i++) {
		if (p2cMapping.compoundId == compound_ids[i])
			compoundlocalid_in_bridge = i;
	}
	if (compoundlocalid_in_bridge == -1) {
		throw std::runtime_error("Could not find compoundlocalid_in_bridge");
	}

	particle_refs[n_particles] = ParticleReference{
		p2cMapping.compoundId,
		p2cMapping.localIdInCompound,
		static_cast<uint8_t>(compoundlocalid_in_bridge)
	};

	n_particles++;
}

template <int N>
bool Belongs(const std::array<int, N> bondIds, int idOfFirstAtomInMolecule, int idOfLastAtomInMolecule) {
	for (const auto& id : bondIds) {
		if (id<idOfFirstAtomInMolecule || id>idOfLastAtomInMolecule) {
			return false;
		}
	}
	return true;
}
template <int N>
void CorrectBondIds(std::array<int, N>& bondIds, int idOfFirstAtomInMolecule) {
	for (auto& id : bondIds) {
		id -= idOfFirstAtomInMolecule;
	}
}
template <typename BondType>
void AddBondsThatBelongsInMolecule(std::vector<BondType>& dst, const std::vector<BondType>& src, int idOfFirstAtomInMolecule, int idOfLastAtomInMolecule) {
	for (const auto& bond : src) {
		if (Belongs<BondType::n>(bond.ids, idOfFirstAtomInMolecule, idOfLastAtomInMolecule)) {
			dst.push_back(bond);
			CorrectBondIds<BondType::n>(dst.back().ids, idOfFirstAtomInMolecule);
		}
	}
}


struct MoleculeTopology {
	// When entire topology file belongs to the molecule
	MoleculeTopology(const TopologyFile::Moleculetype& moltype) 
		: moleculeTypeName(moltype.name)
	{
		atoms = moltype.atoms;
		singlebonds = moltype.singlebonds;
		anglebonds = moltype.anglebonds;
		dihedralbonds = moltype.dihedralbonds;
		improperdihedralbonds = moltype.improperdihedralbonds;
	}

	// When only the provided id's belong to the molecule
	MoleculeTopology(const TopologyFile::Moleculetype& moltype, std::vector<int> idsInMolecule) {
		assert(std::is_sorted(idsInMolecule.begin(), idsInMolecule.end()));

		const int idOfFirstAtomInMolecule = idsInMolecule.front();
		const int idOfLastAtomInMolecule = idsInMolecule.back();

		for (const auto& atom : moltype.atoms) {
			if (atom.id >= idOfFirstAtomInMolecule && atom.id <= idOfLastAtomInMolecule) {
				atoms.push_back(atom);
				atoms.back().id -= idOfFirstAtomInMolecule;
			}
		}
		AddBondsThatBelongsInMolecule(singlebonds, moltype.singlebonds, idOfFirstAtomInMolecule, idOfLastAtomInMolecule);
		AddBondsThatBelongsInMolecule(anglebonds, moltype.anglebonds, idOfFirstAtomInMolecule, idOfLastAtomInMolecule);
		AddBondsThatBelongsInMolecule(dihedralbonds, moltype.dihedralbonds, idOfFirstAtomInMolecule, idOfLastAtomInMolecule);
		AddBondsThatBelongsInMolecule(improperdihedralbonds, moltype.improperdihedralbonds, idOfFirstAtomInMolecule, idOfLastAtomInMolecule);
	}

	std::string moleculeTypeName;
	std::vector<TopologyFile::AtomsEntry> atoms;
	std::vector<TopologyFile::SingleBond> singlebonds;
	std::vector<TopologyFile::AngleBond> anglebonds;
	std::vector<TopologyFile::DihedralBond> dihedralbonds;
	std::vector<TopologyFile::ImproperDihedralBond> improperdihedralbonds;
};




struct TinyMolRef {
	int atomsOffsetInTinymols = 0;
	int atomsOffsetInGrofile = 0;
	//const TopologyFile::MoleculeEntry& molecule;
	MoleculeTopology topology;
};
struct MoleculeRef {
	int atomsOffsetInMolecules = 0;
	int atomsOffsetInGrofile = 0;
	//const TopologyFile::MoleculeEntry& molecule;
	MoleculeTopology topology;
	std::shared_ptr<MoleculeGraph> graph;
};
std::pair<const std::vector<MoleculeRef>, const std::vector<TinyMolRef>> PrepareMolecules(const TopologyFile& topol) {
	std::vector<MoleculeRef> molecules;
	std::vector<TinyMolRef> tinyMolecules;

	int atomsOffsetInMolecules = 0;
	int atomsOffsetInTinyMolecules = 0;
	int atomsOffsetInGrofile = 0;
	for (const auto& molecule : topol.GetSystem().molecules) {
		if (molecule.moleculetype->atoms.empty())
			throw std::runtime_error("Molecule has no atoms");


		auto superGraph = std::make_shared<MoleculeGraph> (*molecule.moleculetype);
		if (!superGraph->GraphIsDisconnected()) {

			if (molecule.moleculetype->atoms.size() <= 3 && molecule.moleculetype->atoms[0].residue != "lxx") // lxx is a lima code that forces it to be a normal compound
			{
				tinyMolecules.emplace_back(TinyMolRef{ atomsOffsetInMolecules, atomsOffsetInGrofile, MoleculeTopology(*molecule.moleculetype) });
				atomsOffsetInTinyMolecules += molecule.moleculetype->atoms.size();
			}
			else {
				molecules.emplace_back(MoleculeRef{ atomsOffsetInMolecules, atomsOffsetInGrofile, MoleculeTopology(*molecule.moleculetype), superGraph });
				atomsOffsetInMolecules += molecule.moleculetype->atoms.size();
			}
			atomsOffsetInGrofile += molecule.moleculetype->atoms.size();
		}
		else {
			auto subGraphs = superGraph->GetListOfListsofConnectedNodeids();
			
			for (auto& graphNodeIds : subGraphs) {
				std::sort(graphNodeIds.begin(), graphNodeIds.end());
				if (graphNodeIds.size() <= 3 && molecule.moleculetype->atoms[0].residue != "lxx") {
					tinyMolecules.emplace_back(TinyMolRef{ atomsOffsetInMolecules, atomsOffsetInGrofile, MoleculeTopology(*molecule.moleculetype, graphNodeIds) });
					atomsOffsetInTinyMolecules += graphNodeIds.size();
				}
				else {
					
					std::unordered_set<int> nodesSet(graphNodeIds.begin(), graphNodeIds.end());
					molecules.emplace_back(MoleculeRef{ atomsOffsetInMolecules, atomsOffsetInGrofile, MoleculeTopology(*molecule.moleculetype, graphNodeIds), std::make_shared<MoleculeGraph>(*molecule.moleculetype, nodesSet)});
					atomsOffsetInMolecules += graphNodeIds.size();
					//throw std::runtime_error("Disconnected molecule with more than 3 atoms is not supported yet.");
				}
				atomsOffsetInGrofile += graphNodeIds.size();
			}
		}




	}
	return { molecules, tinyMolecules };
}

std::vector<TinyMolFactory> LoadTinyMols(const GroFile& grofile, const std::vector<TinyMolRef>& tinyMolsRefs, LIMAForcefield& forcefieldManager) {
	std::vector<TinyMolFactory> tinyMols;
	tinyMols.reserve(grofile.atoms.size());

	for (const TinyMolRef& tinyMolRef : tinyMolsRefs) {
		// FOr now ignore that there are multiple atoms in a tinymol. LOOONG TODO
		tinyMols.emplace_back(TinyMolFactory{ 
			grofile.atoms[tinyMolRef.atomsOffsetInGrofile].position,  
			//forcefieldManager.GetActiveTinymoltypeIndex(tinyMolRef.molecule.moleculetype->atoms[0].type) 
            forcefieldManager.GetActiveTinymoltypeIndex(tinyMolRef.topology.atoms[0].type),
            tinyMolRef.topology.atoms[0].type
			});
	}

	return tinyMols;
}


template <typename BondType, typename BondtypeFactory, typename BondTypeTopologyfile>
void LoadBondIntoTopology(const std::vector<BondTypeTopologyfile>& bondsInTopfile,	int atomIdOffset, LIMAForcefield& forcefield,
	const std::vector<ParticleFactory>& particles, std::vector<BondtypeFactory>& topology)
{
	for (const auto& bondTopol : bondsInTopfile) {
		std::array<uint32_t, BondType::nAtoms> globalIds;
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

const Topology LoadTopology(const std::vector<MoleculeRef>& molecules, LIMAForcefield& forcefield, const GroFile& grofile) {
	Topology topology;	
	if (molecules.empty())
		return topology;

	const int numParticles = molecules.back().atomsOffsetInMolecules + molecules.back().topology.atoms.size();

	topology.particles.reserve(numParticles);
	int uniqueResId = -1;
	for (const auto& mol : molecules) {
		int currentResId = -1;
		for (int atomIndex = 0; atomIndex < mol.topology.atoms.size(); atomIndex++) {
			const TopologyFile::AtomsEntry& atom = mol.topology.atoms[atomIndex];

			if (atom.resnr != currentResId) {
				currentResId = atom.resnr;
				uniqueResId++;
			}

			const int activeLJParamIndex = forcefield.GetActiveLjParameterIndex(atom.type);
			const int indexInGrofile = mol.atomsOffsetInGrofile + atomIndex;
			topology.particles.emplace_back(ParticleFactory{ activeLJParamIndex, &atom, grofile.atoms[indexInGrofile].position, uniqueResId, indexInGrofile });
		}
	}

	topology.singlebonds.reserve(std::accumulate(molecules.begin(), molecules.end(), 0, [](int sum, const MoleculeRef& mol) { return sum + mol.topology.singlebonds.size(); }));
	topology.anglebonds.reserve(std::accumulate(molecules.begin(), molecules.end(), 0, [](int sum, const MoleculeRef& mol) { return sum + mol.topology.anglebonds.size(); }));
	topology.dihedralbonds.reserve(std::accumulate(molecules.begin(), molecules.end(), 0, [](int sum, const MoleculeRef& mol) { return sum + mol.topology.dihedralbonds.size(); }));
	topology.improperdihedralbonds.reserve(std::accumulate(molecules.begin(), molecules.end(), 0, [](int sum, const MoleculeRef& mol) { return sum + mol.topology.improperdihedralbonds.size(); }));
	for (const auto& mol : molecules) 
	{
		LoadBondIntoTopology<SingleBond, SingleBondFactory, TopologyFile::SingleBond>(
			mol.topology.singlebonds, mol.atomsOffsetInMolecules, forcefield, topology.particles, topology.singlebonds);

		LoadBondIntoTopology<AngleUreyBradleyBond, AngleBondFactory, TopologyFile::AngleBond>(
			mol.topology.anglebonds, mol.atomsOffsetInMolecules, forcefield, topology.particles, topology.anglebonds);

		LoadBondIntoTopology<DihedralBond, DihedralBondFactory, TopologyFile::DihedralBond>(
			mol.topology.dihedralbonds, mol.atomsOffsetInMolecules, forcefield, topology.particles, topology.dihedralbonds);

		LoadBondIntoTopology<ImproperDihedralBond, ImproperDihedralBondFactory, TopologyFile::ImproperDihedralBond>(
			mol.topology.improperdihedralbonds, mol.atomsOffsetInMolecules, forcefield, topology.particles, topology.improperdihedralbonds);
	}

	return topology;
}



std::vector<std::vector<int>> MapAtomsToSinglebonds(const Topology& topology) {
	std::vector<std::vector<int>> atomIdToSinglebondsMap(topology.particles.size());

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
	std::set<int> idsOfBondedAtomgroups;
};

bool AreBonded(const AtomGroup& left, const AtomGroup& right, const std::vector<std::vector<int>>& atomIdToSinglebondsMap) {
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

const std::vector<AtomGroup> GroupAtoms(const Topology& topology) {
	const int nResidues = topology.particles.empty() ? 0 : topology.particles.back().uniqueResId + 1;
	std::vector<AtomGroup> atomGroups;
	atomGroups.reserve(nResidues);
	
	// Create groups of atoms
	for (int particleId = 0; particleId < topology.particles.size(); particleId++) {
		const auto& atom = topology.particles[particleId];
		if (atomGroups.empty() || atom.topAtom->section_name.has_value() || atom.uniqueResId != topology.particles[particleId - 1].uniqueResId)
			atomGroups.push_back(AtomGroup{});
		atomGroups.back().atomIds.push_back(particleId);
	}

	// Now figure out which groups are bonded to others
	const auto particleToSinglebondMap = MapAtomsToSinglebonds(topology);
	for (int i = 1; i < atomGroups.size(); i++) {
		if (AreBonded(atomGroups[i-1], atomGroups[i], particleToSinglebondMap)) {
			atomGroups[i].idsOfBondedAtomgroups.insert(i - 1);
			atomGroups[i - 1].idsOfBondedAtomgroups.insert(i);
		}
	}

	return atomGroups;
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

const std::vector<AtomGroup> GroupAtoms1(const std::vector<MoleculeRef>& molecules, const Topology& topology) {
	std::vector<AtomGroup> atomGroups;

	for (const auto& mol : molecules) {
		const std::unordered_map<int, int> nodeIdNumDownstreamNodes = mol.graph->ComputeNumDownstreamNodes();
		const MoleculeTree moleculeTree = mol.graph->ConstructMoleculeTree();

		std::stack<const MoleculeGraph::Node*> nodeStack;
		nodeStack.push(mol.graph->root);

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
					[&nodeIdNumDownstreamNodes](const int& a, const int& b) { return nodeIdNumDownstreamNodes.at(a) < nodeIdNumDownstreamNodes.at(b);}
				) - nodeChildren.begin();
				nodeStack.push(&mol.graph->nodes[nodeChildren[indexOfLongestChain]]);
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
	const auto particleToSinglebondMap = MapAtomsToSinglebonds(topology);
	for (int i = 1; i < atomGroups.size(); i++) {
		if (AreBonded(atomGroups[i - 1], atomGroups[i], particleToSinglebondMap)) {
			atomGroups[i].idsOfBondedAtomgroups.insert(i - 1);
			atomGroups[i - 1].idsOfBondedAtomgroups.insert(i);
		}
	}

	return atomGroups;
}



std::vector<CompoundFactory> CreateCompounds(const Topology& topology, float boxlen_nm, const std::vector<AtomGroup>& atomGroups, BoundaryConditionSelect bc_select)
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

bool ContainsSubset(const std::set<int>& set, const std::set<int>& subset) {
	return std::all_of(subset.begin(), subset.end(), [&set](const int& elem) {
		return set.find(elem) != set.end();
		});
}
 
template<typename BondType>
void AddToSet(const std::vector<BondType>& bonds, std::vector<std::set<int>>& compoundidsInteractedWithByAtoms, const std::vector<ParticleToCompoundMapping>& particleToCompoundidMap) {
	for (const auto bond : bonds) {
		std::set<int> compoundIdsInThisBond;
		for (const auto globalAtomId : bond.global_atom_indexes) {
			compoundIdsInThisBond.insert(particleToCompoundidMap[globalAtomId].compoundId);
		}

		for (const auto globalAtomId : bond.global_atom_indexes) {
			compoundidsInteractedWithByAtoms[globalAtomId].insert(compoundIdsInThisBond.begin(), compoundIdsInThisBond.end());
		}
	}
}

std::vector<std::set<int>> TrimSetsForSubsets(const std::set<std::set<int>>& compoundsBridgeCandidates) {
	std::vector<std::set<int>> compoundsBridgeCandidatesSorted(compoundsBridgeCandidates.begin(), compoundsBridgeCandidates.end());
	std::sort(compoundsBridgeCandidatesSorted.begin(), compoundsBridgeCandidatesSorted.end(), [](const std::set<int>& a, const std::set<int>& b) { return a.size() > b.size(); }); // Largest set first

	std::vector<std::set<int>> compoundsBridgeCandidatesTrimmed;
	for (const auto& candidate : compoundsBridgeCandidatesSorted) {
		bool isSubset = false;
		for (const auto& compound : compoundsBridgeCandidatesTrimmed) {
			if (ContainsSubset(compound, candidate)) {
				isSubset = true;
				break;
			}
		}
		if (!isSubset) {
			compoundsBridgeCandidatesTrimmed.emplace_back(candidate);
		}
	}

	// Sort the bridges by first compoundid, for debugging purposes
	std::sort(compoundsBridgeCandidatesTrimmed.begin(), compoundsBridgeCandidatesTrimmed.end(), [](const std::set<int>& a, const std::set<int>& b) { return *a.begin() < *b.begin(); });

	return compoundsBridgeCandidatesTrimmed;
}




std::vector<BridgeFactory> CreateBridges(const Topology& topology, const std::vector<CompoundFactory>& compounds, 
	const std::vector<ParticleToCompoundMapping>& particleToCompoundidMap)
{
	const int nParticlesTotal = particleToCompoundidMap.size();

	// First figure out which compounds each particle interacts with via bonds
	std::vector<std::set<int>> compoundidsInteractedWithByAtoms(nParticlesTotal);
	AddToSet(topology.singlebonds, compoundidsInteractedWithByAtoms, particleToCompoundidMap);
	AddToSet(topology.anglebonds, compoundidsInteractedWithByAtoms, particleToCompoundidMap);
	AddToSet(topology.dihedralbonds, compoundidsInteractedWithByAtoms, particleToCompoundidMap);

	// Any particle interacting with > 1 compound, is a bridge candidate
	std::set<std::set<int>> compoundsBridgeCandidates;
	for (int i = 0; i < nParticlesTotal; i++) {
		if (compoundidsInteractedWithByAtoms[i].size() > 1)
			compoundsBridgeCandidates.insert(compoundidsInteractedWithByAtoms[i]);
	}

	// Now trim the candidates
	const std::vector<std::set<int>> compoundsBridgeCandidatesTrimmed = TrimSetsForSubsets(compoundsBridgeCandidates);

	std::vector<BridgeFactory> bridges;
	bridges.reserve(compoundsBridgeCandidatesTrimmed.size());
	for (const std::set<int>& compoundIdSet : compoundsBridgeCandidatesTrimmed) {
		if (compoundIdSet.size() < 2)
			throw std::runtime_error("Bridge must contain at least 2 compounds");
		
		const std::vector<int> compoundIdsInBridge{ compoundIdSet.begin(), compoundIdSet.end() };
		bridges.push_back(BridgeFactory{ static_cast<int>(bridges.size()), compoundIdsInBridge });
	}

	for (int i = 0; i < bridges.size(); i++)
		if (bridges[i].bridge_id != i)
			throw std::runtime_error("Bridge id mismatch");

	return bridges;
}


template <int n_ids>
bool SpansTwoCompounds(const std::array<uint32_t, n_ids> bond_globalids, const std::vector<ParticleToCompoundMapping>& p2cMap) {
	const int compoundid_left = p2cMap[bond_globalids[0]].compoundId;

	for (auto globalId : bond_globalids) {
		if (compoundid_left != p2cMap[globalId].compoundId) {
			return true;
		}
	}
	return false;
}

template <int n_ids>
void DistributeLJIgnores(BondedParticlesLUTManagerFactory* bplut_man, const std::vector<ParticleToCompoundMapping>& particleToCompoundMap, const std::array<uint32_t, n_ids>& global_ids) {
	for (auto gid_self : global_ids) {
		for (auto gid_other : global_ids) {
			if (gid_self == gid_other) { continue; }
			
			const ParticleToCompoundMapping& mappingSelf = particleToCompoundMap[gid_self];
			const ParticleToCompoundMapping& mappingOther = particleToCompoundMap[gid_other];

			BondedParticlesLUT& lut = bplut_man->get(mappingSelf.compoundId, mappingOther.compoundId);
			lut.set(mappingSelf.localIdInCompound, mappingOther.localIdInCompound, true);
		}
	}
}

// Returns two compound id's of a bond in a bridge. The id's are sorted with lowest first
template <int N>
std::array<int, 2> getTheTwoDifferentIds(const std::array<uint32_t, N>& particle_global_ids, const std::vector<ParticleToCompoundMapping>& p2cMap) {
	std::array<int, 2> compound_ids = { p2cMap[particle_global_ids[0]].compoundId, -1 };

	for (auto globalId : particle_global_ids) {
		int compoundid_of_particle = p2cMap[globalId].compoundId;
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
BridgeFactory& getBridge(std::vector<BridgeFactory>& bridges, const std::array<uint32_t, N>& particle_global_ids, const std::vector<ParticleToCompoundMapping>& p2cMap
)
{
	std::array<int, 2> compound_ids = getTheTwoDifferentIds<N>(particle_global_ids, p2cMap);

	// If not found in the lookup or no matching bridge, search through all bridges
	for (BridgeFactory& bridge : bridges) {
		if (bridgeContainsTheseTwoCompounds(bridge, compound_ids)) {
			return bridge;
		}
	}

	throw std::runtime_error(std::format("Failed to find the bridge ({})", N).c_str());
}

template <typename BondTypeFactory>
void DistributeBondsToCompoundsAndBridges(const std::vector<BondTypeFactory>& bonds, 
	const std::vector<ParticleToCompoundMapping>& particlesToCompoundIdMap, ParticleToBridgeMap& p2bMap,
	std::vector<CompoundFactory>& compounds, std::vector<BridgeFactory>& bridges,
	BondedParticlesLUTManagerFactory& bpLutManager) 
{
	

	for (const BondTypeFactory& bond : bonds) {
		
		DistributeLJIgnores<BondTypeFactory::nAtoms>(&bpLutManager, particlesToCompoundIdMap, bond.global_atom_indexes);

		// TODO Remove this once we switch to a graph based LJignores search
		if (bond.params.HasZeroParam())
			continue;
		//


		if (SpansTwoCompounds<BondTypeFactory::nAtoms>(bond.global_atom_indexes, particlesToCompoundIdMap)) {
			BridgeFactory& bridge = getBridge<BondTypeFactory::nAtoms>(bridges, bond.global_atom_indexes, particlesToCompoundIdMap);
			bridge.AddBond(particlesToCompoundIdMap, p2bMap, bond);
		}
		else {
			const int compound_id = particlesToCompoundIdMap[bond.global_atom_indexes[0]].compoundId;	// Pick compound using first particle in bond
			CompoundFactory& compound = compounds.at(compound_id);
			compound.AddBond(particlesToCompoundIdMap, bond);
		}

	}
}

void DistributeBondsToCompoundsAndBridges(const Topology& topology, float boxlen_nm, const std::vector<ParticleToCompoundMapping>& particlesToCompoundIdMap, BoundaryConditionSelect bc_select,
	std::vector<CompoundFactory>& compounds, std::vector<BridgeFactory>& bridges, BondedParticlesLUTManagerFactory& bpLutManager)
{
	ParticleToBridgeMap particlesToBridgeMap(particlesToCompoundIdMap.size());

	DistributeBondsToCompoundsAndBridges(topology.singlebonds, particlesToCompoundIdMap, particlesToBridgeMap, compounds, bridges, bpLutManager);
	DistributeBondsToCompoundsAndBridges(topology.anglebonds, particlesToCompoundIdMap, particlesToBridgeMap, compounds, bridges, bpLutManager);
	DistributeBondsToCompoundsAndBridges(topology.dihedralbonds, particlesToCompoundIdMap, particlesToBridgeMap, compounds, bridges, bpLutManager);
	DistributeBondsToCompoundsAndBridges(topology.improperdihedralbonds, particlesToCompoundIdMap, particlesToBridgeMap, compounds, bridges, bpLutManager);

	// While are particle is never bonded to itself, we are not allowed to calc LJ with itself
	// so we can doubly use this lut to avoid doing that

	for (int com_id = 0; com_id < compounds.size(); com_id++) {
		auto& lut = bpLutManager.get(com_id, com_id);
		for (int pid = 0; pid < MAX_COMPOUND_PARTICLES; pid++) {
			lut.set(pid, pid, true);
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
//#pragma omp parallel for // TODO Add OMP here, 
	//for (CompoundFactory& compound : compounds) {
	for (int cid = 0; cid < compounds.size(); cid++) {
		CompoundFactory& compound = compounds[cid];

		//for (auto& p : compound.atomLetters) { p = 'C'; }

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

		//char letter[3] = { 'O', 'N', 'P' };

		//for (auto index : key_indices)
		//	compound.atomLetters[index] = letter[cid];
	}

	// Calc absolute ids of particles in compounds 
	if (compounds.empty()) { return; }
	compounds[0].absoluteIndexOfFirstParticle = 0;
	for (int i = 1; i < compounds.size(); i++) {
		compounds[i].absoluteIndexOfFirstParticle = compounds[i - 1].absoluteIndexOfFirstParticle + compounds[i - 1].n_particles;
	}
}


// Check that we dont have any unrealistic bonds, and warn immediately.
void VerifyBondsAreStable(const Topology& topology, float boxlen_nm, BoundaryConditionSelect bc_select, bool energyMinimizationMode) {	
	const float allowedScalar = energyMinimizationMode ? 7.f : 3.f;//1.9999f;

	for (const auto& bond : topology.singlebonds) 
	{		
		const Float3 pos1 = topology.particles[bond.global_atom_indexes[0]].position;
		const Float3 pos2 = topology.particles[bond.global_atom_indexes[1]].position;
		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos1, pos2, boxlen_nm, bc_select);
		const float bondRelaxedDist = bond.params.b0 * LIMA_TO_NANO;

		if (hyper_dist > bondRelaxedDist * allowedScalar) {
			throw std::runtime_error(std::format("Loading singlebond with illegally large dist ({}). b0: {}", hyper_dist, bond.params.b0 * LIMA_TO_NANO));
		}
		if (hyper_dist < bondRelaxedDist * 0.001)
			throw std::runtime_error(std::format("Loading singlebond with illegally small dist ({}). b0: {}", hyper_dist, bond.params.b0 * LIMA_TO_NANO));
	}
	for (const auto& bond : topology.anglebonds)
	{
		const Float3 pos1 = topology.particles[bond.global_atom_indexes[0]].position;
		const Float3 pos2 = topology.particles[bond.global_atom_indexes[1]].position;
		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(pos1, pos2, boxlen_nm, bc_select);
		if (hyper_dist < 0.001)
			throw std::runtime_error(std::format("Loading singlebond with illegally small dist ({}). b0: {}", hyper_dist, bond.params.ub0 * LIMA_TO_NANO));
	}
}
void VerifyBondsAreLegal(const std::vector<AngleBondFactory>& anglebonds, const std::vector<DihedralBondFactory>& dihedralbonds) {
	for (const auto& bond : anglebonds) {
		if (bond.global_atom_indexes[0] == bond.global_atom_indexes[1] || bond.global_atom_indexes[1] == bond.global_atom_indexes[2] || bond.global_atom_indexes[0] == bond.global_atom_indexes[2])
			throw ("Bond contains the same index multiple times: %d %d %d");
	}
	for (const auto& bond : dihedralbonds) {
		std::unordered_set<int> ids;
		for (int i = 0; i < 4; i++) {
			if (ids.contains(bond.global_atom_indexes[i]))
				throw ("Dihedral bond contains the same index multiple times");
			ids.insert(bond.global_atom_indexes[i]);
		}
	}
}


std::unique_ptr<BoxImage> LIMA_MOLECULEBUILD::buildMolecules(
	const GroFile& grofile,
	const TopologyFile& topol_file,
	VerbosityLevel vl,
	std::unique_ptr<LimaLogger> logger,
	bool ignore_hydrogens,
	const SimParams& simparams
) 
{
	auto [molecules, tinyMolecules] = PrepareMolecules(topol_file);



	LIMAForcefield forcefield{topol_file.forcefieldInclude->contents};

	const Topology topology = LoadTopology(molecules, forcefield, grofile);
	VerifyBondsAreStable(topology, grofile.box_size.x, simparams.bc_select, simparams.em_variant);
	VerifyBondsAreLegal(topology.anglebonds, topology.dihedralbonds);

	// Separate molecules into Compound-sized groups. Each group tracks which other groups it is bonded to
	const std::vector<AtomGroup> atomGroups = GroupAtoms(topology);
	//const std::vector<AtomGroup> atomGroups = GroupAtoms1(molecules, topology);	
	std::vector<CompoundFactory> compounds = CreateCompounds(topology, grofile.box_size.x, atomGroups, simparams.bc_select);

	/*Display d;
	std::vector<std::array<Float3, MAX_COMPOUND_PARTICLES>> positions;
	for (const auto& compound : compounds) {
		positions.push_back({});
		for (int i = 0; i < MAX_COMPOUND_PARTICLES; i++) {
			positions.back()[i] = compound.positions[i];
		}
	}
	std::vector<Compound> compounds2;
	for (const auto& compound : compounds) {
		compounds2.push_back(compound);
	}	
	d.Render(std::make_unique<Rendering::CompoundsTask>(compounds2, positions, grofile.box_size), true);*/

	//printf("%d compounds\n", compounds.size());

	const std::vector<ParticleToCompoundMapping> particleToCompoundidMap = MakeParticleToCompoundidMap(compounds, topology.particles.size());
	std::vector<BridgeFactory> bridges = CreateBridges(topology, compounds, particleToCompoundidMap);



	auto bpLutManager = std::make_unique<BondedParticlesLUTManagerFactory>(compounds.size());
	DistributeBondsToCompoundsAndBridges(topology, grofile.box_size.x, particleToCompoundidMap, simparams.bc_select, compounds, bridges, *bpLutManager);


	for (const auto& bridge : bridges) {
		for (int bid = 0; bid < bridge.n_anglebonds; bid++) {
			//for (const auto& bond : bridge.anglebonds) {
			auto bond = bridge.anglebonds[bid];
			std::unordered_set<int> ids;
			for (auto& lid : bond.atom_indexes) {
				if (ids.contains(lid))
					throw ("Angle bond contains the same index multiple times");
				ids.insert(lid);
			}
		}
	}

	//bpLutManager->get(0, 0)->printMatrix(compounds.begin()->n_particles);

	CalcCompoundMetaInfo(grofile.box_size.x, compounds, simparams.bc_select);


	const std::vector<TinyMolFactory> tinyMols = LoadTinyMols(grofile, tinyMolecules, forcefield);
	const int totalCompoundParticles = std::accumulate(compounds.begin(), compounds.end(), 0, [](int sum, const auto& compound) { return sum + compound.n_particles; });


	std::vector<CompoundBridge> compoundBridges(bridges.size());
	for (int i = 0; i < bridges.size(); i++) {
		compoundBridges[i] = bridges[i];
	}

	return std::make_unique<BoxImage>(
		std::move(compounds),
		static_cast<int>(totalCompoundParticles),
		bpLutManager->Finish(),
		std::move(compoundBridges),
		std::move(tinyMols),
		grofile,	// TODO: wierd ass copy here. Probably make the input a sharedPtr?
		forcefield.GetActiveLjParameters(),
		forcefield.GetTinymolTypes(),
		topology,
		forcefield.GetNonbondedInteractionParams()
	);

}
