#include "CompoundBuilder.h"
#include "Printer.h"
#include "Forcefield.cuh"
#include "ForcefieldMaker.h"

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <format>
#include <span>
#include <array>

using namespace LIMA_Print;

namespace lfs = Filehandler;

// we need this struct because only the topology has a distinct ordering of which atoms (purely position based) belong to which residue
// This is so fucked up insane, but its the standard so :(
struct TopologyAtom {
	const std::string type;			// As it is read in the top/itp file
	const int global_residue_id;	// Given by limi
};

struct Topology {
	std::vector<SingleBondFactory> singlebonds;
	std::vector<AngleBondFactory> anglebonds;
	std::vector<DihedralBondFactory> dihedralbonds;
	std::vector<ImproperDihedralBondFactory> improperdihedralbonds;
};


class MoleculeBuilder {
public:
	MoleculeBuilder(std::unique_ptr<LimaLogger> logger, BoundaryConditionSelect bc, VerbosityLevel vl = SILENT) :
		m_logger(std::move(logger)),
		bc_select(bc),
		verbosity_level(vl)
	{}

	// Upon creation of the BoxImage, the MoleculeBuilder object is no longer valid,
	// as it move's it's data instead of copying!
	std::unique_ptr<BoxImage> buildMolecules(const ParsedGroFile& gro_file, const string& topol_path, bool ignore_hydrogens = true);

private:

	// Members, seen in the order they are filled in
	std::vector<Residue> residues;
	std::vector<Float3> nonsolvent_positions;
	std::vector<Float3> solvent_positions;
	//int32_t n_particles_in_residues = 0;

	ParticleInfoTable particleinfotable;


	std::vector<CompoundFactory> compounds;
	std::vector<BridgeFactory> compound_bridges;

	std::unique_ptr<BondedParticlesLUTManager> bp_lut_manager;

	const BoundaryConditionSelect bc_select;

	std::unique_ptr<LimaLogger> m_logger;
	const VerbosityLevel verbosity_level;

	std::unordered_map<int, std::vector<BridgeFactory*>> compoundToBridgesMap;

	// ------------------------------------ HELPER FUNCTIONS ------------------------------------ //



	ParticleInfoTable loadAtomInfo(const std::string& molecule_dir);

	Topology loadTopology(const std::string& molecule_dir);


	// Only works for fixed positions gro files!! Meaning max, 99.999 particles
	void loadAtomPositions( const ParsedGroFile& grofile);

	/// <summary>
	/// Goes through all the residues, and fills the bondedresidue_ids vector. 
	/// The particle_bonds_lut is used to determine whether particles of two different 
	/// share a bond.
	/// </summary>
	void matchBondedResidues(const std::vector<SingleBondFactory>& singlebonds);

	void createCompounds(const std::vector<SingleBondFactory>& singlebonds, float boxlen_nm);
	void createBridges(const std::vector<SingleBondFactory>& singlebonds);

	void distributeBondsToCompoundsAndBridges(const Topology& topology, float boxlen_nm);

	template <typename BondType>
	void distributeBondsToCompoundsAndBridges(const std::vector<BondType>& bonds);

	// Find index of centermost particle, find "radius" of compound
	void calcCompoundMetaInfo(float boxlen_nm);

	void insertSolventAtom(const GroRecord& record);
};


std::unique_ptr<BoxImage> MoleculeBuilder::buildMolecules(const ParsedGroFile& grofile, const string& molecule_dir, bool ignore_hydrogens)
{
	m_logger->startSection("Building Molecules");
	particleinfotable = loadAtomInfo(molecule_dir);


	const float boxlen_nm = grofile.box_size.x;


	loadAtomPositions(grofile);
	const Topology& topology = loadTopology(molecule_dir);

	// We need each particle to ref singlebonds to determine which residues are bonded
	for (int singlebond_index = 0; singlebond_index < topology.singlebonds.size(); singlebond_index++) {
		for (uint32_t global_id : topology.singlebonds[singlebond_index].global_atom_indexes) {
			particleinfotable[global_id].singlebonds_indices.emplace_back(singlebond_index);
		}
	}

	createCompounds(topology.singlebonds, boxlen_nm);
	createBridges(topology.singlebonds);

	distributeBondsToCompoundsAndBridges(topology, boxlen_nm);

	calcCompoundMetaInfo(boxlen_nm);


	auto bridges_compact = std::make_unique<CompoundBridgeBundleCompact>(
		std::vector<CompoundBridge>(compound_bridges.begin(), compound_bridges.end())
	);

	m_logger->finishSection("Finished building molecules");


	for (const auto& c : compounds) {
		for (int i = 0; i < c.n_singlebonds; i++) {
			const SingleBond& b = c.singlebonds[i];
			if (b.atom_indexes[0] == b.atom_indexes[1]) {
				throw std::runtime_error("CompoundBuilder failed\n");
			}
		}
	}

	for (int j = 0; j < bridges_compact->n_bridges; j++) {
		auto& bridge = bridges_compact->compound_bridges[j];
		for (int i = 0; i < bridge.n_singlebonds; i++) {
			const SingleBond& b = bridge.singlebonds[i];
			if (b.atom_indexes[0] == b.atom_indexes[1]) {
				throw std::runtime_error("CompoundBuilder failed\n");
			}
		}
	}

	return std::make_unique<BoxImage>(
		std::move(compounds),
		static_cast<int>(nonsolvent_positions.size()),
		std::move(bp_lut_manager),
		std::move(bridges_compact),
		std::move(solvent_positions),
		boxlen_nm,	// TODO: Find a better way..
		std::move(particleinfotable),
		ParsedGroFile{ grofile }
	);
}

void MoleculeBuilder::loadAtomPositions(const ParsedGroFile& grofile)
{
	int current_res_id = -1;	// unique

	for (const GroRecord& atom : grofile.atoms) {
		if (atom.atom_name[0] == 'H' && IGNORE_HYDROGEN) {
			continue;
		}

		const bool entry_is_solvent = atom.residue_name == "WATER" || atom.residue_name == "SOL" || atom.residue_name == "HOH";
		if (entry_is_solvent) {
			if (atom.atom_name[0] != 'O') {
				continue;	// Treat all three solvents atoms as a single point, so no hydrogens here
			}

			// Ensure all nonsolvents have been added before adding solvents. Dunno if necessary?
			if (nonsolvent_positions.size() != particleinfotable.size()) {
				throw std::runtime_error(std::format("Trying to add solvents, but nonsolvents added ({}) does not equal the expect amount of .lff file ({})",
					nonsolvent_positions.size(), particleinfotable.size()).c_str());
			}
			solvent_positions.push_back(atom.position);
		}
		else {
			const int assumed_global_id = nonsolvent_positions.size();
			if (assumed_global_id >= particleinfotable.size()) {
				throw std::runtime_error(std::format("Trying to add more noncompounds from gro file than .lff file expects ({})", particleinfotable.size()).c_str());
			}

			// Sadly the topol file will give new gro_ids to the atoms of chains above chain_A
			/*if (record.gro_id != particleinfotable[assumed_global_id].gro_id) {
				throw std::runtime_error("gro_id of .gro file does not match that of .lff file");
			}*/

			if (atom.atom_name != particleinfotable[assumed_global_id].atomname) {
				throw std::runtime_error("atom_name of .gro file does not match that of .lff file");
			}

			const bool is_new_res = particleinfotable[assumed_global_id].unique_res_id != current_res_id
				|| atom.residue_number != residues.back().gro_id;
			if (is_new_res) {
				residues.push_back(Residue{ atom.residue_number, static_cast<int>(residues.size()), atom.residue_name, particleinfotable[assumed_global_id].chain_id });
				current_res_id = particleinfotable[assumed_global_id].unique_res_id;
			}

			residues.back().atoms_globalid.emplace_back(assumed_global_id);
			nonsolvent_positions.emplace_back(atom.position);

			if (residues.back().atoms_globalid.size() > MAX_COMPOUND_PARTICLES) {
				throw std::runtime_error(std::format("Cannot handle residue with {} particles\n", residues.back().atoms_globalid.size()).c_str());
			}
		}
	}

	if (nonsolvent_positions.size() != particleinfotable.size()) {
		throw std::runtime_error("Didn't read the expected amount of nonsolvent particles");
	}
}


Topology MoleculeBuilder::loadTopology(const std::string& molecule_dir)
{
	const string bonded_path = Filehandler::fileExists(Filehandler::pathJoin(molecule_dir, "custom_ffbonded.lff"))
		? Filehandler::pathJoin(molecule_dir, "custom_ffbonded.lff")
		: Filehandler::pathJoin(molecule_dir, "ffbonded.lff");

	MDFiles::ParsedLffFile bondedff{ bonded_path };

	Topology topology{};


	for (const auto& bond : bondedff.singlebonds.entries) {
		topology.singlebonds.emplace_back(SingleBondFactory{ bond.global_ids, bond.b0 * NANO_TO_LIMA, bond.kb / (NANO_TO_LIMA * NANO_TO_LIMA) });
	}
	for (const auto& bond : bondedff.anglebonds.entries) {
		topology.anglebonds.emplace_back(AngleBondFactory{ bond.global_ids, bond.theta0, bond.ktheta});
	}
	for (const auto& bond : bondedff.dihedralbonds.entries) {
		topology.dihedralbonds.emplace_back(DihedralBondFactory{ bond.global_ids, bond.phi0, bond.kphi, bond.n });
	}
	for (const auto& bond : bondedff.improperdihedralbonds.entries) {
		topology.improperdihedralbonds.emplace_back(ImproperDihedralBondFactory{ bond.global_ids, bond.psi0, bond.kpsi});
	}

	return topology;
}


std::vector<ParticleInfo> MoleculeBuilder::loadAtomInfo(const std::string& molecule_dir) {
	const string nonbonded_path = Filehandler::fileExists(Filehandler::pathJoin(molecule_dir, "custom_ffnonbonded.lff"))
		? Filehandler::pathJoin(molecule_dir, "custom_ffnonbonded.lff")
		: Filehandler::pathJoin(molecule_dir, "ffnonbonded.lff");

	const SimpleParsedFile nonbonded_parsed = Filehandler::parseLffFile(nonbonded_path, false);

	ParticleInfoTable atominfotable{};

	for (auto& row : nonbonded_parsed.rows) {
		if (row.section == "atomtype_map") {
			const int global_id = std::stoi(row.words[0]);
			const int gro_id = std::stoi(row.words[1]);
			const int chain_id = std::stoi(row.words[2]);
			const int residue_groid = std::stoi(row.words[3]);
			const int atomtype_id = std::stoi(row.words[4]);
			const auto& atomname = row.words[5];
			const int unique_res_id = std::stoi(row.words[6]);
			atominfotable.emplace_back(ParticleInfo{ global_id, gro_id, chain_id, atomtype_id, atomname, residue_groid, unique_res_id });
		}
	}
	return atominfotable;
}


bool areBonded(const Residue& left, const Residue& right, const ParticleInfoTable& particleinfolut, const std::vector<SingleBondFactory>& singlebonds) {
	assert(left.global_id != right.global_id);

	for (auto& atomleft_gid : left.atoms_globalid) {
		for (int bond_index : particleinfolut[atomleft_gid].singlebonds_indices) {
			const auto& atom_globalids = singlebonds[bond_index].global_atom_indexes;

			for (auto& atomright_gid : right.atoms_globalid) {
				if ((atom_globalids[0] == atomleft_gid && atom_globalids[1] == atomright_gid) ||
					(atom_globalids[1] == atomleft_gid && atom_globalids[0] == atomright_gid))
				{
					return true;
				}
			}
		}
	}
	return false;
}


void MoleculeBuilder::createCompounds(const std::vector<SingleBondFactory>& singlebonds, float boxlen_nm) {
	// Nothing to do if we have no residues
	if (residues.empty()) { return; }
	
	compounds.push_back(CompoundFactory{ 0 });

	int current_residue_id = residues[0].global_id;

	for (int residue_index = 0; residue_index < residues.size(); residue_index++) {
		const Residue& residue = residues[residue_index];

		const bool new_residue = residue.global_id != current_residue_id;
		current_residue_id = residue.global_id;

		// If necessary start new compound
		if (new_residue) {
			const bool is_bonded_with_previous_residue = residue_index == 0 || areBonded(residues[residue_index - 1], residue, particleinfotable, singlebonds);
			const bool compound_has_room_for_residue = compounds.back().hasRoomForRes(residue.atoms_globalid.size());

			// If we are either a new molecule, or same molecule but the current compound has no more room, make new compound
			if (!is_bonded_with_previous_residue || !compound_has_room_for_residue) {
				if (compounds.size() >= MAX_COMPOUNDS) {
					throw std::runtime_error(std::format("Cannot handle more than {} compounds", MAX_COMPOUNDS).c_str());
				}

				compounds.push_back(CompoundFactory{ static_cast<int>(compounds.size()) });
			}
		}

		// Add all atoms of residue to current compound
		for (auto& atom_gid : residue.atoms_globalid) {
			// First add the new information to the particle
			
			ParticleInfo& particleinfo = particleinfotable[atom_gid];

			particleinfo.compound_index = compounds.back().id;
			particleinfo.local_id_compound = compounds.back().n_particles;

			// Then add the particle to the compound
			compounds.back().addParticle(
				nonsolvent_positions[atom_gid],
				particleinfo.atomtype_id,
				Forcefield::atomTypeToIndex(particleinfo.atomname[0]),
				atom_gid,
				boxlen_nm, bc_select
				);
		}
	}

	m_logger->print(std::format("Created {} compounds\n", compounds.size()));
}

bool compoundsAreAlreadyConnected(const std::vector<int> compound_ids, const std::vector<BridgeFactory>& compoundbriges) 
{
	for (auto& bridge : compoundbriges) {

		bool allCompoundIdsPresent = true;
		for (auto& compound_id : compound_ids) {
			if (!bridge.containsCompound(compound_id)) {
				allCompoundIdsPresent = false;
			}
		}

		if (allCompoundIdsPresent) {
			return true;
		}
	}

	return false;
}


void MoleculeBuilder::createBridges(const std::vector<SingleBondFactory>& singlebonds)
{
	// First find which ompounds are bonded to other compounds. I assume singlebonds should find all these, 
	// since angles and dihedrals make no sense if the singlebond is not present, as the distances between the atoms can then be extremely large
	std::vector<std::unordered_set<int>> compound_to_bondedcompound_table(compounds.size());
	for (const auto& bond : singlebonds) {
		const int particle_gid_left = bond.global_atom_indexes[0];
		const int particle_gid_right = bond.global_atom_indexes[1];

		const int compound_gid_left = particleinfotable[particle_gid_left].compound_index;
		const int compound_gid_right = particleinfotable[particle_gid_right].compound_index;

		if (compound_gid_left != compound_gid_right) {
			compound_to_bondedcompound_table[compound_gid_left].insert(compound_gid_right);
			compound_to_bondedcompound_table[compound_gid_right].insert(compound_gid_left);
		}
	}

	// Knowing which compounds are bonded we can create the bridges.
	for (int compound_id_self = 0; compound_id_self < compound_to_bondedcompound_table.size(); compound_id_self++) {
		const std::unordered_set<int>& bonded_compounds_set = compound_to_bondedcompound_table[compound_id_self];

		// if a compound is the entire molecule, this set will be empty, and that is alright
		if (bonded_compounds_set.empty()) {
			continue;
		}

		// A compound can not be bonded to itself
		if (compound_id_self == *std::min_element(bonded_compounds_set.begin(), bonded_compounds_set.end())) {
			throw std::runtime_error("Something went wrong in the createBridges algorithm");
		}


		// Each compound is only allowed to create a bridge to compounds with a larger id
		std::vector<int> bridge_compound_ids{ compound_id_self };
		for (auto compound_id_other : bonded_compounds_set) {
			if (compound_id_other > compound_id_self) {
				bridge_compound_ids.push_back(compound_id_other);
			}			
		}

		// Note that this bridge might contain more compounds than this compound proposes to bridge, as it was made by a compound with a lower id, that is compound does not include!
		if (!compoundsAreAlreadyConnected(bridge_compound_ids, compound_bridges)) {
			compound_bridges.push_back(BridgeFactory{ static_cast<int>(compound_bridges.size()), bridge_compound_ids });

		}

		//const bool thiscompound_is_smallest_id = compound_id_self < *std::min_element(bonded_compounds_set.begin(), bonded_compounds_set.end());
		//if (thiscompound_is_smallest_id) {
		//	std::vector<int> bridge_compound_ids{ compound_id };
		//	bridge_compound_ids.insert(bridge_compound_ids.end(), bonded_compounds_set.begin(), bonded_compounds_set.end());

		//	compound_bridges.push_back(BridgeFactory{ static_cast<int>(compound_bridges.size()), bridge_compound_ids });
		//}
	}
	m_logger->print(std::format("Created {} compound bridges\n", compound_bridges.size()));
}




template <int n_ids>
bool spansTwoCompounds(const uint32_t* bond_globalids, const ParticleInfoTable& particleinfolut) {
	const int compoundid_left = particleinfolut[bond_globalids[0]].compound_index;

	for (int i = 1; i < n_ids; i++) {
		if (compoundid_left != particleinfolut[bond_globalids[i]].compound_index) {
			return true;
		}
	}
	return false;
}


// Returns two compound id's of a bond in a bridge. The id's are sorted with lowest first
template <int n_ids>
std::array<int, 2> getTheTwoDifferentIds(const uint32_t* particle_global_ids, const ParticleInfoTable& particleinfolut) {
	std::array<int, 2> compound_ids = { particleinfolut[particle_global_ids[0]].compound_index, -1 };

	for (int i = 1; i < n_ids; i++) {
		int compoundid_of_particle = particleinfolut[particle_global_ids[i]].compound_index;
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


bool bridgeContainsTheseTwoCompounds(const BridgeFactory& bridge, const std::array<int, 2> compound_ids) {
	for (const int id : compound_ids) {
		bool found = false;

		for (int i = 0; i < bridge.n_compounds; i++) {
			if (bridge.compound_ids[i] == id) {
				found = true;
				break;
			}
		}

		if (!found)
			return false;
	}

	return true;
}



template <int n_ids>
BridgeFactory& getBridge(std::vector<BridgeFactory>& bridges, const uint32_t* particle_global_ids, const ParticleInfoTable& particleinfolut, 
	std::unordered_map<int, std::vector<BridgeFactory*>>& compoundToBridges)
{
	std::array<int, 2> compound_ids = getTheTwoDifferentIds<n_ids>(particle_global_ids, particleinfolut);

	// Check if the compound_id is already associated with bridges
	auto it = compoundToBridges.find(compound_ids[0]);
	if (it != compoundToBridges.end()) {
		for (BridgeFactory* bridge : it->second) {
			if (bridgeContainsTheseTwoCompounds(*bridge, compound_ids)) {
				return *bridge;
			}
		}
	}


	// If not found in the lookup or no matching bridge, search through all bridges
	for (BridgeFactory& bridge : bridges) {
		if (bridgeContainsTheseTwoCompounds(bridge, compound_ids)) {
			// Cache the bridge reference for the first compound_id
			compoundToBridges[compound_ids[0]].push_back(&bridge);
			return bridge;
		}
	}

	throw std::runtime_error(std::format("Failed to find the bridge ({})", n_ids).c_str());
}

template <int n_ids>
void distributeLJIgnores(BondedParticlesLUTManager* bplut_man, const ParticleInfoTable& pinfo, const uint32_t* global_ids) {
	for (int i = 0; i < n_ids; i++) {
		auto gid_self = global_ids[i];
		for (int j = 0; j < n_ids; j++) {
			auto gid_other = global_ids[j];

			if (gid_self == gid_other) { continue; }

			const auto& pinfo_self = pinfo[gid_self];
			const auto pinfo_other = pinfo[gid_other];

			bplut_man->addNewConnectedCompoundIfNotAlreadyConnected(pinfo_self.compound_index, pinfo_other.compound_index);
			BondedParticlesLUT* lut = bplut_man->get(pinfo_self.compound_index, pinfo_other.compound_index);
			lut->set(pinfo_self.local_id_compound, pinfo_other.local_id_compound, true);
		}
	}
}


template <typename BondTypeFactory>
void MoleculeBuilder::distributeBondsToCompoundsAndBridges(const std::vector<BondTypeFactory>& bonds) {

	constexpr int atoms_in_bond = BondTypeFactory::n_atoms;

	for (auto& bond : bonds) {

		if (spansTwoCompounds<atoms_in_bond>(bond.global_atom_indexes, particleinfotable)) {
			BridgeFactory& bridge = getBridge<atoms_in_bond>(compound_bridges, bond.global_atom_indexes, particleinfotable, compoundToBridgesMap);
			bridge.addBond(particleinfotable, bond);
		}
		else {
			const int compound_id = particleinfotable[bond.global_atom_indexes[0]].compound_index;	// Pick compound using first particle in bond
			CompoundFactory& compound = compounds.at(compound_id);
			compound.addBond(particleinfotable, bond);
		}

		distributeLJIgnores<atoms_in_bond>(bp_lut_manager.get(), particleinfotable, bond.global_atom_indexes);
	}
}


void MoleculeBuilder::distributeBondsToCompoundsAndBridges(const Topology& topology, float boxlen_nm) {
	bp_lut_manager = std::make_unique<BondedParticlesLUTManager>();

	// First check that we dont have any unrealistic bonds, and warn immediately.
	for (const auto& bond : topology.singlebonds) {
		int gid1 = bond.global_atom_indexes[0];
		int gid2 = bond.global_atom_indexes[1];

		const Float3 pos1 = nonsolvent_positions[gid1];
		const Float3 pos2 = nonsolvent_positions[gid2];
		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(&pos1, &pos2, boxlen_nm, bc_select);
		if (hyper_dist > bond.b0 * LIMA_TO_NANO * 2.f) {
			throw std::runtime_error(std::format("Loading singlebond with illegally large dist ({}). b0: {}", hyper_dist, bond.b0 * LIMA_TO_NANO).c_str());
		}
	}

	distributeBondsToCompoundsAndBridges(topology.singlebonds);
	m_logger->print(std::format("Added {} singlebonds to molecule\n", topology.singlebonds.size()));
	distributeBondsToCompoundsAndBridges(topology.anglebonds);
	m_logger->print(std::format("Added {} anglebonds to molecule\n", topology.anglebonds.size()));
	distributeBondsToCompoundsAndBridges(topology.dihedralbonds);
	m_logger->print(std::format("Added {} dihedralbonds to molecule\n", topology.dihedralbonds.size()));
	distributeBondsToCompoundsAndBridges(topology.improperdihedralbonds);
	m_logger->print(std::format("Added {} improper dihedralbonds to molecule\n", topology.improperdihedralbonds.size()));


	// While are particle is never bonded to itself, we are not allowed to calc LJ with itself
	// so we can doubly use this lut to avoid doing that

	for (int com_id = 0; com_id < compounds.size(); com_id++) {
		auto* lut = bp_lut_manager->get(com_id, com_id);
		for (int pid = 0; pid < MAX_COMPOUND_PARTICLES; pid++) {
			lut->set(pid, pid, true);
		}
	}

	// Finally tell each compound which other compounds they are bonded to for faster lookup of LUTs
	for (const auto& bridge : compound_bridges) {
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
				float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(&positions[i], &positions[center_indices[j]], boxlen_nm, bc);
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
				float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(&positions[i], &new_centers[j], boxlen_nm, bc);
				if (dist < min_dist) {
					min_dist = dist;
					center_indices[j] = i; // Update the center index to this particle
				}
			}
		}
	}

	return center_indices; // Return the indices of particles that are final centers
}


std::array<float, CompoundInteractionBoundary::k> clusterRadii(Float3* positions, int n_particles, const std::array<Float3, CompoundInteractionBoundary::k>& key_positions, 
	float boxlen_nm, BoundaryConditionSelect bc) {
	std::array<float, CompoundInteractionBoundary::k> radii;
	std::vector<int> labels(n_particles);  // Holds the index of the closest key particle for each particle

	// Assign each particle to the closest key particle
	for (int i = 0; i < n_particles; ++i) {
		float min_dist = std::numeric_limits<float>::infinity();
		for (size_t j = 0; j < CompoundInteractionBoundary::k; ++j) {
			float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(&positions[i], &key_positions[j], boxlen_nm, bc);
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
				float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(&positions[i], &key_positions[j], boxlen_nm, bc);
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
	for (int i = 0; i < n_elems; i++) {
		Float3 pos = positions[i];
		BoundaryConditionPublic::applyHyperposNM(&positions[i], &pos, boxlen_nm, bc);	// Hyperpos around particle 0, since we dont know key position yet
		com += pos;
	}
	return com / static_cast<float>(n_elems);
}


int indexOfParticleClosestToCom(const Float3* positions, int n_elems, const Float3& com, float boxlen_nm, BoundaryConditionSelect bc) {
	int closest_particle_index = 0;
	float closest_particle_distance = std::numeric_limits<float>::infinity();
	for (int i = 0; i < n_elems; i++) {
		//const float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(&positions[i], &com);
		const float dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(&positions[i], &com, boxlen_nm, bc);
		if (dist < closest_particle_distance) {
			closest_particle_distance = dist;
			closest_particle_index = i;
		}
	}
	return closest_particle_index;
}


void MoleculeBuilder::calcCompoundMetaInfo(float boxlen_nm) {
	for (CompoundFactory& compound : compounds) {

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
	}
}


// --------------------------------------------------------------- Factory Functions --------------------------------------------------------------- //
void CompoundFactory::addBond(const ParticleInfoTable& particle_info, const SingleBondFactory& bondtype) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_COMPOUND) { 
		throw std::runtime_error("Failed to add singlebond to compound"); }
	singlebonds[n_singlebonds++] = SingleBond(
		{
			particle_info[bondtype.global_atom_indexes[0]].local_id_compound,
			particle_info[bondtype.global_atom_indexes[1]].local_id_compound,
		},
		bondtype.b0,
		bondtype.kb
	);
}

void CompoundFactory::addBond(const ParticleInfoTable& particle_info, const AngleBondFactory& bondtype) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_COMPOUND) { throw std::runtime_error("Failed to add anglebond to compound"); }
	anglebonds[n_anglebonds++] = AngleBond(
		{
			particle_info[bondtype.global_atom_indexes[0]].local_id_compound,
			particle_info[bondtype.global_atom_indexes[1]].local_id_compound,
			particle_info[bondtype.global_atom_indexes[2]].local_id_compound,
		},		
		bondtype.theta_0,
		bondtype.k_theta
	);
}

void CompoundFactory::addBond(const ParticleInfoTable& particle_info, const DihedralBondFactory& bondtype) {
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_COMPOUND) { 
		throw std::runtime_error("Failed to add dihedralbond to compound"); }
	dihedrals[n_dihedrals++] = DihedralBond(
		{
			particle_info[bondtype.global_atom_indexes[0]].local_id_compound,
			particle_info[bondtype.global_atom_indexes[1]].local_id_compound,
			particle_info[bondtype.global_atom_indexes[2]].local_id_compound,
			particle_info[bondtype.global_atom_indexes[3]].local_id_compound,
		},
		bondtype.phi_0,
		bondtype.k_phi,
		bondtype.n
	);
}

void CompoundFactory::addBond(const ParticleInfoTable& particle_info, const ImproperDihedralBondFactory& bondtype) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND) { throw std::runtime_error("Failed to add improperdihedralbond to compound"); }
	impropers[n_improperdihedrals++] = ImproperDihedralBond(
		{
			particle_info[bondtype.global_atom_indexes[0]].local_id_compound,
			particle_info[bondtype.global_atom_indexes[1]].local_id_compound,
			particle_info[bondtype.global_atom_indexes[2]].local_id_compound,
			particle_info[bondtype.global_atom_indexes[3]].local_id_compound,
		},
		bondtype.psi_0,
		bondtype.k_psi
	);
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






void BridgeFactory::addBond(ParticleInfoTable& particle_info, const SingleBondFactory& bondtype) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add singlebond to bridge"); }
	singlebonds[n_singlebonds++] = SingleBond{
		{
			getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[0]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[1]]),
		},
		bondtype.b0,
		bondtype.kb
	};
}

void BridgeFactory::addBond(ParticleInfoTable& particle_info, const AngleBondFactory& bondtype) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add anglebond to bridge"); }
	anglebonds[n_anglebonds++] = AngleBond{
		{
			getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[0]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[1]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[2]]),
		},		
		bondtype.theta_0,
		bondtype.k_theta
	};
}

void BridgeFactory::addBond(ParticleInfoTable& particle_info, const DihedralBondFactory& bondtype) {

	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_BRIDGE) { 
		throw std::runtime_error("Failed to add dihedralbond to bridge"); }
	dihedrals[n_dihedrals++] = DihedralBond{
		{
		getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[0]]),
		getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[1]]),
		getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[2]]),
		getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[3]]),
		},		
		bondtype.phi_0,
		bondtype.k_phi,
		bondtype.n
	};


}

void BridgeFactory::addBond(ParticleInfoTable& particle_info, const ImproperDihedralBondFactory& bondtype) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add improperdihedralbond to bridge"); }
	impropers[n_improperdihedrals++] = ImproperDihedralBond{
		{
			getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[0]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[1]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[2]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.global_atom_indexes[3]]),
		},
		bondtype.psi_0,
		bondtype.k_psi,
	};
}


uint8_t BridgeFactory::getBridgelocalIdOfParticle(ParticleInfo& particle_info) {

	// Assign bridge id to particle if it doesnt have one
	if (particle_info.bridge_id == -1) {
		particle_info.bridge_id = this->bridge_id;
	}
	// If it has another id, fail as we dont support that for now
	else if (particle_info.bridge_id != this->bridge_id) {
		printf("Particle global id %d\n", particle_info.global_id);
		throw std::runtime_error(std::format("Cannot add particle to this bridge ({}) as it already has another bridge ({})", bridge_id, particle_info.bridge_id).c_str());
	}

	// Particle already has a local id in the bridge
	if (particle_info.local_id_bridge == 255) {	// This limits the amount of particles in a bridge
		addParticle(particle_info);
	}

	return particle_info.local_id_bridge;
}

void BridgeFactory::addParticle(ParticleInfo& particle_info) {
	if (n_particles == MAX_PARTICLES_IN_BRIDGE) { throw std::runtime_error("Failed to add particle to bridge"); }

	particle_info.local_id_bridge = n_particles;

	int compoundlocalid_in_bridge = -1;
	for (int i = 0; i < this->n_compounds; i++) {
		if (particle_info.compound_index == compound_ids[i])
			compoundlocalid_in_bridge = i;
	}
	if (compoundlocalid_in_bridge == -1) {
		throw std::runtime_error("Could not find compoundlocalid_in_bridge");
	}

	particle_refs[n_particles] = ParticleReference{
		particle_info.compound_index,
		particle_info.local_id_compound,
#ifdef LIMAKERNELDEBUGMODE
		particle_info.global_id,
#endif
		static_cast<uint8_t>(compoundlocalid_in_bridge)
	};

	n_particles++;
}


std::unique_ptr<BoxImage> LIMA_MOLECULEBUILD::buildMolecules(
	const std::string& molecule_dir,
	//const std::string& gro_name,
	const ParsedGroFile& gro_file,
	const ParsedTopologyFile& topol_file,
	VerbosityLevel vl,
	std::unique_ptr<LimaLogger> logger,
	bool ignore_hydrogens,
	BoundaryConditionSelect bc_select
) 
{
	LimaForcefieldBuilder::buildForcefield(molecule_dir, molecule_dir, topol_file, EnvMode::Headless);	// Doesnt actually use .gro now..

	MoleculeBuilder mb{ std::move(logger), bc_select, vl };
	return mb.buildMolecules(gro_file, molecule_dir, ignore_hydrogens);
}
