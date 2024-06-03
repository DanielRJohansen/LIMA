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
//struct TopologyAtom {
//	const std::string type;			// As it is read in the top/itp file
//	const int global_residue_id;	// Given by limi
//};

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

	void createCompounds(const Topology& topology, float boxlen_nm);
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

	createCompounds(topology, boxlen_nm);
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
		//std::move(particleinfotable),
		std::vector<AtomRef>{},
		ParsedGroFile{ grofile },
		ForceField_NB()
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
			assert(row.words.size() >= 8);
			Atom atom;
			atom.global_id = std::stoi(row.words[0]);
			atom.gro_id = std::stoi(row.words[1]);
			atom.chain_id = std::stoi(row.words[2]);
			atom.res_id = std::stoi(row.words[3]);
			atom.atomtype_id = std::stoi(row.words[4]);
			atom.atomname = row.words[5];
			atom.unique_res_id = std::stoi(row.words[6]);
			atom.charge = std::stof(row.words[7]);

			// We no longer have atomtype as a string, so jsut use name twice, it doesnt matter from here on
			//Atom atom{ global_id, gro_id, chain_id, residue_groid, atomname, atomname, unique_res_id, charge };
			//atom.atomtype_id = atomtype_id;
			atominfotable.emplace_back(ParticleInfo{ atom });
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


void MoleculeBuilder::createCompounds(const Topology& topology, float boxlen_nm) {
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
			const bool is_bonded_with_previous_residue = residue_index == 0 || areBonded(residues[residue_index - 1], residue, particleinfotable, topology.singlebonds);
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
				particleinfo.atomname[0],
				//Forcefield::atomTypeToIndex(particleinfo.atomname[0]),
				atom_gid,
				boxlen_nm, bc_select,
				particleinfotable[atom_gid].charge
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
			BondedParticlesLUT* lut = bplut_man->get(pinfo_self.compound_index, pinfo_other.compound_index, true);
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
void CompoundFactory::AddBond(const std::vector<AtomRef>& atoms, const SingleBondFactory& bond) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add singlebond to compound");
	}
	singlebonds[n_singlebonds++] = SingleBond(
		{
			atoms[bond.global_atom_indexes[0]].localIdInCompound,
			atoms[bond.global_atom_indexes[1]].localIdInCompound
		},
		bond.b0,
		bond.kb
		);
}
void CompoundFactory::AddBond(const std::vector<AtomRef>& atoms, const AngleBondFactory& bond) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add anglebond to compound");
	}
	anglebonds[n_anglebonds++] = AngleBond(
		{
			atoms[bond.global_atom_indexes[0]].localIdInCompound,
			atoms[bond.global_atom_indexes[1]].localIdInCompound,
			atoms[bond.global_atom_indexes[2]].localIdInCompound
		},
		bond.theta_0,
		bond.k_theta
		);
}
void CompoundFactory::AddBond(const std::vector<AtomRef>& atoms, const DihedralBondFactory& bond) {
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add dihedralbond to compound");
	}
	dihedrals[n_dihedrals++] = DihedralBond(
		{
			atoms[bond.global_atom_indexes[0]].localIdInCompound,
			atoms[bond.global_atom_indexes[1]].localIdInCompound,
			atoms[bond.global_atom_indexes[2]].localIdInCompound,
			atoms[bond.global_atom_indexes[3]].localIdInCompound
		},
		bond.phi_0,
		bond.k_phi,
		bond.n
		);
}
void CompoundFactory::AddBond(const std::vector<AtomRef>& atoms, const ImproperDihedralBondFactory& bond) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND) {
		throw std::runtime_error("Failed to add improperdihedralbond to compound");
	}
	impropers[n_improperdihedrals++] = ImproperDihedralBond(
		{
			atoms[bond.global_atom_indexes[0]].localIdInCompound,
			atoms[bond.global_atom_indexes[1]].localIdInCompound,
			atoms[bond.global_atom_indexes[2]].localIdInCompound,
			atoms[bond.global_atom_indexes[3]].localIdInCompound
		},
		bond.psi_0,
		bond.k_psi
		);
}


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



void BridgeFactory::AddBond(std::vector<AtomRef>& particle_info, const SingleBondFactory& bond) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add singlebond to bridge"); }
	singlebonds[n_singlebonds++] = SingleBond{
		{
			getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[0]]),
			getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[1]]),
		},
		bond.b0,
		bond.kb
	};
}
void BridgeFactory::AddBond(std::vector<AtomRef>& particle_info, const AngleBondFactory& bond) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add anglebond to bridge"); }
	anglebonds[n_anglebonds++] = AngleBond{
		{
			getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[0]]),
			getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[1]]),
			getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[2]]),
		},		
		bond.theta_0,
		bond.k_theta
	};
}
void BridgeFactory::AddBond(std::vector<AtomRef>& particle_info, const DihedralBondFactory& bond) {
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_BRIDGE) {
		throw std::runtime_error("Failed to add dihedralbond to bridge");
	}
	dihedrals[n_dihedrals++] = DihedralBond{
		{
		getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[0]]),
		getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[1]]),
		getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[2]]),
		getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[3]]),
		},
		bond.phi_0,
		bond.k_phi,
		bond.n
	};
}
void BridgeFactory::AddBond(std::vector<AtomRef>& particle_info, const ImproperDihedralBondFactory& bond) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add improperdihedralbond to bridge"); }
	impropers[n_improperdihedrals++] = ImproperDihedralBond{
		{
			getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[0]]),
			getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[1]]),
			getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[2]]),
			getBridgelocalIdOfParticle(particle_info[bond.global_atom_indexes[3]]),
		},
		bond.psi_0,
		bond.k_psi,
	};
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
uint8_t BridgeFactory::getBridgelocalIdOfParticle(AtomRef& particle_info) {

	// Assign bridge id to particle if it doesnt have one
	if (particle_info.bridgeId == -1) {
		particle_info.bridgeId = this->bridge_id;
	}
	// If it has another id, fail as we dont support that for now
	else if (particle_info.bridgeId != this->bridge_id) {
		//printf("Particle global id %d\n", particle_info.globalId);
		throw std::runtime_error(std::format("Cannot add particle to this bridge ({}) as it already has another bridge ({})", bridge_id, particle_info.bridgeId).c_str());
	}

	// Particle already has a local id in the bridge
	if (particle_info.localIdInBridge == 255) {	// This limits the amount of particles in a bridge
		addParticle(particle_info);
	}

	return particle_info.localIdInBridge;
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
void BridgeFactory::addParticle(AtomRef& particle_info) {	// This can be made const, and return it's iD in bridge
	if (n_particles == MAX_PARTICLES_IN_BRIDGE) { throw std::runtime_error("Failed to add particle to bridge"); }

	particle_info.localIdInBridge = n_particles;

	int compoundlocalid_in_bridge = -1;
	for (int i = 0; i < this->n_compounds; i++) {
		if (particle_info.compoundId == compound_ids[i])
			compoundlocalid_in_bridge = i;
	}
	if (compoundlocalid_in_bridge == -1) {
		throw std::runtime_error("Could not find compoundlocalid_in_bridge");
	}

	particle_refs[n_particles] = ParticleReference{
		particle_info.compoundId,
		particle_info.localIdInCompound,
#ifdef LIMAKERNELDEBUGMODE
		particle_info.global_id,
#endif
		static_cast<uint8_t>(compoundlocalid_in_bridge)
	};

	n_particles++;
}












































struct TopologyFileRef {
	TopologyFileRef(int offset, const ParsedTopologyFile& topol) : atomsOffset(offset), topology(topol) {}
	int atomsOffset = 0;
	const ParsedTopologyFile& topology;
};

std::vector<TopologyFileRef> PrepareTopologies(const ParsedTopologyFile& topol) {
	std::vector<TopologyFileRef> topologies;
	int offset = 0;

	topologies.push_back({ offset, topol });
	offset += topol.GetLocalAtoms().size();
	for (const auto& molecule : topol.GetAllSubMolecules()) {
		topologies.emplace_back(TopologyFileRef{offset, *molecule.includeTopologyFile });
		offset += molecule.includeTopologyFile->GetLocalAtoms().size();
	}
	return topologies;
}



std::vector<AtomRef> PrepareAtoms(const std::vector<TopologyFileRef>& topologyFiles, const ParsedGroFile& grofile, LIMAForcefield& forcefield, int nNonsolventAtoms) {
	std::vector<AtomRef> atomRefs(nNonsolventAtoms);
	
	int globalIndex = 0;
	int uniqueResId = 0;
	for (const auto& topology : topologyFiles) {
		bool atomsAddedThisTop = false;
		for (const auto& atom : topology.topology.GetLocalAtoms()) {
			atomsAddedThisTop = true;
			// TODO: Handle ignore atoms (H) somehow?			
			// TODO: When accessing gro file, first check that it is not a solvent

			if (globalIndex > 0 && atomRefs[globalIndex - 1].topAtom->resnr != atom.resnr) {
				uniqueResId++;
			}

			atomRefs[globalIndex] = AtomRef{ &grofile.atoms[globalIndex], &atom, forcefield.GetActiveLjParameterIndex(atom.type), uniqueResId };

			
			globalIndex++;		
		}

		if (atomsAddedThisTop) { uniqueResId++; }
	}

	// TODO: At this point we still need to load the solvents, they should come now

	return atomRefs;
}
std::vector<Float3> LoadSolventPositions(const ParsedGroFile& grofile) {
	std::vector<Float3> solventPositions;// (grofile.atoms.size() - nNonsolventAtoms);
	solventPositions.reserve(grofile.atoms.size());

	for (const auto& atom : grofile.atoms)
		if (atom.residue_name == "SOL" && atom.atom_name[0] == 'O')
			solventPositions.emplace_back(atom.position);
	
	return solventPositions;
}


template <typename BondtypeTopology, typename Bondtype, typename BondtypeFactory>
void LoadBondIntoTopology(const BondtypeTopology& bondTopol, int atomIdOffset, LIMAForcefield& forcefield, const std::vector<AtomRef>& atomRefs, std::vector<BondtypeFactory>& topology)
{
	std::array<uint32_t, bondTopol.n> globalIds{};
	for (int i = 0; i < bondTopol.n; i++) {
		globalIds[i] = bondTopol.ids[i] + atomIdOffset;
		if (globalIds[i] >= atomRefs.size())
			return;
			//throw std::runtime_error("Bond refers to atom that has not been loaded");
	}

	std::array<std::string, bondTopol.n> atomTypenames{};
	for (int i = 0; i < bondTopol.n; i++) {
		atomTypenames[i] = atomRefs[globalIds[i]].topAtom->type;
	}

	const Bondtype& bondparameter = forcefield.GetBondParameters<Bondtype>(atomTypenames); //forcefield.GetSinglebondParameters(atomTypenames);
	topology.emplace_back(BondtypeFactory(globalIds, bondparameter.ToStandardBondRepresentation()));
	//return BondtypeFactory( globalIds, bondparameter.ToStandardBondRepresentation());
}

Topology LoadTopology(const std::vector<TopologyFileRef>& topologyFiles, LIMAForcefield& forcefield, const std::vector<AtomRef>& atomRefs) {
	Topology topology;

	for (const auto& topologyFile : topologyFiles) {
		for (const auto& bondTopol : topologyFile.topology.GetLocalSinglebonds()) {			
			LoadBondIntoTopology<ParsedTopologyFile::SingleBond, Singlebondtype, SingleBondFactory>(
				bondTopol, topologyFile.atomsOffset, forcefield, atomRefs, topology.singlebonds);
		}

		for (const auto& bondTopol : topologyFile.topology.GetLocalAnglebonds()) {			
			LoadBondIntoTopology<ParsedTopologyFile::AngleBond, Anglebondtype, AngleBondFactory>(
				bondTopol, topologyFile.atomsOffset, forcefield, atomRefs, topology.anglebonds);
		}

		for (const auto& bondTopol : topologyFile.topology.GetLocalDihedralbonds()) {			
			LoadBondIntoTopology<ParsedTopologyFile::DihedralBond, Dihedralbondtype, DihedralBondFactory>(
				bondTopol, topologyFile.atomsOffset, forcefield, atomRefs, topology.dihedralbonds);
		}

		for (const auto& bondTopol : topologyFile.topology.GetLocalImproperDihedralbonds()) {
			LoadBondIntoTopology<ParsedTopologyFile::ImproperDihedralBond, Improperdihedralbondtype, ImproperDihedralBondFactory>(
				bondTopol, topologyFile.atomsOffset, forcefield, atomRefs, topology.improperdihedralbonds);
		}
	}

	return topology;
}



std::vector<std::vector<int>> MapAtomsToSinglebonds(const std::vector<AtomRef>& preparedAtoms, const Topology& topology) {
	std::vector<std::vector<int>> atomIdToSinglebondsMap(preparedAtoms.size());


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
};

std::vector<AtomGroup> GroupAtoms(const std::vector<AtomRef>& preparedAtoms, const Topology& topology) {
	const int nResidues = preparedAtoms.empty() ? 0 : preparedAtoms.back().uniqueResId + 1;
	std::vector<AtomGroup> atomGroups;
	atomGroups.reserve(nResidues);

	int currentResidueId = -1;

	for (int particleId = 0; particleId < preparedAtoms.size(); particleId++) {
		const auto atom = preparedAtoms[particleId];
		if (atomGroups.empty() || atom.topAtom->section_name.has_value() || atom.uniqueResId != preparedAtoms[particleId-1].uniqueResId)
			atomGroups.push_back(AtomGroup{});
		atomGroups.back().atomIds.push_back(particleId);
		//atomGroups[atom.uniqueResId].atomIds.push_back(particleId);
	}

	return atomGroups;
}	


bool areBonded(const AtomGroup& left, const AtomGroup& right, const std::vector<AtomRef>& preparedAtoms, const std::vector<std::vector<int>>& atomIdToSinglebondsMap) {
	//assert(left != right);

	for (auto& atomleft_gid : left.atomIds) {
		for (auto& atomright_gid : right.atomIds) {
			const std::vector<int>& atomleft_singlesbonds = atomIdToSinglebondsMap[atomleft_gid];
			const std::vector<int>& atomright_singlesbonds = atomIdToSinglebondsMap[atomright_gid];
			

			for (int i = 0; i < atomleft_singlesbonds.size(); i++) {
				for (int j = 0; j < atomright_singlesbonds.size(); j++) {
					if (atomleft_singlesbonds[i] == atomright_singlesbonds[j]) {
						return true;
					}
				}
			}
			//// If any singlebond id's match, return true
			//// Check if any singlebond id from left matches any in right
			//if (std::any_of(atomleft_singlesbonds.begin(), atomleft_singlesbonds.end(), [&](int id) {
			//	return std::find(atomright_singlesbonds.begin(), atomright_singlesbonds.end(), id) != atomright_singlesbonds.end();
			//	})) {
			//	return true; // Return true if any singlebond id's match
			//}
		}
	}
	return false;
}

std::vector<CompoundFactory> CreateCompounds(const Topology& topology, float boxlen_nm, const std::vector<AtomGroup>& residues,
	std::vector<AtomRef>& preparedAtoms, const std::vector<std::vector<int>>& atomIdToSinglebondsMap, BoundaryConditionSelect bc_select)
{
	std::vector<CompoundFactory> compounds;
	int current_residue_id = -1;

	for (int residue_index = 0; residue_index < residues.size(); residue_index++) {
		const AtomGroup& residue = residues[residue_index];

		const bool is_bonded_with_previous_residue = residue_index > 0 && areBonded(residues[residue_index - 1], residue, preparedAtoms, atomIdToSinglebondsMap);
		const bool compound_has_room_for_residue = residue_index > 0 && compounds.back().hasRoomForRes(residue.atomIds.size());

		// If we are either a new molecule, or same molecule but the current compound has no more room, make new compound
		if (residue_index == 0 || !is_bonded_with_previous_residue || !compound_has_room_for_residue) {
			if (compounds.size() >= MAX_COMPOUNDS) {
				throw std::runtime_error(std::format("Cannot handle more than {} compounds", MAX_COMPOUNDS).c_str());
			}

			compounds.push_back(CompoundFactory{ static_cast<int>(compounds.size()) });

			if (is_bonded_with_previous_residue) {
				compounds.back().indicesOfBondconnectedCompounds.push_back(compounds.size() - 2);
				compounds[compounds.size() - 2].indicesOfBondconnectedCompounds.push_back(compounds.size() - 1);
			}
		}
		

		// Add all atoms of residue to current compound
		for (int i = 0; i < residue.atomIds.size(); i++) {
			const int atom_gid = residue.atomIds[i];
			compounds.back().addParticle(
				preparedAtoms[atom_gid].groAtom->position,
				preparedAtoms[atom_gid].activeLjtypeParameterIndex,
				preparedAtoms[atom_gid].topAtom->type[0],
				atom_gid,
				boxlen_nm, 
				bc_select,
				preparedAtoms[atom_gid].topAtom->charge * elementaryChargeToCoulombPerMole
			);

			preparedAtoms[atom_gid].compoundId = compounds.size() - 1;
			preparedAtoms[atom_gid].localIdInCompound = compounds.back().n_particles - 1;
		}
	}

	return compounds;
}


std::vector<BridgeFactory> CreateBridges(const std::vector<SingleBondFactory>& singlebonds, const std::vector<CompoundFactory>& compounds, const std::vector<AtomRef>& atoms)
{
	std::vector<BridgeFactory> compoundBridges;


	// First find which ompounds are bonded to other compounds. I assume singlebonds should find all these, 
	// since angles and dihedrals make no sense if the singlebond is not present, as the distances between the atoms can then be extremely large
	std::vector<std::unordered_set<int>> compound_to_bondedcompound_table(compounds.size());
	for (const auto& bond : singlebonds) {
		const int particle_gid_left = bond.global_atom_indexes[0];
		const int particle_gid_right = bond.global_atom_indexes[1];

		const int compound_gid_left = atoms[particle_gid_left].compoundId;
		const int compound_gid_right = atoms[particle_gid_right].compoundId;

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
		if (!compoundsAreAlreadyConnected(bridge_compound_ids, compoundBridges)) {
			compoundBridges.push_back(BridgeFactory{ static_cast<int>(compoundBridges.size()), bridge_compound_ids });

		}
	}

	return compoundBridges;
}









template <int n_ids>
bool SpansTwoCompounds(const uint32_t* bond_globalids, const std::vector<AtomRef>& atoms) {
	const int compoundid_left = atoms[bond_globalids[0]].compoundId;

	for (int i = 1; i < n_ids; i++) {
		if (compoundid_left != atoms[bond_globalids[i]].compoundId) {
			return true;
		}
	}
	return false;
}

template <int n_ids>
void DistributeLJIgnores(BondedParticlesLUTManager* bplut_man, const std::vector<AtomRef>& atoms, const uint32_t* global_ids) {
	for (int i = 0; i < n_ids; i++) {
		auto gid_self = global_ids[i];
		for (int j = 0; j < n_ids; j++) {
			auto gid_other = global_ids[j];

			if (gid_self == gid_other) { continue; }
		
			const int compoundIndexSelf = atoms[gid_self].compoundId;
			const int compoundIndexOther = atoms[gid_other].compoundId;


			bplut_man->addNewConnectedCompoundIfNotAlreadyConnected(compoundIndexSelf, compoundIndexOther);
			BondedParticlesLUT* lut = bplut_man->get(compoundIndexSelf, compoundIndexOther, true);
			lut->set(atoms[gid_self].localIdInCompound, atoms[gid_other].localIdInCompound, true);
		}
	}
}

// Returns two compound id's of a bond in a bridge. The id's are sorted with lowest first
template <int N>
std::array<int, 2> getTheTwoDifferentIds(const uint32_t particle_global_ids[N], const std::vector<AtomRef>& atoms) {
	std::array<int, 2> compound_ids = { atoms[particle_global_ids[0]].compoundId, -1 };

	for (int i = 1; i < N; i++) {
		int compoundid_of_particle = atoms[particle_global_ids[i]].compoundId;
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
BridgeFactory& getBridge(std::vector<BridgeFactory>& bridges, const uint32_t particle_global_ids[N], const std::vector<AtomRef>& atoms,
	std::unordered_map<int, std::vector<BridgeFactory*>>& compoundToBridges
)
{
	std::array<int, 2> compound_ids = getTheTwoDifferentIds<N>(particle_global_ids, atoms);

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

	throw std::runtime_error(std::format("Failed to find the bridge ({})", N).c_str());
}

template <typename BondTypeFactory>
void DistributeBondsToCompoundsAndBridges(const std::vector<BondTypeFactory>& bonds, std::vector<AtomRef>& particleinfotable, std::vector<CompoundFactory>& compounds, std::vector<BridgeFactory>& bridges,
	BondedParticlesLUTManager& bpLutManager, std::unordered_map<int, std::vector<BridgeFactory*>>& compoundToBridges) {
	
	constexpr int atoms_in_bond = BondTypeFactory::n_atoms;

	for (const BondTypeFactory& bond : bonds) {

		if (SpansTwoCompounds<atoms_in_bond>(bond.global_atom_indexes, particleinfotable)) {
			BridgeFactory& bridge = getBridge<atoms_in_bond>(bridges, bond.global_atom_indexes, particleinfotable, compoundToBridges);
			bridge.AddBond(particleinfotable, bond);
		}
		else {
			const int compound_id = particleinfotable[bond.global_atom_indexes[0]].compoundId;	// Pick compound using first particle in bond
			CompoundFactory& compound = compounds.at(compound_id);
			compound.AddBond(particleinfotable, bond);
		}

		DistributeLJIgnores<atoms_in_bond>(&bpLutManager, particleinfotable, bond.global_atom_indexes);
	}
}

void DistributeBondsToCompoundsAndBridges(const Topology& topology, float boxlen_nm, std::vector<AtomRef>& atoms, BoundaryConditionSelect bc_select, 
	std::vector<CompoundFactory>& compounds, std::vector<BridgeFactory>& bridges, BondedParticlesLUTManager& bpLutManager)
{
	std::unordered_map<int, std::vector<BridgeFactory*>> compoundToBridges{};

	// First check that we dont have any unrealistic bonds, and warn immediately.
	for (const auto& bond : topology.singlebonds) {
		int gid1 = bond.global_atom_indexes[0];
		int gid2 = bond.global_atom_indexes[1];

		const Float3 pos1 = atoms[gid1].groAtom->position;
		const Float3 pos2 = atoms[gid2].groAtom->position;
		const float hyper_dist = LIMAPOSITIONSYSTEM::calcHyperDistNM(&pos1, &pos2, boxlen_nm, bc_select);
		if (hyper_dist > bond.b0 * LIMA_TO_NANO * 2.f) {
			throw std::runtime_error(std::format("Loading singlebond with illegally large dist ({}). b0: {}", hyper_dist, bond.b0 * LIMA_TO_NANO).c_str());
		}
	}

	DistributeBondsToCompoundsAndBridges(topology.singlebonds, atoms, compounds, bridges, bpLutManager, compoundToBridges);
	//m_logger->print(std::format("Added {} singlebonds to molecule\n", topology.singlebonds.size()));
	DistributeBondsToCompoundsAndBridges(topology.anglebonds, atoms, compounds, bridges, bpLutManager, compoundToBridges);
	//m_logger->print(std::format("Added {} anglebonds to molecule\n", topology.anglebonds.size()));
	DistributeBondsToCompoundsAndBridges(topology.dihedralbonds, atoms, compounds, bridges, bpLutManager, compoundToBridges);
	//m_logger->print(std::format("Added {} dihedralbonds to molecule\n", topology.dihedralbonds.size()));
	DistributeBondsToCompoundsAndBridges(topology.improperdihedralbonds, atoms, compounds, bridges, bpLutManager, compoundToBridges);
	//m_logger->print(std::format("Added {} improper dihedralbonds to molecule\n", topology.improperdihedralbonds.size()));


	// While are particle is never bonded to itself, we are not allowed to calc LJ with itself
	// so we can doubly use this lut to avoid doing that

	for (int com_id = 0; com_id < compounds.size(); com_id++) {
		auto* lut = bpLutManager.get(com_id, com_id);
		for (int pid = 0; pid < MAX_COMPOUND_PARTICLES; pid++) {
			lut->set(pid, pid, true);
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
	// Solvents are not present in top file, so we can't require these to match
	const int nNonsolventAtoms = topol_file.GetElementCount<ParsedTopologyFile::AtomsEntry>();
	assert(gro_file.atoms.size() >= nNonsolventAtoms);
	const std::vector<TopologyFileRef> preparedTopologyFiles = PrepareTopologies(topol_file);

	LIMAForcefield forcefield{};

	std::vector<AtomRef> preparedAtoms = PrepareAtoms(preparedTopologyFiles, gro_file, forcefield, nNonsolventAtoms);

	Topology topology = LoadTopology(preparedTopologyFiles, forcefield, preparedAtoms);

	std::vector<std::vector<int>> atomIdToSinglebondsMap = MapAtomsToSinglebonds(preparedAtoms, topology);

	std::vector<AtomGroup> atomGroups = GroupAtoms(preparedAtoms, topology);


	std::vector<CompoundFactory> compounds = CreateCompounds(topology, gro_file.box_size.x, atomGroups, preparedAtoms, atomIdToSinglebondsMap, bc_select);
	std::vector<BridgeFactory> bridges = CreateBridges(topology.singlebonds, compounds, preparedAtoms);

	auto bpLutManager = std::make_unique<BondedParticlesLUTManager>();
	DistributeBondsToCompoundsAndBridges(topology, gro_file.box_size.x, preparedAtoms, bc_select, compounds, bridges, *bpLutManager);

	CalcCompoundMetaInfo(gro_file.box_size.x, compounds, bc_select);





	LimaForcefieldBuilder::buildForcefield(molecule_dir, molecule_dir, topol_file, EnvMode::Headless);	// Doesnt actually use .gro now..


	auto bridges_compact = std::make_unique<CompoundBridgeBundleCompact>(
		std::vector<CompoundBridge>(bridges.begin(), bridges.end())
	);



	std::vector<Float3> solventPositions = LoadSolventPositions(gro_file);











	//Forcefield ffold(VerbosityLevel::V1, molecule_dir);
	//ForceField_NB nbffold = ffold.getNBForcefield();
	//ForceField_NB nbffnew = forcefield.GetActiveLjParameters();

	//bool matchFound = false;

	//for (int i = 0; i < MAX_ATOM_TYPES; i++) {
	//	matchFound = false; // Reset for each i
	//	for (int j = 0; j < MAX_ATOM_TYPES; j++) {
	//		if (nbff.particle_parameters[i].mass == nbff2.particle_parameters[j].mass &&
	//			nbff.particle_parameters[i].sigma == nbff2.particle_parameters[j].sigma &&
	//			nbff.particle_parameters[i].epsilon == nbff2.particle_parameters[j].epsilon) {
	//			matchFound = true;
	//			break; // Break inner loop if a match is found
	//		}
	//	}
	//	if (!matchFound) {
	//		int c = 0; // Set c only if no match found
	//		// Additional actions if no match is found can be added here
	//	}
	//}

	//MoleculeBuilder mb{ std::move(logger), bc_select, vl };
	//auto boximage = mb.buildMolecules(gro_file, molecule_dir, ignore_hydrogens);



	//for (int cid = 0; cid < compounds.size(); cid++) {
	//	for (int pid = 0; pid < compounds[0].n_particles; pid++) {

	//		if (nbffnew.particle_parameters[compounds[cid].atom_types[pid]].mass != nbffold.particle_parameters[boximage->compounds[cid].atom_types[pid]].mass
	//			|| nbffnew.particle_parameters[compounds[cid].atom_types[pid]].sigma != nbffold.particle_parameters[boximage->compounds[cid].atom_types[pid]].sigma
	//			|| nbffnew.particle_parameters[compounds[cid].atom_types[pid]].epsilon != nbffold.particle_parameters[boximage->compounds[cid].atom_types[pid]].epsilon
	//			)
	//			int c = 0;

	//	}
	//}




	auto a = std::make_unique<BoxImage>(
		std::move(compounds),
		static_cast<int>(preparedAtoms.size()),
		std::move(bpLutManager),
		std::move(bridges_compact),
		std::move(solventPositions),
		gro_file.box_size.x,	// TODO: Find a better way..
		std::move(preparedAtoms),
		ParsedGroFile{ gro_file },	// TODO: wierd ass copy here
		forcefield.GetActiveLjParameters()
	);




	return a;


	//assert(compounds.size() == boximage->compounds.size());
	//for (int i = 0; i < compounds.size(); i++) {
	//	if (compounds[i].n_particles != boximage->compounds[i].n_particles) {
	//		auto a = compounds[i].n_particles;
	//		auto b = boximage->compounds[i].n_particles;
	//		int c = 0;
	//	}
	//}

	//assert(bridges.size() == boximage->bridgebundle->n_bridges);
	//for (int i = 0; i < bridges.size(); i++) {
	//	assert(bridges[i].n_particles == boximage->bridgebundle->compound_bridges[i].n_particles);
	//	for (int j = 0; j < bridges[i].n_anglebonds; j++) {
	//		assert(bridges[i].anglebonds[j].atom_indexes[2] == boximage->bridgebundle->compound_bridges[i].anglebonds[j].atom_indexes[2]);
	//	}
	//}

	//for (int i = 0; i < preparedAtoms.size(); i++) {
	//	if (boximage->particleinfotable[i].compound_index != preparedAtoms[i].compoundId
	//		|| boximage->particleinfotable[i].local_id_compound != preparedAtoms[i].localIdInCompound) {

	//		auto a = boximage->particleinfotable[i];
	//		auto b = preparedAtoms[i];
	//		int c = 0;
	//	}
	//}

	//return boximage;
}
