#include "CompoundBuilder.h"
#include "Printer.h"
#include "EngineUtils.cuh"





#include <algorithm>
#include <format>
#include <span>

using namespace LIMA_Print;


struct GroRecord {
	int residue_number{};
	std::string residue_name{};
	std::string atom_name{};
	int gro_id{};
	Float3 position{};
	Float3 velocity{};
};

struct Topology {
	std::vector<SingleBond> singlebonds;
	std::vector<AngleBond> anglebonds;
	std::vector<DihedralBond> dihedralbonds;
	std::vector<ImproperDihedralBond> improperdihedralbonds;
};


class MoleculeBuilder {
public:
	MoleculeBuilder(Forcefield* ff, std::unique_ptr<LimaLogger>, const std::string& work_dir = "", VerbosityLevel vl = SILENT);

	// Upon creation of the CompoundCollection, the MoleculeBuilder object is no longer valid,
	// as it move's it's data instead of copying!
	CompoundCollection buildMolecules(const string& gro_path, const string& topol_path, bool ignore_hydrogens = true);

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

	std::unique_ptr<LimaLogger> m_logger;
	const VerbosityLevel verbosity_level;
	const Forcefield* forcefield;

	const std::string work_dir;

	// ------------------------------------ HELPER FUNCTIONS ------------------------------------ //



	ParticleInfoTable loadAtomInfo(const std::string& molecule_dir);

	Topology loadTopology(const std::string& molecule_dir);


	// Only works for fixed positions gro files!! Meaning max, 99.999 particles
	void loadAtomPositions(const std::string& gro_path);

	/// <summary>
	/// Goes through all the residues, and fills the bondedresidue_ids vector. 
	/// The particle_bonds_lut is used to determine whether particles of two different 
	/// share a bond.
	/// </summary>
	void matchBondedResidues(const std::vector<SingleBond>& singlebonds);

	void createCompounds(const std::vector<SingleBond>& singlebonds);
	void createBridges(const std::vector<SingleBond>& singlebonds);

	void distributeBondsToCompoundsAndBridges(const Topology& topology);

	template <typename BondType>
	void distributeBondsToCompoundsAndBridges(const std::vector<BondType>& bonds);

	// Find index of centermost particle, find "radius" of compound
	void calcCompoundMetaInfo();

	void insertSolventAtom(const GroRecord& record);
};

MoleculeBuilder::MoleculeBuilder(Forcefield* ff, std::unique_ptr<LimaLogger> logger, const std::string& work_dir, VerbosityLevel vl) :
	m_logger(std::move(logger)),
	verbosity_level(vl), 
	forcefield(ff) 
{}

CompoundCollection MoleculeBuilder::buildMolecules(const string& gro_path, const string& molecule_dir, bool ignore_hydrogens) {

	m_logger->startSection("Building Molecules");
	particleinfotable = loadAtomInfo(molecule_dir);
	loadAtomPositions(gro_path);

	const Topology& topology = loadTopology(molecule_dir);

	// We need each particle to ref singlebonds to determine which residues are bonded
	for (int singlebond_index = 0; singlebond_index < topology.singlebonds.size(); singlebond_index++) {
		for (uint32_t global_id : topology.singlebonds[singlebond_index].atom_indexes) {
			particleinfotable[global_id].singlebonds_indices.emplace_back(singlebond_index);
		}
	}

	createCompounds(topology.singlebonds);
	createBridges(topology.singlebonds);

	distributeBondsToCompoundsAndBridges(topology);

	calcCompoundMetaInfo();


	auto bridges_compact = std::make_unique<CompoundBridgeBundleCompact>(
		std::vector<CompoundBridge>(compound_bridges.begin(), compound_bridges.end())
	);

	m_logger->finishSection("Finished building molecules");

	return CompoundCollection{ 
		std::move(compounds), 
		static_cast<int>(nonsolvent_positions.size()),
		std::move(bp_lut_manager), 
		std::move(bridges_compact),
		std::move(solvent_positions) };
}








GroRecord parseGroLine(const std::string& line) {
	GroRecord record;

	// Parse residue number (5 positions, integer)
	std::istringstream(line.substr(0, 5)) >> record.residue_number;

	// Parse residue name (5 characters)
	record.residue_name = line.substr(5, 5);
	record.residue_name.erase(std::remove(record.residue_name.begin(), record.residue_name.end(), ' '), record.residue_name.end());

	// Parse atom name (5 characters)
	record.atom_name = line.substr(10, 5);
	record.atom_name.erase(std::remove(record.atom_name.begin(), record.atom_name.end(), ' '), record.atom_name.end());

	// Parse atom number (5 positions, integer)
	std::istringstream(line.substr(15, 5)) >> record.gro_id;

	// Parse position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
	std::istringstream(line.substr(20, 8)) >> record.position.x;
	std::istringstream(line.substr(28, 8)) >> record.position.y;
	std::istringstream(line.substr(36, 8)) >> record.position.z;

	if (line.size() >= 68) {
		// Parse velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
		std::istringstream(line.substr(44, 8)) >> record.velocity.x;
		std::istringstream(line.substr(52, 8)) >> record.velocity.y;
		std::istringstream(line.substr(60, 8)) >> record.velocity.z;
	}

	return record;
}

void MoleculeBuilder::loadAtomPositions(const std::string& gro_path) {	// could be const if not for the temporary atomname for rendering

	const SimpleParsedFile parsedfile = Filehandler::parseGroFile(gro_path, false);

	// First check box has correct size
	if (parsedfile.rows.back().section != "box_size") {
		throw std::runtime_error("Parsed gro file did not contain a box-size at the expected line (2nd last line)");
	}

	for (auto& length_str : parsedfile.rows.back().words) {
		if (std::stof(length_str) != BOX_LEN_NM) {
			throw std::runtime_error(".gro file box size does not match the compiler defined box size");
		}
	}


	// Find atom count to resize particle_info
	if (parsedfile.rows[0].section != "n_atoms") {
		throw std::runtime_error("Expected first line of parsed gro file to contain atom count");
	}
	const int n_atoms_total = std::stoi(parsedfile.rows[0].words[0]);

	for (const auto& row : parsedfile.rows) {
		if (row.section != "atoms") { continue; }

		GroRecord record = parseGroLine(row.words[0]);

		if (record.atom_name[0] == 'H' && IGNORE_HYDROGEN) {
			continue; 
		}

		const bool entry_is_solvent = record.residue_name == "WATER" || record.residue_name == "SOL" || record.residue_name == "HOH";
		if (entry_is_solvent) {
			if (record.atom_name[0] != 'O') {
				continue;	// Treat all three solvents atoms as a single point, so no hydrogens here
			}

			if (nonsolvent_positions.size() != particleinfotable.size()) {
				throw std::runtime_error(std::format("Trying to add solvents, but nonsolvents added ({}) does not equal the expect amount of .lff file ({})", 
					nonsolvent_positions.size(), particleinfotable.size()).c_str());
			}
			solvent_positions.push_back(record.position);
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
			if (record.atom_name != particleinfotable[assumed_global_id].atomname) {
				throw std::runtime_error("atom_name of .gro file does not match that of .lff file");
			}

			const bool is_new_res = residues.empty() || residues.back().gro_id != record.residue_number;
			if (is_new_res) {
				residues.push_back(Residue{ record.residue_number, static_cast<int>(residues.size()), record.residue_name, particleinfotable[assumed_global_id].chain_id });
			}
			
			residues.back().atoms_globalid.emplace_back(assumed_global_id);
			nonsolvent_positions.emplace_back(record.position);
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

	const SimpleParsedFile bonded_parsed = Filehandler::parseLffFile(bonded_path, false);

	Topology topology{};

	for (auto& row : bonded_parsed.rows) {
		if (row.section == "singlebonds") {
			assert(row.words.size() == 6);

			std::array<uint32_t, 2> global_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 2; i++) {
				global_ids[i] = std::stoi(row.words[i]);
			}
			const float b0 = std::stof(row.words[4]) * NANO_TO_LIMA;							// convert [nm] to [lm]
			// Units of kb is [J/mol/nm^2]. One nm is for the error in forcecalc, and one distance in integration
			//const float kb = stof(row.words[5]) / (NANO_TO_PICO * NANO_TO_LIMA);		// convert [J/(mol * nm^2)] to [J/(mol *  * lm)
			const float kb = std::stof(row.words[5]) / (NANO_TO_LIMA * NANO_TO_LIMA);

			topology.singlebonds.emplace_back(SingleBond{ global_ids, b0, kb });
		}
		else if (row.section == "anglebonds") {
			assert(row.words.size() == 8);

			std::array<uint32_t, 3> global_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 3; i++) {
				global_ids[i] = std::stoi(row.words[i]);
			}
			const float theta0 = std::stof(row.words[6]);
			const float ktheta = std::stof(row.words[7]);
			topology.anglebonds.emplace_back(AngleBond{ global_ids, theta0, ktheta });
		}
		else if (row.section == "dihedralbonds") {
			assert(row.words.size() == 11);

			std::array<uint32_t, 4> global_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 4; i++) {
				global_ids[i] = std::stoi(row.words[i]);
			}
			const float phi0 = std::stof(row.words[8]);
			const float kphi = std::stof(row.words[9]);
			const float multiplicity = std::stof(row.words[10]);
			topology.dihedralbonds.emplace_back(DihedralBond{ global_ids, phi0, kphi, multiplicity });
		}
		else if (row.section == "improperdihedralbonds") {
			assert(row.words.size() == 10);

			std::array<uint32_t, 4> global_ids; //{ stoi(row.words[0]), stoi(row.words[1]) };
			for (int i = 0; i < 4; i++) {
				global_ids[i] = std::stoi(row.words[i]);
			}
			const float psi0 = std::stof(row.words[8]);
			const float kpsi = std::stof(row.words[9]);
			topology.improperdihedralbonds.emplace_back(ImproperDihedralBond{ global_ids, psi0, kpsi });
		}
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

			atominfotable.emplace_back(ParticleInfo{ global_id, gro_id, chain_id, atomtype_id, atomname, residue_groid });
		}
	}
	return atominfotable;
}


bool areBonded(const Residue& left, const Residue& right, const ParticleInfoTable& particleinfolut, const std::vector<SingleBond>& singlebonds) {
	assert(left.global_id != right.global_id);

	for (auto& atomleft_gid : left.atoms_globalid) {
		for (int bond_index : particleinfolut[atomleft_gid].singlebonds_indices) {
			const auto& atom_globalids = singlebonds[bond_index].atom_indexes;

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

void MoleculeBuilder::createCompounds(const std::vector<SingleBond>& singlebonds) {
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
				forcefield->atomTypeToIndex(particleinfo.atomname[0]),
				atom_gid
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

#include <unordered_set>
void MoleculeBuilder::createBridges(const std::vector<SingleBond>& singlebonds) 
{
	// First find which ompounds are bonded to other compounds. I assume singlebonds should find all these, 
	// since angles and dihedrals make no sense if the singlebond is not present, as the distances between the atoms can then be extremely large
	std::vector<std::unordered_set<int>> compound_to_bondedcompound_table(compounds.size());
	for (const auto& bond : singlebonds) {
		const int particle_gid_left = bond.atom_indexes[0];
		const int particle_gid_right = bond.atom_indexes[1];

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
std::array<int, 2> getTheTwoDifferentIds(const uint32_t* particle_ids, const ParticleInfoTable& particleinfolut) {
	std::array<int, 2> out = { particleinfolut[particle_ids[0]].compound_index, -1 };

	for (int i = 1; i < n_ids; i++) {
		int compoundid_of_particle = particleinfolut[particle_ids[i]].compound_index;
		if (compoundid_of_particle != out[0]) {
			out[1] = compoundid_of_particle;
			break;
		}
	}

	if (out[1] == -1) { 
		throw std::runtime_error("Failed to find the second compound of bridge"); 
	}

	if (out[0] > out[1]) {
		std::swap(out[0], out[1]);
	}

	return out;
}

template <int n_ids>
BridgeFactory& getBridge(std::vector<BridgeFactory>& bridges, const uint32_t* ids, const ParticleInfoTable& particleinfolut) {
	auto compound_ids = getTheTwoDifferentIds<n_ids>(ids, particleinfolut);

	for (BridgeFactory& bridge : bridges) {
		std::vector<int> compoundids_of_bridge;
		for (int i = 0; i < bridge.n_compounds; i++) {
			compoundids_of_bridge.emplace_back(bridge.compound_ids[i]);
		}

		// Iterate over both the compound_ids of the particles, and the compound_ids in the bridge, 
		// and see if both compound_ids are present in the bridge
		bool bothValuesPresent = true;
		for (auto& compound_id : compound_ids) {
			if (std::find(compoundids_of_bridge.begin(), compoundids_of_bridge.end(), compound_id) == compoundids_of_bridge.end()) {
				bothValuesPresent = false;
				break;
			}
		}

		if (bothValuesPresent) {
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

template <typename BondType>
void MoleculeBuilder::distributeBondsToCompoundsAndBridges(const std::vector<BondType>& bonds) {

	constexpr int atoms_in_bond = BondType::n_atoms;

	for (auto& bond : bonds) {

		if (spansTwoCompounds<atoms_in_bond>(bond.atom_indexes, particleinfotable)) {
			BridgeFactory& bridge = getBridge<atoms_in_bond>(compound_bridges, bond.atom_indexes, particleinfotable);
			bridge.addBond(particleinfotable, bond);
		}
		else {
			const int compound_id = particleinfotable[bond.atom_indexes[0]].compound_index;	// Pick compound using first particle in bond
			CompoundFactory& compound = compounds.at(compound_id);
			compound.addBond(particleinfotable, bond);
		}

		distributeLJIgnores<atoms_in_bond>(bp_lut_manager.get(), particleinfotable, bond.atom_indexes);
	}
}

void MoleculeBuilder::distributeBondsToCompoundsAndBridges(const Topology& topology) {
	bp_lut_manager = std::make_unique<BondedParticlesLUTManager>();

	// First check that we dont have any unrealistic bonds, and ward immediately.
	for (const auto& bond : topology.singlebonds) {
		int gid1 = bond.atom_indexes[0];
		int gid2 = bond.atom_indexes[1];

		const Float3 pos1 = nonsolvent_positions[gid1];
		const Float3 pos2 = nonsolvent_positions[gid2];
		const float hyper_dist = EngineUtils::calcHyperDistNM(&pos1, &pos2);
		if (hyper_dist > bond.b0 * LIMA_TO_NANO * 2.f) {
			throw std::runtime_error(std::format("Loading singlebond with illegally large dist ({}). b0: {}", hyper_dist, bond.b0 * LIMA_TO_NANO).c_str());
		}
	}

	distributeBondsToCompoundsAndBridges(topology.singlebonds);
	distributeBondsToCompoundsAndBridges(topology.anglebonds);
	distributeBondsToCompoundsAndBridges(topology.dihedralbonds);
	distributeBondsToCompoundsAndBridges(topology.improperdihedralbonds);

	m_logger->print(std::format("Added {} singlebonds to molecule\n", topology.singlebonds.size()));
	m_logger->print(std::format("Added {} anglebonds to molecule\n", topology.anglebonds.size()));
	m_logger->print(std::format("Added {} dihedralbonds to molecule\n", topology.dihedralbonds.size()));
	m_logger->print(std::format("Added {} improper dihedralbonds to molecule\n", topology.improperdihedralbonds.size()));
}






Float3 calcCOM(const Float3* positions, int n_elems) {
	Float3 com{};
	for (int i = 0; i < n_elems; i++) {
		Float3 pos = positions[i];
		EngineUtils::applyHyperposNM(&positions[i], &pos);	// Hyperpos around particle 0, since we dont know key position yet
		com += pos;
	}
	return com / static_cast<float>(n_elems);
}

int indexOfParticleClosestToCom(const Float3* positions, int n_elems, const Float3& com) {
	int closest_particle_index = 0;
	float closest_particle_distance = std::numeric_limits<float>::infinity();
	for (int i = 0; i < n_elems; i++) {
		const float dist = EngineUtils::calcHyperDistNM(&positions[i], &com);
		if (dist < closest_particle_distance) {
			closest_particle_distance = dist;
			closest_particle_index = i;
		}
	}
	return closest_particle_index;
}

int indexOfParticleFurthestFromCom(Float3* positions, int n_elems, const Float3& com) {
	int furthest_particle_index = 0;
	float furthest_particle_distance = 0.0f;
	for (int i = 0; i < n_elems; i++) {
		const float dist = EngineUtils::calcHyperDistNM(&positions[i], &com);
		if (dist > furthest_particle_distance) {
			furthest_particle_distance = dist;
			furthest_particle_index = i;
		}
	}
	return furthest_particle_index;
}
void MoleculeBuilder::calcCompoundMetaInfo() {
	for (auto& compound : compounds) {
		const Float3 com = calcCOM(compound.positions, compound.n_particles);	
		compound.key_particle_index = indexOfParticleClosestToCom(compound.positions, compound.n_particles, com);


		const Float3& key_pos = compound.positions[compound.key_particle_index];
		const int indexOfFurthestParticle = indexOfParticleFurthestFromCom(compound.positions, compound.n_particles, key_pos);
		const Float3& furthestParticle = compound.positions[indexOfFurthestParticle];
		const float radius = EngineUtils::calcHyperDistNM(&furthestParticle, &key_pos);

		compound.radius = radius * 1.2f;	// Add 20% leeway

		if (compound.radius > MAX_COMPOUND_RADIUS) { 
			throw std::runtime_error(std::format("Compound radius {} spans more than {} nm!", compound.radius, MAX_COMPOUND_RADIUS).c_str()); }
	}
}


// --------------------------------------------------------------- Factory Functions --------------------------------------------------------------- //

void CompoundFactory::addParticle(const Float3& position, int atomtype_id, int atomtype_color_id, int global_id) {
	if (n_particles >= MAX_COMPOUND_PARTICLES) {
		throw std::runtime_error("Failed to add particle to compound");
	}
	
	// Hyperposition each compound relative to particle 0, so we can find key_particle and radius later
	Float3 hyperpos = position;
	if (n_particles > 0) {
		EngineUtils::applyHyperposNM(&positions[0], &hyperpos);
	}

	// Variables only present in factory
	positions[n_particles] = hyperpos;
	global_ids[n_particles] = global_id;

	// Variables present in Compound
	atom_types[n_particles] = atomtype_id;
	atom_color_types[n_particles] = atomtype_color_id;	// wtf is this

	n_particles++;
}

void CompoundFactory::addBond(const ParticleInfoTable& particle_info, const SingleBond& bondtype) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_COMPOUND) { 
		throw std::runtime_error("Failed to add singlebond to compound"); }
	singlebonds[n_singlebonds++] = SingleBond(
		{
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[0]].local_id_compound),
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[1]].local_id_compound),
		},
		bondtype.b0,
		bondtype.kb
	);
}

void CompoundFactory::addBond(const ParticleInfoTable& particle_info, const AngleBond& bondtype) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_COMPOUND) { throw std::runtime_error("Failed to add anglebond to compound"); }
	anglebonds[n_anglebonds++] = AngleBond(
		{
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[0]].local_id_compound),
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[1]].local_id_compound),
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[2]].local_id_compound),
		},		
		bondtype.theta_0,
		bondtype.k_theta
	);
}

void CompoundFactory::addBond(const ParticleInfoTable& particle_info, const DihedralBond& bondtype) {
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_COMPOUND) { 
		throw std::runtime_error("Failed to add dihedralbond to compound"); }
	dihedrals[n_dihedrals++] = DihedralBond(
		{
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[0]].local_id_compound),
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[1]].local_id_compound),
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[2]].local_id_compound),
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[3]].local_id_compound),
		},
		bondtype.phi_0,
		bondtype.k_phi,
		bondtype.n
	);
}

void CompoundFactory::addBond(const ParticleInfoTable& particle_info, const ImproperDihedralBond& bondtype) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND) { throw std::runtime_error("Failed to add improperdihedralbond to compound"); }
	impropers[n_improperdihedrals++] = ImproperDihedralBond(
		std::array<uint32_t, 4>{
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[0]].local_id_compound),
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[1]].local_id_compound),
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[2]].local_id_compound),
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[3]].local_id_compound),
		},
		bondtype.psi_0,
		bondtype.k_psi
	);
}



void BridgeFactory::addBond(ParticleInfoTable& particle_info, const SingleBond& bondtype) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add singlebond to bridge"); }
	singlebonds[n_singlebonds++] = SingleBond{
		{
			getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[0]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[1]]),
		},
		bondtype.b0,
		bondtype.kb
	};
}

void BridgeFactory::addBond(ParticleInfoTable& particle_info, const AngleBond& bondtype) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add anglebond to bridge"); }
	anglebonds[n_anglebonds++] = AngleBond{
		{
			getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[0]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[1]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[2]]),
		},		
		bondtype.theta_0,
		bondtype.k_theta
	};
}

void BridgeFactory::addBond(ParticleInfoTable& particle_info, const DihedralBond& bondtype) {

	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add dihedralbond to bridge"); }
	dihedrals[n_dihedrals++] = DihedralBond{
		{
		getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[0]]),
		getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[1]]),
		getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[2]]),
		getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[3]]),
		},		
		bondtype.phi_0,
		bondtype.k_phi,
		bondtype.n
	};


}

void BridgeFactory::addBond(ParticleInfoTable& particle_info, const ImproperDihedralBond& bondtype) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_BRIDGE) { throw std::runtime_error("Failed to add dihedralbond to bridge"); }
	impropers[n_improperdihedrals++] = ImproperDihedralBond{
		std::array<uint32_t,4>{
			getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[0]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[1]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[2]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[3]]),
		},
		bondtype.psi_0,
		bondtype.k_psi,
	};
}


uint32_t BridgeFactory::getBridgelocalIdOfParticle(ParticleInfo& particle_info) {

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
	if (particle_info.local_id_bridge == -1) {
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









CompoundCollection LIMA_MOLECULEBUILD::buildMolecules(
	Forcefield* ff,
	const std::string& work_dir,
	VerbosityLevel vl,
	const string& gro_path,
	const string& molecule_dir,
	std::unique_ptr<LimaLogger> logger,
	bool ignore_hydrogens
) 
{
	MoleculeBuilder mb{ ff, std::move(logger), work_dir, vl};
	return mb.buildMolecules(gro_path, molecule_dir, ignore_hydrogens);
}
