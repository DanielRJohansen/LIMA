#include "CompoundBuilder.h"

#include "Printer.h"
#include "EngineUtils.cuh"


#include <algorithm>
#include <format>
using namespace LIMA_Print;

class MoleculeBuilder {
public:
	MoleculeBuilder(Forcefield* ff, const std::string& work_dir = "", VerbosityLevel vl = SILENT);

	// Upon creation of the CompoundCollection, the MoleculeBuilder object is no longer valid,
	// as it move's it's data instead of copying!
	CompoundCollection buildMolecules(const string& gro_path, const string& topol_path, bool ignore_hydrogens = true);

private:
	// Members, seen in the order they are filled in
	std::vector<Residue> residues;
	std::vector<Float3> solvent_positions;
	int64_t n_particles_in_residues = 0;

	std::vector<std::array<int, 2>> singlebonds;	// gro id's of particles
	std::vector<std::array<int, 3>> anglebonds;		// gro id's of particles
	std::vector<std::array<int, 4>> dihedralbonds;	// gro id's of particles
	std::vector<ParticleInfo> particle_info;		// Uses gro 1-indexed 

	std::vector<CompoundFactory> compounds;
	std::vector<BridgeFactory> compound_bridges;

	BondedParticlesLUTManager* bp_lut_manager = nullptr;

	LimaLogger logger;
	const VerbosityLevel verbosity_level;
	const Forcefield* forcefield;


	// ------------------------------------ HELPER FUNCTIONS ------------------------------------ //

	// Only works for fixed positions gro files!! Meaning max, 99.999 particles
	void loadResiduesAndSolvents(const std::string gro_path);

	// Loads all bonds, and sets references to these for each particle - HEAVY!
	void loadTopology(const std::string& topol_path);

	/// <summary>
	/// Goes through all the residues, and fills the bondedresidue_ids vector. 
	/// The particle_bonds_lut is used to determine whether particles of two different 
	/// share a bond.
	/// </summary>
	void matchBondedResidues();

	// TODO: Put correct references in the bridges too!
	void createCompoundsAndBridges();

	void distributeBondsToCompoundsAndBridges();

	// Find index of centermost particle, find "radius" of compound
	void calcCompoundMetaInfo();

	template <int n_ids>
	bool spansTwoCompounds(std::array<int, n_ids> bond_ids);
};




bool isLineCommented(const std::string& input, const std::vector<char>& comment_markers = { ';', '#' }) {
	// Iterate through each character in the string
	for (char c : input) {
		// If the character is not a space, check if it's a comment marker
		if (!std::isspace(c)) {
			// Iterate through each comment marker in the vector
			for (char comment_marker : comment_markers) {
				if (c == comment_marker) {
					return true; // The first non-space character is a comment marker
				}
			}
			return false; // The first non-space character is not a comment marker
		}
	}
	return false; // The string is empty or contains only spaces
}

std::vector<std::string> readFile(const std::string& file_path) {
	// Check if file exists
	std::ifstream file(file_path);
	if (!file.good()) {
		std::cerr << "Error: File " << file_path << " does not exist." << std::endl;
		return {};
	}

	// Read file lines
	std::vector<std::string> lines;
	std::string line;
	while (std::getline(file, line)) {
		if (isLineCommented(line)) {continue;} // Skip comments
		lines.push_back(line);
	}

	// Close file
	file.close();

	return lines;
}






MoleculeBuilder::MoleculeBuilder(Forcefield* ff, const std::string& work_dir, VerbosityLevel vl) :
	logger(LimaLogger::LogMode::compact, "moleculebuilder", work_dir),
	verbosity_level(vl), 
	forcefield(ff) 
{}

CompoundCollection MoleculeBuilder::buildMolecules(const string& gro_path, const string& topol_path, bool ignore_hydrogens) {

	printH2("Building molecules", true, false);

	loadResiduesAndSolvents(gro_path);

	loadTopology(topol_path);

	matchBondedResidues();

	createCompoundsAndBridges();

	distributeBondsToCompoundsAndBridges();

	calcCompoundMetaInfo();

	CompoundBridgeBundleCompact bridges_compact(
		std::vector<CompoundBridge>(compound_bridges.begin(), compound_bridges.end())
	);

	// Danger: After this point, the contents of this object is no longer valid.
	return CompoundCollection{ 
		std::move(compounds), 
		n_particles_in_residues,
		bp_lut_manager, 
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
	std::istringstream(line.substr(15, 5)) >> record.atom_number;

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


void MoleculeBuilder::loadResiduesAndSolvents(const std::string gro_path) {

	const bool ignore_hydrogens = true;

	int64_t line_cnt = -1;

	std::vector<string> lines = readFile(gro_path);

	// First setup the lut. +1 because the extern particle indices are 1-indexed
	// We need to set this up now, because we will index it with the gro ids, of which there is almost 1 per line in the .gro file
	particle_info.resize(lines.size() + 1);	// 

	for (const auto& line : lines) {
		line_cnt++;
		
		if (line_cnt == 0) { continue; }	// Title is irrelevant
		if (line_cnt == 1) {				// 2nd line is particle count
			const int64_t atoms_total = stoll(line);	// Only a ballpark number!

			// Reserve enough vectorspace to handle all particles being either in residues or solvents
			residues.reserve(atoms_total / Residue::max_size_soft);	// Misses the mark when some residues are smaller..
			solvent_positions.reserve(atoms_total);
			continue;
		}

		if (line.size() < 44) { continue; }	// Should only happen at the end. TODO: Extract box size here. Make sure it is the same as the constant var!

		const GroRecord record = parseGroLine(line);

		if (record.atom_name[0] == 'H' && ignore_hydrogens) { 			
			continue; 
		}

		// Is Solvent
		if (record.residue_name == "WATER" || record.residue_name == "SOL") {
			solvent_positions.push_back(record.position);
		}
		else {	// Is in residue
			// Belongs to current residue?
			if (!residues.empty() && residues.back().gro_id == record.residue_number) {				
				// Do nothing

				// SAFETY CHECK: This is likely because we have molecules consisting of only 1 residue, so we can figure out which belongs to the same molecule!
				if (residues.back().name != record.residue_name) { throw std::exception("Something went seriously wrong when handling the .gro file!"); }	
			}
			else {	// Create new residue
				int unique_res_id = residues.size();
				residues.push_back(Residue{ record.residue_number, unique_res_id, record.residue_name });
			}



			residues.back().atoms.push_back(AtomEntry{ record.atom_number, record.position, record.atom_name });
			n_particles_in_residues++;
		}

		// In some cases the atoms aren't numbered.- Then we need to expand the lut
		// IF this happens once, it may happen many times more, so we expand to twice the necessary size.
		if (record.atom_number >= particle_info.size()) {
			particle_info.resize(record.atom_number * 2);
		}

		particle_info[record.atom_number].isUsed = true;
		particle_info[record.atom_number].gro_id = record.atom_number;
	}
}








bool setMode(const std::vector<string>& row, TopologyMode& current_mode, int& dihedral_sections_count)
{
	// If empty row, then we are finished with a section	- THIS might be sorted away in the file-read, so...?
	if (row.size() == 0) { return INACTIVE; }

	// Since the above is sorted away in the file read, we probably go from one section directly to the the next section header
	if (row.size() == 1) {
		if (row[0] == "bonds") {
			current_mode = BOND;
		}
		else if (row[0] == "angles") {
			current_mode = ANGLE;
		}
		else if (row[0] == "dihedrals") {
			dihedral_sections_count++;
			if (dihedral_sections_count == 1) { current_mode = DIHEDRAL; }
			else { current_mode = INACTIVE; }
		}
		else {
			current_mode = INACTIVE;
		}
		return true;
	}

	// By default we stay in the current mode.
	return false;
}

template <int n_ids> 
bool allParticlesExist(std::array<int, n_ids> gro_ids, const std::vector<ParticleInfo>& particle_info) {
	for (auto id : gro_ids) {
		if (id >= particle_info.size() || !particle_info[id].isUsed) {
			return false;
		}
	}
	return true;
}

void MoleculeBuilder::loadTopology(const std::string& topol_path) {
	logger.print(std::format("Reading topology from file {}\n", topol_path));

	int line_cnt = 0;
	TopologyMode mode = INACTIVE;
	int dihedral_sections_count = 0;	// bad fix-...

	auto lines = Filehandler::readFile(topol_path, { ';', '#' }, { "[", "]" });

	for (const auto& line : lines) {

		const bool new_section = setMode(line, mode, dihedral_sections_count);
		if (new_section) { continue; }

		switch (mode)
		{
		case INACTIVE:
			break;
		case ATOMS:
			break;
		case BOND:
		{
			std::array extern_indices = { stoi(line[0]), stoi(line[1]) };
			if (!allParticlesExist(extern_indices, particle_info)) { break; }

			singlebonds.push_back(extern_indices);							// First create the bond
			for (int index : extern_indices) {
				particle_info[index].singlebonds_indices.push_back(singlebonds.size()-1);		// Then make all involved particles reference this bond
			}


			break;
		}			
		case ANGLE:
		{
			std::array extern_indices = { stoi(line[0]), stoi(line[1]), stoi(line[2]) };
			if (!allParticlesExist(extern_indices, particle_info)) { break; }

			anglebonds.push_back( extern_indices );
			for (int index : extern_indices) {
				particle_info[index].anglebonds_indices.push_back(anglebonds.size()-1);
			}
			break;
		}			
		case DIHEDRAL:
		{
			std::array extern_indices = { stoi(line[0]), stoi(line[1]), stoi(line[2]), stoi(line[3]) };
			if (!allParticlesExist(extern_indices, particle_info)) { break; }

			dihedralbonds.push_back( extern_indices );
			for (int index : extern_indices) {
				particle_info[index].dihedralbonds_indices.push_back(dihedralbonds.size()-1);
			}
			break;
		}			
		default:
			break;
		}
	}
}

bool areBonded(const Residue& left, const Residue& right, std::vector<ParticleInfo>& particle_bonds_lut, const std::vector<std::array<int, 2>>& all_singlebonds) {
	for (auto& atom_left : left.atoms) {
		for (int bond_index : particle_bonds_lut[atom_left.id].singlebonds_indices) {
			const auto& bond = all_singlebonds[bond_index];

			for (auto& atom_right : right.atoms) {
				if ((bond[0] == atom_left.id && bond[1] == atom_right.id)
					|| (bond[1] == atom_left.id && bond[0] == atom_right.id)) 
				{
					return true;
				}
			}
		}
	}
	return false;
}

void MoleculeBuilder::matchBondedResidues() {
	for (int i = 0; i < static_cast<int>(residues.size()) - 1; i++) {
		Residue& residue_left = residues[i];
		Residue& residue_right = residues[i + 1];

		if (areBonded(residue_left, residue_right, particle_info, singlebonds)) {
			residue_left.bondedresidue_ids.push_back(i + 1);
			residue_right.bondedresidue_ids.push_back(i);
		}
	}
}





#include <algorithm>
#include <iterator>
void MoleculeBuilder::createCompoundsAndBridges() {
	// Nothing to do if we have no residues
	if (residues.empty()) { return; }
	
	compounds.push_back(CompoundFactory{ 0 });

	int current_residue_id = residues[0].id;

	for (int residue_index = 0; residue_index < residues.size(); residue_index++) {
		const Residue& residue = residues[residue_index];

		const bool new_residue = residue.id != current_residue_id;
		bool bridges_to_previous_compound = false;

		// If necessary start new compound
		if (new_residue) {
			const auto& bonded_ids = residue.bondedresidue_ids;
			const bool is_bonded_with_previous_residue = std::find(bonded_ids.begin(), bonded_ids.end(), current_residue_id) != bonded_ids.end();

			
			// If the new residue is a different molecule, start a new compound
			if (!is_bonded_with_previous_residue) {
				compounds.push_back(compounds.size());
			}
			// If we are bonded, but no more room, start new compound and make bridge
			else if (!compounds.back().hasRoomForRes(residue.atoms.size())) {
				compounds.push_back(compounds.size());

				int id_left = compounds.size() - 2;
				int id_right = compounds.size() - 1;

				compound_bridges.push_back(BridgeFactory{ { id_left, id_right} });

				bridges_to_previous_compound = true;
			}

			current_residue_id = residue.id;
		}

		// Add all atoms of residue to current compound
		for (auto& atom : residue.atoms) {
			// First add the new information to the particle
			particle_info[atom.id].compound_index = compounds.back().id;
			particle_info[atom.id].local_id_compound = compounds.back().n_particles;

			// Then add the particle to the compound
			compounds.back().addParticle(
				atom.position,
				forcefield->getAtomtypeID(atom.id),
				forcefield->atomTypeToIndex(atom.name[0]),
				particle_info[atom.id].gro_id
				);
		}
	}
}







template <int n_ids>
bool MoleculeBuilder::spansTwoCompounds(std::array<int, n_ids> bond_ids) {
	const int compound_left = particle_info[bond_ids[0]].compound_index;

	for (int i = 1; i < n_ids; i++) {
		if (compound_left != particle_info[bond_ids[i]].compound_index) {
			return true;
		}
	}
	return false;
}

template <int n_ids>
std::array<int, 2> getTheTwoDifferentIds(std::array<int, n_ids> particle_ids, const std::vector<ParticleInfo>& particle_info) {
	std::array<int, 2> out = { particle_info[particle_ids[0]].compound_index, -1 };

	for (int i = 1; i < n_ids; i++) {
		int compoundid_of_particle = particle_info[particle_ids[i]].compound_index;
		if (compoundid_of_particle != out[0]) {
			out[1] = compoundid_of_particle;
			break;
		}
	}

	if (out[1] == -1) { throw std::exception("Failed to find the second compound of bridge"); }

	if (out[0] > out[1]) {
		std::swap(out[0], out[1]);
	}

	return out;
}

template <int n_ids>
BridgeFactory& getBridge(std::vector<BridgeFactory>& bridges, const std::array<int, n_ids>& ids, const std::vector<ParticleInfo>& particle_info) {
	auto compound_ids = getTheTwoDifferentIds(ids, particle_info);

	for (BridgeFactory& bridge : bridges) {
		if (bridge.compound_id_left == compound_ids[0] && bridge.compound_id_right == compound_ids[1]) {	// I assume they are ordered by size here
			return bridge;
		}
	}

	throw std::exception("Failed to find the bridge");
}


template <int n_ids>
void distributeLJIgnores(BondedParticlesLUTManager* bplut_man, const std::vector<ParticleInfo>& particle_info, const std::array<int, n_ids>& gro_ids) {

	for (auto id_self : gro_ids) {
		for (auto id_other : gro_ids) {
			if (id_self == id_other) { continue; }

			BondedParticlesLUT* lut = bplut_man->get(particle_info[id_self].compound_index, particle_info[id_other].compound_index);
			lut->set(particle_info[id_self].local_id_compound, particle_info[id_other].local_id_compound, true);
		}
	}
}

void MoleculeBuilder::distributeBondsToCompoundsAndBridges() {

	bp_lut_manager = new BondedParticlesLUTManager{0};

	// Distribute single bonds
	for (auto& bond_ids : singlebonds) {	// gro 1-indexed ids
		const SingleBond& bondtype = *forcefield->getBondType({ bond_ids[0], bond_ids[1]});

		// First add the bond information to either bridge or compound
		if (spansTwoCompounds(bond_ids)) {
			BridgeFactory& bridge = getBridge(compound_bridges, bond_ids, particle_info);
			bridge.addSingleBond({ &particle_info[bond_ids[0]], &particle_info[bond_ids[1]]}, bondtype);
		}
		else {
			const int compound_id = particle_info[bond_ids[0]].compound_index;
			CompoundFactory& compound = compounds.at(compound_id);
			compound.addSingleBond({ particle_info[bond_ids[0]], particle_info[bond_ids[1]] }, bondtype);
		}

		distributeLJIgnores(bp_lut_manager, particle_info, bond_ids);
	}

	// Distribute angle bonds
	for (auto& bond_ids : anglebonds) {    // gro 1-indexed ids
		const AngleBond& bondtype = *forcefield->getAngleType({ bond_ids[0], bond_ids[1], bond_ids[2] });

		// First add the bond information to either bridge or compound
		if (spansTwoCompounds(bond_ids)) {
			BridgeFactory& bridge = getBridge(compound_bridges, bond_ids, particle_info);
			bridge.addAngleBond({ &particle_info[bond_ids[0]], &particle_info[bond_ids[1]], &particle_info[bond_ids[2]] }, bondtype);
		}
		else {
			const int compound_id = particle_info[bond_ids[0]].compound_index;
			CompoundFactory& compound = compounds.at(compound_id);
			compound.addAngleBond({ particle_info[bond_ids[0]],	particle_info[bond_ids[1]],	particle_info[bond_ids[2]] }, bondtype);
		}

		distributeLJIgnores(bp_lut_manager, particle_info, bond_ids);
	}

	// Distribute dihedral bonds
	for (auto& bond_ids : dihedralbonds) {    // gro 1-indexed ids
		const DihedralBond& bondtype = *forcefield->getDihedralType({ bond_ids[0], bond_ids[1], bond_ids[2], bond_ids[3] });

		// First add the bond information to either bridge or compound
		if (spansTwoCompounds(bond_ids)) {
			BridgeFactory& bridge = getBridge(compound_bridges, bond_ids, particle_info);
			bridge.addDihedralBond({ &particle_info[bond_ids[0]], &particle_info[bond_ids[1]], &particle_info[bond_ids[2]], &particle_info[bond_ids[3]] }, bondtype);
		}
		else {
			const int compound_id = particle_info[bond_ids[0]].compound_index;
			CompoundFactory& compound = compounds.at(compound_id);
			compound.addDihedralBond({ particle_info[bond_ids[0]], particle_info[bond_ids[1]], particle_info[bond_ids[2]], particle_info[bond_ids[3]] }, bondtype);
		}

		distributeLJIgnores(bp_lut_manager, particle_info, bond_ids);
	}
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
		//float particle_distance = (positions[i] - com).len();
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
		//float particle_distance = (positions[i] - com).len();
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

		if (compound.radius > 1.5f) { 
			throw std::exception("Compound radius spans more than 1.8 nm!"); }
	}
}


// --------------------------------------------------------------- Factory Functions --------------------------------------------------------------- //

void CompoundFactory::addParticle(const Float3& position, int atomtype_id, int atomtype_color_id, int gro_id) {
	if (n_particles >= MAX_COMPOUND_PARTICLES) {
		throw std::exception("Failed to add particle to compound");
	}
	
	// Hyperposition each compound relative to particle 0, so we can find key_particle and radius later
	Float3 hyperpos = position;
	if (n_particles > 0) {
		EngineUtils::applyHyperposNM(&positions[0], &hyperpos);
	}

	// Variables only present in factory
	positions[n_particles] = hyperpos;
	gro_ids[n_particles] = gro_id;

	// Variables present in Compound
	atom_types[n_particles] = atomtype_id;
	atom_color_types[n_particles] = atomtype_color_id;	// wtf is this

	n_particles++;
}

void CompoundFactory::addSingleBond(const std::array<ParticleInfo, 2>& particle_info, const SingleBond& bondtype) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_COMPOUND) { throw std::exception("Failed to add singlebond to compound"); }
	singlebonds[n_singlebonds++] = SingleBond(
		particle_info[0].local_id_compound,
		particle_info[1].local_id_compound,
		bondtype.b0,
		bondtype.kb
	);
}

void CompoundFactory::addAngleBond(const std::array<ParticleInfo, 3>& particle_info, const AngleBond& bondtype) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_COMPOUND) { throw std::exception("Failed to add anglebond to compound"); }
	anglebonds[n_anglebonds++] = AngleBond(
		particle_info[0].local_id_compound,
		particle_info[1].local_id_compound,
		particle_info[2].local_id_compound,
		bondtype.theta_0,
		bondtype.k_theta
	);
}

void CompoundFactory::addDihedralBond(const std::array<ParticleInfo, 4>& particle_info, const DihedralBond& bondtype) {
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_COMPOUND) { throw std::exception("Failed to add dihedralbond to compound"); }
	dihedrals[n_dihedrals++] = DihedralBond(
		particle_info[0].local_id_compound,
		particle_info[1].local_id_compound,
		particle_info[2].local_id_compound,
		particle_info[3].local_id_compound,
		bondtype.phi_0,
		bondtype.k_phi,
		bondtype.n
	);
}


void BridgeFactory::addSingleBond(std::array<ParticleInfo*, 2> particle_info, const SingleBond& bondtype) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_BRIDGE) { throw std::exception("Failed to add singlebond to bridge"); }
	singlebonds[n_singlebonds++] = SingleBond{
		getBridgelocalIdOfParticle(*particle_info[0]),
		getBridgelocalIdOfParticle(*particle_info[1]),
		bondtype.b0,
		bondtype.kb
	};
}

void BridgeFactory::addAngleBond(std::array<ParticleInfo*, 3> particle_info, const AngleBond& bondtype) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_BRIDGE) { throw std::exception("Failed to add anglebond to bridge"); }
	anglebonds[n_anglebonds++] = AngleBond{
		getBridgelocalIdOfParticle(*particle_info[0]),
		getBridgelocalIdOfParticle(*particle_info[1]),
		getBridgelocalIdOfParticle(*particle_info[2]),
		bondtype.theta_0,
		bondtype.k_theta
	};
}

void BridgeFactory::addDihedralBond(std::array<ParticleInfo*, 4> particle_info, const DihedralBond& bondtype) {
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_BRIDGE) { throw std::exception("Failed to add dihedralbond to bridge"); }
	dihedrals[n_dihedrals++] = DihedralBond{
		getBridgelocalIdOfParticle(*particle_info[0]),
		getBridgelocalIdOfParticle(*particle_info[1]),
		getBridgelocalIdOfParticle(*particle_info[2]),
		getBridgelocalIdOfParticle(*particle_info[3]),
		bondtype.phi_0,
		bondtype.k_phi,
		bondtype.n
	};
}

int BridgeFactory::getBridgelocalIdOfParticle(ParticleInfo& particle_info) {
	if (particle_info.local_id_bridge == -1) {
		if (n_particles == MAX_PARTICLES_IN_BRIDGE) { throw std::exception("Failed to add particle to bridge"); }
		particle_info.local_id_bridge = n_particles;
		particle_refs[n_particles++] = ParticleReference{ 
			particle_info.compound_index, 
			particle_info.local_id_compound,
#ifdef LIMADEBUGMODE
			particle_info.gro_id
#endif
		};
	}
	return particle_info.local_id_bridge;
}









CompoundCollection LIMA_MOLECULEBUILD::buildMolecules(
	Forcefield* ff,
	const std::string& work_dir,
	VerbosityLevel vl,
	const string& gro_path,
	const string& topol_path,
	bool ignore_hydrogens) 
{
	MoleculeBuilder mb{ ff, work_dir, vl };
	return mb.buildMolecules(gro_path, topol_path, ignore_hydrogens);
}