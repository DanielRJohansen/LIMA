#include "LIMA_MD/include/CompoundBuilder.h"
#include "LIMA_BASE/include/Printer.h"
#include "LIMA_ENGINE/include/EngineUtils.cuh"





#include <algorithm>
#include <format>
#include <span>

using namespace LIMA_Print;

class MoleculeBuilder {
public:
	MoleculeBuilder(Forcefield* ff, std::unique_ptr<LimaLogger>, const std::string& work_dir = "", VerbosityLevel vl = SILENT);

	// Upon creation of the CompoundCollection, the MoleculeBuilder object is no longer valid,
	// as it move's it's data instead of copying!
	CompoundCollection buildMolecules(const string& gro_path, const string& topol_path, bool ignore_hydrogens = true);

private:

	// Members, seen in the order they are filled in
	std::vector<Residue> residues;
	std::vector<Float3> solvent_positions;
	int32_t n_particles_in_residues = 0;

	std::vector<ParticleInfo> particle_info;		// Uses gro 1-indexed 

	std::vector<CompoundFactory> compounds;
	std::vector<BridgeFactory> compound_bridges;

	//BondedParticlesLUTManager* bp_lut_manager = nullptr;
	std::unique_ptr<BondedParticlesLUTManager> bp_lut_manager;

	std::unique_ptr<LimaLogger> m_logger;
	const VerbosityLevel verbosity_level;
	const Forcefield* forcefield;


	// ------------------------------------ HELPER FUNCTIONS ------------------------------------ //

	// Only works for fixed positions gro files!! Meaning max, 99.999 particles
	void loadResiduesAndSolvents(const std::string gro_path);

	/// <summary>
	/// Goes through all the residues, and fills the bondedresidue_ids vector. 
	/// The particle_bonds_lut is used to determine whether particles of two different 
	/// share a bond.
	/// </summary>
	void matchBondedResidues(const std::vector<SingleBond>& singlebonds);

	// TODO: Put correct references in the bridges too!
	void createCompoundsAndBridges();

	void distributeBondsToCompoundsAndBridges(const Forcefield::Topology& topology);

	template <typename BondType>
	void distributeBondsToCompoundsAndBridges(const std::vector<BondType>& bonds);

	// Find index of centermost particle, find "radius" of compound
	void calcCompoundMetaInfo();
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






MoleculeBuilder::MoleculeBuilder(Forcefield* ff, std::unique_ptr<LimaLogger> logger, const std::string& work_dir, VerbosityLevel vl) :
	m_logger(std::move(logger)),
	verbosity_level(vl), 
	forcefield(ff) 
{}

CompoundCollection MoleculeBuilder::buildMolecules(const string& gro_path, const string& topol_path, bool ignore_hydrogens) {

	m_logger->startSection("Building Molecules");

	loadResiduesAndSolvents(gro_path);

	const Forcefield::Topology& topology = forcefield->getTopology();

	matchBondedResidues(topology.singlebonds);

	createCompoundsAndBridges();

	distributeBondsToCompoundsAndBridges(topology);

	calcCompoundMetaInfo();


	auto bridges_compact = std::make_unique<CompoundBridgeBundleCompact>(
		std::vector<CompoundBridge>(compound_bridges.begin(), compound_bridges.end())
	);

	m_logger->finishSection("Finished building molecules");

	return CompoundCollection{ 
		std::move(compounds), 
		n_particles_in_residues,
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

	const SimpleParsedFile parsedfile = Filehandler::parseGroFile(gro_path, false);


	// First check box has correct size
	if (parsedfile.rows.back().section != "box_size") {
		throw std::exception("Parsed gro file did not contain a box-size at the expected line (2nd last line)");
	}
	for (auto& length_str : parsedfile.rows.back().words) {
		if (std::stof(length_str) != BOX_LEN_NM) {
			throw std::exception(".gro file box size does not match the compiler defined box size");
		}
	}


	// Find atom count to resize particle_info
	if (parsedfile.rows[0].section != "n_atoms") {
		throw std::exception("Expected first line of parsed gro file to contain atom count");
	}
	const int64_t n_atoms = std::stoll(parsedfile.rows[0].words[0]);
	particle_info.resize(n_atoms);



	for (const auto& row : parsedfile.rows) {
		if (row.section != "atoms") { continue; }

		GroRecord record = parseGroLine(row.words[0]);

		if (record.atom_name[0] == 'H' && ignore_hydrogens) {
			continue; 
		}

		// Is Solvent
		if (record.residue_name == "WATER" || record.residue_name == "SOL" || record.residue_name == "HOH") {
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
				const int unique_res_id = static_cast<int>(residues.size());
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


bool areBonded(const Residue& left, const Residue& right, std::vector<ParticleInfo>& particle_bonds_lut, const std::vector<SingleBond>& singlebonds) {
	for (auto& atom_left : left.atoms) {
		for (int bond_index : particle_bonds_lut[atom_left.id].singlebonds_indices) {
			const auto& atom_groids = singlebonds[bond_index].atom_indexes;

			for (auto& atom_right : right.atoms) {
				if ((atom_groids[0] == atom_left.id && atom_groids[1] == atom_right.id)
					|| (atom_groids[1] == atom_left.id && atom_groids[0] == atom_right.id))
				{
					return true;
				}
			}
		}
	}
	return false;
}

void MoleculeBuilder::matchBondedResidues(const std::vector<SingleBond>& singlebonds) {
	for (int singlebond_index = 0; singlebond_index < singlebonds.size(); singlebond_index++) {
		for (uint32_t gro_id : singlebonds[singlebond_index].atom_indexes) {
			particle_info[gro_id].singlebonds_indices.push_back(singlebond_index);
		}
	}


	for (int i = 0; i < static_cast<int>(residues.size()) - 1; i++) {
		Residue& residue_left = residues[i];
		Residue& residue_right = residues[i + 1];

		if (areBonded(residue_left, residue_right, particle_info, singlebonds)) {
		//if (areBonded(residue_left, residue_right, particle_info, singlebond_references)) {
			residue_left.bondedresidue_ids.push_back(i + 1);
			residue_right.bondedresidue_ids.push_back(i);
		}
	}
}

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
				if (compounds.size() >= MAX_COMPOUNDS) {
					throw std::exception(std::format("Cannot handle more than {} compounds", MAX_COMPOUNDS).c_str());
				}

				compounds.push_back(CompoundFactory{ static_cast<int>(compounds.size()) });
			}
			// If we are bonded, but no more room, start new compound and make bridge
			else if (!compounds.back().hasRoomForRes(residue.atoms.size())) {

				if (compounds.size() >= MAX_COMPOUNDS) {
					throw std::exception(std::format("Cannot handle more than {} compounds", MAX_COMPOUNDS).c_str());
				}
				compounds.push_back(CompoundFactory{ static_cast<int>(compounds.size()) });

				const int id_left =  static_cast<int>(compounds.size()) - 2;
				const int id_right = static_cast<int>(compounds.size()) - 1;

				compound_bridges.push_back(BridgeFactory{(int)compound_bridges.size(), { id_left, id_right} });

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

	m_logger->print(std::format("Created {} compounds\n", compounds.size()));
	m_logger->print(std::format("Created {} compound bridges\n", compound_bridges.size()));
}






template <int n_ids>
bool spansTwoCompounds(const uint32_t* bond_groids, const std::vector<ParticleInfo>& particle_info) {
	const int compound_left = particle_info[bond_groids[0]].compound_index;

	for (int i = 1; i < n_ids; i++) {
		if (compound_left != particle_info[bond_groids[i]].compound_index) {
			return true;
		}
	}
	return false;
}


// Returns two compound id's of a bond in a bridge. The id's are sorted with lowest first
template <int n_ids>
std::array<int, 2> getTheTwoDifferentIds(const uint32_t* particle_ids, const std::vector<ParticleInfo>& particle_info) {
	std::array<int, 2> out = { particle_info[particle_ids[0]].compound_index, -1 };

	for (int i = 1; i < n_ids; i++) {
		int compoundid_of_particle = particle_info[particle_ids[i]].compound_index;
		if (compoundid_of_particle != out[0]) {
			out[1] = compoundid_of_particle;
			break;
		}
	}

	if (out[1] == -1) { 
		throw std::exception("Failed to find the second compound of bridge"); 
	}

	if (out[0] > out[1]) {
		std::swap(out[0], out[1]);
	}

	return out;
}

template <int n_ids>
BridgeFactory& getBridge(std::vector<BridgeFactory>& bridges, const uint32_t* ids, const std::vector<ParticleInfo>& particle_info) {
	auto compound_ids = getTheTwoDifferentIds<n_ids>(ids, particle_info);

	for (BridgeFactory& bridge : bridges) {
		if (bridge.compound_id_left == compound_ids[0] && bridge.compound_id_right == compound_ids[1]) {	// I assume they are ordered by size here
			return bridge;
		}
	}

	throw std::exception(std::format("Failed to find the bridge ({})", n_ids).c_str());
}

template <int n_ids>
void distributeLJIgnores(BondedParticlesLUTManager* bplut_man, const std::vector<ParticleInfo>& particle_info, const uint32_t* gro_ids) {
	for (int i = 0; i < n_ids; i++) {
		auto id_self = gro_ids[i];
		for (int j = 0; j < n_ids; j++) {
			auto id_other = gro_ids[j];

			if (id_self == id_other) { continue; }

			bplut_man->addNewConnectedCompoundIfNotAlreadyConnected(particle_info[id_self].compound_index, particle_info[id_other].compound_index);
			BondedParticlesLUT* lut = bplut_man->get(particle_info[id_self].compound_index, particle_info[id_other].compound_index);
			lut->set(particle_info[id_self].local_id_compound, particle_info[id_other].local_id_compound, true);
		}
	}
}

template <typename BondType>
void MoleculeBuilder::distributeBondsToCompoundsAndBridges(const std::vector<BondType>& bonds) {

	constexpr int atoms_in_bond = BondType::n_atoms;

	for (auto& bond : bonds) {

		if (spansTwoCompounds<atoms_in_bond>(bond.atom_indexes, particle_info)) {
			BridgeFactory* bridge;
			try {
				bridge = &getBridge<atoms_in_bond>(compound_bridges, bond.atom_indexes, particle_info);
				
			}
			catch (const std::exception& ex){
				int a = 0;
				std::cout << ex.what() << "\n";
				// TODO: fix this
				// In this case, a bond is made between two non-adjecent residues, to form a plate or a helix
				// For now, simply ignore these bonds
				continue;
			}
			bridge->addBond(particle_info, bond);

		}
		else {
			const int compound_id = particle_info[bond.atom_indexes[0]].compound_index;	// Pick compound using first particle in bond
			CompoundFactory& compound = compounds.at(compound_id);
			compound.addBond(particle_info, bond);
		}

		distributeLJIgnores<atoms_in_bond>(bp_lut_manager.get(), particle_info, bond.atom_indexes);
	}
}

void MoleculeBuilder::distributeBondsToCompoundsAndBridges(const Forcefield::Topology& topology) {
	bp_lut_manager = std::make_unique<BondedParticlesLUTManager>();

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

		if (compound.radius > MAX_COMPOUND_RADIUS) { 
			throw std::exception(std::format("Compound radius {} spans more than {} nm!", compound.radius, MAX_COMPOUND_RADIUS).c_str()); }
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

template <> void CompoundFactory::addBond(const std::vector<ParticleInfo>& particle_info, const SingleBond& bondtype) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_COMPOUND) { 
		throw std::exception("Failed to add singlebond to compound"); }
	singlebonds[n_singlebonds++] = SingleBond(
		{
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[0]].local_id_compound),
			static_cast<uint32_t>(particle_info[bondtype.atom_indexes[1]].local_id_compound),
		},
		bondtype.b0,
		bondtype.kb
	);
}

template <> void CompoundFactory::addBond(const std::vector<ParticleInfo>& particle_info, const AngleBond& bondtype) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_COMPOUND) { throw std::exception("Failed to add anglebond to compound"); }
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

template <> void CompoundFactory::addBond(const std::vector<ParticleInfo>& particle_info, const DihedralBond& bondtype) {
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_COMPOUND) { 
		throw std::exception("Failed to add dihedralbond to compound"); }
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

template <> void CompoundFactory::addBond(const std::vector<ParticleInfo>& particle_info, const ImproperDihedralBond& bondtype) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_COMPOUND) { throw std::exception("Failed to add improperdihedralbond to compound"); }
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



template <> void BridgeFactory::addBond(std::vector<ParticleInfo>& particle_info, const SingleBond& bondtype) {
	if (n_singlebonds >= MAX_SINGLEBONDS_IN_BRIDGE) { throw std::exception("Failed to add singlebond to bridge"); }
	singlebonds[n_singlebonds++] = SingleBond{
		{
			getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[0]]),
			getBridgelocalIdOfParticle(particle_info[bondtype.atom_indexes[1]]),
		},
		bondtype.b0,
		bondtype.kb
	};
}

template <> void BridgeFactory::addBond(std::vector<ParticleInfo>& particle_info, const AngleBond& bondtype) {
	if (n_anglebonds >= MAX_ANGLEBONDS_IN_BRIDGE) { throw std::exception("Failed to add anglebond to bridge"); }
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

template <> void BridgeFactory::addBond(std::vector<ParticleInfo>& particle_info, const DihedralBond& bondtype) {

	if (bridge_id == 377) {
		printf("%d %d %d %d\n", particle_info[bondtype.atom_indexes[0]].gro_id, particle_info[bondtype.atom_indexes[1]].gro_id, particle_info[bondtype.atom_indexes[2]].gro_id, particle_info[bondtype.atom_indexes[3]].gro_id);
	}
	if (n_dihedrals >= MAX_DIHEDRALBONDS_IN_BRIDGE) { throw std::exception("Failed to add dihedralbond to bridge"); }
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

template <> void BridgeFactory::addBond(std::vector<ParticleInfo>& particle_info, const ImproperDihedralBond& bondtype) {
	if (n_improperdihedrals >= MAX_IMPROPERDIHEDRALBONDS_IN_BRIDGE) { throw std::exception("Failed to add dihedralbond to bridge"); }
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
		printf("Particle gro id %d\n", particle_info.gro_id);
		throw std::exception(std::format("Cannot add particle to this bridge ({}) as it already has another bridge ({})", bridge_id, particle_info.bridge_id).c_str());
	}

	// Particle already has a local id in the bridge
	if (particle_info.local_id_bridge != -1) {
		return particle_info.local_id_bridge;
	}
	// Else, make particle reference
	else {
		if (n_particles == MAX_PARTICLES_IN_BRIDGE) { throw std::exception("Failed to add particle to bridge"); }

		particle_info.local_id_bridge = n_particles;

		particle_refs[n_particles++] = ParticleReference{
			particle_info.compound_index,
			particle_info.local_id_compound,
	#ifdef LIMAKERNELDEBUGMODE
			particle_info.gro_id
	#endif
		};

		return particle_info.local_id_bridge;
	}	
}









CompoundCollection LIMA_MOLECULEBUILD::buildMolecules(
	Forcefield* ff,
	const std::string& work_dir,
	VerbosityLevel vl,
	const string& gro_path,
	const string& topol_path,
	std::unique_ptr<LimaLogger> logger,
	bool ignore_hydrogens
) 
{
	MoleculeBuilder mb{ ff, std::move(logger), work_dir, vl};
	return mb.buildMolecules(gro_path, topol_path, ignore_hydrogens);
}