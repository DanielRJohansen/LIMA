#include "CompoundBuilder.h"

#include "Printer.h"

#include <algorithm>
#include <format>
using namespace LIMA_Print;

CompoundBuilder::CompoundBuilder(Forcefield* ff, VerbosityLevel vl) : verbosity_level(vl), forcefield{ ff } {}

CompoundCollection CompoundBuilder::buildCompoundCollection(string gro_path, string top_path, uint32_t max_residue_id, uint32_t min_residue_id, bool ignore_hydrogens) {


	printH2("Building molecule", true, false);
	compound_bridge_bundle = new CompoundBridgeBundle;
	particle_id_maps = new ParticleRef[MAX_ATOMS];

	atom_data = parseGRO(gro_path);
	vector<vector<string>> top_data = parseTOP(top_path);



	// Input:	residue_seq_id
	// Output:	0-indexed index of first atom belonging to that residue
	residueId_to_firstatomindex = makeResidueIdToAtomindexMap();
	bonded_atoms = makeBondedatomsLUT(top_path);

	assignMoleculeIDs();

	CompoundCollection molecule;




	// Load particles into molecule
	loadParticles(&molecule, &atom_data, max_residue_id, min_residue_id, ignore_hydrogens);
	//bonded_interactions_list = new LJ_Ignores[molecule.n_atoms_total * 10];		// DANGER - could get big. We need to ref the lists with particles global id, which comes directly from gro files, thus includes hydrogens and is 1-indexed!

	//printf("%d particles added\n", compound.n_particles);
	loadTopology(&molecule, &top_data);

	//molecule.compound_bridge_bundle
	//molecule.compound_bridge_bundle = new CompoundBridgeBundleCompact;// (&compound_bridge_bundle);		
	molecule.compound_bridge_bundle = new CompoundBridgeBundleCompact{ compound_bridge_bundle, verbosity_level>=V2 };	// Convert the host version to a compact device version, belonging to the molecule

	countElements(&molecule);
	if (verbosity_level >= V1) { printf("CompoundCollection built with %d compounds and %d bridges\n", molecule.n_compounds, molecule.compound_bridge_bundle->n_bridges); }

	for (int i = 0; i < molecule.n_compounds; i++) {
		molecule.compounds[i].calcParticleSphere();
	}

	
	if (verbosity_level >= V3) {
		for (int i = 1; i < MAX_ATOMS; i++) {
			if (particle_id_maps[i].global_id == -1)
				continue;
			printf("compound %d    local %d    global %d\n", particle_id_maps[i].compound_id, particle_id_maps[i].local_id_compound, particle_id_maps[i].global_id);
		}

		if (1) {
			auto lut = molecule.bonded_particles_lut_manager->get(0, 0);
			for (int i = 0; i < MAX_COMPOUND_PARTICLES; i++) {
				for (int ii = 0; ii < MAX_COMPOUND_PARTICLES; ii++) {
					printf("%d ", *lut->get(i, ii));
				}
				printf("\n");
			}
		}
		

		for (int i = 0; i < molecule.n_compounds; i++) {
			auto& compound = molecule.compounds[i];
			printf("´\nCompound %d\n", i);
			for (int p = 0; p < compound.n_particles; p++) {
				printf("\t p%d atomtype %d\n", p, compound.atom_types[p]);
			}
		}
	}


	

	delete[] particle_id_maps;
	delete compound_bridge_bundle;
	//delete[] bonded_interactions_list;

	printH2("Finished building molecule", false, true);

	return molecule;
}


vector<Float3> CompoundBuilder::getSolventPositions(string gro_path) {
	vector<Record_ATOM> atom_data = parseGRO(gro_path);
	
	vector<Float3> solvent_positions;
	for (Record_ATOM record : atom_data) {
		if (record.residue_name == "SOL" && record.atom_name[0] == 'O') {	// Easy solution, just say Oxygen is the center of the solvent. Ignore the hydrogens
			solvent_positions.push_back(record.coordinate);
		
		}
	}
	return solvent_positions;
}

void CompoundBuilder::loadParticles(CompoundCollection* compound_collection, vector<CompoundBuilder::Record_ATOM>* pdb_data, uint32_t max_residue_id, uint32_t min_residue_id, bool ignore_protons) {
	int current_res_id = -1;
	int current_compound_id = -1;
	int current_molecule_id = -1;
	CompoundCarrier* current_compound = nullptr;

	//ignore_protons = false;	// temp..

	for (Record_ATOM record : *pdb_data) {

		if (record.residue_seq_number < min_residue_id) { continue; }

		if (record.residue_seq_number > max_residue_id) { break; }

		if (ignore_protons && record.atom_name[0] == 'H') { continue; }

		if (record.residue_name == "SOL") { continue; }		

		const bool new_molecule = record.moleculeID != current_molecule_id;
		current_molecule_id = record.moleculeID;

		const bool new_residue = record.residue_seq_number != current_res_id;

		if (new_residue) {
			if (current_compound == nullptr || !current_compound->hasRoomForRes() || new_molecule) {	// TODO: Make this better............. probe how many particles beforehand
				//molecule->compound_bridge_bundle.addBridge(current_compound_id, current_compound_id + 1);

				if (current_compound_id != -1 && !new_molecule)		// Dont add bridge to first compound
					compound_bridge_bundle->addBridge(current_compound_id, current_compound_id + 1);

				current_compound_id++;
				compound_collection->compounds.emplace_back(CompoundCarrier{ current_compound_id });
				current_compound = &compound_collection->compounds.back();
				//current_compound = &molecule->compounds[current_compound_id];
				compound_collection->n_compounds++;

			}
			current_res_id = record.residue_seq_number;
		}


		//particle_id_maps[record.atom_serial_number] = IDMap(record.atom_serial_number, current_compound_id, molecule->compounds[current_compound_id].n_particles);
		if (particle_id_maps[record.atom_serial_number].compound_id != -1) {
			printf("WARNING: Added particle is NOT unique");
			exit(1);
		}
			
		particle_id_maps[record.atom_serial_number] = ParticleRef(record.atom_serial_number, current_compound_id, compound_collection->compounds[current_compound_id].n_particles);

		current_compound->addParticle(
			forcefield->getAtomtypeID(record.atom_serial_number),
			record.coordinate,
			forcefield->atomTypeToIndex(record.atom_name[0]),
			record.atom_serial_number);

		compound_collection->n_atoms_total++;
	}
}

void CompoundBuilder::loadTopology(CompoundCollection* compound_collection, vector<vector<string>>* top_data)
{
	//compound_collection->bonded_particles_lut_manager->get(12, 2)->set(2, 3, true);


	dihedral_sections_count = 0;	// bad fix-...
	TopologyMode mode = INACTIVE;
	for (vector<string> record : *top_data) {
		if (record.size() == 0) {
			mode = INACTIVE;
			continue;
		}
		bool new_section = setMode(record, mode);
		if (new_section) { continue; }

		//if (mode == INACTIVE) {
		//	mode = setMode(record[0]);

		//	if (mode == DIHEDRAL)			// Bad fix, but for now we ignore the bottom dihedral bonds, as i think they are IMPROPER DIHEDRALS
		//		dihedral_cnt++;

		//	continue;
		//}
		//if (mode == DIHEDRAL && dihedral_cnt > 1)
		//	continue;

		addGeneric(compound_collection, &record, mode);	// CANNOT ADD DIHEDRALS RIGHT NOW; DUE TO THE SECOND DIHEDRAL SECTION ALREADY BEING NOTED IN ANOTHER FUNCTION!! DANGER DANGER!
	}
}


bool CompoundBuilder::setMode(vector<string>& row, TopologyMode& current_mode)
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

		return true;
	}
	
	// By default we stay in the current mode.
	return false;
}

void CompoundBuilder::loadMaps(ParticleRef* maps, vector<string>* record, int n) {
	for (int i = 0; i < n; i++) {
		maps[i] = particle_id_maps[stoi((*record)[i])];
	}		
}

void CompoundBuilder::addGeneric(CompoundCollection* compound_collection, vector<string>* record, TopologyMode mode) {
	ParticleRef maps[4];

	switch (mode)
	{
	case INACTIVE:
		break;
	case BOND:
		addBond(compound_collection, maps, record);
		break;
	case ANGLE:
		addAngle(compound_collection, maps, record);
		break;

	case DIHEDRAL:
		addDihedral(compound_collection, maps, record);

		break;
	default:
		printf("Default case???!!?\n");
		return;
	}
}

void CompoundBuilder::addBond(CompoundCollection* compound_collection, ParticleRef* maps, vector<string>*record) {
	loadMaps(maps, record, 2);
	GenericBond g_bond = GenericBond(maps, 2);
	if (!g_bond.allParticlesExist()) {
		return;
	}
		

	const SingleBond* bondtype = forcefield->getBondType({ maps[0].global_id, maps[1].global_id });

	distributeLJIgnores(compound_collection, maps, 2);				// DANGER

	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &compound_collection->compounds[maps[0].compound_id];
		if (compound->n_singlebonds == MAX_PAIRBONDS) {
			printf("Too many bonds in compound\n");
			exit(0);
		}

		compound->singlebonds[compound->n_singlebonds++] = SingleBond(maps[0].local_id_compound, maps[1].local_id_compound, bondtype->b0, bondtype->kb);


	}
	else {

		// First, we need to make sure all bond particles are added to the bridge.
		// To create the single-bond we need to access bridge_local_indexes			
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, compound_collection);
		//bridge->addSinglebond(SingleBond(bondtype->b0, bondtype->kb, maps[0].global_id, maps[1].global_id));
		bridge->addGenericBond(SingleBond(bondtype->b0, bondtype->kb, maps[0].global_id, maps[1].global_id));
	}
}

void CompoundBuilder::addAngle(CompoundCollection* compound_collection, ParticleRef* maps, vector<string>* record) {
	loadMaps(maps, record, 3);
	GenericBond g_bond = GenericBond(maps, 3);
	if (!g_bond.allParticlesExist())
		return;





	const AngleBond* angletype = forcefield->getAngleType({ maps[0].global_id, maps[1].global_id, maps[2].global_id });

	distributeLJIgnores(compound_collection, maps, 3);



	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &compound_collection->compounds[maps[0].compound_id];
		if (compound->n_anglebonds == MAX_ANGLEBONDS) {
			printf("Too many angles in compound\n");
			exit(0);
		}

	//	printf("		Adding anglebond %d %d\n", maps[0].local_id_compound, maps[1].local_id_compound, maps[2].local_id_compound);
		compound->anglebonds[compound->n_anglebonds++] = AngleBond(maps[0].local_id_compound, maps[1].local_id_compound, maps[2].local_id_compound, angletype->theta_0, angletype->k_theta);
	}
	else {
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, compound_collection);
		bridge->addGenericBond(AngleBond(maps[0].global_id, maps[1].global_id, maps[2].global_id, angletype->theta_0, angletype->k_theta));
	}
}

void CompoundBuilder::addDihedral(CompoundCollection* molecule, ParticleRef* maps, vector<string>* record) {
	loadMaps(maps, record, 4);

	GenericBond g_bond = GenericBond(maps, 4);
	if (!g_bond.allParticlesExist())
		return;


	const DihedralBond* dihedraltype = forcefield->getDihedralType({ maps[0].global_id, maps[1].global_id, maps[2].global_id, maps[3].global_id });

	distributeLJIgnores(molecule, maps, 4);


	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &molecule->compounds[maps[0].compound_id];
		if (compound->n_dihedrals == MAX_DIHEDRALS) {
			printf("Too many dihedrals in compound\n");
			exit(0);
		}
		compound->dihedrals[compound->n_dihedrals++] = DihedralBond(maps[0].local_id_compound, maps[1].local_id_compound, maps[2].local_id_compound, maps[3].local_id_compound, dihedraltype->phi_0, dihedraltype->k_phi, dihedraltype->n);
	}
	else {
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, molecule);
		bridge->addGenericBond(DihedralBond(maps[0].global_id, maps[1].global_id, maps[2].global_id, maps[3].global_id, dihedraltype->phi_0, dihedraltype->k_phi, dihedraltype->n));
	}
}

void CompoundBuilder::distributeLJIgnores(CompoundCollection* molecule, ParticleRef* particle_refs, int n) {	// Works whether particle spans multiple compounds or not.
	for (int id_self = 0; id_self < n; id_self++) {
		for (int id_other = 0; id_other < n; id_other++) {
			if (id_self == id_other) { continue; }

			BondedParticlesLUTManager* lut_man = molecule->bonded_particles_lut_manager;
			BondedParticlesLUT* lut = lut_man->get(particle_refs[id_self].compound_id, particle_refs[id_other].compound_id);
			lut->set(particle_refs[id_self].local_id_compound, particle_refs[id_other].local_id_compound, true);
		}
	}
}

void CompoundBuilder::assignMoleculeIDs()
{
	printf("Assigning IDs to molecules\n");
	// Worst case scenario where each residue is a separate atom and a separate molecule (so solvents)
	vector<uint32_t> residueID_to_moleculeID(atom_data.size() + 1);	

	uint32_t current_molecule_id = 0;

	// Residue 1 (GMX 1-indexed) is always the first molecule (LIMA 0-indexed)
	residueID_to_moleculeID[1] = current_molecule_id;
	

	uint32_t residue_id_left = 1;
	uint32_t residue_id_right = 2;

	while (residue_id_right <= atom_data.back().residue_seq_number) {
		if (!areResiduesBonded(residue_id_left, residue_id_right)) { 
			current_molecule_id++; 
		}

		residueID_to_moleculeID[residue_id_right] = current_molecule_id;
		residue_id_left++;
		residue_id_right++;
	}

	for (auto& record : atom_data) {
		record.moleculeID = residueID_to_moleculeID[record.residue_seq_number];
	}
}

// This function assumes that the bonds in the topol.top file are always smallest on left side.
bool CompoundBuilder::areResiduesBonded(uint32_t res1, uint32_t res2) const 
{
	const uint32_t res1_firstatom_index = residueId_to_firstatomindex[res1];
	for (int i = res1_firstatom_index; i < atom_data.size(); i++) {
		if (atom_data[i].residue_seq_number != res1) { return false; }	// Finished, no match


		const uint32_t atom_id = atom_data[i].atom_serial_number;	// GMX 1-indexed
		auto& atom_bonds_to = bonded_atoms[atom_id];

		const uint32_t res2_firstatom_index = residueId_to_firstatomindex[res2];
		for (int j = res2_firstatom_index; j < atom_data.size(); j++) {

			if (atom_data[j].residue_seq_number != res2) { return false; }	// Finished, no match

			const uint32_t atom_id_query = atom_data[j].atom_serial_number;
			if (find(atom_bonds_to.begin(), atom_bonds_to.end(), atom_id_query) != atom_bonds_to.end()) {	// If found
				return true;
			}
		}
	}

	return false;
}


vector<uint32_t> CompoundBuilder::makeResidueIdToAtomindexMap()
{
	vector<uint32_t> residueIdToAtomindexMap;
	// Worst case where every residue contains 1 atom (So solvent e.g.). +1 for GMX 1-indexing
	residueIdToAtomindexMap.resize(atom_data.size() + 1);	

	int64_t currentResidueID = -1;

	for (uint32_t i = 0; i < atom_data.size(); i++) {
		auto& elem = atom_data[i];

		if (elem.residue_seq_number != currentResidueID) {
			currentResidueID = elem.residue_seq_number;
			residueIdToAtomindexMap[elem.residue_seq_number] = i;
		}
	}
	return residueIdToAtomindexMap;
}

// This function is BAAAD and it needs to go.
vector<vector<uint32_t>> CompoundBuilder::makeBondedatomsLUT(const string& topol_path)
{
	Topology topology = parseTop1(topol_path);
	vector<vector<uint32_t >> bondedatoms_lut;
	//bondedatoms_lut.resize(atom_data.size() + 1);	// +1 to account for GMX 1-indexing
	bondedatoms_lut.resize(std::max(atom_data.size() + 1, topology.bonds_data.size())*4);	// This is not right... but it works for a while

	for (auto& elem : topology.bonds_data) {
		bondedatoms_lut[elem.ai].push_back(elem.aj);
	}	
	return bondedatoms_lut;
}


Topology CompoundBuilder::parseTop1(string path) {		// Naive file segmentation, DANGEROUS!
	auto rows = Filehandler::readFile(path, {';', '#'}, {"[", "]"});
	Topology topology;



	int dihedral_cnt = 0;
	TopologyMode mode = INACTIVE;
	for (auto& row : rows) {
		bool new_section = setMode(row, mode);
		if (new_section) { continue; }

		switch (mode)
		{
		case INACTIVE:
			break;
		case ATOMS:
			//topology.addAtomsdataEntry(atol(row[0].c_str()), row[1], atol(row[2].c_str()));
			break;
		case BOND:
			topology.addBondsdataEntry(atol(row[0].c_str()), atol(row[1].c_str()), atoi(row[2].c_str()));
			break;
		case ANGLE:
			break;
		case DIHEDRAL:
			break;
		default:
			break;
		}
	}
	return topology;
}

vector<vector<string>> CompoundBuilder::parseTOP(string path)		// Naive file segmentation, DANGEROUS!
{
	if (verbosity_level >= V1) { cout << "Reading topology from file " << path << "\n"; }
	fstream file;
	file.open(path);
	int line_cnt = 0;
	string line;
	string space_delimiter = " ";

	vector<vector<string>> records;
	printf("done\n");
	auto rows = Filehandler::readFile(path, { ';', '#' }, { "[", "]" });
	printf("done\n");

	while (getline(file, line)) {
		vector<string> record;


		stringstream ss(line);
		string word;
		while (getline(ss, word, ' ')) {
			if (word == "" || word == "[" || word == "]") {
				continue;
			}
			record.push_back(word);
		}
		if (record.size() > 0) {
			if (record[0] == ";")
				continue;
		}
		records.push_back(std::move(record));
	}
	return records;
}

vector<CompoundBuilder::Record_ATOM> CompoundBuilder::parseGRO(string path)
{
	//cout << "Reading particles from file" << path << endl;
	if (verbosity_level >= V1) { cout <<"Reading particles from file " << path << "\n"; }
	fstream file;
	file.open(path);
	int line_cnt = 0;
	string line;
	string space_delimiter = " ";

	vector<Record_ATOM> records;

	while (getline(file, line)) {
		vector<string> words;
		stringstream ss(line);
		string word;
		while (getline(ss, word, ' ')) {
			if (word == "" || word == "[" || word == "]") {
				continue;
			}
			words.push_back(word);
		}


		ResidueComboId res = parseResidueID(words.at(0));

		if (words.size() == 5)
			words = splitAtomnameFromId(words);									// Disgusting function, and i hate GROMACS for making me do this!

		if (!res.valid || words.size() != 6) continue;

		if (words[0][0] == ';') continue;
			

		Record_ATOM record(
			stoi(words.at(2)),
			words.at(1),
			' ',
			res.name,
			' ',
			res.id,
			' ',
			Float3(stof(words.at(3)), stof(words.at(4)), stof(words.at(5)))
		);
		records.push_back(record);
	}

	file.close();

	return records;
}

CompoundBuilder::ResidueComboId CompoundBuilder::parseResidueID(string s)
{
	string id = "";
	string name = "";
	for (char c : s) {
		if (isAsciiNumber(c))
			id = id + c;
		else
			name = name + c;
	}
	
	if (id == "" || name == "")
		return ResidueComboId();
	return ResidueComboId(stoi(id), name);
}

void CompoundBuilder::countElements(CompoundCollection* molecule) {
	int counters[4] = { 0 };
	for (int c = 0; c < molecule->n_compounds; c++) {
		Compound* C = &molecule->compounds[c];
		counters[0] += C->n_particles;
		counters[1] += C->n_singlebonds;
		counters[2] += C->n_anglebonds;
		counters[3] += C->n_dihedrals;
	}

	if (verbosity_level >= 1) {
		printf("CompoundCollection created with\n");
		LIMA_Printer::printNameValuePairs("Particles", counters[0], "Bonds", counters[1], "Angles", counters[2], "Dihedrals", counters[3]);
	}
	/*
	printf("%d particles added\n", counters[0]);
	printf("%d singlebonds added\n", counters[1]);
	printf("%d anglebonds added\n", counters[2]);
	printf("%d dihedrals added\n", counters[3]);*/
}

vector<string> CompoundBuilder::splitAtomnameFromId(vector<string> words) {
	vector<string> words_;
	words_.push_back(words.at(0));
	if (words.at(1)[0] == 'O') {	// This is the most painful code i've been forced to write my entire life..
		words_.push_back("OW");
		words_.push_back(&words.at(1)[2]);
	}
	else {
		words_.push_back("HW");		// should be followed by 1 or 2, too lazy to implement
		words_.push_back(&words.at(1)[3]);
	}
	for (int i = 2; i < 5; i++)
		words_.push_back(words.at(i));

	return words_;

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

void MoleculeBuilder::buildMolecules(const string& gro_path, const string& topol_path, bool ignore_hydrogens) {

	printH2("Building molecules", true, false);

	loadResiduesAndSolvents(gro_path);

	loadTopology(topol_path);

	matchBondedResidues();

	createCompoundsAndBridges();

	distributeBondsToCompoundsAndBridges();

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
			const int64_t atoms_total = stoll(line);

			// Reserve enough vectorspace to handle all particles being either in residues or solvents
			residues.reserve(atoms_total / Residue::max_size_soft);	// Misses the mark when some residues are smaller..
			solvent_positions.reserve(atoms_total);
		}

		if (line.size() < 44) { continue; }	// Should only happen at the end. TODO: Extract box size here. Make sure it is the same as the constant var!

		GroRecord record = parseGroLine(line);

		if (record.atom_name[0] == 'H' && ignore_hydrogens) { 
			particle_info[record.atom_number].ignore = true;
			continue; 
		}

		// Is Solvent
		if (record.residue_name == "WATER" || record.residue_name == "SOL") {
			solvent_positions.push_back(record.position);
		}
		else {	// Is in residue
			// Belongs to current residue?
			if (!residues.empty() && residues.back().gro_id == record.residue_number) {				

				// SAFETY CHECK: This is likely because we have molecules consisting of only 1 residue, so we can figure out which belongs to the same molecule!
				if (residues.back().name != record.residue_name) { throw std::exception("Something went seriously wrong when handling the .gro file!"); }	
				// Do nothing
			}
			else {	// Create new residue
				int unique_res_id = residues.size();
				residues.push_back(Residue{ record.residue_number, unique_res_id, record.residue_name });
			}

			residues.back().atoms.push_back(AtomEntry{ record.atom_number, record.position, record.atom_name });
			n_particles_in_residues++;
		}
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
		if (particle_info[id].ignore) {
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
	for (size_t i = 0; i < residues.size() - 1; i++) {
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
				forcefield->atomTypeToIndex(atom.name[0])
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
		swap(out[0], out[1]);
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
void distributeLJIgnores(BondedParticlesLUTManager* bplut_man, const std::vector<ParticleInfo>& particle_info, const std::array<int, n_ids> gro_ids) {

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
			bridge.addSingleBond({ particle_info[bond_ids[0]], particle_info[bond_ids[1]]}, bondtype);
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
			bridge.addAngleBond({ particle_info[bond_ids[0]], particle_info[bond_ids[1]], particle_info[bond_ids[2]] }, bondtype);
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
			bridge.addDihedralBond({ particle_info[bond_ids[0]], particle_info[bond_ids[1]], particle_info[bond_ids[2]], particle_info[bond_ids[3]] }, bondtype);
		}
		else {
			const int compound_id = particle_info[bond_ids[0]].compound_index;
			CompoundFactory& compound = compounds.at(compound_id);
			compound.addDihedralBond({ particle_info[bond_ids[0]], particle_info[bond_ids[1]], particle_info[bond_ids[2]], particle_info[bond_ids[3]] }, bondtype);
		}

		distributeLJIgnores(bp_lut_manager, particle_info, bond_ids);
	}
}

/*
* 	SingleBond* bondtype = forcefield->getBondType({ maps[0].global_id, maps[1].global_id });

	distributeLJIgnores(compound_collection, maps, 2);				// DANGER

	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &compound_collection->compounds[maps[0].compound_id];
		if (compound->n_singlebonds == MAX_PAIRBONDS) {
			printf("Too many bonds in compound\n");
			exit(0);
		}

		compound->singlebonds[compound->n_singlebonds++] = SingleBond(maps[0].local_id_compound, maps[1].local_id_compound, bondtype->b0, bondtype->kb);


	}
	else {

		// First, we need to make sure all bond particles are added to the bridge.
		// To create the single-bond we need to access bridge_local_indexes			
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, compound_collection);
		//bridge->addSinglebond(SingleBond(bondtype->b0, bondtype->kb, maps[0].global_id, maps[1].global_id));
		bridge->addGenericBond(SingleBond(bondtype->b0, bondtype->kb, maps[0].global_id, maps[1].global_id));
*/








void CompoundFactory::addParticle(const Float3& position, int atomtype_id, int atomtype_color_id) {
	if (n_particles >= MAX_COMPOUND_PARTICLES) {
		printf("ERROR: Cannot add particle to compound!\n");
		exit(1);
	}

	positions[n_particles] = position;

	atom_types[n_particles] = atomtype_id;
	atom_color_types[n_particles] = atomtype_color_id;	// wtf is this


	n_particles++;
}

void CompoundFactory::addSingleBond(std::array<ParticleInfo, 2> particle_info, const SingleBond& bondtype)
{
	singlebonds[n_singlebonds++] = SingleBond(
		particle_info[0].local_id_compound,
		particle_info[1].local_id_compound,
		bondtype.b0,
		bondtype.kb
	);
}

void CompoundFactory::addAngleBond(std::array<ParticleInfo, 3> particle_info, const AngleBond& bondtype)
{
	anglebonds[n_anglebonds++] = AngleBond(
		particle_info[0].local_id_compound,
		particle_info[1].local_id_compound,
		particle_info[2].local_id_compound,
		bondtype.theta_0,
		bondtype.k_theta
	);
}

void CompoundFactory::addDihedralBond(std::array<ParticleInfo, 4> particle_info, const DihedralBond& bondtype)
{
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


void BridgeFactory::addSingleBond(std::array<ParticleInfo, 2> particle_info, const SingleBond& bondtype) {
	singlebonds[n_singlebonds++] = SingleBond{
		getBridgelocalIdOfParticle(particle_info[0]),
		getBridgelocalIdOfParticle(particle_info[1]),
		bondtype.b0,
		bondtype.kb
	};
}

void BridgeFactory::addAngleBond(std::array<ParticleInfo, 3> particle_info, const AngleBond& bondtype) {
	anglebonds[n_anglebonds++] = AngleBond{
		getBridgelocalIdOfParticle(particle_info[0]),
		getBridgelocalIdOfParticle(particle_info[1]),
		getBridgelocalIdOfParticle(particle_info[2]),
		bondtype.theta_0,
		bondtype.k_theta
	};
}

void BridgeFactory::addDihedralBond(std::array<ParticleInfo, 4> particle_info, const DihedralBond& bondtype) {
	dihedrals[n_dihedrals++] = DihedralBond{
		getBridgelocalIdOfParticle(particle_info[0]),
		getBridgelocalIdOfParticle(particle_info[1]),
		getBridgelocalIdOfParticle(particle_info[2]),
		getBridgelocalIdOfParticle(particle_info[3]),
		bondtype.phi_0,
		bondtype.k_phi,
		bondtype.n
	};
}

int BridgeFactory::getBridgelocalIdOfParticle(ParticleInfo& particle_info) {
	if (particle_info.local_id_bridge == -1) {
		particle_info.local_id_bridge = n_particles;
		particle_refs[n_particles++] = ParticleRefCompact{ particle_info.compound_index, particle_info.local_id_compound };
	}
	return particle_info.local_id_bridge;
}