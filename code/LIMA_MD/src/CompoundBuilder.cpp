#include "CompoundBuilder.h"

#include "Printer.h"

#include <algorithm>

using namespace LIMA_Print;

CompoundBuilder::CompoundBuilder(Forcefield* ff, VerbosityLevel vl) : verbosity_level(vl), forcefield{ ff } {}

Molecule CompoundBuilder::buildMolecule(string gro_path, string top_path, int max_residue_id, int min_residue_id, bool ignore_hydrogens) {


	printH2("Building molecule", true, false);
	compound_bridge_bundle = new CompoundBridgeBundle;
	particle_id_maps = new ParticleRef[MAX_ATOMS];

	atom_data = parseGRO(gro_path);
	vector<vector<string>> top_data = parseTOP(top_path);
	topology = parseTop1(top_path);


	// Input:	residue_seq_id
	// Output:	0-indexed index of first atom belonging to that residue
	residueId_to_firstatomindex = makeResidueIdToAtomindexMap();
	bonded_atoms = makeBondedatomsLUT();

	assignMoleculeIDs();

	Molecule molecule;




	// Load particles into molecule
	loadParticles(&molecule, &atom_data, max_residue_id, min_residue_id, ignore_hydrogens);
	//bonded_interactions_list = new LJ_Ignores[molecule.n_atoms_total * 10];		// DANGER - could get big. We need to ref the lists with particles global id, which comes directly from gro files, thus includes hydrogens and is 1-indexed!



	//printf("%d particles added\n", compound.n_particles);
	loadTopology(&molecule, &top_data);

	//molecule.compound_bridge_bundle
	//molecule.compound_bridge_bundle = new CompoundBridgeBundleCompact;// (&compound_bridge_bundle);		
	molecule.compound_bridge_bundle = new CompoundBridgeBundleCompact{ compound_bridge_bundle, verbosity_level>=V2 };	// Convert the host version to a compact device version, belonging to the molecule

	countElements(&molecule);
	if (verbosity_level >= V1) { printf("Molecule built with %d compounds and %d bridges\n", molecule.n_compounds, molecule.compound_bridge_bundle->n_bridges); }

	for (int i = 0; i < molecule.n_compounds; i++) {
		molecule.compounds[i].calcParticleSphere();
	}

	
	if (verbosity_level >= V3) {
		for (int i = 1; i < MAX_ATOMS; i++) {
			if (particle_id_maps[i].global_id == -1)
				continue;
			printf("compound %d    local %d    global %d\n", particle_id_maps[i].compound_id, particle_id_maps[i].local_id_compound, particle_id_maps[i].global_id);
		}

		if (0) {
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
			solvent_positions.push_back(record.coordinate / LIMA_SCALE);
		
		}
	}
	return solvent_positions;
}

void CompoundBuilder::loadParticles(Molecule* molecule, vector<CompoundBuilder::Record_ATOM>* pdb_data, int max_residue_id, int min_residue_id, bool ignore_protons) {
	int current_res_id = -1;
	int current_compound_id = -1;
	int current_molecule_id = -1;
	Compound* current_compound = nullptr;



	for (Record_ATOM record : *pdb_data) {

		if (record.residue_seq_number < min_residue_id) { continue; }

		if (record.residue_seq_number > max_residue_id) { break; }

		if (ignore_protons && record.atom_name[0] == 'H') { continue; }

		if (record.residue_name == "SOL") { continue; }		

		bool new_molecule = record.moleculeID != current_molecule_id;
		current_molecule_id = record.moleculeID;

		if (record.residue_seq_number != current_res_id) {
			if (current_compound == nullptr || !current_compound->hasRoomForRes() || new_molecule) {	// TODO: Make this better............. probe how many particles beforehand
				//molecule->compound_bridge_bundle.addBridge(current_compound_id, current_compound_id + 1);

				if (current_compound_id != -1 && !new_molecule)		// Dont add bridge to first compound
					compound_bridge_bundle->addBridge(current_compound_id, current_compound_id + 1);

				current_compound_id++;
				current_compound = &molecule->compounds[current_compound_id];
				molecule->n_compounds++;

			}
			current_res_id = record.residue_seq_number;
		}


		//particle_id_maps[record.atom_serial_number] = IDMap(record.atom_serial_number, current_compound_id, molecule->compounds[current_compound_id].n_particles);
		if (particle_id_maps[record.atom_serial_number].compound_id != -1) {
			printf("WARNING: Added particle is NOT unique");
			exit(1);
		}
			
		particle_id_maps[record.atom_serial_number] = ParticleRef(record.atom_serial_number, current_compound_id, molecule->compounds[current_compound_id].n_particles);

		current_compound->addParticle(
			forcefield->getAtomtypeID(record.atom_serial_number),
			record.coordinate / LIMA_SCALE,
			forcefield->atomTypeToIndex(record.atom_name[0]),
			record.atom_serial_number);

		molecule->n_atoms_total++;
	}
}

void CompoundBuilder::loadTopology(Molecule* molecule, vector<vector<string>>* top_data)
{
	molecule->bonded_particles_lut_manager->get(12, 2)->set(2, 3, true);


	int dihedral_cnt = 0;
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

		addGeneric(molecule, &record, mode);	// CANNOT ADD DIHEDRALS RIGHT NOW; DUE TO THE SECOND DIHEDRAL SECTION ALREADY BEING NOTED IN ANOTHER FUNCTION!! DANGER DANGER!
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

void CompoundBuilder::addGeneric(Molecule* molecule, vector<string>* record, TopologyMode mode) {
	ParticleRef maps[4];

	switch (mode)
	{
	case CompoundBuilder::INACTIVE:
		break;
	case CompoundBuilder::BOND:
		addBond(molecule, maps, record);
		break;
	case CompoundBuilder::ANGLE:
		addAngle(molecule, maps, record);
		break;

	case CompoundBuilder::DIHEDRAL:
		addDihedral(molecule, maps, record);

		break;
	default:
		printf("Default case???!!?\n");
		return;
	}
}

void CompoundBuilder::addBond(Molecule* molecule, ParticleRef* maps, vector<string>*record) {
	loadMaps(maps, record, 2);
	GenericBond g_bond = GenericBond(maps, 2);
	if (!g_bond.allParticlesExist()) {
		return;
	}
		

	PairBond* bondtype = forcefield->getBondType(maps[0].global_id, maps[1].global_id);

	distributeLJIgnores(molecule, maps, 2);				// DANGER

	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &molecule->compounds[maps[0].compound_id];
		if (compound->n_singlebonds == MAX_PAIRBONDS) {
			printf("Too many bonds in compound\n");
			exit(0);
		}

		compound->singlebonds[compound->n_singlebonds++] = PairBond(maps[0].local_id_compound, maps[1].local_id_compound, bondtype->b0, bondtype->kb);


	}
	else {

		// First, we need to make sure all bond particles are added to the bridge.
		// To create the single-bond we need to access bridge_local_indexes			
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, molecule);
		//bridge->addSinglebond(PairBond(bondtype->b0, bondtype->kb, maps[0].global_id, maps[1].global_id));
		bridge->addGenericBond(PairBond(bondtype->b0, bondtype->kb, maps[0].global_id, maps[1].global_id));
	}
}

void CompoundBuilder::addAngle(Molecule* molecule, ParticleRef* maps, vector<string>* record) {
	loadMaps(maps, record, 3);
	GenericBond g_bond = GenericBond(maps, 3);
	if (!g_bond.allParticlesExist())
		return;





	AngleBond* angletype = forcefield->getAngleType(maps[0].global_id, maps[1].global_id, maps[2].global_id);

	distributeLJIgnores(molecule, maps, 3);



	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &molecule->compounds[maps[0].compound_id];
		if (compound->n_anglebonds == MAX_ANGLEBONDS) {
			printf("Too many angles in compound\n");
			exit(0);
		}

	//	printf("		Adding anglebond %d %d\n", maps[0].local_id_compound, maps[1].local_id_compound, maps[2].local_id_compound);
		compound->anglebonds[compound->n_anglebonds++] = AngleBond(maps[0].local_id_compound, maps[1].local_id_compound, maps[2].local_id_compound, angletype->theta_0, angletype->k_theta);
	}
	else {
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, molecule);
		bridge->addGenericBond(AngleBond(maps[0].global_id, maps[1].global_id, maps[2].global_id, angletype->theta_0, angletype->k_theta));
	}
}

void CompoundBuilder::addDihedral(Molecule* molecule, ParticleRef* maps, vector<string>* record) {
	loadMaps(maps, record, 4);

	GenericBond g_bond = GenericBond(maps, 4);
	if (!g_bond.allParticlesExist())
		return;


	DihedralBond* dihedraltype = forcefield->getDihedralType(maps[0].global_id, maps[1].global_id, maps[2].global_id, maps[3].global_id);

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

void CompoundBuilder::distributeLJIgnores(Molecule* molecule, ParticleRef* particle_refs, int n) {	// Works whether particle spans multiple compounds or not.
	for (int id_self = 0; id_self < n; id_self++) {
		for (int id_other = 0; id_other < n; id_other++) {
			
			// Not sure if i should not include these too??
			if (id_self == id_other)
				continue;


			BondedParticlesLUTManager* lut_man = molecule->bonded_particles_lut_manager;
			BondedParticlesLUT* lut = lut_man->get(particle_refs[id_self].compound_id, particle_refs[id_other].compound_id);
			lut->set(particle_refs[id_self].local_id_compound, particle_refs[id_other].local_id_compound, true);
		}
	}
}

void CompoundBuilder::assignMoleculeIDs()
{
	vector<uint32_t> residueID_to_moleculeID;
	// Worst case scenario where each residue is a separate atom and a separate molecule (so solvents)
	residueID_to_moleculeID.resize(atom_data.size() + 1);	

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
bool CompoundBuilder::areResiduesBonded(uint32_t res1, uint32_t res2)
{
	uint32_t res1_firstatom_index = residueId_to_firstatomindex[res1];
	for (int i = res1_firstatom_index; i < atom_data.size(); i++) {
		if (atom_data[i].residue_seq_number != res1) { return false; }	// Finished, no match


		uint32_t atom_id = atom_data[i].atom_serial_number;	// GMX 1-indexed
		auto& atom_bonds_to = bonded_atoms[atom_id];

		uint32_t res2_firstatom_index = residueId_to_firstatomindex[res2];
		for (int j = res2_firstatom_index; j < atom_data.size(); j++) {

			if (atom_data[j].residue_seq_number != res2) { return false; }	// Finished, no match

			uint32_t atom_id_query = atom_data[j].atom_serial_number;
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

vector<vector<uint32_t>> CompoundBuilder::makeBondedatomsLUT()
{
	vector<vector<uint32_t >> bondedatoms_lut;
	bondedatoms_lut.resize(atom_data.size() + 1);	// +1 to account for GMX 1-indexing

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
		case CompoundBuilder::INACTIVE:
			break;
		case CompoundBuilder::ATOMS:
			//topology.addAtomsdataEntry(atol(row[0].c_str()), row[1], atol(row[2].c_str()));
			break;
		case CompoundBuilder::BOND:
			topology.addBondsdataEntry(atol(row[0].c_str()), atol(row[1].c_str()), atoi(row[2].c_str()));
			break;
		case CompoundBuilder::ANGLE:
			break;
		case CompoundBuilder::DIHEDRAL:
			break;
		default:
			break;
		}
	}
	return topology;
}

vector<vector<string>> CompoundBuilder::parseTOP(string path)		// Naive file segmentation, DANGEROUS!
{
	fstream file;
	file.open(path);
	int line_cnt = 0;
	string line;
	string space_delimiter = " ";

	vector<vector<string>> records;

	auto rows = Filehandler::readFile(path, { ';', '#' }, { "[", "]" });


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
		records.push_back(record);
	}
	return records;
}

vector<CompoundBuilder::Record_ATOM> CompoundBuilder::parsePDB(string path)
{
	fstream file;
	file.open(path);
	int line_cnt = 0;

	int endpoints[] = { 4, 11, 16 , 17, 20, 22, 26, 27, 38, 46, 54 };


	vector<Record_ATOM> records;

	string line;
	while (getline(file, line)) {
		stringstream ss(line);
		string row_type;
		getline(ss, row_type, ' ');
		if (row_type != "ATOM")
			continue;

		vector<string> data_buffer;

		int ptr = 0;
		for (int stop : endpoints) {
			string word = "";
			
			while (ptr < stop) {
				if (line[ptr] != ' ')
					word = word + line[ptr];
				ptr++;
			}
			//cout << "Word:" << word << endl;
			data_buffer.push_back(word);
		}

		cout << data_buffer[6] << endl;
		//printf("%d\n", stoi(data_buffer[6]));
		
		records.push_back(Record_ATOM(
			stoi(data_buffer[1]),
			data_buffer[2],
			data_buffer[3][0],			// alt loc
			data_buffer[4],
			data_buffer[5][0],		// chain id
			stoi(data_buffer[6]),
			data_buffer[7][0],
			Float3(stof(data_buffer[8]), stof(data_buffer[9]), stof(data_buffer[10])) * 0.1f	// Convert A to nm right off the bat!
		));
		
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

void CompoundBuilder::countElements(Molecule* molecule) {
	int counters[4] = { 0 };
	for (int c = 0; c < molecule->n_compounds; c++) {
		Compound* C = &molecule->compounds[c];
		counters[0] += C->n_particles;
		counters[1] += C->n_singlebonds;
		counters[2] += C->n_anglebonds;
		counters[3] += C->n_dihedrals;
	}

	if (verbosity_level >= 1) {
		printf("Molecule created with\n");
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






