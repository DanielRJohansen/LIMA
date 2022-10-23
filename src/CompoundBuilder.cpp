#include "CompoundBuilder.h"





Molecule CompoundBuilder::buildMolecule(string gro_path, string itp_path, int max_residue_id, int min_residue_id, bool ignore_hydrogens) {
	//vector<vector<string>> pdb_data = readFile(pdb_path);
	//vector<Record_ATOM> atom_data = parsePDB(pdb_path);
	printf("\n\n############################# BUILDING MOLECULE #############################\n\n\n");

	compound_bridge_bundle = new CompoundBridgeBundle;
	particle_id_maps = new ParticleRef[MAX_ATOMS];


	


	vector<Record_ATOM> atom_data = parseGRO(gro_path);
	vector<vector<string>> top_data = parseTOP(itp_path);


	//Compound compound;													// position of particles stored in compund.particles.pos_tsub1
	Molecule molecule;

	FFM->buildForcefield();



	//loadParticles(&compound, &atom_data, max_residue_id, true);
	loadParticles(&molecule, &atom_data, max_residue_id, min_residue_id, ignore_hydrogens);
	bonded_interactions_list = new LJ_Ignores[molecule.n_atoms_total * 10];		// DANGER - could get big. We need to ref the lists with particles global id, which comes directly from gro files, thus includes hydrogens and is 1-indexed!



	//printf("%d particles added\n", compound.n_particles);
	loadTopology(&molecule, &top_data);

	//molecule.compound_bridge_bundle
	molecule.compound_bridge_bundle = new CompoundBridgeBundleCompact;// (&compound_bridge_bundle);		// Convert the host version to a compact device version, belonging to the molecule
	*molecule.compound_bridge_bundle = CompoundBridgeBundleCompact(compound_bridge_bundle);

	countElements(&molecule);
	printf("Molecule built with %d compounds and %d bridges\n\n\n", molecule.n_compounds, molecule.compound_bridge_bundle->n_bridges);

	for (int i = 0; i < molecule.n_compounds; i++) {
		molecule.compounds[i].calcParticleSphere();
	}

	
	if (0) {
		for (int i = 1; i < MAX_ATOMS; i++) {
			if (particle_id_maps[i].global_id == -1)
				continue;
			printf("compound %d    local %d    global %d\n", particle_id_maps[i].compound_id, particle_id_maps[i].local_id_compound, particle_id_maps[i].global_id);
		}	
	}
	


	delete[] particle_id_maps;
	delete compound_bridge_bundle;
	delete[] bonded_interactions_list;

	printf("\n\n############################# FINISHED BUILDING MOLECULE #############################\n\n\n");

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




















void CompoundBuilder::loadParticles(Molecule* molecule, vector<CompoundBuilder::Record_ATOM>* pdb_data, int max_residue_id, int min_residue_id, bool ignore_protons) {
	int current_res_id = -1;
	int current_compound_id = -1;
	Compound* current_compound = nullptr;
	//int prev_atom_id;


	for (Record_ATOM record : *pdb_data) {

			/*
			int particle_id = stoi(record[1]);
			string atom_name = record[2];
			string monomer_name = record[3];
			cout << record[3] <<endl;
			cout << record[4] << endl;
			int monomer_id = stoi(record[4]);
			Float3 coord(stod(record[5]), stod(record[6]), stod(record[7]));
			coord *= 0.1f;	// convert A to nm
			*/
			//printf("res %d   max res %d\n", record.residue_seq_number, max_residue_id);

		if (record.residue_seq_number < min_residue_id)
			continue;

		if (record.residue_seq_number > max_residue_id)
			break;

		if (ignore_protons && record.atom_name[0] == 'H')
			continue;

		if (record.residue_name == "SOL")
			continue;
		



		if (record.residue_seq_number != current_res_id) {
			if (current_compound == nullptr || !current_compound->hasRoomForRes()) {
				//molecule->compound_bridge_bundle.addBridge(current_compound_id, current_compound_id + 1);

				if (current_compound_id != -1)		// Dont add bridge to first compound
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
		//current_compound->addParticle(FFM.atomTypeToIndex(record.atom_name[0]), CompactParticle(record.coordinate));


		//ParticleParameters pp1 = FFM.getForcefield().particle_parameters[FFM.atomTypeToIndex(record.atom_name[0])];
		//ParticleParameters pp2 = FFM.getNBForcefield().particle_parameters[FFM.getAtomtypeID(record.atom_serial_number)];
		//printf("Change: %f %f %f\n\n", pp2.mass / pp1.mass, pp2.sigma / pp1.sigma, pp2.epsilon / pp1.epsilon);

		//current_compound->addParticle(FFM->atomTypeToIndex(record.atom_name[0]), record.coordinate);
		//current_compound->addParticle(FFM->getAtomtypeID(record.atom_serial_number), record.coordinate);
		current_compound->addParticle(FFM->getAtomtypeID(record.atom_serial_number), record.coordinate, FFM->atomTypeToIndex(record.atom_name[0]), record.atom_serial_number);
		//record.coordinate.print('p');
		molecule->n_atoms_total++;
	}
}

void CompoundBuilder::loadTopology(Molecule* molecule, vector<vector<string>>* top_data)
{	
	/*bonded_interactions_list = new uint32_t * [molecule->n_atoms_total];			// Just delete this section now
	for (uint32_t i = 0; i < molecule->n_atoms_total; i++) {
		bonded_interactions_list[i] = new uint32_t[MAX_BONDED_INTERACTIONS];
		for (uint32_t ii = 0; ii < MAX_BONDED_INTERACTIONS; ii++) {
			bonded_interactions_list[i][ii] = UINT32_MAX;
		}
	}*/



	int dihedral_cnt = 0;
	TopologyMode mode = INACTIVE;
	for (vector<string> record : *top_data) {
		if (record.size() == 0) {
			mode = INACTIVE;
			continue;
		}

		if (mode == INACTIVE) {
			mode = setMode(record[0]);

			if (mode == DIHEDRAL)			// Bad fix, but for now we ignore the lowest dihedral bonds, as i think they are IMPROPER DIHEDRALS
				dihedral_cnt++;

			continue;
		}
		if (mode == DIHEDRAL && dihedral_cnt > 1)
			continue;

		addGeneric(molecule, &record, mode);
	}
}


CompoundBuilder::TopologyMode CompoundBuilder::setMode(string entry)
{
	if (entry == "bonds")
		return BOND;
	if (entry == "angles")
		return ANGLE;
	if (entry == "dihedrals")
		return DIHEDRAL;
	return INACTIVE;
}

void CompoundBuilder::loadMaps(ParticleRef* maps, vector<string>* record, int n) {
	for (int i = 0; i < n; i++) {
		maps[i] = particle_id_maps[stoi((*record)[i])];
	}
		
}

void CompoundBuilder::addGeneric(Molecule* molecule, vector<string>* record, TopologyMode mode) {
	//GenericBond g_bond;
	//PairBond s_bond;
	//AngleBond a_bond;
	ParticleRef maps[4];

	//PairBond* bondtype;
	//AngleBond* angletype;
	//DihedralBond* dihedraltype;


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
		

	PairBond* bondtype = FFM->getBondType(maps[0].global_id, maps[1].global_id);

	distributeLJIgnores(molecule, maps, 2);				// DANGER
	bool apply_BLJI = checkIfFirstBondedInteraction(molecule, maps, 2);

	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &molecule->compounds[maps[0].compound_id];
		if (compound->n_singlebonds == MAX_PAIRBONDS) {
			printf("Too many bonds in compound\n");
			exit(0);
		}

		compound->singlebonds[compound->n_singlebonds++] = PairBond(bondtype->b0, bondtype->kb, maps[0].local_id_compound, maps[1].local_id_compound, apply_BLJI);		
	}
	else {

		// First, we need to make sure all bond particles are added to the bridge.
		// To create the single-bond we need to access bridge_local_indexes			
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, molecule);
		//bridge->addSinglebond(PairBond(bondtype->b0, bondtype->kb, maps[0].global_id, maps[1].global_id));
		bridge->addGenericBond(PairBond(bondtype->b0, bondtype->kb, maps[0].global_id, maps[1].global_id, apply_BLJI));
	}
}

void CompoundBuilder::addAngle(Molecule* molecule, ParticleRef* maps, vector<string>* record) {
	loadMaps(maps, record, 3);
	GenericBond g_bond = GenericBond(maps, 3);
	if (!g_bond.allParticlesExist())
		return;





	AngleBond* angletype = FFM->getAngleType(maps[0].global_id, maps[1].global_id, maps[2].global_id);

	distributeLJIgnores(molecule, maps, 3);
	bool apply_BLJI = checkIfFirstBondedInteraction(molecule, maps, 3);



	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &molecule->compounds[maps[0].compound_id];
		if (compound->n_anglebonds == MAX_ANGLEBONDS) {
			printf("Too many angles in compound\n");
			exit(0);
		}

	//	printf("		Adding anglebond %d %d\n", maps[0].local_id_compound, maps[1].local_id_compound, maps[2].local_id_compound);
		compound->anglebonds[compound->n_anglebonds++] = AngleBond(maps[0].local_id_compound, maps[1].local_id_compound, maps[2].local_id_compound, angletype->theta_0, angletype->k_theta, apply_BLJI);
	}
	else {
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, molecule);
		bridge->addGenericBond(AngleBond(maps[0].global_id, maps[1].global_id, maps[2].global_id, angletype->theta_0, angletype->k_theta, apply_BLJI));
	}
}

void CompoundBuilder::addDihedral(Molecule* molecule, ParticleRef* maps, vector<string>* record) {
	loadMaps(maps, record, 4);

	GenericBond g_bond = GenericBond(maps, 4);
	if (!g_bond.allParticlesExist())
		return;


	DihedralBond* dihedraltype = FFM->getDihedralType(maps[0].global_id, maps[1].global_id, maps[2].global_id, maps[3].global_id);

	distributeLJIgnores(molecule, maps, 4);
	bool apply_BLJI = checkIfFirstBondedInteraction(molecule, maps, 4);


	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &molecule->compounds[maps[0].compound_id];
		if (compound->n_dihedrals == MAX_DIHEDRALS) {
			printf("Too many dihedrals in compound\n");
			exit(0);
		}
		compound->dihedrals[compound->n_dihedrals++] = DihedralBond(maps[0].local_id_compound, maps[1].local_id_compound, maps[2].local_id_compound, maps[3].local_id_compound, dihedraltype->phi_0, dihedraltype->k_phi, dihedraltype->n, apply_BLJI);
	}
	else {
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, molecule);
		bridge->addGenericBond(DihedralBond(maps[0].global_id, maps[1].global_id, maps[2].global_id, maps[3].global_id, dihedraltype->phi_0, dihedraltype->k_phi, dihedraltype->n, apply_BLJI));
	}
}

void CompoundBuilder::distributeLJIgnores(Molecule* molecule, ParticleRef* particle_refs, int n) {	// Works whether particle spans multiple compounds or not.
	// This way only connects the edges of each bond, angle or dihedral
	// ALWAYS do angles first, or this breaks!


	// System 3, bond object based




	/*
	// System 2, never finished
	if (n == 4)  {	// First check if connection is already made. I guess i could do this irrelevant of n
		Compound* compound = &molecule->compounds[particle_refs[0].compound_id];
		if (compound->lj_ignore_list[particle_refs[0].local_id_compound].checkAlreadyConnected(particle_refs[n-1].global_id)) {
			compound->lj_ignore_list[particle_refs[0].local_id_compound].assignDoublyConnectedID(unique_doublyconnected_id);

			compound = &molecule->compounds[particle_refs[n-1].compound_id];
			compound->lj_ignore_list[particle_refs[n-1].local_id_compound].assignDoublyConnectedID(unique_doublyconnected_id);
			printf("Assigning ID %d - %d %d\n", unique_doublyconnected_id, particle_refs[0].global_id, particle_refs[n-1].global_id);
			unique_doublyconnected_id++;
		}
	}
	// Then make connection
	Compound* compound = &molecule->compounds[particle_refs[0].compound_id];
	compound->lj_ignore_list[particle_refs[0].local_id_compound].addIgnoreTarget(particle_refs[n - 1].global_id);

	compound = &molecule->compounds[particle_refs[n - 1].compound_id];
	compound->lj_ignore_list[particle_refs[n - 1].local_id_compound].addIgnoreTarget(particle_refs[0].global_id);
	*/

	// This approach connects all elements of the bond, angle or dihedral
		// OLD system
	for (int id_self = 0; id_self < n; id_self++) {
		for (int id_other = 0; id_other < n; id_other++) {
	//for (int id_self = 0; id_self < n; id_self+=(n-1)) {
		//for (int id_other = 0; id_other < n; id_other+=(n-1)) {
			if (id_self == id_other)
				continue;

			Compound* compound_self = &molecule->compounds[particle_refs[id_self].compound_id];
			compound_self->lj_ignore_list[particle_refs[id_self].local_id_compound].addIgnoreTarget(particle_refs[id_other].local_id_compound, particle_refs[id_other].compound_id);
		}
	}
}

bool CompoundBuilder::checkIfFirstBondedInteraction(Molecule* molecule, ParticleRef* particle_maps, int n) {	// Works whether particle spans multiple compounds or not.
	uint32_t global_id_left = particle_maps[0].global_id;		// left and right is arbitrary.
	uint32_t global_id_right = particle_maps[n-1].global_id;
	
	bool first_connection = !bonded_interactions_list[global_id_left].isAlreadyConnected(global_id_right);

	if (min(global_id_left, global_id_right) == 512 && max(global_id_left, global_id_right) == 515) {
		//printf("Checking the two, n: %d     ac: %d\n", n, first_connection);
	}

	if (!first_connection) {
		//printf("Not first %d\n", n);
	}

	bonded_interactions_list[global_id_left].addConnection(global_id_right);
	bonded_interactions_list[global_id_right].addConnection(global_id_left);


	//printf("%d %d %d\n", particle_refs[0].compound_id, particle_refs[0].global_id, molecule->n_atoms_total);
	
	return first_connection;
}



vector<vector<string>> CompoundBuilder::parseTOP(string path)		// Naive file segmentation, DANGEROUS!
{
	fstream file;
	file.open(path);
	int line_cnt = 0;
	string line;
	string space_delimiter = " ";

	vector<vector<string>> records;

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

	for (auto record : records) {
		for (auto w : record) {
			//cout << w << " ";
		}
		//printf("\n");
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
			Float3(stof(data_buffer[8]), stof(data_buffer[9]), stof(data_buffer[10])) * 0.1	// Convert A to nm right off the bat!
		));
		
	}
	return records;
}

vector<CompoundBuilder::Record_ATOM> CompoundBuilder::parseGRO(string path)
{
	cout << "Reading particles from file" << path << endl;
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

		if (!res.valid || words.size() != 6)
			continue;

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

	printf("Molecule created with %d compounds\n", molecule->n_compounds);
	printf("%d particles added\n", counters[0]);
	printf("%d singlebonds added\n", counters[1]);
	printf("%d anglebonds added\n", counters[2]);
	printf("%d dihedrals added\n", counters[3]);
}















Molecule::Molecule() {
	compounds = new Compound[MAX_COMPOUNDS];
	//compound_bridge_bundle = new CompoundBridgeBundleCompact;
}


