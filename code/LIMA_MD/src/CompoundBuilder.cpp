#include "CompoundBuilder.h"

#include "Printer.h"
using namespace LIMA_Print;

CompoundBuilder::CompoundBuilder(ForceFieldMaker* ffm, VerbosityLevel vl) : verbosity_level(vl) {
	FFM = ffm;
	FFM->buildForcefield();
}

Molecule CompoundBuilder::buildMolecule(string gro_path, string itp_path, int max_residue_id, int min_residue_id, bool ignore_hydrogens) {


	printH2("Building molecule", true, false);
	compound_bridge_bundle = new CompoundBridgeBundle;
	particle_id_maps = new ParticleRef[MAX_ATOMS];

	vector<Record_ATOM> atom_data = parseGRO(gro_path);
	vector<vector<string>> top_data = parseTOP(itp_path);

	Molecule molecule;




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

		auto lut = molecule.bonded_particles_lut_manager->get(0, 0);
		for (int i = 0; i < MAX_COMPOUND_PARTICLES; i++) {
			for (int ii = 0; ii < MAX_COMPOUND_PARTICLES; ii++) {
				printf("%d ", *lut->get(i, ii));
			}
			printf("\n");
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

void CompoundBuilder::loadParticles(Molecule* molecule, vector<CompoundBuilder::Record_ATOM>* pdb_data, int max_residue_id, int min_residue_id, bool ignore_protons) {
	int current_res_id = -1;
	int current_compound_id = -1;
	Compound* current_compound = nullptr;
	//int prev_atom_id;


	for (Record_ATOM record : *pdb_data) {

		if (record.residue_seq_number < min_residue_id) { continue; }

		if (record.residue_seq_number > max_residue_id) { break; }

		if (ignore_protons && record.atom_name[0] == 'H') { continue; }

		if (record.residue_name == "SOL") { continue; }		


		if (record.residue_seq_number != current_res_id) {
			if (current_compound == nullptr || !current_compound->hasRoomForRes()) {	// TODO: Make this better............. probe how many particles beforehand
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
	molecule->bonded_particles_lut_manager->get(12, 2)->set(2, 3, true);


	int dihedral_cnt = 0;
	TopologyMode mode = INACTIVE;
	for (vector<string> record : *top_data) {
		if (record.size() == 0) {
			mode = INACTIVE;
			continue;
		}

		if (mode == INACTIVE) {
			mode = setMode(record[0]);

			if (mode == DIHEDRAL)			// Bad fix, but for now we ignore the bottom dihedral bonds, as i think they are IMPROPER DIHEDRALS
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
		

	PairBond* bondtype = FFM->getBondType(maps[0].global_id, maps[1].global_id);

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





	AngleBond* angletype = FFM->getAngleType(maps[0].global_id, maps[1].global_id, maps[2].global_id);

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


	DihedralBond* dihedraltype = FFM->getDihedralType(maps[0].global_id, maps[1].global_id, maps[2].global_id, maps[3].global_id);

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

	for (auto& record : records) {
		for (auto& w : record) {
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








