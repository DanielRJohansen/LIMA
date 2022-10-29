#include "Environment.h"


Environment::Environment() {
}
Environment::Environment(string conf_filename, string topol_filename) {
	simulation = new Simulation();
	verifySimulationParameters();


	display = new DisplayV2();


	ForceFieldMaker* forcefieldmaker = new ForceFieldMaker();
	compoundbuilder = new CompoundBuilder(forcefieldmaker);


	
	int min_res_id = 0;
	int max_res_id = 20;
	bool ignore_hydrogens = true;
	Molecule mol_6lzm_10 = compoundbuilder->buildMolecule(MOL_FOLDER + conf_filename, MOL_FOLDER + topol_filename, max_res_id, min_res_id, ignore_hydrogens);
	//Molecule mol_6lzm_10 = compoundbuilder->buildMolecule(MOL_FOLDER + "conf.gro", MOL_FOLDER + "topol.top", max_res_id, min_res_id, ignore_hydrogens);
	//Molecule mol_6lzm_10 = compoundbuilder->buildMolecule(MOL_FOLDER + "conf_test.gro", MOL_FOLDER + "topol_test.top", max_res_id, min_res_id, ignore_hydrogens);
	//vector<Float3> solvent_positions = compoundbuilder->getSolventPositions(MOL_FOLDER + "box7.gro");
	vector<Float3> solvent_positions = compoundbuilder->getSolventPositions(MOL_FOLDER + "conf.gro");


	boxbuilder.buildBox(simulation);
	boxbuilder.addSingleMolecule(simulation, &mol_6lzm_10);

#ifdef ENABLE_SOLVENTS
	boxbuilder.solvateBox(simulation, &solvent_positions);
#endif
	//exit(1);

	//boxbuilder.addScatteredMolecules(simulation, &mol_dpc, N_LIPID_COPIES);
	//boxbuilder.solvateBox(simulation);	// Always do after placing compounds
	delete[] mol_6lzm_10.compounds;
	boxbuilder.finishBox(simulation);


	simulation->moveToDevice();	// Only moves the Box to the device
	verifyBox();


 

	if (print_compound_positions) {	
		for (int c = 0; c < simulation->n_compounds; c++) {
			Compound* comp = &simulation->compounds_host[c];
			for (int p = 0; p < comp->n_particles; p++) {
				printf("%d   ", comp->particle_global_ids[p]);
				simulation->box->compound_state_array[c].positions[p].print();
			}
		}
	}

	engine = new Engine(simulation, forcefieldmaker->getNBForcefield());
}


void Environment::verifySimulationParameters() {	// Not yet implemented
	assert(THREADS_PER_COMPOUNDBLOCK >= MAX_COMPOUND_PARTICLES);
	//assert(THREADS_PER_SOLVENTBLOCK >= N_SOLVATE_MOLECULES);
	assert(BOX_LEN > 3.f);
	//assert(BOX_LEN >= CUTOFF + 0.5f);
	//assert(simulation->n_compounds <= 1);	// Otherwise data_GAN goes haywire

	assert(simulation->n_steps % STEPS_PER_LOGTRANSFER == 0);
	assert(simulation->n_steps % STEPS_PER_THERMOSTAT == 0);
	assert(simulation->n_steps % STEPS_PER_TRAINDATATRANSFER == 0);

	assert(STEPS_PER_THERMOSTAT % STEPS_PER_LOGTRANSFER == 0);		// Change to trajtransfer later
	//assert(STEPS_PER_THERMOSTAT >= STEPS_PER_LOGTRANSFER);
	assert(THREADS_PER_SOLVENTBLOCK >= MAX_COMPOUND_PARTICLES);

	assert(STEPS_PER_THERMOSTAT >= STEPS_PER_LOGTRANSFER);


	printf("Simulation parameters verified\n");
}

void Environment::verifyBox() {
	for (int c = 0; c < simulation->n_compounds; c++) {
		//printf("Compound radius: %f\t center: %f %f %f\n", simulation->compounds_host[c].confining_particle_sphere, simulation->compounds_host[c].center_of_mass.x, simulation->compounds_host[c].center_of_mass.y, simulation->compounds_host[c].center_of_mass.z);
		if ((simulation->compounds_host[c].confining_particle_sphere * 1.1) > BOX_LEN_HALF) {
			printf("Compound %d too large for simulation-box\n", c);
			exit(1);
		}
	}
	printf("Environment::verifyBox success\n");
}







void Environment::run() {
	printf("Simulation started\n\n");

	time0 = std::chrono::high_resolution_clock::now();

	while (display->checkWindowStatus()) {
		
		engine->deviceMaster();		// Device first, otherwise offloading data always needs the last datapoint!
		engine->hostMaster();

		handleStatus(simulation);
		handleDisplay(simulation);

		if (handleTermination(simulation)) {
			break;
		}
	}
	printf("\n\n\n########################## SIMULATION FINISHED ##########################\n\n\n\n");

	if (simulation->finished || simulation->box->critical_error_encountered) {
		postRunEvents();
	}
}



void Environment::postRunEvents() {
	simulation->out_dir += "\\Steps_" + to_string(simulation->getStep()) + "\\";

#ifndef __linux__
	int check = _mkdir(&(simulation->out_dir[0]));		// Not recognized on linux
#else
	// WEhat is the linux equivalent?? TODO: !!
#endif

	dumpToFile(simulation->logging_data, 10 * simulation->getStep(), simulation->out_dir + "logdata.bin");


	for (int i = 0; i < simulation->getStep(); i++) {
		//simulation->traindata_buffer[0 + i * N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->n_compounds].print();
	}

	if (simulation->box->critical_error_encountered) {
		dumpToFile(simulation->traindata_buffer,
			(uint64_t) N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->n_compounds * simulation->getStep(),
			simulation->out_dir + "sim_traindata.bin");
	}
	

	dumpToFile(simulation->traj_buffer, simulation->getStep() * simulation->total_particles_upperbound, simulation->out_dir + "trajectory.bin");

	printf("temp %f \n", simulation->temperature_buffer[0]);

	Analyzer::AnalyzedPackage analyzed_package = analyzer.analyzeEnergy(simulation);
	dumpToFile(analyzed_package.energy_data, analyzed_package.n_energy_values, simulation->out_dir + "energy.bin");
	dumpToFile(analyzed_package.temperature_data, analyzed_package.n_temperature_values, simulation->out_dir + "temperature.bin");
	printf("temp %f \n", simulation->temperature_buffer[0]);
	printf("temp %f %f\n", analyzed_package.temperature_data[0], analyzed_package.temperature_data[analyzed_package.n_temperature_values - 1]);

	dumpToFile(simulation->potE_buffer, simulation->getStep() * simulation->total_particles_upperbound, simulation->out_dir + "potE.bin");

#ifndef __linux__
	if (!simulation->box->critical_error_encountered && 0) {	// Skipping for now
		string data_processing_command = "C:\\Users\\Daniel\\git_repo\\Quantom\\LIMA_services\\x64\\Release\\LIMA_services.exe "
			+ simulation->out_dir + " "
			+ to_string(simulation->getStep()) + " "
			+ "0" + " "											// do_shuffle
			+ to_string(simulation->n_compounds) + " "
			+ to_string(MAX_COMPOUND_PARTICLES)
			;

		cout << data_processing_command << "\n\n";
		system(&data_processing_command[0]);
	}

	
#endif
}



void Environment::handleStatus(Simulation* simulation) {
	if (!(simulation->getStep() % simulation->steps_per_render)) {
		printf("\r\tStep #%06d", simulation->box->step);
		double duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time0).count();
		int remaining_minutes = (int)(1.f / 1000 * duration / simulation->steps_per_render * (simulation->n_steps - simulation->box->step) / 60);
		printf("\tAvg. step time: %.2fms (%05d/%05d/%05d) \tRemaining: %04d min", duration / simulation->steps_per_render, engine->timings.x / simulation->steps_per_render, engine->timings.y / simulation->steps_per_render, engine->timings.z/simulation->steps_per_render, remaining_minutes);
		engine->timings = Int3(0, 0, 0);


		time0 = std::chrono::high_resolution_clock::now();
	}
}



void Environment::handleDisplay(Simulation* simulation) {	
	if (!(simulation->getStep() % simulation->steps_per_render)) {
		display->render(simulation);
	}
}

bool Environment::handleTermination(Simulation* simulation)
{
	if (simulation->finished)
		return true;
	if (simulation->getStep() >= simulation->n_steps) {
		simulation->finished = true;
		return true;
	}		
	if (simulation->box->critical_error_encountered)
		return true;

	return false;
}

void Environment::prepFF(string conf_path, string topol_path) {
	string program_command = "C:\\Users\\Daniel\\git_repo\\Quantom\\LIMA_ForcefieldMaker\\Release\\LIMA_ForcefieldMaker.exe "
		+ (string) "prepsim" + " "
		+ conf_path + " "
		+ topol_path + " "
		;

	cout << program_command << "\n\n";
	system(&program_command[0]);
}

void Environment::renderTrajectory(string trj_path)
{
	/*
	Trajectory* trj = new Trajectory(trj_path);
	for (int i = 0; i < trj->n_particles; i++) {
		trj->particle_type[i] = 0;
	}
	trj->particle_type[0] = 1;

	display->animate(trj);
	*/
}

void Environment::makeVirtualTrajectory(string trj_path, string waterforce_path) {
	Trajectory* trj = new Trajectory(trj_path);
	Trajectory* force_buffer = new Trajectory(waterforce_path);
	int n_steps = trj->n_steps;

	printf(" part: %d\n", trj->n_particles);


	Float3* particle_position = new Float3[n_steps];
	for (int step = 0; step < n_steps; step++)
		particle_position[step] = trj->positions[0 + step * trj->n_particles];
	Float3* forces = force_buffer->positions;
	

	VirtualPathMaker VPM;
	Float3* vp_path = VPM.makeVirtualPath(particle_position, forces, n_steps);

	std::ofstream myfile("D:\\Quantom\\virtrj.csv");
	for (int step = 0; step < n_steps; step++) {

		for (int k = 0; k < 3; k++) {
			myfile << particle_position[step].at(k) << ";";
		}
		for (int k = 0; k < 3; k++) {
			myfile << vp_path[step].at(k) << ";";
		}

		myfile << "\n";
	}
	myfile.close();

}





template <typename T>
void Environment::dumpToFile(T* data, uint64_t n_datapoints, string file_path_s) {	
	char* file_path;
	file_path = &file_path_s[0];
	printf("Writing %.03Lf MB to binary file ", (long double) sizeof(T) * n_datapoints * 1e-6);
	cout << file_path << endl;

	FILE* file;
#ifndef __linux__
	fopen_s(&file, file_path, "wb");
#else
	file = fopen(file_path, "wb");
#endif

	std::printf("Check %d %lld\n", static_cast<int>(sizeof(T)), n_datapoints);

	fwrite(data, sizeof(T), n_datapoints, file);
	fclose(file);
}
