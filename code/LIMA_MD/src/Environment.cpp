#include "Environment.h"

#include "Printer.h"

#include <filesystem>

using namespace LIMA_Print;
using std::string;

Environment::Environment() {
	display = new DisplayV2();
}

void Environment::CreateSimulation(string conf_path, string topol_path, string work_folder) {
	simulation = std::make_unique<Simulation>(sim_params);

	verifySimulationParameters();

	this->work_folder = work_folder;

	prepFF(conf_path, topol_path);									// TODO: Make check in here whether we can skip this!
	forcefield.loadForcefield(work_folder + "/molecule");

	CompoundBuilder compoundbuilder(&forcefield, V1);
	CompoundCollection mol_6lzm_10 = compoundbuilder.buildCompoundCollection(conf_path, topol_path, 5);

	boxbuilder.buildBox(simulation.get());
	boxbuilder.addCompoundCollection(simulation.get(), &mol_6lzm_10);

#ifdef ENABLE_SOLVENTS
	//boxbuilder.solvateBox(simulation);
	vector<Float3> solvent_positions = compoundbuilder.getSolventPositions(conf_path);
	boxbuilder.solvateBox(simulation.get(), &solvent_positions);
#endif

	//delete[] mol_6lzm_10.compounds;
}


void Environment::verifySimulationParameters() {	// Not yet implemented
	static_assert(THREADS_PER_COMPOUNDBLOCK >= MAX_COMPOUND_PARTICLES, "Illegal kernel parameter");
	static_assert(THREADS_PER_SOLVENTBLOCK >= THREADS_PER_COMPOUNDBLOCK, "Illegal kernel parameter");
	//assert(THREADS_PER_SOLVENTBLOCK >= N_SOLVATE_MOLECULES);
	static_assert(BOX_LEN > 3.f, "Box too small");
	//assert(BOX_LEN >= CUTOFF + 0.5f);
	//assert(simulation->n_compounds <= 1);	// Otherwise data_GAN goes haywire

	assert(simulation->n_steps % STEPS_PER_LOGTRANSFER == 0);
	//assert(simulation->n_steps % STEPS_PER_THERMOSTAT == 0);
	//assert(simulation->n_steps % STEPS_PER_TRAINDATATRANSFER == 0);

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


	if (std::abs(SOLVENT_MASS - engine->getForcefield().particle_parameters[0].mass) > 1e-3f) {
		printf("Error in solvent mass");
		exit(0);
	}


	if (print_compound_positions) {
		for (int c = 0; c < simulation->n_compounds; c++) {
			Compound* comp = &simulation->compounds_host[c];
			for (int p = 0; p < comp->n_particles; p++) {
				printf("%d   ", comp->particle_global_ids[p]);
			}
		}
	}

	printf("Environment::verifyBox success\n");
}

void Environment::prepareForRun() {
	boxbuilder.finishBox(simulation.get(), forcefield.getNBForcefield());

	simulation->moveToDevice();	// Only moves the Box to the device
	engine = std::make_unique<Engine>(simulation.get(), forcefield.getNBForcefield());

	verifyBox();
	ready_to_run = true;
}





void Environment::run() {
	if (!ready_to_run) { prepareForRun(); }
	printH1("Simulation started", true, false);

	time0 = std::chrono::high_resolution_clock::now();

	while (display->checkWindowStatus()) {

		engine->deviceMaster();		// Device first, otherwise offloading data always needs the last datapoint!
		engine->hostMaster();

		handleStatus(simulation.get());
		handleDisplay(simulation.get());

		if (handleTermination(simulation.get())) {
			break;
		}

	}
	engine->terminateSimulation();
	printH1("Simulation Finished", true, true);

	if (simulation->finished || simulation->box->critical_error_encountered) {
		postRunEvents();
	}
}



void Environment::postRunEvents() {
	
	const std::string out_dir = work_folder + "/Steps_" + to_string(simulation->getStep()) + "/";

	std::filesystem::current_path(work_folder);
	std::filesystem::create_directories(out_dir);

	// Nice to have for matlab stuff	
	printH2();
	LIMA_Printer::printNameValuePairs("n steps", static_cast<int>(simulation->getStep()), "n solvents", simulation->n_solvents, "max comp particles", MAX_COMPOUND_PARTICLES, "n compounds", simulation->n_compounds);
	printH2();
	
	if (0) {
		dumpToFile(simulation->logging_data, 10 * simulation->getStep(), out_dir + "logdata.bin");
	}

	if (simulation->box->critical_error_encountered) {
		dumpToFile(simulation->traindata_buffer,
			(uint64_t) N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->n_compounds * simulation->getStep(),
			out_dir + "sim_traindata.bin");
	}
	
	if (DUMP_TRAJ) {
		dumpToFile(simulation->traj_buffer, simulation->getStep() * simulation->total_particles_upperbound, out_dir + "trajectory.bin");
	}

	if (POSTSIM_ANAL) {
		postsim_anal_package = analyzer.analyzeEnergy(simulation.get());
		dumpToFile(
			postsim_anal_package.energy_data.data(),
			postsim_anal_package.energy_data.size(),
			out_dir + "energy.bin"
		);
		//dumpToFile(analyzed_package.temperature_data, analyzed_package.n_temperature_values, out_dir + "temperature.bin");
	}

	if (DUMP_POTE) {
		dumpToFile(simulation->potE_buffer, simulation->getStep() * simulation->total_particles_upperbound, out_dir + "potE.bin");
	}

#ifndef __linux__
	if (!simulation->box->critical_error_encountered && 0) {	// Skipping for now
		string data_processing_command = "C:\\Users\\Daniel\\git_repo\\Quantom\\LIMA_services\\x64\\Release\\LIMA_services.exe "
			+ out_dir + " "
			+ to_string(simulation->getStep()) + " "
			+ "0" + " "											// do_shuffle
			+ to_string(simulation->n_compounds) + " "
			+ to_string(MAX_COMPOUND_PARTICLES)
			;

		cout << data_processing_command << "\n\n";
		system(&data_processing_command[0]);
	}
#endif

	// KILL simulation
	simulation.release();
	engine.release();
	ready_to_run = false;

	printH2("Post-run events finished Finished", true, true);
	printf("\n\n\n");
}



void Environment::handleStatus(Simulation* simulation) {
	if (!(simulation->getStep() % simulation->steps_per_render)) {
		printf("\r\tStep #%06d", simulation->box->step);
		double duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time0).count();
		int remaining_minutes = (int)(1.f / 1000 * duration / simulation->steps_per_render * (simulation->n_steps - simulation->box->step) / 60);
		printf("\tAvg. step time: %.2fms (%05d/%05d/%05d) \tRemaining: %04d min", duration / simulation->steps_per_render, engine->timings.x / simulation->steps_per_render, engine->timings.y / simulation->steps_per_render, engine->timings.z/simulation->steps_per_render, remaining_minutes);
		engine->timings = Int3(0, 0, 0);


		// Deadspin to slow down rendering for visual debugging :)
		while (simulation->getStep() >= FIRST_INTERRENDER_WAIT && (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time0).count() < FORCED_INTERRENDER_TIME) {}

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
	//if (simulation->box->critical_error_encountered)
		//return true;

	return false;
}

void Environment::prepFF(string conf_path, string topol_path) {
	ForcefieldMaker FFM(work_folder);	// Not to confuse with the engine FFM!!!!=!?!
	FFM.prepSimulationForcefield();
	//const std::string conf_name = "conf.gro";;
	//const std::string topol_name = "topol.top";
	//string program_command = "C:\\Users\\Daniel\\git_repo\\Quantom\\LIMA_ForcefieldMaker\\Release\\LIMA_ForcefieldMaker.exe "
	//	+ (string) "prepsim" + " "
	//	+ conf_name + " "
	//	+ topol_name + " "
	//	;

	//cout << program_command << "\n\n";
	//system(&program_command[0]);
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
	if (!fopen_s(&file, file_path, "wb")) {

		std::printf("Check %d %lld\n", static_cast<int>(sizeof(T)), n_datapoints);

		fwrite(data, sizeof(T), n_datapoints, file);
		fclose(file);
	}
#else
	file = fopen(file_path, "wb");
#endif
}



void Environment::loadSimParams(const std::string& path) {	
	auto param_dict = Filehandler::parseINIFile(path);
	sim_params.overloadParams(param_dict);
}

SimulationParams* Environment::getSimparamRef() {
	return &sim_params;
}

Simulation* Environment::getSim() {
	return simulation.get();
}

Analyzer::AnalyzedPackage* Environment::getAnalyzedPackage()
{
	return &postsim_anal_package;
}

CompoundCoords* Environment::getCoordarrayPtr(std::string selector)
{
	if (selector == "current") return boxbuilder.coordarray;
	if (selector == "prev") return boxbuilder.coordarray_prev;
	assert(false);
	return NULL;
}
