#include "LIMA_MD/include/Environment.h"

#include "LIMA_BASE/include/Printer.h"

#include <filesystem>

using namespace LIMA_Print;
using std::string;
using std::cout;
using std::printf;

Environment::Environment() {
	sayHello();
	display = new Display();
}

void Environment::CreateSimulation(string gro_path, string topol_path, string work_folder) {
	simulation = std::make_unique<Simulation>(sim_params);
	boxbuilder = new BoxBuilder(work_folder);

	verifySimulationParameters();

	this->work_folder = work_folder;

	prepFF(gro_path, topol_path);									// TODO: Make check in here whether we can skip this!
	forcefield.loadForcefield(work_folder + "/molecule");

	CompoundCollection collection = LIMA_MOLECULEBUILD::buildMolecules(&forcefield, work_folder, SILENT, gro_path, topol_path, true);

	boxbuilder->buildBox(simulation.get());
	//boxbuilder->addCompoundCollection(simulation.get(), &mol_6lzm_10);
	boxbuilder->addCompoundCollection(simulation.get(), collection);

#ifdef ENABLE_SOLVENTS
	boxbuilder->solvateBox(simulation.get(), collection.solvent_positions);
#endif
}

void Environment::CreateSimulation(const Box& boxorigin) {

}

void Environment::verifySimulationParameters() {	// Not yet implemented
	static_assert(THREADS_PER_COMPOUNDBLOCK >= MAX_COMPOUND_PARTICLES, "Illegal kernel parameter");
	static_assert(THREADS_PER_SOLVENTBLOCK >= THREADS_PER_COMPOUNDBLOCK, "Illegal kernel parameter");
	//assert(THREADS_PER_SOLVENTBLOCK >= N_SOLVATE_MOLECULES);
	static_assert(BOX_LEN > 3.f, "Box too small");
	//assert(BOX_LEN >= CUTOFF + 0.5f);
	//assert(simulation->n_compounds <= 1);	// Otherwise data_GAN goes haywire

	//assert(simulation->n_steps % STEPS_PER_LOGTRANSFER == 0);
	//assert(simulation->n_steps % STEPS_PER_THERMOSTAT == 0);
	//assert(simulation->n_steps % STEPS_PER_TRAINDATATRANSFER == 0);

	static_assert(STEPS_PER_THERMOSTAT % STEPS_PER_LOGTRANSFER == 0);		// Change to trajtransfer later
	//assert(STEPS_PER_THERMOSTAT >= STEPS_PER_LOGTRANSFER);
	static_assert(THREADS_PER_SOLVENTBLOCK >= MAX_COMPOUND_PARTICLES);

	static_assert(STEPS_PER_THERMOSTAT >= STEPS_PER_LOGTRANSFER);
	
	//auto a = std::roundf(std::abs(BOX_LEN / SolventBlockGrid::node_len)) * SolventBlockGrid::node_len;// -BOX_LEN_NM;

	auto a = static_cast<int>(static_cast<double>(_BOX_LEN_PM) * 1000);
	// Assert that boxlen is a multiple of nodelen
	assert((static_cast<int>(static_cast<double>(_BOX_LEN_PM)*1000) % BOXGRID_NODE_LEN_pico) == 0 && "BOXLEN must be a multiple of nodelen (1.2 nm)");

	printf("Simulation parameters verified\n");
}

void Environment::verifyBox() {
	for (int c = 0; c < simulation->n_compounds; c++) {
		//printf("Compound radius: %f\t center: %f %f %f\n", simulation->compounds_host[c].confining_particle_sphere, simulation->compounds_host[c].center_of_mass.x, simulation->compounds_host[c].center_of_mass.y, simulation->compounds_host[c].center_of_mass.z);
		if ((simulation->compounds_host[c].radius * 1.1) > BOX_LEN_HALF) {
			printf("Compound %d too large for simulation-box\n", c);
			exit(1);
		}
	}


	//if (std::abs(SOLVENT_MASS - engine->getForcefield().particle_parameters[0].mass) > 1e-3f) {
	if (std::abs(SOLVENT_MASS - forcefield.getNBForcefield().particle_parameters[0].mass) > 1e-3f) {
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

bool Environment::prepareForRun() {
	if (simulation->finished) { 
		printf("Cannot prepare run, since simulation has already finished");
		return false; 
	}

	boxbuilder->finishBox(simulation.get(), forcefield.getNBForcefield());

	simulation->moveToDevice();	// Only moves the Box to the device
	verifyBox();
	simulation->ready_to_run = true;

	engine = std::make_unique<Engine>(simulation.get(), forcefield.getNBForcefield());

	

	

	return simulation->ready_to_run;
}


void Environment::sayHello() {
	std::ifstream file(main_dir+"resources/logo_ascii.txt");
	if (!file) {
		throw std::runtime_error("Failed to open logo file");
	}

	std::string file_contents((std::istreambuf_iterator<char>(file)),
		std::istreambuf_iterator<char>());

	cout << "\033[1;32m"; // set text color to green
	cout << file_contents;
	printf(" \t\t<< Welcome to LIMA Molecular Dynamics >>\n\n");
	cout << "\033[0m"; // reset text color
}


void Environment::run(bool em_variant) {
	if (!prepareForRun()) { return; }

	printH1("Simulation started", true, false);

	time0 = std::chrono::high_resolution_clock::now();

	while (display->checkWindowStatus()) {

		if (em_variant)
			engine->step<true>();
		else
			engine->step<false>();

		handleStatus(simulation.get());
		handleDisplay(simulation.get());

		// Deadspin to slow down rendering for visual debugging :)
		while ((double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time0).count() < FORCED_INTERRENDER_TIME) {}

		if (handleTermination(simulation.get())) {
			break;
		}

	}

	simulation->finished = true;
	simulation->ready_to_run = false;

	engine->terminateSimulation();
	printH1("Simulation Finished", true, true);

	if (simulation->finished || simulation->simparams_device->critical_error_encountered) {
		postRunEvents();
	}
}

void Environment::postRunEvents() {
	
	const std::string out_dir = work_folder + "/Steps_" + std::to_string(simulation->getStep()) + "/";

	std::filesystem::current_path(work_folder);
	std::filesystem::create_directories(out_dir);

	// Nice to have for matlab stuff	
	printH2();
	LIMA_Printer::printNameValuePairs("n steps", static_cast<int>(simulation->getStep()), "n solvents", simulation->n_solvents, "max comp particles", MAX_COMPOUND_PARTICLES, "n compounds", simulation->n_compounds);
	printH2();
	
	if (0) {
		dumpToFile(simulation->logging_data, 10 * simulation->getStep(), out_dir + "logdata.bin");
	}

	if (simulation->simparams_device->critical_error_encountered) {
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

#ifdef USEDEBUGF3
	dumpToFile(simulation->box->debugdataf3, simulation->getStep() * simulation->total_particles_upperbound * DEBUGDATAF3_NVARS, out_dir + "debugf3.bin");
#endif 

#ifndef __linux__
	if (!simulation->simparams_device->critical_error_encountered && 0) {	// Skipping for now
		string data_processing_command = "C:\\Users\\Daniel\\git_repo\\Quantom\\LIMA_services\\x64\\Release\\LIMA_services.exe "
			+ out_dir + " "
			+ std::to_string(simulation->getStep()) + " "
			+ "0" + " "											// do_shuffle
			+ std::to_string(simulation->n_compounds) + " "
			+ std::to_string(MAX_COMPOUND_PARTICLES)
			;

		cout << data_processing_command << "\n\n";
		system(&data_processing_command[0]);
	}
#endif

	// KILL simulation
	//simulation.release();
	//engine.release();
	simulation->ready_to_run = false;

	printH2("Post-run events finished Finished", true, true);
	printf("\n\n\n");
}



void Environment::handleStatus(Simulation* simulation) {
	if (!(simulation->getStep() % simulation->steps_per_render)) {
		printf("\r\tStep #%06d", simulation->getStep());
		double duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time0).count();
		int remaining_minutes = (int)(1.f / 1000 * duration / simulation->steps_per_render * (simulation->simparams_device->constparams.n_steps - simulation->simparams_host.step) / 60);
		printf("\tAvg. step time: %.2fms (%05d/%05d/%05d) \tRemaining: %04d min", duration / simulation->steps_per_render, engine->timings.x / simulation->steps_per_render, engine->timings.y / simulation->steps_per_render, engine->timings.z/simulation->steps_per_render, remaining_minutes);
		//engine->timings = Int3(0, 0, 0);
		engine->timings.x = 0;
		engine->timings.y = 0;

		
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
	if (simulation->getStep() >= simulation->simparams_device->constparams.n_steps) {
		simulation->finished = true;
		return true;
	}		
	if (simulation->simparams_device->critical_error_encountered)
		return true;

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




// Todo: move this to the utilities.h file
template <typename T>
void Environment::dumpToFile(T* data, uint64_t n_datapoints, string file_path_s) {	
	char* file_path;
	file_path = &file_path_s[0];
	printf("Writing %.03Lf MB to binary file ", (long double) sizeof(T) * n_datapoints * 1e-6);
	cout << file_path << '\n';

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

InputSimParams* Environment::getSimparamRef() {
	return &sim_params;
}

Simulation* Environment::getSim() {
	return simulation.get();
}

Analyzer::AnalyzedPackage* Environment::getAnalyzedPackage()
{
	return &postsim_anal_package;
}

//std::array<CompoundCoords, MAX_COMPOUNDS>& Environment::getCoordarrayRef(std::string selector)
//{
//	//if (selector == "current") return boxbuilder->coordarray;
//	//if (selector == "prev") return boxbuilder->coordarray_prev;
//	assert(false);
//	return std::array<CompoundCoords, MAX_COMPOUNDS>{};
//}

SolventBlockGrid* Environment::getAllSolventBlocksPrev()
{
	assert(!simulation->ready_to_run);	// Only valid before simulation is locked
	return boxbuilder->solventblocks_prev;
}

std::unique_ptr<SolventBlockGrid> Environment::getCurrentSolventblockGrid()
{
	auto gridptr_device = CoordArrayQueueHelpers::getSolventblockGridPtr(simulation->box->solventblockgrid_circular_queue, simulation->getStep());

	auto grid_host = std::make_unique<SolventBlockGrid>();
	cudaMemcpy(grid_host.get(), gridptr_device, sizeof(SolventBlockGrid), cudaMemcpyDeviceToHost);

	return grid_host;
}
