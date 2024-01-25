#include "Environment.h"
#include "Printer.h"
#include "MDFiles.h"
#include "SimulationBuilder.h"
#include "EngineUtils.cuh"

#include <filesystem>
#include <stdexcept>

using namespace LIMA_Print;
using std::string;
using std::cout;
using std::printf;
namespace lfs = Filehandler;

Environment::Environment(const string& wf, EnvMode mode, bool save_output)
	: work_dir(wf)
	, m_mode(mode)
	, m_logger{ LimaLogger::compact, m_mode, "environment", wf }
	, save_output(save_output)
{
	switch (mode)
	{
	case EnvMode::Full:
		display = std::make_unique<Display>();
		[[fallthrough]];
	case EnvMode::ConsoleOnly:
		sayHello();
		[[fallthrough]];
	case EnvMode::Headless:
		break;
	}

	boxbuilder = std::make_unique<BoxBuilder>(std::make_unique<LimaLogger>(LimaLogger::normal, m_mode, "boxbuilder", work_dir));
}

void Environment::CreateSimulation(float boxsize_nm) {
	SimParams simparams{};
	setupEmptySimulation(simparams);
	boxbuilder->buildBox(simulation.get(), boxsize_nm);
	simulation->box_host->boxparams.dims = Float3{ boxsize_nm };
}

void Environment::CreateSimulation(string gro_path, string topol_path, const SimParams params) {
	const ParsedGroFile gro_file = MDFiles::loadGroFile(gro_path);
	const std::unique_ptr<ParsedTopologyFile> topol_file = MDFiles::loadTopologyFile(topol_path);

	CreateSimulation(gro_file, *topol_file, params);
}

void Environment::CreateSimulation(const ParsedGroFile& grofile, const ParsedTopologyFile& topolfile, const SimParams& params) 
{
	setupEmptySimulation(params);

	boximage = LIMA_MOLECULEBUILD::buildMolecules(
		lfs::pathJoin(work_dir, "molecule"),
		grofile,
		topolfile,
		V1,
		std::make_unique<LimaLogger>(LimaLogger::normal, m_mode, "moleculebuilder", work_dir),
		IGNORE_HYDROGEN,
		simulation->simparams_host.bc_select
		);

	//TODO Find a better place for this
	simulation->forcefield = std::make_unique<Forcefield>(m_mode == Headless ? SILENT : V1, lfs::pathJoin(work_dir, "molecule"));

	boxbuilder->buildBox(simulation.get(), boximage->box_size);

	boxbuilder->addBoxImage(simulation.get(), *boximage);

#ifdef ENABLE_SOLVENTS
	boxbuilder->solvateBox(simulation.get(), boximage->solvent_positions);
#endif
}

void Environment::CreateSimulation(Simulation& simulation_src, const SimParams params) {
	// If we already have a box, we must have a forcefield too, no?

	//SimParams simparams{ ip };
	setupEmptySimulation(params);
	//boxbuilder->copyBoxState(simulation.get(), simulation_src.sim_dev->box, simulation_src.simparams_host, simulation_src.simparams_host.step);	// TODO: Fix this again
	boxbuilder->copyBoxState(simulation.get(), std::move(simulation_src.box_host), simulation_src.simsignals_host, simulation_src.simsignals_host.step);	// TODO: Fix this again
	simulation->extraparams = simulation_src.extraparams;

	//TODO Find a better place for this
	simulation->forcefield = std::make_unique<Forcefield>(m_mode == Headless ? SILENT : V1, lfs::pathJoin(work_dir, "molecule"));
}

void Environment::createMembrane(bool carryout_em) {
	// Load in the lipid types, for now just POPC
	const std::string POPC_PATH = main_dir + "/resources/Lipids/POPC/";
	const ParsedGroFile inputgrofile = MDFiles::loadGroFile(POPC_PATH + "popc.gro");
	std::unique_ptr<ParsedTopologyFile> inputtopologyfile = MDFiles::loadTopologyFile(POPC_PATH + "POPC.itp");

	// Insert the x lipids with plenty of distance in a non-pbc box
	auto monolayerfiles = SimulationBuilder::buildMembrane({ inputgrofile, *inputtopologyfile }, simulation->box_host->boxparams.dims);

	// Create simulation and run on the newly created files in the workfolder
	monolayerfiles.first.printToFile(lfs::pathJoin(work_dir, "/molecule/monolayer.gro"));
	monolayerfiles.second.printToFile(lfs::pathJoin(work_dir, "/molecule/monolayer.top"));

	// Monolayer energy Minimization NoBC
	{
		SimParams ip{};
		ip.bc_select = NoBC;
		ip.n_steps = carryout_em ? 10000 : 0;
		ip.snf_select = HorizontalSqueeze;
		ip.em_variant = true;
		//CreateSimulation(lfs::pathJoin(work_dir, "/molecule/membrane.gro"), lfs::pathJoin(work_dir, "/molecule/membrane.top"), ip);
		CreateSimulation(monolayerfiles.first, monolayerfiles.second, ip);

		// Draw each lipid towards the center - no pbc
		run();

		if (carryout_em)
			if (!boxbuilder->verifyAllParticlesIsInsideBox(*simulation, 0.06f)) { return; }	// FAIL
	}
	
	const ParsedGroFile monolayer_grofile_em = carryout_em ? writeBoxCoordinatesToFile() : monolayerfiles.first;
	//auto em_monolayer_grofile = MDFiles::loadGroFile(work_dir + "/molecule/membrane.gro");
	

	// Copy each particle, and flip them around the xy plane, so the monolayer becomes a bilayer
	auto bilayerfiles = SimulationBuilder::makeBilayerFromMonolayer({ monolayer_grofile_em, monolayerfiles.second }, simulation->box_host->boxparams.dims);

	// Run EM for a while - with pbc
	{
		SimParams ip{};
		ip.n_steps = carryout_em ? 30000 : 0;
		ip.dt = 50.f;
		ip.bc_select = carryout_em ? PBC : NoBC;	// Cannot insert compounds with PBC, if they are not in box
		CreateSimulation(bilayerfiles.first, bilayerfiles.second, ip);

		// Draw each lipid towards the center - no pbc
		run();
	}
	const ParsedGroFile bilayer_grofile_em = carryout_em ? writeBoxCoordinatesToFile() : bilayerfiles.first;

	
	// Save box to .gro and .top file
	bilayer_grofile_em.printToFile(lfs::pathJoin(work_dir, "/molecule/membrane.gro"));
	//bilayerfiles.first.printToFile(lfs::pathJoin(work_dir, "/molecule/membrane.gro"));	// TEMP
	bilayerfiles.second.printToFile(lfs::pathJoin(work_dir, "/molecule/membrane.top"));
}

void Environment::setupEmptySimulation(const SimParams& simparams) {
	//assert(forcefield.forcefield_loaded && "Forcefield was not loaded before creating simulation!");


	simulation = std::make_unique<Simulation>(simparams, work_dir + "/molecule/", m_mode);

	verifySimulationParameters();
}

void Environment::verifySimulationParameters() {	// Not yet implemented
	static_assert(THREADS_PER_COMPOUNDBLOCK >= MAX_COMPOUND_PARTICLES, "Illegal kernel parameter");
	//static_assert(BOX_LEN > 3.f, "Box too small");
	//static_assert(BOX_LEN > CUTOFF_NM *2.f, "CUTOFF too large relative to BOXLEN");

	static_assert(STEPS_PER_THERMOSTAT % STEPS_PER_LOGTRANSFER == 0);		// Change to trajtransfer later
	static_assert(STEPS_PER_THERMOSTAT >= STEPS_PER_LOGTRANSFER);
	
	//auto a = std::roundf(std::abs(BOX_LEN / SolventBlockGrid::node_len)) * SolventBlockGrid::node_len;// -BOX_LEN_NM;

	auto a = static_cast<int>(static_cast<double>(_BOX_LEN_PM) * 1000);
	// Assert that boxlen is a multiple of nodelen
	//assert((static_cast<int>(static_cast<double>(_BOX_LEN_PM)*1000) % BOXGRID_NODE_LEN_pico) == 0 && "BOXLEN must be a multiple of nodelen (1.2 nm)");
}

void Environment::verifyBox() {
	for (int c = 0; c < simulation->boxparams_host.n_compounds; c++) {
		//printf("Compound radius: %f\t center: %f %f %f\n", simulation->compounds_host[c].confining_particle_sphere, simulation->compounds_host[c].center_of_mass.x, simulation->compounds_host[c].center_of_mass.y, simulation->compounds_host[c].center_of_mass.z);
		/*if ((simulation->compounds_host[c].radius * 1.1) > BOX_LEN_HALF) {
			throw std::runtime_error(std::format("Compound {} too large for simulation-box", c).c_str());
		}*/
		for (int i = 0; i < CompoundInteractionBoundary::k; i++) {
			/*if ((simulation->compounds_host[c].interaction_boundary.radii[i] * 1.1) > BOX_LEN_HALF) {
				throw std::runtime_error(std::format("Compound {} too large for simulation-box", c).c_str());
			}*/
		}
		
	}

	if (simulation->simparams_host.bc_select == NoBC && simulation->boxparams_host.n_solvents != 0) {
		throw std::runtime_error("A simulation with no Boundary Condition may not contain solvents, since they may try to acess a solventblock outside the box causing a crash");
	}

	

	if (std::abs(SOLVENT_MASS - simulation->forcefield->getNBForcefield().particle_parameters[0].mass) > 1e-3f) {
		throw std::runtime_error("Error: Solvent mass is unreasonably large");
	}


#ifdef LIMAKERNELDEBUGMODE
	if (print_compound_positions) {
		for (int c = 0; c < simulation->boxparams_host.n_compounds; c++) {
			Compound* comp = &simulation->compounds_host[c];
			for (int p = 0; p < comp->n_particles; p++) {
				printf("%d   ", comp->particle_global_ids[p]);
			}
		}
	}
#endif
}

bool Environment::prepareForRun() {
	if (simulation->finished) { 
		printf("Cannot prepare run, since simulation has already finished");
		assert(false);
		return false; 
	}

	m_logger.startSection("Simulation started");

	step_at_last_render = 0;
	step_at_last_update = 0;
	time0 = std::chrono::high_resolution_clock::now();

	if (simulation->ready_to_run) { return true; }

	boxbuilder->finishBox(simulation.get());
	//simulation->moveToDevice();	// Only moves the Box to the device


	
	verifyBox();
	simulation->ready_to_run = true;

	// TEMP, this is a bad solution
	compounds = &simulation->compounds_host;
	boxparams = simulation->boxparams_host;

	engine = std::make_unique<Engine>(
		std::move(simulation),
		simulation->simparams_host.bc_select,
		std::make_unique<LimaLogger>(LimaLogger::compact, m_mode, "engine", work_dir));

	return true;
}


void Environment::sayHello() {
	std::ifstream file(main_dir+"/resources/logo_ascii.txt");
	if (!file) {
		throw std::runtime_error("Failed to open logo file");
	}

	std::string file_contents((std::istreambuf_iterator<char>(file)),
		std::istreambuf_iterator<char>());

	//cout << "\033[1;32m"; // set text color to green
	cout << file_contents;
	//printf(" \t\t<< Welcome to LIMA Molecular Dynamics >>\n\n");
	//cout << "\033[0m"; // reset text color
}


void Environment::run(bool em_variant) {
	if (!prepareForRun()) { return; }


	while (true) {
		if (engine->runstatus.simulation_finished) { break; }
		if (em_variant)
			engine->step();
		else
			engine->step();



		handleStatus(engine->runstatus.current_step, 0);	// TODO fix the 0

		if (!handleDisplay(*compounds, boxparams)) { break; }

		// Deadspin to slow down rendering for visual debugging :)
		while ((double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time0).count() < FORCED_INTERRENDER_TIME) {}

	}

	// Transfers the remaining traj data and more
	engine->terminateSimulation();

	simulation = engine->takeBackSim();

	simulation->finished = true;
	simulation->ready_to_run = false;

	
	m_logger.finishSection("Simulation Finished");

	//if (simulation->finished || simulation->sim_dev->params->critical_error_encountered) {
	if (simulation->finished) {
		postRunEvents();
	}
}

ParsedGroFile Environment::writeBoxCoordinatesToFile() {
	ParsedGroFile outputfile{ boximage->grofile };

	if (simulation->boxparams_host.n_solvents != 0) {
		throw std::runtime_error("This function is now designed to handle solvents");
	}

	for (int i = 0; i < outputfile.atoms.size(); i++) {
		const int cid = boximage->particleinfotable[i].compound_index;
		const int pid = boximage->particleinfotable[i].local_id_compound;

		int index = LIMALOGSYSTEM::getMostRecentDataentryIndex(simulation->simsignals_host.step-1);
		const Float3 new_position = simulation->traj_buffer->getCompoundparticleDatapointAtIndex(cid, pid, index);
		
		outputfile.atoms[i].position = new_position;
	}

	// Handle solvents somehow

	//outputfile.printToFile(work_dir + "/molecule/" + filename + ".gro");
	return outputfile;
}


void Environment::postRunEvents() {
	if (simulation->getStep() == 0) { return; }

	if (POSTSIM_ANAL) {
		Analyzer analyzer(std::make_unique<LimaLogger>(LimaLogger::compact, m_mode, "analyzer", work_dir));
		postsim_anal_package = analyzer.analyzeEnergy(simulation.get());
	}

	if (!save_output) { return; }


	const std::string out_dir = work_dir + "Steps_" + std::to_string(simulation->getStep()) + "/";

	const std::filesystem::path out_path{ out_dir };
	std::filesystem::create_directories(out_path);
	//std::filesystem::current_path(work_folder);
	//std::filesystem::create_directories(out_dir);

	// Nice to have for matlab stuff
	if (m_mode != Headless) {
		printH2();
		LIMA_Printer::printNameValuePairs("n steps", static_cast<int>(simulation->getStep()), "n solvents", simulation->boxparams_host.n_solvents, 
			"max comp particles", MAX_COMPOUND_PARTICLES, "n compounds", simulation->boxparams_host.n_compounds, "total p upperbound", simulation->boxparams_host.total_particles_upperbound);
		printH2();
	}

	//if (simulation->sim_dev->params->critical_error_encountered) {
	//	Filehandler::dumpToFile(simulation->trainingdata.data(),
	//		(uint64_t) N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->boxparams_host.n_compounds * simulation->getStep(),
	//		out_dir + "sim_traindata.bin");
	//}
	
	if (DUMP_TRAJ) {
		//dumpToFile(simulation->traj_buffer->data(), simulation->getStep() * simulation->total_particles_upperbound, out_dir + "trajectory.bin");
		MDFiles::TrrFile::dumpToFile(simulation.get(), out_dir + "trajectory.trr");
	}

	if (POSTSIM_ANAL) {
		Filehandler::dumpToFile(
			postsim_anal_package.energy_data.data(),
			postsim_anal_package.energy_data.size(),
			out_dir + "energy.bin"
		);
	}

	if (DUMP_POTE) {
		Filehandler::dumpToFile(simulation->potE_buffer->getBufferAtIndex(0), simulation->getStep() * simulation->boxparams_host.total_particles_upperbound, out_dir + "potE.bin");
	}

#ifdef USEDEBUGF3
	dumpToFile(simulation->box->debugdataf3, simulation->getStep() * simulation->total_particles_upperbound * DEBUGDATAF3_NVARS, out_dir + "debugf3.bin");
#endif 

#ifndef __linux__
	//if (!simulation->sim_dev->params->critical_error_encountered && 0) {	// Skipping for now
	//	string data_processing_command = "C:\\Users\\Daniel\\git_repo\\Quantom\\LIMA_services\\x64\\Release\\LIMA_services.exe "
	//		+ out_dir + " "
	//		+ std::to_string(simulation->getStep()) + " "
	//		+ "0" + " "											// do_shuffle
	//		+ std::to_string(simulation->boxparams_host.n_compounds) + " "
	//		+ std::to_string(MAX_COMPOUND_PARTICLES)
	//		;

	//	cout << data_processing_command << "\n\n";
	//	system(&data_processing_command[0]);
	//}
#endif

	simulation->ready_to_run = false;

	m_logger.finishSection("Post-run events finished Finished");
}



void Environment::handleStatus(const int64_t step, const int64_t n_steps) {
	if (m_mode == Headless) {
		return;
	}

	const int steps_since_update = step - step_at_last_update;
	if (steps_since_update >= STEPS_PER_RENDER) {

		const double duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time0).count();
		const int remaining_minutes = (int)(1.f / 1000 * duration / steps_since_update * (n_steps - step) / 60);

		printf("\r\tStep #%06llu", step);
		printf("\tAvg. step time: %.2fms (%05d/%05d/%05d/%05d) \tRemaining: %04d min", 
			duration / steps_since_update,
			engine->timings.compound_kernels / steps_since_update,
			engine->timings.solvent_kernels / steps_since_update,
			engine->timings.cpu_master/ steps_since_update,
			engine->timings.nlist/ steps_since_update,
			remaining_minutes);

		step_at_last_update = step;
		engine->timings.reset();
		time0 = std::chrono::high_resolution_clock::now();
	}
}



bool Environment::handleDisplay(const std::vector<Compound>& compounds_host, const BoxParams& boxparams) {	
	if (!display) { return true; }	// Headless or ConsoleOnly

	if (engine->runstatus.current_step - step_at_last_render > STEPS_PER_RENDER && engine->runstatus.most_recent_positions != nullptr) {
		
		display->render(engine->runstatus.most_recent_positions, compounds_host, boxparams, engine->runstatus.current_step, engine->runstatus.current_temperature);
		step_at_last_render = engine->runstatus.current_step;
	}

	const bool displayStillExists = display->checkWindowStatus();
	return displayStillExists;
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



//
//// Todo: move this to the utilities.h file
//template <typename T>
//void Environment::dumpToFile(T* data, uint64_t n_datapoints, string file_path_s) {	
//	char* file_path;
//	file_path = &file_path_s[0];
//
//	const std::string str = std::to_string((long double)sizeof(T) * n_datapoints * 1e-6);
//	m_logger.print("Writing " + str + "MB to binary file " + file_path + "\n");
//
//	FILE* file;
//
//#ifndef __linux__
//	if (!fopen_s(&file, file_path, "wb")) {
//
//		assert(sizeof(T));
//		assert(n_datapoints);
//
//		fwrite(data, sizeof(T), n_datapoints, file);
//		fclose(file);
//	}
//#else
//	file = fopen(file_path, "wb");
//#endif
//}



SimParams Environment::loadSimParams(const std::string& path) {
	SimParams simparams{};
	auto param_dict = Filehandler::parseINIFile(path);
	simparams.overloadParams(param_dict);
	return simparams;
}

std::unique_ptr<Simulation> Environment::getSim() {
	// Should we delete the forcefield here?
	//boxbuilder.reset();
	engine.reset();
	return std::move(simulation);
}

Simulation* Environment::getSimPtr() {
	if (simulation) { 
		return simulation.get(); 
	}
	return nullptr;
}

Analyzer::AnalyzedPackage* Environment::getAnalyzedPackage()
{
	return &postsim_anal_package;
}
//
//SolventBlockGrid* Environment::getSolventBlocksPrevRef()
//{
//	if (simulation->sim_dev != nullptr) { throw "Can't ref solventblocks when box is on device"; }
//	auto grid = CoordArrayQueueHelpers::getSolventblockGridPtr(simulation->box_host->solventblockgrid_circular_queue, SolventBlockGrid::first_step_prev);
//	return grid;
//}
//
//const std::unique_ptr<SolventBlockGrid> Environment::getCurrentSolventblockGrid()
//{
//	auto gridptr_device = CoordArrayQueueHelpers::getSolventblockGridPtr(simulation->sim_dev->box->solventblockgrid_circular_queue, simulation->getStep());
//
//	auto grid_host = std::make_unique<SolventBlockGrid>();
//	cudaMemcpy(grid_host.get(), gridptr_device, sizeof(SolventBlockGrid), cudaMemcpyDeviceToHost);
//
//	return grid_host;
//}

void Environment::resetEnvironment() {

}
