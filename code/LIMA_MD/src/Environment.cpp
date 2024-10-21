#include <chrono>
#include <filesystem>
#include <string>
#include <assert.h>  


#include "Environment.h"
#include "Printer.h"
#include "MDFiles.h"
#include "CompoundBuilder.h"
#include "VirtualPathMaker.h"
#include "Display.h"
#include "BoxBuilder.cuh"
#include "Engine.cuh"
#include "Forcefield.h"

using namespace LIMA_Print;

namespace lfs = FileUtils;
namespace fs = std::filesystem;


// ------------------------------------------------ Display Parameters ------------------------------------------ //
const int STEPS_PER_UPDATE = 50;
constexpr float MIN_STEP_TIME = 0.f;		// [ms] Set to 0 for full speed sim
// -------------------------------------------------------------------------------------------------------------- //

Environment::Environment(const fs::path& workdir, EnvMode mode)
	: work_dir(workdir)
	, m_mode(mode)
	, m_logger{ LimaLogger::compact, m_mode, "environment", workdir.string()}	// .string() is temp
{
	switch (mode)
	{
	case EnvMode::Full:
		[[fallthrough]];
	case EnvMode::ConsoleOnly:
		sayHello();
		[[fallthrough]];
	case EnvMode::Headless:
		break;
	}
}

Environment::~Environment() {}

void Environment::CreateSimulation(float boxsize_nm) {
	SimParams simparams{};
	simulation = std::make_unique<Simulation>(simparams, std::make_unique<Box>(boxsize_nm));
	simulation->box_host->boxparams.boxSize = static_cast<int>(boxsize_nm);
}

void Environment::CreateSimulation(std::string gro_path, std::string topol_path, const SimParams params) {
	const auto groFile = std::make_unique<GroFile>(gro_path);
	const auto topFile = std::make_unique<TopologyFile>(topol_path);
	CreateSimulation(*groFile, *topFile, params);
}

void Environment::CreateSimulation(const GroFile& grofile, const TopologyFile& topolfile, const SimParams& params) 
{

	boximage = LIMA_MOLECULEBUILD::buildMolecules(
		grofile,
		topolfile,
		V1,
		std::make_unique<LimaLogger>(LimaLogger::normal, m_mode, "moleculebuilder", work_dir),
		IGNORE_HYDROGEN,
		params
		);

	simulation = std::make_unique<Simulation>(params, BoxBuilder::BuildBox(params, *boximage));
	simulation->forcefield = boximage->forcefield;
}

void Environment::CreateSimulation(Simulation& simulation_src, const SimParams params) {

	simulation.reset(new Simulation(params));
	BoxBuilder::copyBoxState(*simulation, std::move(simulation_src.box_host), simulation_src.getStep());

	simulation->forcefield = simulation_src.forcefield;
}


void Environment::createSimulationFiles(float boxlen) {
	GroFile grofile{};
	grofile.m_path = work_dir / "molecule/conf.gro";
	grofile.box_size = Float3{ boxlen };
	grofile.printToFile();

	TopologyFile topfile{};
	topfile.path = work_dir / "molecule/topol.top";
	topfile.printToFile();

	SimParams simparams{};
	simparams.dumpToFile(work_dir / "sim_params.txt");
}

void constexpr Environment::verifySimulationParameters() {	// Not yet implemented
	static_assert(THREADS_PER_COMPOUNDBLOCK >= MAX_COMPOUND_PARTICLES, "Illegal kernel parameter");
}

void Environment::verifyBox() {
	for (int c = 0; c < simulation->box_host->boxparams.n_compounds; c++) {
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

	

	if (simulation->simparams_host.bc_select == NoBC && simulation->box_host->boxparams.n_solvents != 0) {
		throw std::runtime_error("A simulation with no Boundary Condition may not contain solvents, since they may try to acess a solventblock outside the box causing a crash");
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
	if (simulation == nullptr)// TEMP, ENv should never give sim to engine
		return true;

	if (simulation->finished) { 
		printf("Cannot prepare run, since simulation has already finished");
		assert(false);
		return false; 
	}

	m_logger.startSection("Simulation started");
	time0 = std::chrono::steady_clock::now();

	if (simulation->ready_to_run) { return true; }

	simulation->PrepareDataBuffers();
	
	verifyBox();
	simulation->ready_to_run = true;

	// TEMP, this is a bad solution ?? TODO NOW
	this->compounds = simulation->box_host->compounds;

	boxparams = simulation->box_host->boxparams;
	coloringMethod = simulation->simparams_host.coloring_method;


	engine = std::make_unique<Engine>(
		std::move(simulation),
		simulation->simparams_host.bc_select,
		std::make_unique<LimaLogger>(LimaLogger::compact, m_mode, "engine", work_dir));

	return true;
}


void Environment::sayHello() {
	std::ifstream file(FileUtils::GetLimaDir() / "resources/logo/logo_ascii.txt");
	if (!file) {
		throw std::runtime_error("Failed to open logo file");
	}

	std::string file_contents((std::istreambuf_iterator<char>(file)),
		std::istreambuf_iterator<char>());

	std::cout << file_contents;
}
#include <optional>

void Environment::run(bool doPostRunEvents) {
	const bool emVariant = simulation->simparams_host.em_variant;
	if (!prepareForRun()) { return; }

	std::unique_ptr<Display> display = nullptr;

	if (m_mode == Full) {
		display = std::make_unique<Display>();
		display->WaitForDisplayReady();
	}

	simulationTimer.emplace(TimeIt{ "Simulation" });
	while (true) {
		auto stepStartTime = std::chrono::steady_clock::now();

		if (engine->runstatus.simulation_finished) { 
			break; 
		}
		
		engine->step();

		handleStatus(engine->runstatus.current_step, 0);	// TODO fix the 0

		if (!handleDisplay(compounds, boxparams, *display, emVariant)) {
			break; 
		}
		
		// Deadspin to slow down rendering for visual debugging :)
		while ((double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - stepStartTime).count() < MIN_STEP_TIME) {}
	}
	simulationTimer->stop();

	// Transfers the remaining traj data and more
	engine->terminateSimulation();

	simulation = engine->takeBackSim();

	simulation->finished = true;
	simulation->ready_to_run = false ;

	
	m_logger.finishSection("Simulation Finished");

	if (simulation->finished && doPostRunEvents) {
		postRunEvents();
	}
}


void Environment::WriteBoxCoordinatesToFile(GroFile& grofile, std::optional<int64_t> _step) {
	if (boximage->total_compound_particles + boximage->solvent_positions.size() != grofile.atoms.size()) {
		throw std::runtime_error("Number of particles in grofile does not match the number of particles in the simulation");
	}

	const int64_t stepToLoadFrom = _step.value_or(simulation->getStep() - 1);

	for (int i = 0; i < boximage->total_compound_particles; i++) {
		const int cid = boximage->particleinfos[i].compoundId;
		const int pid = boximage->particleinfos[i].localIdInCompound;
		const Float3 new_position = simulation->traj_buffer->GetMostRecentCompoundparticleDatapoint(cid, pid, stepToLoadFrom);
		grofile.atoms[i].position = new_position;
	}

	// Handle solvents 
	const int firstSolventIndex = boximage->total_compound_particles;
	for (int solventId = 0; solventId < simulation->box_host->boxparams.n_solvents; solventId++) {
		const Float3 new_position = simulation->traj_buffer->GetMostRecentSolventparticleDatapointAtIndex(solventId, stepToLoadFrom);
		grofile.atoms[firstSolventIndex + solventId].position = new_position;
	}
}
GroFile Environment::WriteBoxCoordinatesToFile(const std::optional<std::string> filename) {
	GroFile outputfile{ boximage->grofile };

	if (filename.has_value()) {
		outputfile.m_path = work_dir / "molecule" / (filename.value() + ".gro");
	}

	WriteBoxCoordinatesToFile(outputfile);

	return outputfile;
}


void Environment::postRunEvents() {
	if (simulation->getStep() == 0) { return; }

	const fs::path out_dir = work_dir / ("Steps_" + std::to_string(simulation->getStep()));
	std::filesystem::create_directories(out_dir);

	// Nice to have for matlab stuff
	if (m_mode != Headless) {
		//printH2();
		//LIMA_Printer::printNameValuePairs("n steps", static_cast<int>(simulation->getStep()), "n solvents", simulation->box_host->boxparams.n_solvents,
		//	"max comp particles", MAX_COMPOUND_PARTICLES, "n compounds", simulation->box_host->boxparams.n_compounds, "total p upperbound", simulation->box_host->boxparams.total_particles_upperbound);
		//printH2();
	}
	
	if (simulation->simparams_host.save_trajectory)
		simulation->ToTracjectoryFile()->Dump(out_dir / "trajectory.trr");

	/*if (simulation->simparams_host.save_energy)
		FileUtils::dumpToFile(postsim_anal_package.energy_data.data(), postsim_anal_package.energy_data.size(), out_dir.string() + "energy.bin");*/

	simulation->ready_to_run = false;
	m_logger.finishSection("Post-run events finished Finished");
}



void Environment::handleStatus(const int64_t step, const int64_t n_steps) {
	if (m_mode == Headless) {
		return;
	}

	if (step % STEPS_PER_UPDATE == STEPS_PER_UPDATE-1) {

		const double duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - time0).count();

		// First clear the current line
		printf("\r\033[K");

		printf("\rStep #%06llu", step);
		printf("\tAvg. time: %.2fms (%05d/%05d/%05d/%05d/%05d) \tRemaining: %04d min         ", 
			duration / STEPS_PER_UPDATE,
			engine->timings.compound_kernels / STEPS_PER_UPDATE,
			engine->timings.solvent_kernels / STEPS_PER_UPDATE,
			engine->timings.cpu_master/ STEPS_PER_UPDATE,
			engine->timings.nlist/ STEPS_PER_UPDATE,
			engine->timings.electrostatics / STEPS_PER_UPDATE,
			0);

		engine->timings.reset();
		time0 = std::chrono::steady_clock::now();
	}
}



bool Environment::handleDisplay(const std::vector<Compound>& compounds_host, const BoxParams& boxparams, Display& display, bool emVariant) {
	if (m_mode != Full) {
		return true;
	}

	auto displayException = display.displayThreadException;
	if (displayException) {
		std::rethrow_exception(displayException);
	}

	if (engine->runstatus.stepForMostRecentData != step_at_last_render && engine->runstatus.most_recent_positions != nullptr) {

		const std::string info = emVariant
			? std::format("Step {:d} MaxForce {:.02f}", static_cast<int>(engine->runstatus.current_step), static_cast<float>(engine->runstatus.greatestForce))
			: std::format("Step {:d} Temp {:.02f}", static_cast<int>(engine->runstatus.current_step), static_cast<float>(engine->runstatus.current_temperature));

		display.Render(std::make_unique<Rendering::SimulationTask>( 
			engine->runstatus.most_recent_positions, compounds_host, boxparams, info, coloringMethod
		));
		step_at_last_render = engine->runstatus.current_step;
	}

	return !display.DisplaySelfTerminated();
}






void Environment::renderTrajectory(std::string trj_path)
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

void Environment::makeVirtualTrajectory(std::string trj_path, std::string waterforce_path) {
	//Trajectory* trj = new Trajectory(trj_path);
	//Trajectory* force_buffer = new Trajectory(waterforce_path);
	//int n_steps = trj->n_steps;

	//printf(" part: %d\n", trj->n_particles);


	//Float3* particle_position = new Float3[n_steps];
	//for (int64_t step = 0; step < n_steps; step++)
	//	particle_position[step] = trj->positions[0 + step * trj->n_particles];
	//Float3* forces = force_buffer->positions;
	//

	//VirtualPathMaker VPM;
	//Float3* vp_path = VPM.makeVirtualPath(particle_position, forces, n_steps);

	//std::ofstream myfile("D:\\Quantom\\virtrj.csv");
	//for (int64_t step = 0; step < n_steps; step++) {

	//	for (int k = 0; k < 3; k++) {
	//		myfile << particle_position[step].at(k) << ";";
	//	}
	//	for (int k = 0; k < 3; k++) {
	//		myfile << vp_path[step].at(k) << ";";
	//	}

	//	myfile << "\n";
	//}
	//myfile.close();
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

const SimAnalysis::AnalyzedPackage& Environment::getAnalyzedPackage()
{
	if (simulation == nullptr)
		throw std::runtime_error("Env has no simulation");
	// TODO: make some check here that the simulation has finished
	if (!postsim_anal_package.has_value())
		postsim_anal_package = SimAnalysis::analyzeEnergy(simulation.get());
	return postsim_anal_package.value();
}
