#include <chrono>
#include <filesystem>
#include <string>
#include <optional>

#include "Environment.h"
#include "MDFiles.h"
#include "CompoundBuilder.h"
#include "Display.h"
#include "BoxBuilder.cuh"
#include "Engine.cuh"

namespace lfs = FileUtils;
namespace fs = std::filesystem;


// ------------------------------------------------ Display Parameters ------------------------------------------ //
const int STEPS_PER_UPDATE = 100;
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
	simulation = std::make_unique<Simulation>(simparams, std::make_unique<Box>(Float3(boxsize_nm)));
	simulation->box_host->boxparams.boxSize = static_cast<int>(boxsize_nm);
}

void Environment::CreateSimulation(const std::string& gro_path, const std::string& topol_path, const SimParams& params) {
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
	simulation->forcefieldTinymol = boximage->tinymolTypes;
	simulation->forcefieldTest = boximage->nonbondedInteractionParams;
}

void Environment::CreateSimulation(Simulation& simulation_src, const SimParams params) {

	simulation.reset(new Simulation(params));
	BoxBuilder::copyBoxState(*simulation, std::move(simulation_src.box_host), simulation_src.getStep());

	simulation->forcefield = simulation_src.forcefield;
	simulation->forcefieldTinymol = simulation_src.forcefieldTinymol;
	simulation->forcefieldTest = simulation_src.forcefieldTest;
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

	if (simulation->ready_to_run) { return true; }

	simulation->PrepareDataBuffers();
	
	verifyBox();
	simulation->ready_to_run = true;

	avgStepTimes.reserve((simulation->simparams_host.n_steps + 1) / STEPS_PER_UPDATE);

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

std::chrono::duration<double> Environment::run() {
	const bool emVariant = simulation->simparams_host.em_variant;
    if (!prepareForRun()) { return {}; }

	std::unique_ptr<Display> display = nullptr;

	if (m_mode == Full) {
		display = std::make_unique<Display>();
		display->WaitForDisplayReady();
	}

	simulationTimer.emplace(TimeIt{ "Simulation" });
	time0 = std::chrono::steady_clock::now();
    auto t0 = std::chrono::steady_clock::now();
	while (true) {

		if (!handleDisplay(compounds, boxparams, display.get(), emVariant)) {
			break;
		}

		auto stepStartTime = std::chrono::steady_clock::now();
		
		engine->step();

		handleStatus(engine->runstatus.current_step);
		
		if (engine->runstatus.simulation_finished) {
			break;
		}

		// Deadspin to slow down rendering for visual debugging :)
		while ((double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - stepStartTime).count() < MIN_STEP_TIME) {}
	}
    auto t1 = std::chrono::steady_clock::now();
	simulationTimer->stop();

	// Transfers the remaining traj data and more
	engine->terminateSimulation();

	simulation = engine->takeBackSim();

	simulation->finished = true;
	simulation->ready_to_run = false ;

	
	m_logger.finishSection("Simulation Finished");

    return t1-t0;
}

void Environment::WriteBoxCoordinatesToFile(GroFile& grofile, std::optional<int64_t> _step) {	 	 
	int particlesUpdated = 0;

	const int64_t stepToLoadFrom = _step.value_or(simulation->getStep())-1;

	for (int cid = 0; cid < boximage->compounds.size(); cid++) {
		for (int pid = 0; pid < boximage->compounds[cid].n_particles; pid++) {
			const Float3 newPos = simulation->traj_buffer->GetMostRecentCompoundparticleDatapoint(cid, pid, stepToLoadFrom);
			grofile.atoms[boximage->compounds[cid].indicesInGrofile[pid]].position = newPos;
			particlesUpdated++;
		}
	}

	for (int tinymolId = 0; tinymolId < simulation->box_host->boxparams.n_solvents; tinymolId++) {

		const TinyMolFactory tinymol = boximage->solvent_positions[tinymolId];
		const int nAtomsInTinymol = tinymol.nParticles;

		const Float3 new_position = simulation->traj_buffer->GetMostRecentSolventparticleDatapointAtIndex(tinymolId, stepToLoadFrom);
		const Float3 deltaPos = new_position - grofile.atoms[tinymol.firstParticleIdInGrofile].position;

		assert(grofile.atoms[tinymol.firstParticleIdInGrofile].atomName[0] == tinymol.atomType[0]);

		for (int i = 0; i < nAtomsInTinymol; i++) {
			grofile.atoms[tinymol.firstParticleIdInGrofile + i].position += deltaPos;
			particlesUpdated++;
		}		
	}

	if (particlesUpdated != grofile.atoms.size()) {
		throw std::runtime_error(std::format("Only {} out of {} particles were updated", particlesUpdated, grofile.atoms.size()));
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
std::vector<Float3> Environment::GetForces(int64_t step) const {
	int particlesUpdated = 0;

	std::vector<Float3> forces(boximage->grofile.atoms.size()); // [kJ/mol/nm]



		for (int cid = 0; cid < boximage->compounds.size(); cid++) {
			for (int pid = 0; pid < boximage->compounds[cid].n_particles; pid++) {
				forces[boximage->compounds[cid].indicesInGrofile[pid]] = simulation->forceBuffer->GetMostRecentCompoundparticleDatapoint(cid, pid, step) / KILO;
				particlesUpdated++;
			}
		}

		for (int tinymolId = 0; tinymolId < simulation->box_host->boxparams.n_solvents; tinymolId++) {
			const TinyMolFactory tinymol = boximage->solvent_positions[tinymolId];
			const int nAtomsInTinymol = tinymol.nParticles;

			for (int i = 0; i < nAtomsInTinymol; i++) {
				forces[tinymol.firstParticleIdInGrofile + i] = simulation->forceBuffer->GetMostRecentSolventparticleDatapointAtIndex(tinymolId, step);
				particlesUpdated++;
			}
		}

		if (particlesUpdated != boximage->grofile.atoms.size()) {
			throw std::runtime_error(std::format("Only {} out of {} particles were updated", particlesUpdated, boximage->grofile.atoms.size()));
		}
	

	return forces;
}

Trajectory Environment::WriteSimToTrajectory() const {

	const int nSteps = simulation->getStep();
	const int nAtoms = boximage->grofile.atoms.size();

	Trajectory trajectory(nSteps, nAtoms, boximage->grofile.box_size, simulation->simparams_host.dt);


	for (int step = 0; step < nSteps; step += simulation->simparams_host.data_logging_interval) {

		for (int cid = 0; cid < boximage->compounds.size(); cid++) {
			for (int pid = 0; pid < boximage->compounds[cid].n_particles; pid++) {
				const int atomIndex = boximage->compounds[cid].indicesInGrofile[pid];
				trajectory.Set(step, atomIndex, simulation->traj_buffer->GetMostRecentCompoundparticleDatapoint(cid, pid, step));
			}
		}

		for (int tinymolId = 0; tinymolId < simulation->box_host->boxparams.n_solvents; tinymolId++) {
			const TinyMolFactory tinymol = boximage->solvent_positions[tinymolId];
			const int nAtomsInTinymol = tinymol.nParticles;

			const Float3 new_position = simulation->traj_buffer->GetMostRecentSolventparticleDatapointAtIndex(tinymolId, step);
			const Float3 deltaPos = new_position - boximage->grofile.atoms[tinymol.firstParticleIdInGrofile].position;

			for (int i = 0; i < nAtomsInTinymol; i++) {
				const int atomId = tinymol.firstParticleIdInGrofile + i;
				const Float3 newPos = boximage->grofile.atoms[atomId].position + deltaPos;
				trajectory.Set(step, atomId, newPos);
			}
		}
	}

	return trajectory;
}


void Environment::handleStatus(const int64_t step) {
	if (m_mode == Headless) {
		return;
	}

	if (step % STEPS_PER_UPDATE == STEPS_PER_UPDATE-1) {
		const std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - time0);
		const double duration_ms = duration.count();

		//// First clear the current line
		//printf("\r\033[K");
		// Move cursor to the beginning of the line and clear it
		printf("\033[1000D\033[K");

		printf("Step #%06llu", step);
		printf("\tAvg. time: %.2fms", duration_ms / STEPS_PER_UPDATE);

		time0 = std::chrono::steady_clock::now();
		avgStepTimes.emplace_back(duration_ms / STEPS_PER_UPDATE);
	}
}



bool Environment::handleDisplay(const std::vector<Compound>& compounds_host, const BoxParams& boxparams, Display* const display, bool emVariant) {
	if (m_mode != Full) {
		return true;
	}

	auto displayException = display->displayThreadException;
	if (displayException) {
		std::rethrow_exception(displayException);
	}

	if (engine->runstatus.stepForMostRecentData != step_at_last_render && engine->runstatus.most_recent_positions != nullptr) {
		
		const std::string info = emVariant
			? std::format("Step {:d} MaxForce {:.02f}", static_cast<int>(engine->runstatus.current_step), static_cast<float>(engine->runstatus.greatestForce))
			: std::format("Step {:d} Temp {:.02f}", static_cast<int>(engine->runstatus.current_step), static_cast<float>(engine->runstatus.current_temperature));

		display->Render(std::make_unique<Rendering::SimulationTask>(
			engine->runstatus.most_recent_positions, compounds_host, boxparams, info, coloringMethod
		));
		step_at_last_render = engine->runstatus.current_step;
		engine->runstatus.most_recent_positions = nullptr;
	}

	return !display->DisplaySelfTerminated();
}

std::unique_ptr<Simulation> Environment::getSim() {
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
