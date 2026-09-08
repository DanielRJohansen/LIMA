#include <chrono>
#include <filesystem>
#include <string>
#include <optional>
#include <numeric>

#include "Environment.h"
#include "MDFiles.h"
#include "BoxImageBuilder.h"
#include "Display.h"
#include "BoxBuilder.cuh"
#include "Engine.cuh"
#include "UpgradeableFileFormat.h"
#include "SimulationBuilder.h"
#include "MoleculeUtils.h"

namespace lfs = FileUtils;
namespace fs = std::filesystem;


// ------------------------------------------------ Display Parameters ------------------------------------------ //
const int STEPS_PER_UPDATE = 100;
constexpr float MIN_STEP_TIME = 0.f;		// [ms] Set to 0 for full speed sim
// -------------------------------------------------------------------------------------------------------------- //

Environment::SimulationSession::SimulationSession(std::unique_ptr<Simulation> simulation, EnvMode mode, const fs::path& workDir)
	: simulation(std::move(simulation))
	, logger(LimaLogger::compact, mode, "environment", workDir.string()) {}

Environment::SimulationSession::~SimulationSession() = default;

Environment::SimulationSession& Environment::Session() {
	if (!simulationSession)
		throw std::runtime_error("Environment has no simulation session");
	return *simulationSession;
}

const Environment::SimulationSession& Environment::Session() const {
	if (!simulationSession)
		throw std::runtime_error("Environment has no simulation session");
	return *simulationSession;
}

void Environment::SetSimulation(std::unique_ptr<Simulation> simulation) {
	if (!simulationSession) {
		simulationSession = std::make_unique<SimulationSession>(std::move(simulation), m_mode, workDir);
		return;
	}

	simulationSession->simulation = std::move(simulation);
}

Environment::Environment(const fs::path& workdir, EnvMode mode)
	: workDir(workdir)
	, m_mode(mode)
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

void Environment::CreateSimulation(const Float3& boxsize_nm) {
	SimParams simparams{};
	SetSimulation(std::make_unique<Simulation>(simparams, std::make_unique<Box>(Float3(boxsize_nm))));
}

void Environment::CreateSimulation(const std::string& gro_path, const std::string& topol_path, const SimParams& params) {
	const auto groFile = std::make_unique<GroFile>(gro_path);
	const auto topFile = std::make_unique<TopologyFile>(topol_path);
	CreateSimulation(*groFile, *topFile, params);
}

void Environment::CreateSimulation(const GroFile& grofile, const TopologyFile& topolfile, const SimParams& params) 
{
	auto boxImage = LIMA_MOLECULEBUILD::buildMolecules(
		grofile,
		topolfile,
		V1,
		std::make_unique<LimaLogger>(LimaLogger::normal, m_mode, "moleculebuilder", workDir),
		IGNORE_HYDROGEN,
		params
		);

	auto simulation = std::make_unique<Simulation>(params, BoxBuilder::BuildBox(params, *boxImage));
	simulation->boxImage = std::shared_ptr<BoxImage>(std::move(boxImage));
	SetSimulation(std::move(simulation));
	SimulationSession& session = Session();

	if (display) {
		display->Render(std::make_unique<Rendering::AtomRenderTask>(
			session.simulation->box->persistentClusters, session.simulation->box->persistentClustersMetadata,
			session.simulation->box->boxparams, session.simStatus, session.simulation->box->backboneChains
		));
	}
}

void Environment::CreateSimulation(Simulation& simulation_src, const SimParams params) {

	auto simulation = std::make_unique<Simulation>(params);
	BoxBuilder::copyBoxState(*simulation, std::move(simulation_src.box), simulation_src.getStep());
	simulation->boxImage = std::move(simulation_src.boxImage);
	SetSimulation(std::move(simulation));
}


std::tuple<GroFile, TopologyFile, SimParams> Environment::CreateSimulationFiles(Float3 boxlen) {
	GroFile grofile{};
	grofile.m_path = workDir / "conf.gro";
	grofile.box_size = Float3{ boxlen };
	grofile.printToFile();

	TopologyFile topfile{};
	topfile.SetSystem("MySystem");
	topfile.path = workDir / "topol.top";
	//topfile.forcefieldInclude = TopologyFile::ForcefieldInclude("charmm27.ff/forcefield.itp");
	topfile.printToFile();
	topfile = TopologyFile{workDir / "topol.top"}; // Reload the topfile to parse the ffinclude


	SimParams simparams{};
	simparams.DumpToFile(workDir / "sim_params.txt");

	return { grofile, topfile, simparams };
}

void Environment::verifySimulationParameters() {	// Not yet implemented
	const auto& simulation = Session().simulation;
	if (simulation->simParams.cutoff_nm != 1.2f) {// TODO: DANGER
		//throw std::runtime_error("Currently only cutoff 1.2 nm is supported, as that is hardcoded into the Coulumbforce Chebyshev Coefficients"); // TODO: figure out how to support other cutoff's again
	}
}

fs::path Environment::FixPath(const fs::path& path) const {
	if (path.is_absolute())
		return path;
	if (fs::exists(workDir / path))
		return workDir / path;
	if (fs::exists( "./" / path))
		return "./" / path;
	return path;
}

void Environment::verifyBox() {


	

	




//#ifdef LIMAKERNELDEBUGMODE
//	if (print_compound_positions) {
//		for (int c = 0; c < simulation->boxparams_host.n_compounds; c++) {
//			Compound* comp = &simulation->compounds_host[c];
//			for (int p = 0; p < comp->n_particles; p++) {
//				printf("%d   ", comp->particle_global_ids[p]);
//			}
//		}
//	}
//#endif
}

bool Environment::prepareForRun() {
	SimulationSession& session = Session();
	auto& simulation = session.simulation;
	if (simulation == nullptr)// TEMP, ENv should never give sim to engine
		return true;

	if (simulation->finished) { 
		printf("Cannot prepare run, since simulation has already finished");
		assert(false);
		return false; 
	}

	session.logger.startSection("Simulation started");

	if (simulation->ready_to_run) { return true; }

	simulation->PrepareDataBuffers();
	
	verifySimulationParameters();
	verifyBox();
	simulation->ready_to_run = true;

	session.avgStepTimes.reserve((simulation->simParams.n_steps + 1) / STEPS_PER_UPDATE);


	session.engine = std::make_unique<Engine>(
		simulation.get(),
		simulation->simParams.bc_select);

	return true;
}


void Environment::sayHello() {
	static bool hasSaidHello = false;
	if (hasSaidHello)
		return;
	hasSaidHello = true;

	std::ifstream file(FileUtils::GetLimaDir() / "resources/logo/logo_ascii.txt");
	if (!file) {
		throw std::runtime_error("Failed to open logo file");
	}

	std::string file_contents((std::istreambuf_iterator<char>(file)),
		std::istreambuf_iterator<char>());

	std::cout << file_contents;
}

std::chrono::duration<double> Environment::run() {
	SimulationSession& session = Session();
	auto& simulation = session.simulation;
	auto& engine = session.engine;
	auto& simStatus = session.simStatus;
	auto& time0 = session.time0;
	auto& simulationTimer = session.simulationTimer;
	auto& engineTime = session.engineTime;
	const bool emVariant = simulation->simParams.em_variant;
	const bool stepwise = simulation->simParams.stepwise;
	//simparamsCopy = simulation->simParams;

    if (!prepareForRun()) { return {}; }

	std::unique_ptr<Display> display = nullptr;

	if (m_mode == Full) {
		display = std::make_unique<Display>();
		display->WaitForDisplayReady();
		display->Render(std::make_unique<Rendering::AtomRenderTask>(
			simulation->box->persistentClusters, simulation->box->persistentClustersMetadata,
			simulation->box->boxparams, simStatus, simulation->box->backboneChains
		), stepwise);
	}

	simulationTimer.emplace(TimeIt{ "Simulation" });
	time0 = std::chrono::steady_clock::now();
    auto t0 = std::chrono::steady_clock::now();
	while (true) {

		if (!handleDisplay(simulation->box->boxparams, display.get(), emVariant, stepwise)) {
			break;
		}

		auto stepStartTime = std::chrono::steady_clock::now();
		
		engine->step();

		UpdateSimstatus(true, true);
		
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

	//simulation = engine->takeBackSim();
	//simparamsCopy.reset();

	simulation->finished = true;
	simulation->ready_to_run = false ;

	
	session.logger.finishSection("Simulation Finished");

	engineTime = t1 - t0;
    return t1-t0;
}





void Environment::WriteBoxCoordinatesToFile(GroFile& grofile, std::optional<int64_t> _step) {	 	 
	SimulationSession& session = Session();
	auto& simulation = session.simulation;
	auto& engine = session.engine;
	int particlesUpdated = 0;

	
	// First offload the current state from engine to host - if there is no engine, the state in the current boxhost IS the current state
	if (engine) {
		CudaBuffer<PersistentCluster>& pcBuffer = engine->OffloadPclusterState();
		simulation->box->persistentClusters = GenericCopyToHost(pcBuffer.Get(), simulation->box->persistentClusters.size()); // TODO: Reuse mem here somehow, this'll be slow..
	}

	for (int pcId = 0; pcId < simulation->box->persistentClusters.size(); pcId++) {
		const PersistentClusterMeta& pcMeta = simulation->box->persistentClustersMetadata[pcId];
		const PersistentCluster& pcData = simulation->box->persistentClusters[pcId];
		for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
			const int pidGlobal = pcMeta.particleIdsGlobal[pid];
			if (pidGlobal != -1) {
				grofile.atoms[pidGlobal].position = pcData.pqd[pid].position;
				particlesUpdated++;
			}
		}
	}

	if (AllAtom && particlesUpdated != grofile.atoms.size()) {
		throw std::runtime_error(std::format("Only {} out of {} particles were updated", particlesUpdated, grofile.atoms.size()));
	}
}
GroFile Environment::WriteBoxCoordinatesToFile(const std::optional<std::string> filename) {
	const auto& boximage = Session().simulation->boxImage;
	if (!boximage)
		throw std::runtime_error("Cannot write coordinates without a BoxImage");
	GroFile outputfile{ boximage->grofile };

	if (filename.has_value()) {
		outputfile.m_path = workDir / "molecule" / (filename.value() + ".gro");
	}

	WriteBoxCoordinatesToFile(outputfile);

	return outputfile;
}
std::vector<Float3> Environment::GetForces(int64_t step) const {
	const auto& simulation = Session().simulation;
	const auto& boximage = simulation->boxImage;
	if (!boximage)
		throw std::runtime_error("Cannot get forces without a BoxImage");
	std::vector<Float3> forces(boximage->grofile.atoms.size());		// [kJ/mol/nm]
	const auto& forcesBuffer = *simulation->forceBuffer;			// [J/mol/nm]
	for (int pcid = 0; pcid < simulation->box->persistentClusters.size(); pcid++) {
		for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
			const int gpid = simulation->box->persistentClustersMetadata[pcid].particleIdsGlobal[pid]; 			
			if (gpid == -1)
				continue;

			const Float3 force = forcesBuffer.GetDatapointAtStep(pcid, pid, step);
			forces[gpid] = force / KILO;
		}
	}

	return forces;
}

Trajectory Environment::WriteSimToTrajectory() const {
	const auto& simulation = Session().simulation;
	const auto& boximage = simulation->boxImage;
	if (!boximage)
		throw std::runtime_error("Cannot write a trajectory without a BoxImage");

	const int nSteps = simulation->getStep();
	const int nAtoms = boximage->grofile.atoms.size();

	Trajectory trajectory(nSteps, nAtoms, boximage->grofile.box_size, simulation->simParams.dt);

	// TODO!!
	//for (int step = 0; step < nSteps; step += simulation->simParams.data_logging_interval) {

	//	for (int cid = 0; cid < boximage->compounds.size(); cid++) {
	//		for (int pid = 0; pid < boximage->compounds[cid].n_particles; pid++) {
	//			const int atomIndex = boximage->compounds[cid].indicesInGrofile[pid];
	//			trajectory.Set(step, atomIndex, simulation->traj_buffer->GetMostRecentCompoundparticleDatapoint(cid, pid, step));
	//		}
	//	}

	//	for (int tinymolId = 0; tinymolId < simulation->box->boxparams.nTinymols; tinymolId++) {
	//		const TinyMolFactory tinymol = boximage->solvent_positions[tinymolId];
	//		const int nAtomsInTinymol = tinymol.nParticles;

	//		const Float3 new_position = simulation->traj_buffer->GetMostRecentSolventparticleDatapointAtIndex(tinymolId, step);
	//		const Float3 deltaPos = new_position - boximage->grofile.atoms[tinymol.firstParticleIdInGrofile].position;

	//		for (int i = 0; i < nAtomsInTinymol; i++) {
	//			const int atomId = tinymol.firstParticleIdInGrofile + i;
	//			const Float3 newPos = boximage->grofile.atoms[atomId].position + deltaPos;
	//			trajectory.Set(step, atomId, newPos);
	//		}
	//	}
	//}

	return trajectory;
}

void Environment::WriteTrajectoryAsUff(const fs::path& path) const {
	const auto& simulation = Session().simulation;
	const auto& boximage = simulation->boxImage;
	if (!boximage)
		throw std::runtime_error("Cannot write a trajectory without a BoxImage");
	const int nSteps = simulation->getStep();
	const int nAtoms = boximage->grofile.atoms.size();

	UpgradeableFileFormat file(path);

	file.WriteSection("numAtoms", std::vector{ nAtoms });
	file.WriteSection("numFrames", std::vector{nSteps / simulation->simParams.data_logging_interval});
	file.WriteSection("trajectory", simulation->traj_buffer->GetBuffer());
}

void Environment::UpdateSimstatus(bool printToConsole, bool alwaysUpdate) {
	SimulationSession& session = Session();
	auto& simulation = session.simulation;
	auto& engine = session.engine;
	auto& simStatus = session.simStatus;
	auto& time0 = session.time0;
	auto& avgStepTimes = session.avgStepTimes;
	auto& simulationTimer = session.simulationTimer;
	auto& forceWriteSimstatusToDisplay = session.forceWriteSimstatusToDisplay;
	if (!simulation || !engine) {
		return;
	}

	const int64_t step = simulation->getStep();
	if ((step % STEPS_PER_UPDATE == STEPS_PER_UPDATE-1) || forceWriteSimstatusToDisplay) {		
		forceWriteSimstatusToDisplay = false;
		auto duration = std::chrono::steady_clock::now() - time0;		
		const double duration_ms = std::chrono::duration_cast<std::chrono::microseconds>(duration).count() * 1e-3;
		const double avgSteptime = duration_ms / (double) STEPS_PER_UPDATE;

		if (printToConsole && m_mode == Full) {
			//// First clear the current line
			//printf("\r\033[K");
			// Move cursor to the beginning of the line and clear it
			printf("\033[1000D\033[K");

			printf("Step #%06llu", step);
			printf("\tAvg. time: %.2fms", avgSteptime);
		}

		time0 = std::chrono::steady_clock::now();
		avgStepTimes.emplace_back(avgSteptime);




		const int nStepsSinceLast = engine->runstatus.current_step - *simStatus.step;
		const double totalNsSimulated = nStepsSinceLast * simulation->simParams.dt; // [ns]
		const double wall_time_sec = duration_ms * 1e-3;
		const double ns_per_day = totalNsSimulated / (wall_time_sec / 86400.0);  // 86400 seconds in a day
		const double completionFraction = (double)step / (double)simulation->simParams.n_steps;
		const std::optional<std::chrono::duration<double>> expectedTimeToFinish = simulation->simParams.n_steps > 0 && simulationTimer.has_value()
			? std::optional<std::chrono::duration<double>> {simulationTimer->Elapsed()* (1. / completionFraction * (1.-completionFraction))}
			: std::nullopt;

		SimStatus newStatus{};
		newStatus.step = engine->runstatus.current_step;
		newStatus.avgStepTime = avgStepTimes.empty() ? 0.f : avgStepTimes.back();
		newStatus.expectedTimeToFinish = expectedTimeToFinish;
		if (simulation->simParams.em_variant) {
		}
		else {
			newStatus.simulationPerformance = ns_per_day;
		}

		simStatus = newStatus;
	}

	// "Free" updates
	if (simulation->simParams.em_variant) {
		simStatus.maxForce = engine->runstatus.greatestForce;
	}
	else {
		if (!std::isnan(engine->runstatus.current_temperature))
			simStatus.temperature = engine->runstatus.current_temperature;
	}
}



bool Environment::handleDisplay(const BoxParams& boxparams, Display* const display, bool emVariant, bool stepwise) {
	SimulationSession& session = Session();
	auto& engine = session.engine;
	auto& simStatus = session.simStatus;
	auto& step_at_last_render = session.stepAtLastRender;
	if (m_mode != Full) {
		return true;
	}

	auto displayException = display->displayThreadException;
	if (displayException) {
		std::rethrow_exception(displayException);
	}

	int64_t stepForMostRecentData = engine ? engine->runstatus.stepForMostRecentData : -1;
	Float3* renderPositions = engine ? engine->runstatus.most_recent_positions : nullptr;
	std::string info{};

	if (engine) {
		info = emVariant
			? std::format("Step {:d} MaxForce {:.02f}", static_cast<int>(engine->runstatus.current_step), static_cast<float>(engine->runstatus.greatestForce))
			: std::format("Step {:d} Temp {:.02f}", static_cast<int>(engine->runstatus.current_step), static_cast<float>(engine->runstatus.current_temperature));
	}

	if (stepForMostRecentData > step_at_last_render) {
		display->Render(std::make_unique<Rendering::SimulationTaskUpdate>(
			renderPositions, nullptr, simStatus
		), stepwise);
		step_at_last_render = stepForMostRecentData;
		//engine->runstatus.most_recent_positions = nullptr;
	}

	return !display->DisplaySelfTerminated();
}

std::unique_ptr<Simulation> Environment::GetSim() {
	SimulationSession& session = Session();
	auto& engine = session.engine;
	auto& simulation = session.simulation;
	engine.reset();
	return std::move(simulation);
}

void Environment::ReleaseEngine() {
	Session().engine.reset();
}

Simulation* Environment::getSimPtr() {
	const auto& simulation = Session().simulation;
	if (simulation) { 
		return simulation.get(); 
	}
	return nullptr;
}

const SimAnalysis::AnalyzedPackage& Environment::getAnalyzedPackage()
{
	SimulationSession& session = Session();
	auto& simulation = session.simulation;
	auto& postsim_anal_package = session.analyzedPackage;
	if (simulation == nullptr)
		throw std::runtime_error("Env has no simulation");
	// TODO: make some check here that the simulation has finished
	if (!postsim_anal_package.has_value())
		postsim_anal_package = SimAnalysis::analyzeEnergy(simulation.get());
	return postsim_anal_package.value();
}

const std::optional<TimeIt>& Environment::SimulationTimer() const {
	return Session().simulationTimer;
}

const std::vector<float>& Environment::AverageStepTimes() const {
	return Session().avgStepTimes;
}

void Environment::QueueLiveEditCommand(LiveEdit::Command command) {
	Session().liveEditCommandsQueue.push_back(std::move(command));
}

void Environment::PrintTiming() const {	
	const SimulationSession& session = Session();
	const auto& simulation = session.simulation;
	const auto& engineTime = session.engineTime;
	if (!engineTime || !simulation)
		return;

	const double wall_time_sec = engineTime->count();
	const double totalNsSimulated = static_cast<double>(simulation->getStep()) * simulation->simParams.dt;

	// Calculate performance metrics
	const double ns_per_day = totalNsSimulated / (wall_time_sec / 86400.0);  // 86400 seconds in a day
	const double hr_per_ns = (wall_time_sec / totalNsSimulated) / 3600.0;    // convert to hours per ns

	// Print time and performance info in the GROMACS-like format
	printf("\n");
	printf("               Wall t (s)\n");
	printf("       Time:    %10.3f\n", wall_time_sec);
	printf("                 (ns/day)    (hour/ns)\n");
	printf("Performance:    %10.3f     %10.3f\n", ns_per_day, hr_per_ns);

}
