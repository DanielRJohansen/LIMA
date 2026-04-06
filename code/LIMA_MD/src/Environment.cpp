#include <chrono>
#include <filesystem>
#include <string>
#include <optional>
#include <numeric>

#include "Environment.h"
#include "MDFiles.h"
#include "CompoundBuilder.h"
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

void Environment::CreateSimulation(const Float3& boxsize_nm) {
	SimParams simparams{};
	simulation = std::make_unique<Simulation>(simparams, std::make_unique<Box>(Float3(boxsize_nm)));
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
	//simulation->forcefieldTinymol = boximage->tinymolTypes;
	simulation->forcefieldTest = boximage->nonbondedInteractionParams;
}

void Environment::CreateSimulation(Simulation& simulation_src, const SimParams params) {

	simulation.reset(new Simulation(params));
	BoxBuilder::copyBoxState(*simulation, std::move(simulation_src.box_host), simulation_src.getStep());

	simulation->forcefield = simulation_src.forcefield;
	//simulation->forcefieldTinymol = simulation_src.forcefieldTinymol;
	simulation->forcefieldTest = simulation_src.forcefieldTest;
}


std::tuple<GroFile, TopologyFile, SimParams> Environment::CreateSimulationFiles(Float3 boxlen) {
	GroFile grofile{};
	grofile.m_path = work_dir / "conf.gro";
	grofile.box_size = Float3{ boxlen };
	grofile.printToFile();

	TopologyFile topfile{};
	topfile.SetSystem("MySystem");
	topfile.path = work_dir / "topol.top";
	topfile.forcefieldInclude = TopologyFile::ForcefieldInclude("charmm27.ff/forcefield.itp");
	topfile.printToFile();
	topfile = TopologyFile{work_dir / "topol.top"}; // Reload the topfile to parse the ffinclude


	SimParams simparams{};
	simparams.dumpToFile(work_dir / "sim_params.txt");

	return { grofile, topfile, simparams };
}

void constexpr Environment::verifySimulationParameters() {	// Not yet implemented
	if (simulation->simparams_host.cutoff_nm != 1.2f) {// TODO: DANGER
		//throw std::runtime_error("Currently only cutoff 1.2 nm is supported, as that is hardcoded into the Coulumbforce Chebyshev Coefficients"); // TODO: figure out how to support other cutoff's again
	}
}

fs::path Environment::FixPath(const fs::path& path) const {
	if (path.is_absolute())
		return path;
	if (fs::exists(work_dir / path))
		return work_dir / path;
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
	
	verifySimulationParameters();
	verifyBox();
	simulation->ready_to_run = true;

	avgStepTimes.reserve((simulation->simparams_host.n_steps + 1) / STEPS_PER_UPDATE);

	//boxparams = simulation->box_host->boxparams;
	coloringMethod = simulation->simparams_host.coloring_method;


	engine = std::make_unique<Engine>(
		simulation.get(),
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
	const bool stepwise = simulation->simparams_host.stepwise;
	simparamsCopy = simulation->simparams_host;

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

		if (!handleDisplay(simulation->box_host->boxparams, display.get(), emVariant, stepwise)) {
			break;
		}

		auto stepStartTime = std::chrono::steady_clock::now();
		
		engine->step();

		handleStatus(engine->runstatus.current_step, emVariant);
		
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
	simparamsCopy.reset();

	simulation->finished = true;
	simulation->ready_to_run = false ;

	
	m_logger.finishSection("Simulation Finished");

	engineTime = t1 - t0;
    return t1-t0;
}



void Environment::InsertMolecule(GroFile& grofile, TopologyFile& topfile, LiveEdit::InsertMolecule& insertionCmd, SimParams simparams) {
	insertionCmd.groPath = FixPath(insertionCmd.groPath);
	insertionCmd.topPath = FixPath(insertionCmd.topPath);

	// First load the new data
	GroFile newmolGro(insertionCmd.groPath);	
	auto newmolTop = std::make_shared<TopologyFile>(insertionCmd.topPath);
	MoleculeUtils::CenterMolecule(newmolGro, newmolTop->GetMoleculeType());	// Make molecule whole

	BoundingBox newmolBb(newmolGro.atoms | std::views::transform([](auto& a) { return a.position; }));


	// Save the current state to current files
	WriteBoxCoordinatesToFile(grofile);	
	Float3 defaultInsertSite = Float3{ grofile.box_size.x / 2, grofile.box_size.y / 2, grofile.box_size.z - (newmolBb.Dimensions().z / 2.f) }; 
	Float3 insertionPosition = insertionCmd.position.value_or(defaultInsertSite);

	engine.reset(); // Kill the engine, since it holds a pointer to the sim which we will now change under it. We will make a new engine after creating the new sim
	SimulationBuilder::InsertSubmoleculeInSimulation(grofile, topfile, newmolGro, newmolTop, insertionPosition);
	CreateSimulation(grofile, topfile, simparams);
	//SimulationBuilder
	//FileUtils::mer
	//MDFiles::
}

void GatherPositionsIntoVector(std::vector<Float3>& dst, CudaBuffer<PersistentCluster>& pcBuffer, int nPclusters) {
	dst.resize(nPclusters * PersistentCluster::maxParticles);
	std::vector<PersistentCluster> pcHost = GenericCopyToHost(pcBuffer.Get(), nPclusters); // TODO: Reuse mem here somehow, this'll be slow..
	for (int pcId = 0; pcId < nPclusters; pcId++) {
		const PersistentCluster& pc = pcHost[pcId];
		for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
			dst[pcId * PersistentCluster::maxParticles + pid] = pc.pqd[pid].position;
		}
	}
}

void Environment::HandleDragMoleculeCommand(const LiveEdit::DragMolecule& newDragCommand, const LiveEdit::DragMolecule& prevDragCommand, std::vector<int>& affectedParticleIds, std::vector<Float3>&fixedMovements) {

	if (newDragCommand.particleId == prevDragCommand.particleId && newDragCommand.draggingForce == prevDragCommand.draggingForce) {
		return;
	}
	
	affectedParticleIds.clear();
	fixedMovements.clear();

	if (newDragCommand.particleId >= 0 && newDragCommand.draggingForce.len() > 0) {
		fixedMovements.resize(simulation->box_host->boxparams.totalParticles, Float3{0.f});
		for (const auto& node : boximage->systemGraph->BFS(newDragCommand.particleId)) {
			affectedParticleIds.push_back(node.atomid);
			fixedMovements[node.atomid] = newDragCommand.draggingForce * .1f;
			//fixedMovements.push_back(newDragCommand.draggingForce * .1f);
		}
	}


	//bool newId = newDragCommand.particleId != prevDragCommand.particleId;
	//bool newForce = newDragCommand.draggingForce != prevDragCommand.draggingForce;

	//if (newId || newForce) {
	//	for (const auto& id : affectedParticleIds) {
	//		fixedMovements[id] = Float3{ 0 };
	//	}
	//}

	//if (newId) {
	//	affectedParticleIds.clear();
	//	for (const auto& node : boximage->systemGraph->BFS(newDragCommand.particleId)) {
	//		affectedParticleIds.push_back(node.atomid);
	//	}		
	//}

	//if (newForce || newDragCommand.draggingForce.len() > 0) {
	//	float scale = .1f;
	//	for (const auto& id : affectedParticleIds) {
	//		fixedMovements[id] = Float3{ 0.01f, 0.f, 0.f };// newDragCommand.draggingForce* scale;
	//	}
	//}
}

void Environment::LiveEdit(GroFile& grofile, TopologyFile& topfile) {
	std::unique_ptr<Display> display = nullptr;

	simulation->simparams_host.n_steps = 0;
	simulation->simparams_host.data_logging_interval = 0;
	simulation->simparams_host.em_variant = false;

	display = std::make_unique<Display>();
	display->WaitForDisplayReady();	

	bool shouldExit = false;
	std::vector<Float3> positionData;
	bool shouldUpdateRender = true;

	// MoleculeDragging
	LiveEdit::DragMolecule prevDragmoleculeCmd{};
	std::vector<int> affectedParticleIds;
	std::vector<Float3> fixedMovements;

	int remainingStepsCount = 0;

	auto GetNextCommand = [&]() -> std::optional<LiveEdit::Command> {
		if (!liveEditCommandsQueue.empty()) {
			LiveEdit::Command cmd = liveEditCommandsQueue.front();
			liveEditCommandsQueue.pop_front();
			return cmd;
		}
		return display->GetLiveEditCommand();
		};

	while (true) {
		if (shouldExit) {
			break;
		}

		if (display->DisplaySelfTerminated()) {
			break;
		}

		// Poll interface for new commands, and execute if any		
		if (auto newCmd = GetNextCommand()) {
			std::visit(
				[&](auto&& cmd) {
					using T = std::decay_t<decltype(cmd)>;

					if constexpr (std::is_same_v<T, LiveEdit::Invalid>) {
						return;
					}
					else if constexpr (std::is_same_v<T, LiveEdit::InsertMolecule>) {
						InsertMolecule(grofile, topfile, cmd, simulation->simparams_host);
						fixedMovements.resize(simulation->box_host->boxparams.totalParticles, Float3{ 0 });
						remainingStepsCount = 50;
						prevDragmoleculeCmd = LiveEdit::DragMolecule{};
					}
					else if constexpr (std::is_same_v<T, LiveEdit::DragMolecule>) {
						HandleDragMoleculeCommand(cmd, prevDragmoleculeCmd, affectedParticleIds, fixedMovements);
						prevDragmoleculeCmd = cmd;
						engine->SetFixedParticleMovementBuffer(fixedMovements);
						remainingStepsCount = 50;
					}
				},
				*newCmd
			);
		}

		// Decide if engine should run, if so new params?

		// Run engine
		if (!engine && simulation->box_host->boxparams.totalParticles > 0) {
			engine = std::make_unique<Engine>(
				simulation.get(),
				simulation->simparams_host.bc_select,
				std::make_unique<LimaLogger>(LimaLogger::compact, m_mode, "engine", work_dir));

			 auto& pcBuffer = engine->OffloadPclusterState();
			 GatherPositionsIntoVector(positionData, pcBuffer, simulation->box_host->persistentClusters.size());
			 shouldUpdateRender = true;
		}
		//printf("Step count %d\n", remainingStepsCount);
		if (engine && remainingStepsCount > 0) {
			// Add step logic here
			//shouldUpdateRender = true;
			engine->step();
			auto& pcBuffer = engine->OffloadPclusterState();
			GatherPositionsIntoVector(positionData, pcBuffer, simulation->box_host->persistentClusters.size());
			shouldUpdateRender = true;
			remainingStepsCount--;
			if (remainingStepsCount == 0) {
				// check engine if we should continue..
			}
		}

		if (shouldUpdateRender) {
			display->Render(std::make_unique<Rendering::SimulationTask>(
				positionData.data(), simulation->box_host->persistentClusters, simulation->box_host->persistentClustersMetadata, simulation->box_host->boxparams, "", coloringMethod, simStatus
			), false);
			shouldUpdateRender = false;
		}
	}
}

void Environment::WriteBoxCoordinatesToFile(GroFile& grofile, std::optional<int64_t> _step) {	 	 
	int particlesUpdated = 0;

	
	// First offload the current state from engine to host - if there is no engine, the state in the current boxhost IS the current state
	if (engine) {
		CudaBuffer<PersistentCluster>& pcBuffer = engine->OffloadPclusterState();
		simulation->box_host->persistentClusters = GenericCopyToHost(pcBuffer.Get(), simulation->box_host->persistentClusters.size()); // TODO: Reuse mem here somehow, this'll be slow..
	}

	for (int pcId = 0; pcId < simulation->box_host->persistentClusters.size(); pcId++) {
		const PersistentClusterMeta& pcMeta = simulation->box_host->persistentClustersMetadata[pcId];
		const PersistentCluster& pcData = simulation->box_host->persistentClusters[pcId];
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

//TODO!!

		//for (int cid = 0; cid < boximage->compounds.size(); cid++) {
		//	for (int pid = 0; pid < boximage->compounds[cid].n_particles; pid++) {
		//		forces[boximage->compounds[cid].indicesInGrofile[pid]] = simulation->forceBuffer->GetMostRecentCompoundparticleDatapoint(cid, pid, step) / KILO;
		//		particlesUpdated++;
		//	}
		//}

		//for (int tinymolId = 0; tinymolId < simulation->box_host->boxparams.nTinymols; tinymolId++) {
		//	const TinyMolFactory tinymol = boximage->solvent_positions[tinymolId];
		//	const int nAtomsInTinymol = tinymol.nParticles;

		//	for (int i = 0; i < nAtomsInTinymol; i++) {
		//		forces[tinymol.firstParticleIdInGrofile + i] = simulation->forceBuffer->GetMostRecentSolventparticleDatapointAtIndex(tinymolId, step);
		//		particlesUpdated++;
		//	}
		//}

		//if (particlesUpdated != boximage->grofile.atoms.size()) {
		//	throw std::runtime_error(std::format("Only {} out of {} particles were updated", particlesUpdated, boximage->grofile.atoms.size()));
		//}
	

	return forces;
}

Trajectory Environment::WriteSimToTrajectory() const {

	const int nSteps = simulation->getStep();
	const int nAtoms = boximage->grofile.atoms.size();

	Trajectory trajectory(nSteps, nAtoms, boximage->grofile.box_size, simulation->simparams_host.dt);

	// TODO!!
	//for (int step = 0; step < nSteps; step += simulation->simparams_host.data_logging_interval) {

	//	for (int cid = 0; cid < boximage->compounds.size(); cid++) {
	//		for (int pid = 0; pid < boximage->compounds[cid].n_particles; pid++) {
	//			const int atomIndex = boximage->compounds[cid].indicesInGrofile[pid];
	//			trajectory.Set(step, atomIndex, simulation->traj_buffer->GetMostRecentCompoundparticleDatapoint(cid, pid, step));
	//		}
	//	}

	//	for (int tinymolId = 0; tinymolId < simulation->box_host->boxparams.nTinymols; tinymolId++) {
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
	const int nSteps = simulation->getStep();
	const int nAtoms = boximage->grofile.atoms.size();

	UpgradeableFileFormat file(path);

	file.WriteSection("numAtoms", std::vector{ nAtoms });
	file.WriteSection("numFrames", std::vector{nSteps / simulation->simparams_host.data_logging_interval});
	file.WriteSection("trajectory", simulation->traj_buffer->GetBuffer());
}

void Environment::handleStatus(const int64_t step, bool emVariant) {
	if (m_mode == Headless) {
		return;
	}

	if (step % STEPS_PER_UPDATE == STEPS_PER_UPDATE-1) {		
		auto duration = std::chrono::steady_clock::now() - time0;
		const double duration_ms = std::chrono::duration_cast<std::chrono::microseconds>(duration).count() * 1e-3;
		const double avgSteptime = duration_ms / (double) STEPS_PER_UPDATE;
		//// First clear the current line
		//printf("\r\033[K");
		// Move cursor to the beginning of the line and clear it
		printf("\033[1000D\033[K");

		printf("Step #%06llu", step);
		printf("\tAvg. time: %.2fms", avgSteptime);

		time0 = std::chrono::steady_clock::now();
		avgStepTimes.emplace_back(avgSteptime);



		SimStatus newStatus{};
		newStatus.step = engine->runstatus.current_step;
		newStatus.maxForce = emVariant ? std::optional<float>(engine->runstatus.greatestForce) : std::nullopt;
		newStatus.temperature = !emVariant ? std::optional<float>(engine->runstatus.current_temperature) : std::nullopt;
		newStatus.avgStepTime = avgStepTimes.empty() ? 0.f : avgStepTimes.back();
		const int nStepsSinceLast = engine->runstatus.current_step - simStatus.step;
		const double totalNsSimulated = nStepsSinceLast * simparamsCopy->dt; // [ns]
		const double wall_time_sec = duration_ms * 1e-3;
		const double ns_per_day = totalNsSimulated / (wall_time_sec / 86400.0);  // 86400 seconds in a day
		newStatus.simulationPerformance = ns_per_day;
		simStatus = newStatus;
	}
}



bool Environment::handleDisplay(const BoxParams& boxparams, Display* const display, bool emVariant, bool stepwise) {
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
		display->Render(std::make_unique<Rendering::SimulationTask>(
			renderPositions, simulation->box_host->persistentClusters, simulation->box_host->persistentClustersMetadata, boxparams, info, coloringMethod, simStatus
		), stepwise);
		step_at_last_render = stepForMostRecentData;
		//engine->runstatus.most_recent_positions = nullptr;
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

void Environment::PrintTiming() const {	
	if (!engineTime || !simulation)
		return;

	const double wall_time_sec = engineTime->count();
	const double totalNsSimulated = static_cast<double>(simulation->getStep()) * simulation->simparams_host.dt;

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