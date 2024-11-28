#pragma once

#include "Analyzer.cuh"
#include "Bodies.cuh"
#include "TimeIt.h"
#include "MDFiles.h"

#include <memory>
#include <chrono>

class Display;
struct BoxImage;
class Engine;


namespace fs = std::filesystem;


class Environment
{
public:
	Environment() = delete;
	Environment(const Environment&) = delete;
	Environment(const fs::path& workdir, EnvMode mode);

	~Environment();

	/// <summary>
	/// Create a simulation, and create the necessary files in process, if the defaults
	/// (conf.gro and topol.top and simparams.txt) are not available
	/// </summary>
	void CreateSimulation(float boxsize_nm);

	/// <summary>
	/// Create a simulation from existing files
	/// </summary>
	void CreateSimulation(const std::string conf_filename, std::string topol_filename, SimParams);	// TODO make constref

	// The basic createSim
	void CreateSimulation(const GroFile&, const TopologyFile&, const SimParams&);


	/// <summary>
	/// Create a simulation that starts from where boxorigin is currently
	/// </summary>
	void CreateSimulation(Simulation& simulation_src, SimParams);

	/// <summary>
	/// Create .gro .top and simparams.txt files in the current directory
	/// </summary>
	void createSimulationFiles(float boxlen);

	// Run a standard MD sim
	void run(bool doPostRunEvents=true);

	/// <summary>
	/// Intended to be called after a sim run, uses the BoxImage to write new coordinates for the
	/// atoms in the input coordinate file.
	/// </summary>
	/// <param name="filename">New name of the file. Defaults to same name as it was as input</param>
	/// <returns></returns>
	GroFile WriteBoxCoordinatesToFile(const std::optional<std::string> filename= "out");
	void WriteBoxCoordinatesToFile(GroFile& grofile, std::optional<int64_t> step=std::nullopt);

	void RenderSimulation();
	
	
	
	
	void renderTrajectory(std::string trj_path);
	
	void makeVirtualTrajectory(std::string trj_path, std::string waterforce_path);

	// Functions for dev only : TODO move to child whioch inherits all as public
	std::unique_ptr<Simulation> getSim();
	Simulation* getSimPtr();
	const SimAnalysis::AnalyzedPackage& getAnalyzedPackage();
	SolventBlocksCircularQueue* getSolventBlocks();

	std::string getWorkdir() { return work_dir.string(); }

	std::chrono::steady_clock::time_point time0;

	const fs::path work_dir = "";	// Main dir of the current simulation

	std::optional<TimeIt> simulationTimer;
	std::vector<float> avgStepTimes; // [ms] - averaged over STEP_PER_UPDATE

private:

	void constexpr verifySimulationParameters();			// Constants before doing anything
	void verifyBox();							// Checks wheter the box will break
	
	void postRunEvents();
	void handleStatus(int64_t step);

	// Returns false if display has been closed by user
	bool handleDisplay(const std::vector<Compound>& compounds_host, const BoxParams& boxparams, Display* const display, bool emVariant);

	void sayHello();

	bool prepareForRun();


	EnvMode m_mode;

	int64_t step_at_last_render = INT64_MIN;


	//std::unique_ptr<BoxBuilder> boxbuilder;
	LimaLogger m_logger;



	std::unique_ptr<Engine> engine;
	std::unique_ptr<Simulation> simulation;

	ColoringMethod coloringMethod;	// Not ideal to have here..

	// TEMP: Cache some constants here before we give ownership to engine. DO NOT READ VOLATILE VALUES FROM THESE
	std::vector<Compound> compounds;
	BoxParams boxparams;

	std::unique_ptr<BoxImage> boximage;

	std::optional<SimAnalysis::AnalyzedPackage> postsim_anal_package;
};
