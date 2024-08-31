#pragma once

#include <memory>
#include <chrono>

#include "Analyzer.cuh"
#include "SimulationBuilder.h"
#include "Bodies.cuh"
#include "DisplayV2.h"
#include "TimeIt.h"

class BoxBuilder;
class Display;
struct BoxImage;
class Engine;




namespace fs = std::filesystem;


class Environment
{
public:
	Environment() = delete;
	Environment(const Environment&) = delete;
	Environment(const fs::path& workdir, EnvMode mode, bool save_output);

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
	/// Create a lipid bi-layer in the x-y plane.
	/// </summary>
	/// <param name="carryout_em">Carry out an energy minimization with no boundary condition, 
	/// which ensures all particles are inside the box</param>
	//void createMembrane(LipidsSelection& lipidselection, bool carryout_em = true);

	/// <summary>
	/// Create .gro .top and simparams.txt files in the current directory
	/// </summary>
	void createSimulationFiles(float boxlen);



	// Run a standard MD sim
	void run(bool doPostRunEvents=true);

	// Intended to be called after a sim run, uses the BoxImage to write new coordinates for the
	// atoms in the input coordinate file.
	GroFile writeBoxCoordinatesToFile(const std::string& filename="out");

	void RenderSimulation();
	
	
	
	
	void renderTrajectory(std::string trj_path);
	
	void makeVirtualTrajectory(std::string trj_path, std::string waterforce_path);

	// Functions for dev only : TODO move to child whioch inherits all as public
	std::unique_ptr<Simulation> getSim();
	Simulation* getSimPtr();
	Analyzer::AnalyzedPackage* getAnalyzedPackage();
	SolventBlocksCircularQueue* getSolventBlocks();

	std::string getWorkdir() { return work_dir.string(); }

#ifdef __linux__
	std::chrono::system_clock::time_point time0;
	std::string main_dir = "/opt/LIMA";
#else
	std::chrono::steady_clock::time_point time0;

	// Main folder of the lima package, contains code/lib/resources and more
	const std::string main_dir = "C:/Users/Daniel/git_repo/LIMA/";
#endif

	const fs::path mainDir = Filehandler::GetLimaDir();

	const fs::path work_dir = "";	// Main dir of the current simulation

	std::optional<TimeIt> simulationTimer;


private:
	void resetEnvironment();	// TODO: This shouldn't exist, rather we should just destruct the entire object

	void setupEmptySimulation(const SimParams&);
	void constexpr verifySimulationParameters();			// Constants before doing anything
	void verifyBox();							// Checks wheter the box will break
	
	void postRunEvents();
	void handleStatus(int64_t step, int64_t n_steps);

	// Returns false if display has been closed by user
	bool handleDisplay(const std::vector<Compound>& compounds_host, const BoxParams& boxparams, Display& display);

	void sayHello();

	bool prepareForRun();


	EnvMode m_mode;
	const bool save_output;

	//std::unique_ptr<Display> display = nullptr;
	int step_at_last_render = 0;
	int step_at_last_update = 0;


	std::unique_ptr<BoxBuilder> boxbuilder;
	LimaLogger m_logger;



	std::unique_ptr<Engine> engine;
	std::unique_ptr<Simulation> simulation;
	std::unique_ptr<Display> display = nullptr;

	ColoringMethod coloringMethod;	// Not ideal to have here..

	// TEMP: Cache some constants here before we give ownership to engine. DO NOT READ VOLATILE VALUES FROM THESE
	std::vector<Compound>* compounds = nullptr;
	BoxParams boxparams;

	std::unique_ptr<BoxImage> boximage;

	Analyzer::AnalyzedPackage postsim_anal_package;
};
