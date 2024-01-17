#pragma once

#include "Bodies.cuh"
#include "DisplayV2.h"
#include "Interface.h"
#include "Engine.cuh"
#include "Analyzer.cuh"
#include "CompoundBuilder.h"
#include "VirtualPathMaker.h"
#include "BoxBuilder.cuh"


// For logging
#include <fstream>
#include <string>
#include <assert.h>  
#include <stdlib.h>
#include <stdio.h>
#include <memory>


#ifndef __linux__
#include <direct.h>
#endif


#include "ForcefieldMaker.h"






class Environment
{

public:
	Environment(const Environment&) = delete;
	Environment(const std::string& wf, EnvMode mode, bool save_output);

	/// <summary>
	/// Create a simulation, and create the necessary files in process, if the defaults
	/// (conf.gro and topol.top and simparams.txt) are not available
	/// </summary>
	void CreateSimulation(float boxsize_nm);

	/// <summary>
	/// Create a simulation from existing files
	/// </summary>
	void CreateSimulation(std::string conf_filename, std::string topol_filename, SimParams);

	// The basic createSim
	void CreateSimulation(const ParsedGroFile&, const ParsedTopologyFile&, const SimParams&);


	/// <summary>
	/// Create a simulation that starts from where boxorigin is currently
	/// </summary>
	void CreateSimulation(Simulation& simulation_src, SimParams);

	/// <summary>
	/// Create a lipid bi-layer in the x-y plane.
	/// </summary>
	/// <param name="carryout_em">Carry out an energy minimization with no boundary condition, 
	/// which ensures all particles are inside the box</param>
	void createMembrane(bool carryout_em=true);





	// Run a standard MD sim
	void run(bool em_variant=false);

	// Intended to be called after a sim run, uses the BoxImage to write new coordinates for the
	// atoms in the input coordinate file.
	ParsedGroFile writeBoxCoordinatesToFile();




	// Return if cannot run
	bool prepareForRun();

	static SimParams loadSimParams(const std::string& path);
	void renderTrajectory(std::string trj_path);
	void makeVirtualTrajectory(std::string trj_path, std::string waterforce_path);

	// Functions for dev only : TODO move to child whioch inherits all as public
	//const InputSimParams getSimparams();
	std::unique_ptr<Simulation> getSim();
	Simulation* getSimPtr();
	Analyzer::AnalyzedPackage* getAnalyzedPackage();
	SolventBlocksCircularQueue* getSolventBlocks();
	//SolventBlockGrid* getSolventBlocksPrevRef();
	//const std::unique_ptr<SolventBlockGrid> getCurrentSolventblockGrid();
	const std::string& getWorkdir() { return work_dir; }

private:
	void resetEnvironment();

	void setupEmptySimulation(const SimParams&);
	void verifySimulationParameters();			// Constants before doing anything
	void verifyBox();							// Checks wheter the box will break
	
	void postRunEvents();
	void handleStatus(int64_t step, int64_t n_steps);

	// Returns false if display has been closed by user
	bool handleDisplay(const std::vector<Compound>& compounds_host, const BoxParams& boxparams);

	void sayHello();

	EnvMode m_mode;
	const bool save_output;

	std::unique_ptr<Display> display;
	int step_at_last_render = 0;
	int step_at_last_update = 0;

	std::unique_ptr<BoxBuilder> boxbuilder;
	LimaLogger m_logger;

	const std::string work_dir = "";	// Main dir of the current simulation


	std::unique_ptr<Engine> engine;
	std::unique_ptr<Simulation> simulation;

	// TEMP: Cache some constants here before we give ownership to engine. DO NOT READ VOLATILE VALUES FROM THESE
	std::vector<Compound>* compounds = nullptr;
	BoxParams boxparams;

	std::unique_ptr<BoxImage> boximage;

	Analyzer::AnalyzedPackage postsim_anal_package;
#ifdef __linux__
	std::chrono::system_clock::time_point time0;
	std::string main_dir = "/opt/LIMA";
#else
	std::chrono::steady_clock::time_point time0;

	// Main folder of the lima package, contains code/lib/resources and more
	const std::string main_dir = "C:/Users/Daniel/git_repo/LIMA/";
#endif




};
