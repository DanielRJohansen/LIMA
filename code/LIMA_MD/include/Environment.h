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
	Environment(const string& wf, EnvMode mode);

	void CreateSimulation(string conf_filename, string topol_filename, InputSimParams);

	/// <summary>
	/// Create a simulation that starts from where boxorigin is currently
	/// </summary>
	void CreateSimulation(const Simulation& simulation_src, InputSimParams);

	void run(bool em_variant=false);


	// Return if cannot run
	bool prepareForRun();

	static InputSimParams loadInputSimParams(const std::string& path);
	void renderTrajectory(string trj_path);
	void makeVirtualTrajectory(string trj_path, string waterforce_path);

	// Functions for dev only : TODO move to child whioch inherits all as public
	//const InputSimParams getSimparams();
	std::unique_ptr<Simulation> getSim();
	Simulation* getSimPtr();
	Analyzer::AnalyzedPackage* getAnalyzedPackage();
	SolventBlocksCircularQueue* getSolventBlocks();
	//SolventBlockGrid* getSolventBlocksPrevRef();
	//const std::unique_ptr<SolventBlockGrid> getCurrentSolventblockGrid();
	const std::string& getWorkdir() { return work_folder; }

private:
	void resetEnvironment();

	void setupEmptySimulation(const SimParams&);
	void verifySimulationParameters();			// Constants before doing anything
	void verifyBox();							// Checks wheter the box will break
	
	void postRunEvents();
	void handleStatus(int64_t step, int64_t n_steps);

	// Returns false if display has been closed by user
	bool handleDisplay(const std::vector<Compound>& compounds_host, const BoxParams& boxparams);
	void prepFF();

	void sayHello();

	EnvMode m_mode;

	std::unique_ptr<Display> display;
	int step_at_last_render = 0;

	std::unique_ptr<BoxBuilder> boxbuilder;
	LimaLogger m_logger;

	const std::string work_folder = "";


	std::unique_ptr<Engine> engine;
	std::unique_ptr<Simulation> simulation;

	// TEMP: Cache some constants here before we give ownership to engine. DO NOT READ VOLATILE VALUES FROM THESE
	std::vector<Compound>* compounds = nullptr;
	BoxParams boxparams;

	Analyzer::AnalyzedPackage postsim_anal_package;
#ifdef __linux__
	std::chrono::system_clock::time_point time0;
	std::string main_dir = "/home/lima/Desktop/git_repo/LIMA";
#else
	std::chrono::steady_clock::time_point time0;
	std::string main_dir = "C:/Users/Daniel/git_repo/LIMA/";
#endif




};
