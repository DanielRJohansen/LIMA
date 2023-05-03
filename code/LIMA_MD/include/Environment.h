#pragma once

//#include "QuantomTypes.cuh"

#include "LIMA_BASE/include/Bodies.cuh"
//#include "Display.h"
#include "LIMA_MD/include/DisplayV2.h"
#include "LIMA_MD/include/Interface.h"
#include "LIMA_ENGINE/include/Engine.cuh"
#include "LIMA_MD/include/Analyzer.cuh"
#include "LIMA_MD/include/CompoundBuilder.h"
#include "LIMA_MD/include/VirtualPathMaker.h"
#include "LIMA_MD/src/BoxBuilder.cuh"


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


#include "LIMA_FORCEFIELDMAKER/include/ForcefieldMaker.h"







class Environment
{
public:
	Environment();

	void CreateSimulation(string conf_filename, string topol_filename, std::string work_folder);

	void run();
	void postRunEvents();
	void handleStatus(Simulation* simulation);
	void handleDisplay(Simulation* simulation);
	bool handleTermination(Simulation* simulation);
	void prepFF(string conf_path, string topol_path);

	void loadSimParams(const std::string& path);
	void renderTrajectory(string trj_path);
	void makeVirtualTrajectory(string trj_path, string waterforce_path);

	// Functions for dev only : TODO move to child whioch inherits all as public
	SimulationParams* getSimparamRef();
	Simulation* getSim();
	Analyzer::AnalyzedPackage* getAnalyzedPackage();
	std::array<CompoundCoords, MAX_COMPOUNDS>& getCoordarrayRef(std::string selector = "current" /*"current"|"prev"*/);
	SolventBlockGrid* getAllSolventBlocksPrev();


private:
	void verifySimulationParameters();			// Constants before doing anything
	void verifyBox();							// Checks wheter the box will break
	void prepareForRun();

	void sayHello();

	//Display* display;
	Display* display = nullptr;
	//Interface* interface = nullptr;
	
	Forcefield forcefield{VerbosityLevel::V1};
	Analyzer analyzer;
	BoxBuilder* boxbuilder = nullptr;

	std::string work_folder = "";
	SimulationParams sim_params{};

	bool ready_to_run = false;

	// These should be in interface maybe?
	template <typename T>
	void dumpToFile(T* data, uint64_t n_datapoints, string file_path);

	std::unique_ptr<Engine> engine;
	std::unique_ptr<Simulation> simulation;

	Analyzer::AnalyzedPackage postsim_anal_package;
#ifdef __linux__
	std::chrono::system_clock::time_point time0;
	std::string main_dir = "../"
#else
	std::chrono::steady_clock::time_point time0;
	std::string main_dir = "C:/Users/Daniel/git_repo/LIMA/";
#endif




};

