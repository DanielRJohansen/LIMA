#pragma once

//#include "QuantomTypes.cuh"

#include "Bodies.cuh"
//#include "Display.h"
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
	static void prepFF(string conf_path, string topol_path);

	void loadSimParams(std::string path);
	void renderTrajectory(string trj_path);
	void makeVirtualTrajectory(string trj_path, string waterforce_path);

	// Functions for dev only
	SimulationParams* getSimparamRef();
	Simulation* getSim();
	Analyzer::AnalyzedPackage* getAnalyzedPackage();


private:
	void verifySimulationParameters();			// Constants before doing anything
	void verifyBox();							// Checks wheter the box will break
	void prepareForRun();
	//Display* display;
	DisplayV2* display = nullptr;
	//Interface* interface = nullptr;
	
	ForceFieldMaker forcefieldmaker{VerbosityLevel::V1};
	Analyzer analyzer;
	CompoundBuilder* compoundbuilder = nullptr;
	BoxBuilder boxbuilder;

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
#else
	std::chrono::steady_clock::time_point time0;
#endif

};

