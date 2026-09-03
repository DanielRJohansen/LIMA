#pragma once

#include "Analyzer.cuh"
#include "Bodies.cuh"
#include "TimeIt.h"
#include "MDFiles.h"
#include "Trajectory.h"
#include "LiveEditCommands.h"

#include <memory>
#include <chrono>

class Display;
struct BoxImage;
class Engine;
struct LiveEditData;

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
	void CreateSimulation(const Float3& boxsize_nm);

	/// <summary>
	/// Create a simulation from existing files
	/// </summary>
	void CreateSimulation(const std::string& conf_filename, const std::string& topol_filename, const SimParams&);

	// The basic createSim
	void CreateSimulation(const GroFile&, const TopologyFile&, const SimParams&);

	/// <summary>
	/// Create a simulation that starts from where boxorigin is currently
	/// </summary>
	void CreateSimulation(Simulation& simulation_src, SimParams);

	/// <summary>
	/// Create .gro .top and simparams.txt files in the current directory, and returns them in memory for optional use
	/// </summary>
	std::tuple<GroFile, TopologyFile, SimParams> CreateSimulationFiles(Float3 boxlen);

	// Run a standard MD sim
    /// <returns>Elapsed Engine time in seconds</returns>
    std::chrono::duration<double> run();



	////////////////// LIVE EDIT //////////////////

	/// <summary>
	/// A mode where the user can continously give inputs to the program
	/// </summary>
	void LiveEdit(GroFile& grofile, TopologyFile& topfile);
private:
	void InsertMolecule(LiveEditData*, GroFile& grofile, TopologyFile& topfile, LiveEdit::InsertMolecule& insertionCmd, SimParams simparams);
	void BuildMembrane(LiveEditData*, const LiveEdit::BuildMembrane& cmd, GroFile& grofile, TopologyFile& topfile);
	void HandleMoveMoleculeCommand(LiveEditData*, const LiveEdit::MoveMolecule& newMoveCommand);
	void UpdateForcemask(LiveEditData*, const LiveEdit::AddForcemaskToSelection&);
	void UpdateSelection(LiveEditData*, const LiveEdit::AtomSelected&);
	void UpdateSelection(LiveEditData*, const LiveEdit::SelectAtomsBasedOnQualifier&);
	void UpdateElasticPosition(LiveEditData*, const LiveEdit::ElasticPosition&);
	void EM(LiveEditData*);
	////////////////// ////////////////// ////////////////// 
public:



	/// <summary>
	/// Intended to be called after a sim run, uses the BoxImage to write new coordinates for the
	/// atoms in the input coordinate file.
	/// </summary>
	/// <param name="filename">New name of the file. Defaults to same name as it was as input</param>
	/// <returns></returns>
	GroFile WriteBoxCoordinatesToFile(const std::optional<std::string> filename= "out");
	void WriteBoxCoordinatesToFile(GroFile& grofile, std::optional<int64_t> step=std::nullopt);

	// Returns a vector of forces (in kJ/mol/nm) for each particle, in the order they were in the gro file
	std::vector<Float3> GetForces(int64_t step) const;

	Trajectory WriteSimToTrajectory() const;
	void WriteTrajectoryAsUff(const fs::path& path) const;

	void RenderSimulation();
	

	
	
	void renderTrajectory(std::string trj_path);
	
	void makeVirtualTrajectory(std::string trj_path, std::string waterforce_path);

	// Functions for dev only : TODO move to child whioch inherits all as public
	std::unique_ptr<Simulation> GetSim();
	Simulation* getSimPtr();
	const SimAnalysis::AnalyzedPackage& getAnalyzedPackage();

	std::string getWorkdir() { return workDir.string(); }

	void PrintTiming() const;

	std::chrono::steady_clock::time_point time0;

	const fs::path workDir = "";	// Main dir of the current simulation

	std::optional<TimeIt> simulationTimer;
	std::vector<float> avgStepTimes; // [ms] - averaged over STEP_PER_UPDATE
	std::optional<std::chrono::duration<double>> engineTime;

	std::deque<LiveEdit::Command> liveEditCommandsQueue;	

	SimStatus simStatus{};
	bool forceWriteSimstatusToDisplay = false;

	bool prepareForRun();
private:

	fs::path FixPath(const fs::path& path) const;

	void constexpr verifySimulationParameters();			// Constants before doing anything
	void verifyBox();							// Checks wheter the box will break
	
	void UpdateSimstatus(bool printToConsole, bool alwaysUpdate/*Performance hit*/);

	// Returns false if display has been closed by user
	bool handleDisplay(const BoxParams& boxparams, Display* const display, bool emVariant, bool stepwise);

	void sayHello();

	


	EnvMode m_mode;

	int64_t step_at_last_render = INT64_MIN;

	LimaLogger m_logger;


	std::unique_ptr<Display> display = nullptr;
	std::unique_ptr<Engine> engine = nullptr;
	std::unique_ptr<Simulation> simulation = nullptr;
	std::unique_ptr<BoxImage> boximage = nullptr;

	std::optional<SimAnalysis::AnalyzedPackage> postsim_anal_package;
};
