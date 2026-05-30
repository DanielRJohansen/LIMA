#include "Programs.h"
#include "TestUtils.h"
#include "TimeIt.h"
#include "MoleculeUtils.h"
#include "Statistics.h"

namespace Benchmarks {

	using namespace TestUtils;
	namespace fs = std::filesystem;
	
	static void ReadGroFile(EnvMode mode) {
		assert(ENABLE_FILE_CACHING == false);
		TimeIt timer("ReadGroFile", true);
		const fs::path workDir = simulations_dir / "MembraneAndPsome";
		GroFile psomeGrofile{ workDir / "molecule/membrane_with_psome.gro" };
		printf("N atoms: %d\n", psomeGrofile.atoms.size());
	}

	static void MembraneWithPsome(EnvMode envmode) {
		//const fs::path workDir = simulations_dir / "MembraneAndPsome";
		// 
		//bool buildFromScratch = true;
		//if (buildFromScratch) {
		//	Environment env{ workDir.string(), envmode, false };
		//	const float boxlen = 23.f;

		//	env.CreateSimulation(boxlen);
		//	Lipids::Selection lipidselection;
		//	for (const auto& name : Lipids::Select::valid_lipids) {
		//		lipidselection.emplace_back(Lipids::Select{ name, name == "POPC" ? 50 : 10 });	// 10% of each lipid, except 50% POPC
		//	}

		//	auto [membraneGrofile, membraneTopfile] = Programs::CreateMembrane(workDir, lipidselection, Float3{boxlen}, boxlen / 4.f, envmode);


		//	GroFile psomeGrofile{ workDir / "molecule/psome.gro" };
		//	auto psomeTopFile = std::make_shared<TopologyFile>(workDir / "molecule/psome.top");
		//	Programs::SetMoleculeCenter(psomeGrofile, Float3{ boxlen / 2.f, boxlen / 2.f, 17.f });
		//	Programs::EnergyMinimize(env, psomeGrofile, *psomeTopFile, true, boxlen);
		//	
		//	MDFiles::MergeFiles(*membraneGrofile, *membraneTopfile, psomeGrofile, psomeTopFile);

		//	// Now the "membrane" files also has the psome. Print it to file
		//	membraneGrofile->printToFile(workDir / "molecule/membrane_with_psome.gro");
		//	membraneTopfile->printToFile(workDir / "molecule/membrane_with_psome.top");
		//}
		////return;

		//SimParams emparams{ 2000, 20, true, PBC };

		//const GroFile conf{workDir / "molecule" / "membrane_with_psome.gro"};
		//const TopologyFile topol{workDir / "molecule" / "membrane_with_psome.top"};
		//const fs::path simpar = workDir / "sim_params.txt";

		//Environment env{ workDir, envmode, false};

		//const SimParams ip{ simpar.string()};
		//env.CreateSimulation(conf, topol, ip);

		//// Do em
		//env.run(true);
		//

		////// Do sim
		//////InputSimParams simparams{ 100, 2000 };
		////const fs::path work_folder = simulations_dir / folder_name;
		////const std::string simpar_path = work_folder + "/sim_params.txt";
		////SimParams params{ simpar_path };
		////auto sim = env->GetSim();
		////env->CreateSimulation(*sim, params);
		////env->run();
		//////Analyzer::findAndDumpPiecewiseEnergies(*env->getSimPtr(), env->getWorkdir());

		////const auto analytics = env->getAnalyzedPackage();

		////if (envmode != Headless) {
		////	analytics.Print();
		////	LIMA_Print::printMatlabVec("cv", std::vector<float>{ analytics->variance_coefficient});
		////	LIMA_Print::printMatlabVec("energy_gradients", std::vector<float>{ analytics->energy_gradient});
		////}

		////const auto result = evaluateTest({ analytics->variance_coefficient }, max_vc, { analytics->energy_gradient }, max_gradient);
		////const auto status = result.first == true ? true : false;
	}

	static LimaUnittestResult Bench(const fs::path& workDir, const GroFile& grofile, const TopologyFile& topfile, SimParams ip, std::chrono::microseconds allowedTimePerStep, std::string name, std::optional<int> nSteps = std::nullopt) {
		EnvMode envmode = ConsoleOnly;

		ip.data_logging_interval = 20;
		//ip.dt = 0.5f * FEMTO_TO_NANO;
		ip.enable_electrostatics = true;
		if (nSteps)
			ip.n_steps = nSteps.value();
		Environment env{ workDir, envmode };
		env.CreateSimulation(grofile, topfile, ip);
		env.run();

		ASSERT(env.getSimPtr()->getStep() == env.getSimPtr()->simParams.n_steps, "Simulation did not run fully");

		auto duration = env.simulationTimer->GetTiming();
		const std::chrono::microseconds timePerStep = std::chrono::duration_cast<std::chrono::microseconds>(duration / ip.n_steps);
		env.PrintTiming();
		return LimaUnittestResult{ timePerStep < allowedTimePerStep, std::format("{} - Time per step: {} [us] Allowed: {} [us]", name, timePerStep.count(), allowedTimePerStep.count()), envmode != Headless };
	}


	static LimaUnittestResult Psome(std::optional<int> nSteps=std::nullopt) {
		// if (envmode== Full)
		//	 envmode = ConsoleOnly;	// Cant go fast in Full

		//const fs::path workDir = simulations_dir / "psome";
		//
		//bool em = false;
		//if (em) {
		//	GroFile grofile{ workDir / "molecule" / "conf.gro" };
		//	TopologyFile topfile{ workDir / "molecule" / "topol.top" };

		//	MoleculeUtils::CenterMolecule(grofile, topfile.GetMoleculeType());
		//	//SimulationBuilder::SolvateGrofile(grofile);
		//	auto sim = Programs::EnergyMinimize(grofile, topfile, true, workDir, envmode, false);
		//	grofile.printToFile(std::string{ "em.gro" });

		//	SimAnalysis::PlotPotentialEnergyDistribution(*sim, workDir, {0,1000, 2000, 3000, 4000 - 1});
		//}

		//GroFile grofile{ workDir / "molecule" / "em.gro" };
		//TopologyFile topfile{ workDir / "molecule" / "topol.top" };
		//SimParams ip{ workDir / "sim_params.txt" };
		//ip.data_logging_interval = 20;
		//ip.dt = 0.5f * FEMTO_TO_NANO;
		//ip.enable_electrostatics = true;
		//if (nSteps)
		//	ip.n_steps = nSteps.value();
		//Environment env{ workDir, envmode };
		//env.CreateSimulation(grofile, topfile, ip);
		//env.run();

		//ASSERT(env.getSimPtr()->getStep() == env.getSimPtr()->simParams.n_steps, "Simulation did not run fully");

		//auto duration = env.simulationTimer->GetTiming();
		//const std::chrono::microseconds timePerStep = std::chrono::duration_cast<std::chrono::microseconds>(duration / ip.n_steps);
		//const std::chrono::microseconds allowedTimePerStep{ 4000 };

		//return LimaUnittestResult { timePerStep < allowedTimePerStep, std::format("Time per step: {} [us] Allowed: {} [us]", timePerStep.count(), allowedTimePerStep.count()), envmode!=Headless};

		const fs::path workDir = simulations_dir / "psome";
		GroFile grofile{ workDir / "molecule" / "conf.gro"};
		TopologyFile topfile{ workDir / "molecule" / "topol.top"};
		SimParams ip{ workDir / "sim_params.txt" };
		return Bench(workDir, grofile, topfile, ip, std::chrono::microseconds{ 4000 }, "Psome", 30);
	}

	static LimaUnittestResult STMV(int nSteps) {
		const fs::path workDir = simulations_dir / "benchmarking" / "stmv";
		GroFile grofile{ workDir  / "conf.gro" };
		TopologyFile topfile{ workDir  / "topol.top" };
		SimParams ip{ workDir / "sim_params.txt" };
		return Bench(workDir, grofile, topfile, ip, std::chrono::microseconds{ 4500 }, "STMV", nSteps);
	}

	static LimaUnittestResult ManyT4(EnvMode envmode) {
		if (envmode== Full)
		    envmode = ConsoleOnly;	// Cant go fast in Full

		const fs::path workDir  = simulations_dir / "manyt4";		
		TopologyFile topfile(workDir / "t4_many.top");
		GroFile grofile(workDir / "t4_many_em.gro");

		//bool em = false;
		//if (em) {
		//	GroFile grofile{ workDir  / "molecule" / "conf.gro" };
		//	TopologyFile topfile{ workDir  / "molecule" / "topol.top" };

		//	MoleculeUtils::CenterMolecule(grofile, topfile.GetMoleculeType());
		//	//SimulationBuilder::SolvateGrofile(grofile);
		//	auto sim = Programs::EnergyMinimize(grofile, topfile, true, workDir , envmode, false);
		//	grofile.printToFile(std::string{ "em.gro" });

		//	SimAnalysis::PlotPotentialEnergyDistribution(*sim, workDir , { 0,1000, 2000, 3000, 4000 - 1 });
		//}
		SimParams ip{ workDir  / "sim_params.txt" };
		ip.data_logging_interval = 50;
		ip.dt = 1.f * FEMTO_TO_NANO;
		ip.n_steps = 4000;
		Environment env{ workDir , envmode };
		env.CreateSimulation(grofile, topfile, ip);
		env.run();

		ASSERT(env.getSimPtr()->getStep() == env.getSimPtr()->simParams.n_steps, "Simulation did not run fully");

		auto duration = env.simulationTimer->GetTiming();
		const std::chrono::microseconds timePerStep = std::chrono::duration_cast<std::chrono::microseconds>(duration / ip.n_steps);
		const std::chrono::microseconds allowedTimePerStep{ 4000 };

		return LimaUnittestResult{ timePerStep < allowedTimePerStep, std::format("Time per step: {} [us] Allowed: {} [us]", timePerStep.count(), allowedTimePerStep.count()), envmode != Headless };
	}

	// Returns {avg ms/step, stdDev}
	static std::pair<float, float> Benchmark(const fs::path& dir, std::optional<std::string> name = std::nullopt, std::optional<int> nSteps = std::nullopt) {
		
		if (!IS_FAST_MODE) {
			TestUtils::setConsoleTextColorYellow();
			printf("Warning: Benchmarking with debug mode enabled. Results may be significantly slower than expected.\n");
			TestUtils::setConsoleTextColorDefault();
		}

		if (!ALL_PHYSICS_ENABLED) {
			TestUtils::setConsoleTextColorYellow();
			printf("Warning: Benchmarking with some physics disabled. Results may be significantly faster than expected.\n");
			TestUtils::setConsoleTextColorDefault();
		}


		const fs::path workDir = simulations_dir / "benchmarking"/dir;
		fs::path topPath, groPath;

		if (name) {
			groPath = workDir / (*name + ".gro");
			topPath = workDir / (*name + ".top");
		}
		else {
			for (const auto& entry : fs::directory_iterator(workDir)) {
				auto ext = entry.path().extension();
				if (ext == ".top") {
					if (!topPath.empty()) throw std::runtime_error("Multiple .top files found");
					topPath = entry.path();
				}
				else if (ext == ".gro") {
					if (!groPath.empty()) throw std::runtime_error("Multiple .gro files found");
					groPath = entry.path();
				}
			}
		}
		TopologyFile topfile(topPath);
		GroFile grofile(groPath);

		fs::path spPath = fs::exists(workDir / "sim_params.txt") ? workDir / "sim_params.txt" : workDir / ".." / "sim_params.txt";
		SimParams params{ spPath };
		if (nSteps) 
			params.n_steps = *nSteps;
		//params.dt = 1.f * FEMTO_TO_NANO; 		
		Environment env{ workDir , ConsoleOnly };
		//Environment env{ workDir , Full };
		env.CreateSimulation(grofile, topfile, params);
		env.run();

		if (env.getSimPtr()->getStep() != env.getSimPtr()->simParams.n_steps) {
			throw std::runtime_error("Simulation did not run fully");
		}

		const float meanSteptime = Statistics::Mean(env.avgStepTimes);
		const float stdDev = Statistics::StdDev(env.avgStepTimes);
		printf("Env time: %f [ms/step]\n", std::chrono::duration_cast<std::chrono::milliseconds>(env.simulationTimer->GetTiming()).count() / (float)params.n_steps);
		printf("Average step time: %f [ms] StdDev: %f [ms]\n", meanSteptime, stdDev);

		env.PrintTiming();

		return { meanSteptime, stdDev};
	}


	static void Benchmark(const std::vector<fs::path>& dirs) {
		// Header for the output table
		std::cout << std::left << std::setw(20) << "Directory"
			<< std::setw(15) << "Avg Time (ms)"
			<< std::setw(15) << "Std Dev (ms)" << std::endl;
		std::cout << std::string(50, '-') << std::endl;

		for (const auto& dir : dirs) {
			auto [meanTime, stdDev] = Benchmark(dir);

			// Format the output
			std::cout << std::left << std::setw(20) << dir.filename().string()
				<< std::setw(15) << std::fixed << std::setprecision(2) << meanTime
				<< std::setw(15) << std::fixed << std::setprecision(2) << stdDev
				<< std::endl;
		}
	}	

	static LimaUnittestResult PrepareSimulation_stmv(EnvMode envmode) {
		TimeIt timer("Load Sim");
		const fs::path workDir = simulations_dir / "benchmarking"/"stmv";

		GroFile grofile{ workDir / "conf.gro" };
		TopologyFile topfile{ workDir /  "topol.top" };
		SimParams ip{};
		ip.n_steps = 1;
		ip.data_logging_interval = 20;
		ip.enable_electrostatics = true;
		Environment env{ workDir, envmode };
		env.CreateSimulation(grofile, topfile, ip);
		env.prepareForRun();
		const std::chrono::duration<double> elapsedTime = timer.elapsed();
		
		const std::chrono::duration<double> maxTime{ 8. }; // [s]
		
		return LimaUnittestResult{ elapsedTime < maxTime, std::format("Elapsed time: {:.2f} [s] Allowed: {:.2f} [s]", elapsedTime.count(), maxTime.count()), envmode != Headless };
	}


}