#include "Programs.h"
#include "TestUtils.h"
#include "TimeIt.h"
#include "MoleculeUtils.h"
#include "Statistics.h"

namespace Benchmarks {

	using namespace TestUtils;
	namespace fs = std::filesystem;
	constexpr int automatedTestRuns = 3;

	template<typename Duration>
	struct PerformanceBounds {
		Duration min;
		Duration max;
	};
	
	const fs::path TestsDir() {
		return HeavyTestsDir();
	}

	static LimaUnittestResult ToGmxLargeCif(EnvMode envmode) {
		const fs::path input = TestsDir() / "fileconversions" / "3J3Q.cif";
		ASSERT(fs::is_regular_file(input), "Missing ToGmx benchmark input: " + input.string());
		TimeIt timer;
		const auto conversion = Programs::ToGmx(input);
		const auto elapsed = timer.stop();
		const std::chrono::seconds allowedTime{ 10 };

		ASSERT(!conversion.grofile.atoms.empty(), "ToGmx benchmark produced no atoms");
		return LimaUnittestResult{ elapsed < allowedTime,
			std::format("3J3Q.cif elapsed: {:.3f} allowed: {:.3f}",
				std::chrono::duration<double>(elapsed).count(), std::chrono::duration<double>(allowedTime).count()),
			envmode == Full };
	}

	static void ReadGroFile(EnvMode mode) {
		assert(ENABLE_FILE_CACHING == false);
		TimeIt timer("ReadGroFile", true);
		const fs::path workDir = TestsDir() / "MembraneAndPsome";
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

	static LimaUnittestResult Bench(EnvMode envmode, const fs::path& workDir, const GroFile& grofile, const TopologyFile& topfile,
		SimParams ip, PerformanceBounds<std::chrono::microseconds> allowedTimePerStep, std::string name,
		std::optional<int> nSteps = std::nullopt, int nRuns = 1) {

		ip.data_logging_interval = 20;
		//ip.dt = 0.5f * FEMTO_TO_NANO;
		ip.enable_electrostatics = true;
		if (nSteps)
			ip.n_steps = nSteps.value();
		std::vector<std::chrono::microseconds> timesPerStep;
		timesPerStep.reserve(nRuns);
		for (int run = 0; run < nRuns; run++) {
			Environment env{ workDir, EnvMode::Headless /*envmode*/ };
			env.CreateSimulation(grofile, topfile, ip);
			RunOnGpu(env);

			ASSERT(env.getSimPtr()->getStep() == env.getSimPtr()->simParams.n_steps, "Simulation did not run fully");
			const auto duration = env.simulationTimer->GetTiming();
			timesPerStep.push_back(std::chrono::duration_cast<std::chrono::microseconds>(duration / ip.n_steps));
			if (envmode == EnvMode::Full)
				env.PrintTiming();
		}

		const auto [fastest, slowest] = std::minmax_element(timesPerStep.begin(), timesPerStep.end());
		const bool withinBounds = *fastest >= allowedTimePerStep.min && *slowest <= allowedTimePerStep.max;
		return LimaUnittestResult{ withinBounds,
			std::format("({:.3f}-{:.3f}) / ({:.3f}-{:.3f}) [ms/step] ({} runs)",
				fastest->count() / 1000., slowest->count() / 1000., allowedTimePerStep.min.count() / 1000.,
				allowedTimePerStep.max.count() / 1000., nRuns), envmode != Headless };
	}


	static LimaUnittestResult Psome(EnvMode envmode, std::optional<int> nSteps=std::nullopt) {
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

		const fs::path workDir = TestsDir() / "psome";
		GroFile grofile{ workDir / "molecule" / "conf.gro"};
		TopologyFile topfile{ workDir / "molecule" / "topol.top"};
		SimParams ip{ workDir / "sim_params.txt" };
		return Bench(envmode, workDir, grofile, topfile, ip,
			{ std::chrono::microseconds{ 0 }, std::chrono::microseconds{ 4000 } }, "Psome", 30);
	}

	static LimaUnittestResult STMV(EnvMode envmode, int nSteps, int nRuns = 1) {
		const fs::path workDir = TestsDir() / "benchmarking" / "stmv";
		GroFile grofile{ workDir  / "conf.gro" };
		TopologyFile topfile{ workDir  / "topol.top" };
		SimParams ip{ workDir / "sim_params.txt" };
		return Bench(envmode, workDir, grofile, topfile, ip,
			{ std::chrono::microseconds{ 12000 }, std::chrono::microseconds{ 13500 } }, "STMV", nSteps, nRuns);
	}

	static LimaUnittestResult T4(EnvMode envmode, int nSteps=500, int nRuns = 1) {
		const fs::path workDir = TestsDir() / "benchmarking" / "t4";
		GroFile grofile{ workDir  / "conf.gro" };
		TopologyFile topfile{ workDir  / "topol.top" };
		SimParams ip{ workDir / "../sim_params.txt" };
		return Bench(envmode, workDir, grofile, topfile, ip,
			{ std::chrono::microseconds{ 200 }, std::chrono::microseconds{ 300 } }, "T4", nSteps, nRuns);
	}

	static LimaUnittestResult ManyT4(EnvMode envmode) {
		const fs::path workDir  = TestsDir() / "manyt4";		
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
		RunOnGpu(env);

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


		const fs::path workDir = TestsDir() / "benchmarking"/dir;
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
		RunOnGpu(env);

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
		const fs::path workDir = TestsDir() / "benchmarking"/"stmv";

		GroFile grofile{ workDir / "conf.gro" };
		TopologyFile topfile{ workDir /  "topol.top" };
		SimParams ip{};
		ip.n_steps = 1;
		ip.data_logging_interval = 20;
		ip.enable_electrostatics = true;
		Environment env{ workDir, envmode };
		env.CreateSimulation(grofile, topfile, ip);
		WithGpu([&] {
			env.prepareForRun();
			env.ReleaseEngine();
		});
		const std::chrono::duration<double> elapsedTime = timer.elapsed();
		
		const std::chrono::duration<double> maxTime{ 8. }; // [s]
		
		return LimaUnittestResult{ elapsedTime < maxTime, std::format("Elapsed time: {:.2f} [s] Allowed: {:.2f} [s]", elapsedTime.count(), maxTime.count()), envmode != Headless };
	}


}
