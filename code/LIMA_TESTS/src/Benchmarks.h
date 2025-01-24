#include "Programs.h"
#include "TestUtils.h"
#include "TimeIt.h"
#include "MoleculeUtils.h"

namespace Benchmarks {

	using namespace TestUtils;
	namespace fs = std::filesystem;
	
	static void ReadGroFile(EnvMode mode) {
		assert(ENABLE_FILE_CACHING == false);
		TimeIt timer("ReadGroFile", true);
		const fs::path work_dir = simulations_dir / "MembraneAndPsome";
		GroFile psomeGrofile{ work_dir / "molecule/membrane_with_psome.gro" };
		printf("N atoms: %d\n", psomeGrofile.atoms.size());
	}

	static void MembraneWithPsome(EnvMode envmode) {
		//const fs::path work_dir = simulations_dir / "MembraneAndPsome";
		// 
		//bool buildFromScratch = true;
		//if (buildFromScratch) {
		//	Environment env{ work_dir.string(), envmode, false };
		//	const float boxlen = 23.f;

		//	env.CreateSimulation(boxlen);
		//	Lipids::Selection lipidselection;
		//	for (const auto& name : Lipids::Select::valid_lipids) {
		//		lipidselection.emplace_back(Lipids::Select{ name, name == "POPC" ? 50 : 10 });	// 10% of each lipid, except 50% POPC
		//	}

		//	auto [membraneGrofile, membraneTopfile] = Programs::CreateMembrane(work_dir, lipidselection, Float3{boxlen}, boxlen / 4.f, envmode);


		//	GroFile psomeGrofile{ work_dir / "molecule/psome.gro" };
		//	auto psomeTopFile = std::make_shared<TopologyFile>(work_dir / "molecule/psome.top");
		//	Programs::SetMoleculeCenter(psomeGrofile, Float3{ boxlen / 2.f, boxlen / 2.f, 17.f });
		//	Programs::EnergyMinimize(env, psomeGrofile, *psomeTopFile, true, boxlen);
		//	
		//	MDFiles::MergeFiles(*membraneGrofile, *membraneTopfile, psomeGrofile, psomeTopFile);

		//	// Now the "membrane" files also has the psome. Print it to file
		//	membraneGrofile->printToFile(work_dir / "molecule/membrane_with_psome.gro");
		//	membraneTopfile->printToFile(work_dir / "molecule/membrane_with_psome.top");
		//}
		////return;

		//SimParams emparams{ 2000, 20, true, PBC };

		//const GroFile conf{work_dir / "molecule" / "membrane_with_psome.gro"};
		//const TopologyFile topol{work_dir / "molecule" / "membrane_with_psome.top"};
		//const fs::path simpar = work_dir / "sim_params.txt";

		//Environment env{ work_dir, envmode, false};

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
		////auto sim = env->getSim();
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


	static LimaUnittestResult Psome(EnvMode envmode) {
		 //if (envmode== Full)
			// envmode = ConsoleOnly;	// Cant go fast in Full

		const fs::path work_dir = simulations_dir / "psome";
		
		bool em = false;
		if (em) {
			GroFile grofile{ work_dir / "molecule" / "conf.gro" };
			TopologyFile topfile{ work_dir / "molecule" / "topol.top" };

			MoleculeUtils::CenterMolecule(grofile, topfile.GetMoleculeType());
			//SimulationBuilder::SolvateGrofile(grofile);
			auto sim = Programs::EnergyMinimize(grofile, topfile, true, work_dir, envmode, false);
			grofile.printToFile(std::string{ "em.gro" });

			SimAnalysis::PlotPotentialEnergyDistribution(*sim, work_dir, {0,1000, 2000, 3000, 4000 - 1});
		}

		GroFile grofile{ work_dir / "molecule" / "em.gro" };
		TopologyFile topfile{ work_dir / "molecule" / "topol.top" };
		SimParams ip{ work_dir / "sim_params.txt" };
		ip.data_logging_interval = 20;
		ip.dt = 0.5f * FEMTO_TO_NANO;
		ip.enable_electrostatics = true;
		Environment env{ work_dir, envmode };
		env.CreateSimulation(grofile, topfile, ip);
		env.run(false);

		ASSERT(env.getSimPtr()->getStep() == env.getSimPtr()->simparams_host.n_steps, "Simulation did not run fully");

		auto duration = env.simulationTimer->GetTiming();
		const std::chrono::microseconds timePerStep = std::chrono::duration_cast<std::chrono::microseconds>(duration / ip.n_steps);
		const std::chrono::microseconds allowedTimePerStep{ 4000 };

		return LimaUnittestResult { timePerStep < allowedTimePerStep, std::format("Time per step: {} [ys] Allowed: {} [ys]", timePerStep.count(), allowedTimePerStep.count()), envmode!=Headless};
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
		env.run(false);

		ASSERT(env.getSimPtr()->getStep() == env.getSimPtr()->simparams_host.n_steps, "Simulation did not run fully");

		auto duration = env.simulationTimer->GetTiming();
		const std::chrono::microseconds timePerStep = std::chrono::duration_cast<std::chrono::microseconds>(duration / ip.n_steps);
		const std::chrono::microseconds allowedTimePerStep{ 4000 };

		return LimaUnittestResult{ timePerStep < allowedTimePerStep, std::format("Time per step: {} [ys] Allowed: {} [ys]", timePerStep.count(), allowedTimePerStep.count()), envmode != Headless };
	}

	// Returns {avg ms/step, stdDev}
	static std::pair<float, float> Benchmark(const fs::path& dir) {

		const fs::path workDir = simulations_dir / "benchmarking"/dir;
		fs::path topPath, groPath;

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
		TopologyFile topfile(topPath);
		GroFile grofile(groPath);

		SimParams params{ workDir / "../sim_params.txt" };
		params.dt = 1.f * FEMTO_TO_NANO; 		
		Environment env{ workDir , ConsoleOnly };
		//Environment env{ workDir , Full };
		env.CreateSimulation(grofile, topfile, params);
		env.run(false);

		if (env.getSimPtr()->getStep() != env.getSimPtr()->simparams_host.n_steps) {
			throw std::runtime_error("Simulation did not run fully");
		}

		const float meanSteptime = Statistics::Mean(env.avgStepTimes);
		const float stdDev = Statistics::StdDev(env.avgStepTimes);
		printf("Env time: %f [ms/step]\n", std::chrono::duration_cast<std::chrono::milliseconds>(env.simulationTimer->GetTiming()).count() / (float)params.n_steps);
		printf("Average step time: %f [ms] StdDev: %f [ms]\n", meanSteptime, stdDev);

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
		const fs::path work_dir = simulations_dir / "benchmarking"/"stmv";

		GroFile grofile{ work_dir / "conf.gro" };
		TopologyFile topfile{ work_dir /  "topol.top" };
		SimParams ip{};
		ip.n_steps = 1;
		ip.data_logging_interval = 20;
		ip.enable_electrostatics = true;
		Environment env{ work_dir, envmode };
		env.CreateSimulation(grofile, topfile, ip);
		env.prepareForRun();
		const std::chrono::duration<double> elapsedTime = timer.elapsed();
		
		//env.run();

		const std::chrono::duration<double> maxTime{ 30. }; // [s]
		
		return LimaUnittestResult{ elapsedTime < maxTime, std::format("Elapsed time: {:.2f} [s] Allowed: {:.2f} [s]", elapsedTime.count(), maxTime.count()), envmode != Headless };
	}


}