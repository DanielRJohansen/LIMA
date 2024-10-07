#include "Programs.h"
#include "TestUtils.h"
#include "TimeIt.h"

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
		////	Analyzer::printEnergy(analytics);
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
		float boxlen = 23.f;
		Environment env{ work_dir, envmode, false };
		
		bool em = false;
		if (em) {
			GroFile grofile{ work_dir / "molecule" / "conf.gro" };
			grofile.box_size = Float3{ boxlen, boxlen, boxlen };
			TopologyFile topfile{ work_dir / "molecule" / "topol.top" };

			Programs::SetMoleculeCenter(grofile, Float3{ boxlen / 2.f, boxlen / 2.f, boxlen / 2.f });
			SimulationBuilder::SolvateGrofile(grofile);
			Programs::EnergyMinimize(env, grofile, topfile, true, boxlen);

			SimAnalysis::PlotPotentialEnergyDistribution(*env.getSimPtr(), env.work_dir, { 0,1000, 2000, 3000, 4000 - 1 });

			GroFile emGro = env.WriteBoxCoordinatesToFile();
			emGro.printToFile(std::string{ "em.gro" });
		}

		GroFile grofile{ work_dir / "molecule" / "em.gro" };
		TopologyFile topfile{ work_dir / "molecule" / "topol.top" };
		SimParams ip{ work_dir / "sim_params.txt" };
		ip.data_logging_interval = 20;
		ip.dt = 50;
		ip.enable_electrostatics = true;
		env.CreateSimulation(grofile, topfile, ip);
		env.run(false);

		ASSERT(env.getSimPtr()->getStep() == env.getSimPtr()->simparams_host.n_steps, "Simulation did not run fully");

		auto duration = env.simulationTimer->GetTiming();
		const std::chrono::microseconds timePerStep = std::chrono::duration_cast<std::chrono::microseconds>(duration / ip.n_steps);
		const std::chrono::microseconds allowedTimePerStep{ 4000 };

		return LimaUnittestResult { timePerStep < allowedTimePerStep, std::format("Time per step: {} [ys] Allowed: {} [ys]", timePerStep.count(), allowedTimePerStep.count()), envmode!=Headless};
	}
}