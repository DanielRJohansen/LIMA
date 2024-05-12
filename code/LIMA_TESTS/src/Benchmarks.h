

#include "TestUtils.h"
#include "Programs.h"

namespace Benchmarks {

	using namespace TestUtils;
	namespace fs = std::filesystem;


	static void MembraneWithPsome(EnvMode envmode) {
		const fs::path work_dir = simulations_dir + "/MembraneAndPsome";

		bool buildFromScratch = false;
		if (buildFromScratch) {
			Environment env{ work_dir.string(), envmode, false };
			const float boxlen = 20.f;

			env.CreateSimulation(boxlen);
			LipidsSelection lipidselection;
			for (const auto& name : LipidSelect::valid_lipids) {
				lipidselection.emplace_back(LipidSelect{ name, name == "POPC" ? 50 : 10 });	// 10% of each lipid, except 50% POPC
			}

			auto [membraneGrofile, membraneTopfile] = Programs::CreateMembrane(env, lipidselection, true, boxlen / 4.f, false);

			ParsedGroFile psomeGrofile{ work_dir / "molecule/psome.gro" };
			Programs::SetMoleculeCenter(psomeGrofile, Float3{ boxlen / 2.f, boxlen / 2.f, 15.f });

			std::unique_ptr<ParsedTopologyFile> psomeTopFile = std::make_unique<ParsedTopologyFile>(work_dir / "molecule/psome.top");
			MDFiles::MergeFiles(*membraneGrofile, *membraneTopfile, psomeGrofile, std::move(psomeTopFile));

			// Now the "membrane" files also has the psome. Print it to file
			membraneGrofile->printToFile(work_dir / "molecule/membrane_with_psome.gro");
			membraneTopfile->printToFile(work_dir / "molecule/membrane_with_psome.top");
		}


		SimParams emparams{ 2000, 20, true, PBC };

		const fs::path conf = work_dir / "molecule" / "membrane_with_psome.gro";
		const fs::path topol = work_dir / "molecule" / "membrane_with_psome.top";
		const fs::path simpar = work_dir / "sim_params.txt";

		Environment env{ work_dir, envmode, false};

		const SimParams ip{ simpar.string()};
		env.CreateSimulation(conf, topol, ip);

		// Do em
		//env.run(true);
		

		//// Do sim
		////InputSimParams simparams{ 100, 2000 };
		//const fs::path work_folder = simulations_dir / folder_name;
		//const std::string simpar_path = work_folder + "/sim_params.txt";
		//SimParams params{ simpar_path };
		//auto sim = env->getSim();
		//env->CreateSimulation(*sim, params);
		//env->run();
		////Analyzer::findAndDumpPiecewiseEnergies(*env->getSimPtr(), env->getWorkdir());

		//const auto analytics = env->getAnalyzedPackage();

		//if (envmode != Headless) {
		//	Analyzer::printEnergy(analytics);
		//	LIMA_Print::printMatlabVec("cv", std::vector<float>{ analytics->variance_coefficient});
		//	LIMA_Print::printMatlabVec("energy_gradients", std::vector<float>{ analytics->energy_gradient});
		//}

		//const auto result = evaluateTest({ analytics->variance_coefficient }, max_vc, { analytics->energy_gradient }, max_gradient);
		//const auto status = result.first == true ? LimaUnittestResult::SUCCESS : LimaUnittestResult::FAIL;
	}

}