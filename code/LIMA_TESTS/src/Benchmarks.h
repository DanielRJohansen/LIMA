

#include "TestUtils.h"
#include "Programs.h"

namespace Benchmarks {

	using namespace TestUtils;
	namespace fs = std::filesystem;


	static void MembraneWithPsome(EnvMode envmode) {
		const fs::path work_dir = simulations_dir + "/MembraneAndPsome";
		Environment env{ work_dir.string(), envmode, false };

		const float boxlen = 20.f;

		env.CreateSimulation(boxlen);
		LipidsSelection lipidselection;
		for (const auto& name : LipidSelect::valid_lipids) {
			lipidselection.emplace_back(LipidSelect{ name, name == "POPC" ? 50 : 10 });	// 10% of each lipid, except 50% POPC
		}

		auto [membraneGrofile, membraneTopfile] = Programs::CreateMembrane(env, lipidselection, true, boxlen/4.f, false);

		ParsedGroFile psomeGrofile{simulations_dir + "/MembraneAndPsome/molecule/psome.gro"};
		Programs::SetMoleculeCenter(psomeGrofile, Float3{ boxlen / 2.f, boxlen / 2.f, 15.f });

		std::unique_ptr<ParsedTopologyFile> psomeTopFile = std::make_unique<ParsedTopologyFile>(simulations_dir + "/MembraneAndPsome/molecule/psome.top");
		MDFiles::MergeFiles(*membraneGrofile, *membraneTopfile, psomeGrofile, std::move(psomeTopFile) );


		// Now the "membrane" files also has the psome. Print it to file
		membraneGrofile->printToFile(simulations_dir + "/MembraneAndPsome/molecule/membrane_with_psome.gro");
		membraneTopfile->printToFile(simulations_dir + "/MembraneAndPsome/molecule/membrane_with_psome.top");
	}

}