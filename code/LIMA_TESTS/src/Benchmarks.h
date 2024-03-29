

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

		Programs::CreateMembrane(env, lipidselection, true, boxlen/4.f);

	}

}