#include "ForceCorrectness.h"
#include "MDStability.cuh"
#include "EnergyMinimization.cuh"
#include "MembraneBuilder.h"
#include "MinorPrograms.h"


using namespace TestUtils;
using namespace ForceCorrectness;
using namespace TestMDStability;
using namespace StressTesting;
using namespace TestMembraneBuilder;
using namespace TestMinorPrograms;
void runAllUnitTests();

//void makeLipids() {
// Currently working POPC, POPE DMPC, chol, DPPC, DOPC, 
// SM16 almost works, but H11 is only defined in Slipids_2020, which i need to integrate into the forcefieldmaker, but thats for later
//	std::vector<std::string> targets = { "DOPC" };
////std::vector<std::string> targets = { "DAPC", "DPPC", "DSPC", "POPC", "DOPC", "DLPC", "DMPC", "DUPC", "DIPC", "DPPS", "DOPS", "ENET", "ENTF", "ENTW", "PDPC", "PiLPC", "POPE", "POPG", "POPS", "SAPC", "SDPC", "SIET", "SOPC", "choleterol", "SM16" };
////const string target = "DAPC";
//for (auto target : targets) {
//	try {
//		const std::string to_folder = "C:/Users/Daniel/git_repo/LIMA/resources/Lipids/" + target + "/";
//		const std::string from_folder = "C:/Users/Daniel/git_repo/LIMA/resources/Lipids/" + target + "/";
//		LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(to_folder, from_folder, target);
//	}
//	catch (...) {
//		// do nothing
//	}
//}
//}






int main() {
	try {
		constexpr auto envmode = EnvMode::Full;


		//loadAndRunBasicSimulation("DisplayTest", envmode);

		//doPoolBenchmark(envmode);			// Two 1-particle molecules colliding
		//doPoolCompSolBenchmark(envmode);	// One 1-particle molecule colliding with 1 solvent

		//doSinglebondBenchmark(envmode);
		//doAnglebondBenchmark(envmode);
		//doDihedralbondBenchmark(envmode);
		//TestUtils::loadAndRunBasicSimulation("torsion2", envmode, 0.0002);
		//doImproperDihedralBenchmark(envmode);
		//TestUtils::loadAndRunBasicSimulation("improper", envmode, 7e-5, 2.3e-7);

		//doMethionineBenchmark(envmode);
		//doPhenylalanineBenchmark(envmode);
		//TestUtils::loadAndRunBasicSimulation("TenSolvents", envmode, 0.0004, 1.2e-6);

		//doEightResiduesNoSolvent(envmode);
		//loadAndRunBasicSimulation("Solventsonly", envmode, 4e-6f);


		//loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 2.240e-5, 2e-5);
		//loadAndRunBasicSimulation("T4Lysozyme", envmode, 9.16e-5, 2.9e-7);

		//loadAndRunBasicSimulation("manyt4", envmode, 1.6e-3);
		//loadAndRunBasicSimulation("psome", envmode, 7.6e-5, 1.1e-6);
		//doPool50x(EnvMode::Headless);
	


		
		/*const fs::path work_dir = simulations_dir + "/test";
		Environment env{ work_dir.string(), envmode, false };
		env.CreateSimulation(20.f);
		LipidsSelection lipids;
		lipids.emplace_back(LipidSelect{ "POPC", 50 });
		lipids.emplace_back(LipidSelect{ "DMPC", 40 });
		lipids.emplace_back(LipidSelect{ "cholesterol", 10 });
		env.createMembrane(lipids, true);*/



		//testReorderMoleculeParticles(envmode);
		//testBuildmembraneSmall(envmode, false);
		//loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 4.4e-5, 2e-5);		

		runAllUnitTests();
	}
	catch (std::runtime_error ex) {
		std::cerr << "Caught runtime_error: " << ex.what() << std::endl;
	}
	catch (const std::exception& ex) {
		std::cerr << "Caught exception: " << ex.what() << std::endl;
	}
	catch (...) {
		std::cerr << "Caught unnamed exception";
	}

	return 0;
}

#define ADD_TEST(testman, description, execution_function) \
    testman.addTest(std::make_unique<LimaUnittest>(LimaUnittest{ description, [](){ return execution_function;} }))

// Runs all unit tests with the fastest/crucial ones first
void runAllUnitTests() {
	LimaUnittestManager testman;
	constexpr auto envmode = EnvMode::Headless;

	
	// Singled out forces test
	ADD_TEST(testman, "doPoolBenchmark", doPoolBenchmark(envmode));
	ADD_TEST(testman, "doPoolCompSolBenchmark", doPoolCompSolBenchmark(envmode));
	ADD_TEST(testman, "doSinglebondBenchmark", doSinglebondBenchmark(envmode));
	ADD_TEST(testman, "doAnglebondBenchmark", doAnglebondBenchmark(envmode));
	ADD_TEST(testman, "doDihedralbondBenchmark", doDihedralbondBenchmark(envmode));
	ADD_TEST(testman, "Dihedral_exaggerated", TestUtils::loadAndRunBasicSimulation("Dihedralbond2", envmode, 2e-4, 2.2e-7));
	ADD_TEST(testman, "doImproperDihedralBenchmark", doImproperDihedralBenchmark(envmode, 4.3e-3));
	ADD_TEST(testman, "Improper_exaggerated_scaled-up", TestUtils::loadAndRunBasicSimulation("Improperbond2", envmode, 7e-5, 2.9e-7));


	// Smaller compound tests
	ADD_TEST(testman, "doMethionineBenchmark", TestUtils::loadAndRunBasicSimulation("Met", envmode, 4.1e-4, 2e-6));
	ADD_TEST(testman, "doPhenylalanineBenchmark", TestUtils::loadAndRunBasicSimulation("Phe", envmode, 3.77e-4f, 8e-8f););
	ADD_TEST(testman, "TenSolvents", TestUtils::loadAndRunBasicSimulation("TenSolvents", envmode, 7.3e-6, 1.2e-6));
	ADD_TEST(testman, "doEightResiduesNoSolvent", doEightResiduesNoSolvent(envmode));


	// Larger tests
	ADD_TEST(testman, "SolventBenchmark", loadAndRunBasicSimulation("Solventsonly", envmode, 2.85e-6f, 1.01e-7));
	ADD_TEST(testman, "T4Lysozyme", loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 4.9e-5, 2e-5));

	// Programs test
	ADD_TEST(testman, "BuildSmallMembrane", testBuildmembraneSmall(envmode, false));
	ADD_TEST(testman, "ReorderMoleculeParticles", testReorderMoleculeParticles(envmode));



	// Meta tests
	//doPool50x(EnvMode::Headless);

	// Total test status will print as testman is destructed
}
