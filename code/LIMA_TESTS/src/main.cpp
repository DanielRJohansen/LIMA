#include "ForceCorrectness.h"
#include "MDStability.cuh"
#include "MembraneBuilder.h"
#include "MinorPrograms.h"
#include "ElectrostaticsTests.h"
#include "Benchmarks.h"
#include "FileTests.h"
#include "ForcefieldTests.h"

using namespace TestUtils;
using namespace ForceCorrectness;
using namespace TestMDStability;
using namespace StressTesting;
using namespace TestMembraneBuilder;
using namespace TestMinorPrograms;
using namespace ElectrostaticsTests;
using namespace VerletintegrationTesting;
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

#include <ostream>



int main() {
	try {
		constexpr auto envmode = EnvMode::Full;


		//loadAndRunBasicSimulation("DisplayTest", envmode);

		//doPoolBenchmark(envmode);			// Two 1-particle molecules colliding
		//loadAndRunBasicSimulation("PoolElectrostatic", envmode);
		//doPoolCompSolBenchmark(envmode);	// One 1-particle molecule colliding with 1 solvent

		//SinglebondForceAndPotentialSanityCheck(envmode);
		//SinglebondOscillationTest(envmode);
		//doSinglebondBenchmark(envmode);
		//doAnglebondBenchmark(envmode);
		//doDihedralbondBenchmark(envmode);
		//TestUtils::loadAndRunBasicSimulation("Dihedralbond2", envmode, 0.0002);
		//doImproperDihedralBenchmark(envmode);
		//TestUtils::loadAndRunBasicSimulation("improper", envmode, 7e-5, 2.3e-7);

		//TestUtils::loadAndRunBasicSimulation("Met", envmode, 5.6e-4, 2e-6);
		//loadAndEMAndRunBasicSimulation("Met", envmode, 4.1e-4, 2e-6);
		//TestUtils::loadAndRunBasicSimulation("Phe", envmode, 4.1e-4, 2e-6);
		//doPhenylalanineBenchmark(envmode);
		//TestUtils::loadAndRunBasicSimulation("TenSolvents", envmode, 0.0004, 1.2e-6);

		//doEightResiduesNoSolvent(envmode);
		//loadAndRunBasicSimulation("Solventsonly", envmode, 2.85e-6f, 1.1e-7);


		//loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 4.9e-5, 2e-5);
		//loadAndRunBasicSimulation("T4Lysozyme", envmode, 1.15e-4, 2.e-6);

		//loadAndRunBasicSimulation("manyt4", envmode, 1.6e-3);

		//SimParams params;
		////params.enable_electrostatics = true;
		//params.n_steps = 5000;
		//loadAndRunBasicSimulation("psome", envmode, 1.5e-4, 1.1e-6);
		// 
		//doPool50x(EnvMode::Headless);
	

		//TestIntegration(envmode);

		//Benchmarks::ReadGroFile(envmode);
		//Benchmarks::MembraneWithPsome(envmode);



		//const fs::path work_dir = simulations_dir + "/test";
		//Environment env{ work_dir.string(), envmode, false };
		//env.CreateSimulation(25.f);
		//LipidsSelection lipids;
		//lipids.emplace_back(LipidSelect{ "POPC", 50 });
		//lipids.emplace_back(LipidSelect{ "DMPC", 40 });
		//lipids.emplace_back(LipidSelect{ "cholesterol", 10 });		//Programs::CreateMembrane(env, lipids, true, 3.5f, true);


		//FileTests::TestFilesAreCachedAsBinaries(Headless);


		
	/*	GroFile grofile;
		grofile.box_size = Float3{ 5.f };
		TopologyFile topfile;
		Programs::MakeLipidVesicle(grofile, topfile, { {"POPC", 10}, {"cholesterol", 30}, {"DMPC", 60} });*/









		//testReorderMoleculeParticles(envmode);
		//testBuildmembraneSmall(envmode, false);
		//loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 9.333e-5, 2e-5);	
		

		//TestLongrangeEsNoLJ(envmode);
		//MakeChargeParticlesSim();
		//TestChargedParticlesVelocityInUniformElectricField(envmode);
		//CoulombForceSanityCheck(envmode);
		//TestElectrostaticsManyParticles(envmode);
		//doPoolBenchmarkES(envmode);
		//TestAttractiveParticlesInteractingWithESandLJ(envmode);

		//Benchmarks::Psome(envmode);

		
		//TestLimaChosesSameBondparametersAsGromacs(envmode);
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

#define ADD_TEST(description, execution_function) \
    testman.addTest(std::make_unique<LimaUnittest>(LimaUnittest{ description, [](){ return execution_function;} }))

// Runs all unit tests with the fastest/crucial ones first
void runAllUnitTests() {
	LimaUnittestManager testman;
	constexpr auto envmode = EnvMode::Headless;

	
	// Isolated forces sanity checks
	ADD_TEST("SinglebondForceAndPotentialSanityCheck", SinglebondForceAndPotentialSanityCheck(envmode));
	ADD_TEST("SinglebondOscillationTest", SinglebondOscillationTest(envmode));
	ADD_TEST("TestIntegration", TestIntegration(envmode));

	// Stability tests
	ADD_TEST("doPoolBenchmark", doPoolBenchmark(envmode));
	ADD_TEST("doPoolCompSolBenchmark", doPoolCompSolBenchmark(envmode));
	ADD_TEST("doSinglebondBenchmark", doSinglebondBenchmark(envmode));
	ADD_TEST("doAnglebondBenchmark", doAnglebondBenchmark(envmode));
	ADD_TEST("doDihedralbondBenchmark", doDihedralbondBenchmark(envmode));
	ADD_TEST("doImproperDihedralBenchmark", doImproperDihedralBenchmark(envmode));


	// Smaller compound tests
	ADD_TEST("doMethionineBenchmark", TestUtils::loadAndRunBasicSimulation("Met", envmode, 5.6e-4, 2e-6));
	ADD_TEST("TenSolvents", TestUtils::loadAndRunBasicSimulation("TenSolvents", envmode, 7.3e-6, 1.2e-6));
	ADD_TEST("doEightResiduesNoSolvent", doEightResiduesNoSolvent(envmode));


	// Larger tests
	ADD_TEST("SolventBenchmark", loadAndRunBasicSimulation("Solventsonly", envmode, 2.85e-6f, 1.1e-7));
	ADD_TEST("T4Lysozyme", loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 1.4e-4, 2e-5));

	// Electrostatics
	ADD_TEST("CoulombForceSanityCheck", CoulombForceSanityCheck(envmode));
	
	// Test Forcefield and compoundbuilder
	ADD_TEST("TestLimaChosesSameBondparametersAsGromacs", CoulombForceSanityCheck(envmode));


	ADD_TEST("doPoolBenchmarkES", doPoolBenchmarkES(envmode));
	ADD_TEST("TestElectrostaticsManyParticles", TestElectrostaticsManyParticles(envmode));
	ADD_TEST("TestChargedParticlesVelocityInUniformElectricField", TestChargedParticlesVelocityInUniformElectricField(envmode));
	

	// Programs test
	ADD_TEST("BuildSmallMembrane", testBuildmembraneSmall(envmode, false));
	ADD_TEST("ReorderMoleculeParticles", testReorderMoleculeParticles(envmode));
	ADD_TEST("TestFilesAreCachedAsBinaries", FileTests::TestFilesAreCachedAsBinaries(envmode));
	ADD_TEST("TestLimaChosesSameBondparametersAsGromacs", TestLimaChosesSameBondparametersAsGromacs(envmode));

	// Performance test
	//ADD_TEST(testman, "Benchmark Psome", Benchmarks::Psome(envmode));

	// Meta tests
	//doPool50x(EnvMode::Headless);

	// Total test status will print as testman is destructed
}
