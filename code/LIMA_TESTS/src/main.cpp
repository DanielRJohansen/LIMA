#include "ForceCorrectness.h"
#include "MDStability.cuh"
#include "MembraneBuilder.h"
#include "MinorPrograms.h"
#include "ElectrostaticsTests.h"
#include "Benchmarks.h"
#include "FileTests.h"
#include "ForcefieldTests.h"
#include "SetupTests.h"

using namespace TestUtils;
using namespace ForceCorrectness;
using namespace TestMDStability;
using namespace StressTesting;
using namespace TestMembraneBuilder;
using namespace TestMinorPrograms;
using namespace ElectrostaticsTests;
using namespace VerletintegrationTesting;
void runAllUnitTests();




int main() {
	try {
		constexpr auto envmode = EnvMode::Full;

		//Lipids::_MakeLipids(true, false);

		//loadAndRunBasicSimulation("DisplayTest", envmode);
		//Display::TestDisplay();
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

		//loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 1.4e-4, 2e-5);
		//loadAndRunBasicSimulation("T4Lysozyme", envmode, 1.15e-4, 2.e-6);



		//const fs::path work_dir = simulations_dir / "test";
		//Lipids::Selection lipids;
		//lipids.emplace_back(Lipids::Select{ "DPPE", work_dir, 30.5 });
		//lipids.emplace_back(Lipids::Select{ "DMPG", work_dir, 39.5 });
		//lipids.emplace_back(Lipids::Select{ "Cholesterol", work_dir, 10 });
		//lipids.emplace_back(Lipids::Select{ "SM18", work_dir, 20 });
		//auto [grofile, topfile] = SimulationBuilder::CreateMembrane(lipids, Float3{ 10.f }, 5.f);
		////simulationBuilder::CreateMembrane(*grofile, *topfile, lipids, 15.f);
		//auto sim = Programs::EnergyMinimizeMax(*grofile, *topfile, true, work_dir, envmode);
		//grofile->printToFile(work_dir / "membrane.gro");
		//topfile->printToFile(work_dir / "membrane.top");

		/*SimParams params{};
		params.n_steps = 1000;
		params.data_logging_interval = 2;
		params.save_trajectory = true;
		Environment env(work_dir, envmode);
		GroFile grofile{ work_dir / "membrane.gro" };
		TopologyFile topfile{ work_dir / "membrane.top" };
		env.CreateSimulation(grofile, topfile, params);
		env.run(true);*/

		//TestBuildmembraneWithCustomlipidAndCustomForcefield(envmode);
		//TestBuildmembraneSmall(envmode, false);
		//loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 9.333e-5, 2e-5);	
		

		//TestLongrangeEsNoLJ(envmode);
		//MakeChargeParticlesSim();
		//TestChargedParticlesVelocityInUniformElectricField(envmode);
		//CoulombForceSanityCheck(envmode);
		//TestElectrostaticsManyParticles(envmode);
		//doPoolBenchmarkES(envmode);
		//TestAttractiveParticlesInteractingWithESandLJ(envmode);

		//Benchmarks::Psome(envmode);

		//TestAllStockholmlipids(envmode);
		
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
	ADD_TEST("TestLimaChosesSameBondparametersAsGromacs", TestLimaChosesSameBondparametersAsGromacs(envmode));

	// Test Setup
	ADD_TEST("TestBoxIsSavedCorrectlyBetweenSimulations", TestBoxIsSavedCorrectlyBetweenSimulations(envmode));	


	ADD_TEST("doPoolBenchmarkES", doPoolBenchmarkES(envmode));
	ADD_TEST("TestElectrostaticsManyParticles", TestElectrostaticsManyParticles(envmode));
	ADD_TEST("TestChargedParticlesVelocityInUniformElectricField", TestChargedParticlesVelocityInUniformElectricField(envmode));
	

	// Programs test
	ADD_TEST("BuildSmallMembrane", TestBuildmembraneSmall(envmode, false));
	ADD_TEST("TestBuildmembraneWithCustomlipidAndCustomForcefield", TestBuildmembraneWithCustomlipidAndCustomForcefield(envmode));
	ADD_TEST("TestAllStockholmlipids", TestAllStockholmlipids(envmode));

	//ADD_TEST("ReorderMoleculeParticles", testReorderMoleculeParticles(envmode));
	ADD_TEST("TestFilesAreCachedAsBinaries", FileTests::TestFilesAreCachedAsBinaries(envmode));

	// Performance test
	//ADD_TEST(testman, "Benchmark Psome", Benchmarks::Psome(envmode));

	// Meta tests
	//doPool50x(EnvMode::Headless);

	// Total test status will print as testman is destructed
}
