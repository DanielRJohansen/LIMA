#include "ForceCorrectness.h"
#include "MDStability.h"
#include "MembraneBuilder.h"
#include "MinorPrograms.h"
#include "ElectrostaticsTests.h"
#include "Benchmarks.h"
#include "FileTests.h"
#include "ForcefieldTests.h"
#include "SetupTests.h"
#include "Userinterface.h"
#include "Display.h"
#include "ForceComparisons.h"
#include "ProgramsTests.h"
#include "AlgorithmTests.h"


using namespace TestUtils;
using namespace ForceCorrectness;
using namespace TestMDStability;
using namespace StressTesting;
using namespace TestMembraneBuilder;
using namespace TestMinorPrograms;
using namespace ElectrostaticsTests;
using namespace VerletintegrationTesting;

void RunAllUnitTests();

void TestDisplayT4() {
	auto env = TestUtils::basicSetup("T4Lysozyme", std::nullopt, EnvMode::Full);
	Display display{};

	auto& box = env->getSimPtr()->box;

	std::vector<Float3> positions;

	for (auto pc : box->persistentClusters) {
		for (auto pqd : pc.pqd)
			if (pqd.Valid())
				positions.push_back(pqd.position);
	}


	display.Render(std::make_unique<Rendering::SimulationTask>(box->persistentClusters, box->persistentClustersMetadata, box->boxparams), true);
	display.Render(std::make_unique<Rendering::SimulationTaskUpdate>(positions.data(), nullptr, SimStatus{}), true);
}

void LiveEditTest() {
	Environment env({ R"(C:\Users\Daniel\git_repo\LIMA_data\LiveEditTest)" }, EnvMode::Full);
	auto [grofile, topfile, simparams] = env.CreateSimulationFiles(Float3(13.f));
	env.CreateSimulation(grofile, topfile, simparams);

	//// TODO: Being able to set this is like super dangerous, and the ff is then not parsed... Always need a file reset..
	topfile.forcefieldInclude = TopologyFile::ForcefieldInclude("combined/forcefield.itp");
	topfile.printToFile();
	topfile = TopologyFile{ topfile.path };


	//std::vector<std::tuple<std::string, double>> lipids = { {"DPPE", 100.}};
	////std::vector<std::tuple<std::string, double>> lipids = { {"DPPE", 30.5}, {"DMPG", 39.5}, {"cholesterol", 10}, {"SM18", 20} };
	//env.liveEditCommandsQueue.push_back(LiveEdit::BuildMembrane{ lipids, MembraneGeometry::Plane{ 3.f } });
	//env.liveEditCommandsQueue.push_back(LiveEdit::SelectAtomsBasedOnQualifier{ LiveEdit::SelectAtomsBasedOnQualifier::Qualifier::All });
	//env.liveEditCommandsQueue.push_back(LiveEdit::ElasticPosition{ false, false, true });
	//env.liveEditCommandsQueue.push_back(LiveEdit::InsertMolecule{ "t4/conf.gro", "t4/topol_T4.itp" });
	

	env.liveEditCommandsQueue.push_back(LiveEdit::InsertMolecule{ "t4/conf.gro", "t4/topol_T4.itp" });
	env.liveEditCommandsQueue.push_back(LiveEdit::SelectAtomsBasedOnQualifier{ LiveEdit::SelectAtomsBasedOnQualifier::Qualifier::All });
	env.liveEditCommandsQueue.push_back(LiveEdit::ElasticPosition{ true, false, false });

	env.LiveEdit(grofile, topfile);
}

void BuildCellTest() {
	Environment env({ R"(C:\Users\Daniel\git_repo\LIMA_data\LiveEditTest)" }, EnvMode::Full);
	auto [grofile, topfile, simparams] = env.CreateSimulationFiles(Float3(40.f));
	env.CreateSimulation(grofile, topfile, simparams);

	//// TODO: Being able to set this is like super dangerous, and the ff is then not parsed... Always need a file reset..
	topfile.forcefieldInclude = TopologyFile::ForcefieldInclude("combined/forcefield.itp");
	topfile.printToFile();
	topfile = TopologyFile{ topfile.path };


	//std::vector<std::tuple<std::string, double>> lipids = { {"DPPE", 100.}};
	std::vector<std::tuple<std::string, double>> lipids = { {"DPPE", 30.5}, {"DMPG", 39.5}, {"cholesterol", 10}, {"SM18", 20} };
	//env.liveEditCommandsQueue.push_back(LiveEdit::BuildMembrane{ lipids, MembraneGeometry::Plane{ 4.f } });
	env.liveEditCommandsQueue.push_back(LiveEdit::BuildMembrane{ lipids, MembraneGeometry::Sphere{Float3{20,20,20  }, 15.f }});

	env.LiveEdit(grofile, topfile);
}

int main() {
	try {
		constexpr auto envmode = EnvMode::Full;

		//TestDisplayT4();
		//LiveEditTest();
		BuildCellTest();


		//loadAndRunBasicSimulation("Singleatom", envmode);

		//PlotPmePotAsFactorOfDistance(envmode);
		//TestConsistentEnergyWhenGoingFromLresToSres(envmode);
		//TestLongrangeEsNoLJTwoParticles(envmode);
		//TestLongrangeEsNoLJManyParticles(envmode);
		//Lipids::_MakeLipids(true, false);
		//PairbondForceAndPotentialSanityCheck(envmode);
		//loadAndRunBasicSimulation("DisplayTest", envmode);
		//Display::TestDisplay();
		//TestDisplayT4();
		//doPoolBenchmark(envmode);			// Two 1-particle molecules colliding
		//loadAndRunBasicSimulation("PoolElectrostatic", envmode);
		//doPoolCompSolBenchmark(envmode);	// One 1-particle molecule colliding with 1 solvent
		//SinglebondForceAndPotentialSanityCheck(envmode);
		//UreyBradleyForceAndPotentialSanityCheck(envmode);
		//SinglebondOscillationTest(envmode);
		//doSinglebondBenchmark(envmode);
		//doAnglebondBenchmark(envmode);
		//doDihedralbondBenchmark(envmode);
		//loadAndRunBasicSimulation("SinglebondDaisychained", envmode, 0.0002);
		//TestUtils::loadAndRunBasicSimulation("Dihedralbond2", envmode, 0.0002);
		//doImproperDihedralBenchmark(envmode);
		//TestUtils::loadAndRunBasicSimulation("improper", envmode, 7e-5, 2.3e-7);
		//TestUtils::loadAndRunBasicSimulation("Met", envmode, 6.3e-4, 2e-6);
		//loadAndEMAndRunBasicSimulation("Met", envmode, 4.1e-4, 2e-6);
		//TestUtils::loadAndRunBasicSimulation("Phe", envmode, 4.1e-4, 2e-6);
		//doPhenylalanineBenchmark(envmode);
		//doEightResiduesNoSolvent(envmode);
		//loadAndRunBasicSimulation("Solvents", envmode, 5.85e-6f, 1.1e-7);
				//TestLongrangeEsNoLJ(envmode);
		//MakeChargeParticlesSim();
		//TestChargedParticlesVelocityInUniformElectricField(envmode);
		//CoulombForceSanityCheck(envmode);
		//TestElectrostaticsManyParticles(envmode);
		//doPoolBenchmarkES(envmode);
		//TestAttractiveParticlesInteractingWithESandLJ(envmode);
		//TestIntegration(envmode);

		//TestUtils::TestIsDeterministic([]() {return loadAndEMAndRunBasicSimulation("T4Lysozyme", Headless, 2.8e-2, 5e-4); }, 2, envmode);
		//loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 2.8e-2, 5e-4);
		//loadAndRunBasicSimulation("T4Lysozyme", envmode, 1.15e-4, 2.e-6);


		//const fs::path workDir = simulations_dir / "test";
		//Lipids::Selection lipids;
		//lipids.emplace_back(Lipids::Select{ "DPPE", workDir, 30.5 });
		//lipids.emplace_back(Lipids::Select{ "DMPG", workDir, 39.5 });
		//lipids.emplace_back(Lipids::Select{ "Cholesterol", workDir, 10 });
		//lipids.emplace_back(Lipids::Select{ "SM18", workDir, 20 });
		//auto [grofile, topfile] = SimulationBuilder::CreateMembrane(lipids, Float3{ 20.f }, 5.f);
		//SimulationBuilder::CreateMembrane(*grofile, *topfile, lipids, 15.f);'
		//Programs::EnergyMinimize(*grofile, *topfile, true, workDir, envmode, true);
		//grofile->printToFile(workDir / "membrane.gro");
		//topfile->printToFile(workDir / "membrane.top");

		//TestBuildmembraneWithCustomlipidAndCustomForcefield(envmode);
		//TestBuildmembraneSmall(envmode, false);
		//TestAllStockholmlipids(envmode);

		//Lipids::_MakeLipid("cholesterol");

		//TestLimaChosesSameBondparametersAsGromacs(envmode);


		//TestMinorPrograms::InsertMoleculesAndDoStaticbodyEM(envmode);

		//TestForces1To1(envmode);

		//ForceComparisons::T4RmsdAndRmsf();
		//ForceComparisons::DoAllForceComparisons(envmode);

		//Benchmarks::Benchmark({ "t4", "membrane20", "manyt4" });		
		//Benchmarks::Benchmark({ "t4", "manyt4" });
		//Benchmarks::Benchmark("membrane20", "membranesolvated_em");
		//Benchmarks::Benchmark("manyt4", "manyt4sol");
		//Benchmarks::Benchmark("stmv", std::nullopt, 1000);
		/*for (int i = 0; i < 10; i++)
			Benchmarks::PrepareSimulation_stmv(envmode);*/
		//Benchmarks::STMV(500);
		//Benchmarks::Psome(envmode);
		

		//Benchmarks::PrepareSimulation_stmv(envmode);
		//Benchmarks::Psome(envmode);
		//
		// TopologyFile topfile1{ R"(C:\Users\Daniel\git_repo\LIMA_data\Solvents\molecule\topol.top)" };


		//{
		//	GroFile grofile{ R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\manyt4\manyt4.gro)" };
		//	TopologyFile topfile{ R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\manyt4\manyt4.top)" };
		//	SimulationBuilder::SolvateGrofile(grofile, topfile);

		//	grofile.printToFile(fs::path{ R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\manyt4\manyt4_solvated.gro)" });
		//	topfile.printToFile(fs::path{ R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\manyt4\manyt4_solvated.top)" });
		//}

		//{
		//	GroFile grofile{ R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\membrane20\membranesolvated.gro)" };
		//	TopologyFile topfile{ R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\membrane20\membranesolvated.top)" };
		//	Programs::EnergyMinimize(grofile, topfile, true, fs::current_path(), Full, false, 800.f);

		//	grofile.printToFile("membranesolvated_em.gro");
		//	topfile.printToFile("membranesolvated_em.top");
		//}
		/*GroFile grofile{ R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\stmv\em.gro)" };
		TopologyFile topfile{ R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\stmv\topol.top)" };
		Environment env(R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\stmv)", EnvMode::Full);
		env.CreateSimulation(grofile, topfile, SimParams{});
		env.run();*/
		//Programs::EnergyMinimize(grofile, topfile, true, fs::current_path(), Full, false, 800.f);
		
		//ForceComparisons::DoAllForceComparisons(envmode);

		//KernelAlgorithms::WarpSort64_Unittest(envmode);
		//Benchmarks::Psome(envmode);
//Benchmarks::ManyT4(envmode);
//Benchmarks::PrepareSimulation_stmv(envmode);
		//RunAllUnitTests();

	}
	catch (std::runtime_error ex) {
		std::cerr << "\nCaught runtime_error: " << ex.what() << std::endl;
	}
	catch (const std::exception& ex) {
		std::cerr << "\nCaught exception: " << ex.what() << std::endl;
	}
	catch (...) {
		std::cerr << "\nCaught unnamed exception";
	}

	return 0;
}


#define ADD_TEST(description, execution_function) \
    testman.addTest(std::make_unique<LimaUnittest>(LimaUnittest{ description, [](){ return execution_function;} }))

// Runs all unit tests with the fastest/crucial ones first
void RunAllUnitTests() {
	LimaUnittestManager testman;
	constexpr auto envmode = EnvMode::Headless;

	if (!ALL_PHYSICS_ENABLED) {
		TestUtils::setConsoleTextColorRed();
		std::cout << "WARNING: Not all physics modules are enabled, expect tests to fail!" << std::endl;
		TestUtils::setConsoleTextColorDefault();
	}



	// Isolated forces sanity checks
	ADD_TEST("SinglebondForceAndPotentialSanityCheck", SinglebondForceAndPotentialSanityCheck(envmode));
	ADD_TEST("SinglebondOscillationTest", SinglebondOscillationTest(envmode));
	ADD_TEST("UreyBradleyForceAndPotentialSanityCheck", UreyBradleyForceAndPotentialSanityCheck(envmode));
	ADD_TEST("PairbondForceAndPotentialSanityCheck", PairbondForceAndPotentialSanityCheck(envmode));
	ADD_TEST("TestIntegration", TestIntegration(envmode));

	// Stability tests
	ADD_TEST("doPoolBenchmark", doPoolBenchmark(envmode));
	ADD_TEST("doPoolCompSolBenchmark", doPoolCompSolBenchmark(envmode));
	ADD_TEST("doSinglebondBenchmark", doSinglebondBenchmark(envmode));
	ADD_TEST("doAnglebondBenchmark", doAnglebondBenchmark(envmode));
	ADD_TEST("doDihedralbondBenchmark", doDihedralbondBenchmark(envmode));
	ADD_TEST("doImproperDihedralBenchmark", doImproperDihedralBenchmark(envmode));

	// Smaller compound tests
	ADD_TEST("doMethionineBenchmark", loadAndRunBasicSimulation("Met", envmode, 1.199e-3, 2e-6));
	ADD_TEST("doEightResiduesNoSolvent", doEightResiduesNoSolvent(envmode));

	// Larger tests
	ADD_TEST("SolventBenchmark", loadAndRunBasicSimulation("Solvents", envmode, 3.22e-7, 1.1e-7));
	ADD_TEST("T4Lysozyme", loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 1.757e-3, 5e-4));
	ADD_TEST("Deterministic Simulations", TestUtils::TestIsDeterministic([]() {return loadAndEMAndRunBasicSimulation("T4Lysozyme", Headless, 2.8e-2, 5e-4); }, 2, envmode));


	// Electrostatics
	ADD_TEST("CoulombForceSanityCheck", CoulombForceSanityCheck(envmode));
	ADD_TEST("TestLongrangeEsNoLJTwoParticles", TestLongrangeEsNoLJTwoParticles(envmode));
	ADD_TEST("TestLongrangeEsNoLJManyParticles", TestLongrangeEsNoLJManyParticles(envmode));
	//ADD_TEST("TestElectrostaticsManyParticles", TestElectrostaticsManyParticles(envmode));
	ADD_TEST("TestChargedParticlesVelocityInUniformElectricField", TestChargedParticlesVelocityInUniformElectricField(envmode));

	// Test Forcefield and compoundbuilder
	ADD_TEST("TestLimaChosesSameBondparametersAsGromacs", TestLimaChosesSameBondparametersAsGromacs(envmode));

	// Test Setup
	ADD_TEST("TestBoxIsSavedCorrectlyBetweenSimulations", TestBoxIsSavedCorrectlyBetweenSimulations(envmode));

	// Programs test
	ADD_TEST("BuildSmallMembrane", TestBuildmembraneSmall(envmode, false));
	ADD_TEST("BuildSphericalMembrane", TestSphericalMembraneBuilder(envmode));
	ADD_TEST("TestBuildmembraneWithCustomlipidAndCustomForcefield", TestBuildmembraneWithCustomlipidAndCustomForcefield(envmode));
	ADD_TEST("TestAllStockholmlipids", TestAllStockholmlipids(envmode));

	// Gromacs correctness
	ADD_TEST("ForceComparisons", ForceComparisons::DoAllForceComparisons(envmode));

	//ADD_TEST("InsertMoleculesAndDoStaticbodyEM", TestMinorPrograms::InsertMoleculesAndDoStaticbodyEM(envmode));

	//ADD_TEST("ReorderMoleculeParticles", testReorderMoleculeParticles(envmode));
	//ADD_TEST("TestFilesAreCachedAsBinaries", FileTests::TestFilesAreCachedAsBinaries(envmode)); too slow to run...

	// Performance test
	//ADD_TEST(testman, "Benchmark Psome", Benchmarks::Psome(envmode));

	// Meta tests
	//doPool50x(EnvMode::Headless);


	//ADD_TEST("TestBuildmembranesInterface", UserinterfaceTests::TestBuildmembranesInterface(envmode));

	// Total test status will print as testman is destructed
}
