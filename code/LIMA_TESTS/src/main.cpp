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
#include "BatchingTests.h"
#include "EngineBatchTests.h"
#include "DisplayTests.h"


using namespace TestUtils;
using namespace ForceCorrectness;
using namespace TestMDStability;
using namespace StressTesting;
using namespace TestMembraneBuilder;
using namespace TestMinorPrograms;
using namespace ElectrostaticsTests;
using namespace VerletintegrationTesting;
using namespace BatchingTests;

void RunAllUnitTests();

void TestDisplayT4() {
	const fs::path workDir = AutomatedTestsDir() / "T4Lysozyme";
	Environment& environment = Environment::Get();
	SimulationJob job;
	job.workDir = workDir;
	job.run = false;
	auto result = environment.Submit(std::move(job)).Get();
	Display display{};

	auto& box = result.simulation->box;

	std::vector<Float3> positions;

	for (auto pc : box->persistentClusters) {
		for (auto pqd : pc.pqd)
			if (pqd.Valid())
				positions.push_back(pqd.position);
	}


	display.Submit(0, std::make_unique<Rendering::AtomRenderTask>(
		box->persistentClusters, box->persistentClustersMetadata, box->boxparams, SimStatus{}, box->backboneChains), true);
	display.Submit(0, std::make_unique<Rendering::SimulationTaskUpdate>(positions.data(), nullptr, SimStatus{}), true);
}

void LiveEditTest() {
	Environment& env = Environment::Get();
	auto [grofile, topfile, simparams] = env.BeginLiveEdit(
		R"(C:\Users\Daniel\git_repo\LIMA_data\LiveEditTest)", EnvMode::Full, Float3(25.f));

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
	

	env.QueueLiveEditCommand(LiveEdit::InsertMolecule{ "t4/conf.gro", "t4/topol_T4.itp", Float3{8, 10, 10 }});
	env.QueueLiveEditCommand(LiveEdit::InsertMolecule{ "t4/conf.gro", "t4/topol_T4.itp", Float3{16, 10, 10 }});
	env.QueueLiveEditCommand(LiveEdit::SelectAtomsBasedOnQualifier{ LiveEdit::SelectAtomsBasedOnQualifier::Qualifier::All });
	env.QueueLiveEditCommand(LiveEdit::ElasticPosition{ true, false, false });

	env.LiveEdit(grofile, topfile);
}

void BuildCellTest() {
	Environment& env = Environment::Get();
	auto [grofile, topfile, simparams] = env.BeginLiveEdit(
		R"(C:\Users\Daniel\git_repo\LIMA_data\LiveEditTest)", EnvMode::Full, Float3(30, 30, 40));

	//// TODO: Being able to set this is like super dangerous, and the ff is then not parsed... Always need a file reset..
	topfile.forcefieldInclude = TopologyFile::ForcefieldInclude("combined/forcefield.itp");
	topfile.printToFile();
	topfile = TopologyFile{ topfile.path };


	//std::vector<std::tuple<std::string, double>> lipids = { {"DPPE", 100.}};
	std::vector<std::tuple<std::string, double>> lipids = { {"DPPE", 30.5}, {"DMPG", 39.5}, {"cholesterol", 10}, {"SM18", 20} };
	//env.liveEditCommandsQueue.push_back(LiveEdit::BuildMembrane{ lipids, MembraneGeometry::Plane{ 5.f } });
	env.QueueLiveEditCommand(LiveEdit::BuildMembrane{ lipids, MembraneGeometry::Ellipsoid{Float3{15,15,20  }, Float3{12, 12, 18 }} });

	env.LiveEdit(grofile, topfile);
}

void TestDisplayT4Batch() {
	Environment& environment = Environment::Get();
	std::array<SimulationHandle, 4> handles;
	for (auto& handle : handles) {
		auto job = BatchingTests::MakeT4Job(EnvMode::Full, false);
		job.mode = EnvMode::Full;
		job.preprocess = [](GroFile&, TopologyFile&, SimParams& params) {
			params.n_steps = 90000;
			params.data_logging_interval = 200;
		};
		handle = environment.Submit(std::move(job));
	}

	// Submit the entire batch before waiting, otherwise each Get serializes the jobs.
	for (auto& handle : handles) {
		auto result = handle.Get();
		if (!result.simulation || result.execution.batchSize != handles.size())
			throw std::runtime_error("Display T4 simulations did not execute as one batch");
	}
}


// Demonstrates LIMA's parallel simulation workflow.
// Builds two membrane compositions using three independent seeds each, then energy-minimizes all six systems.
// Each system is subsequently simulated at 300 K and 340 K, yielding 12 production simulations for comparing membrane stability across composition and temperature.
void ShowcaseMultisim() {
	const fs::path workDir = TestUtils::HeavyTestsDir() / "etc";
	const fs::path outputDir = workDir / "showcase_multisim";
	constexpr std::array temperatures{ 300.f, 340.f };
	constexpr std::array seeds{ 101, 202, 303 };	
	std::vector<Lipids::Selection> lipidSelections{
		{Lipids::Select{ "DPPC", workDir, 70. },	Lipids::Select{ "DOPC", workDir, 30. }},
		{Lipids::Select{ "DPPC", workDir, 40. },	Lipids::Select{ "DOPC", workDir, 60. }}
	};
	

	struct PreparedMembrane {
		std::string name;
		int seed = 0;
		fs::path initialCoordinates;
		fs::path minimizedCoordinates;
		fs::path topology;
	};
	std::vector<std::shared_ptr<PreparedMembrane>> membranes;
	std::vector<SimulationHandle> minimizations;
	Environment& environment = Environment::Get();

	for (const auto& lipidSelection : lipidSelections) {
		for (int i = 0; i < seeds.size(); i++) {
			auto seed = seeds[i];
			auto membrane = std::make_shared<PreparedMembrane>();
			membrane->name = Lipids::NameSelection(lipidSelection) + std::format("_seed{}", seed);
			membrane->seed = seed;
			const fs::path membraneDir = outputDir / membrane->name;
			membrane->initialCoordinates = membraneDir / "initial.gro";
			membrane->minimizedCoordinates = membraneDir / "minimized.gro";
			membrane->topology = membraneDir / "topol.top";
			membranes.push_back(membrane);

			SimulationJob job;
			job.workDir = workDir;
			job.grofile.emplace();
			job.grofile->box_size = Float3{ 12.f };
			job.topfile.emplace();
			job.simParams = SimParams::BasicEMSimParams(800.f);
			job.simParams->n_steps = 5000;
			job.mode = Full;
			job.preprocess = [workDir, lipidSelection, seed, membrane](
				GroFile& gro, TopologyFile& top, SimParams&) {
				SimulationBuilder::CreateMembrane(gro, top, lipidSelection, 4.f, seed);
				gro.printToFile(membrane->initialCoordinates);
				top.printToFile(membrane->topology);
			};
			job.configureSimulation = [seed](Simulation& simulation) {
				simulation.name += std::format(" seed {}", seed);
			};
			minimizations.push_back(environment.Submit(std::move(job)));
		}
	}

	for (size_t i = 0; i < minimizations.size(); ++i) {
		auto result = minimizations[i].Get();
		if (!result.simulation)
			throw std::runtime_error("Multisim membrane minimization did not return a simulation");
		GroFile gro{ membranes[i]->initialCoordinates };
		result.WriteCoordinatesTo(gro);
		gro.printToFile(membranes[i]->minimizedCoordinates);
	}

	std::vector<SimulationHandle> productions;
	std::vector<std::pair<size_t, float>> productionMetadata;
	for (size_t membraneId = 0; membraneId < membranes.size(); ++membraneId) {
		for (const float temperature : temperatures) {
			SimParams params;
			params.n_steps = 5000;
			params.apply_thermostat = true;
			params.ref_t = temperature;
			params.save_energy = true;
			params.data_logging_interval = 100;
			const auto& membrane = *membranes[membraneId];
			SimulationJob job;
			job.workDir = workDir;
			job.groPath = membrane.minimizedCoordinates;
			job.topPath = membrane.topology;
			job.simParams = params;
			job.mode = Full;
			job.configureSimulation = [seed = membrane.seed, temperature](Simulation& simulation) {
				simulation.name += std::format(" seed:{} temperature:{} K", seed, static_cast<int>(temperature));
			};
			productions.push_back(environment.Submit(std::move(job)));
			productionMetadata.emplace_back(membraneId, temperature);
		}
	}

	for (size_t i = 0; i < productions.size(); ++i) {
		auto result = productions[i].Get();
		if (!result.simulation || result.simulation->simParams.ref_t != productionMetadata[i].second)
			throw std::runtime_error("Multisim production simulation used the wrong temperature");
		const auto& [membraneId, temperature] = productionMetadata[i];
		const auto& membrane = *membranes[membraneId];
		const fs::path runDir = outputDir / membrane.name / std::format("{}K", static_cast<int>(temperature));
		GroFile finalCoordinates{ membrane.minimizedCoordinates };
		result.WriteCoordinatesTo(finalCoordinates);
		finalCoordinates.printToFile(runDir / "out.gro");
	}

	std::cout << "Multisim showcase completed: 6 minimized membranes and 12 production simulations\n";
}



int main(int argc, char** argv) {
	if (argc == 2 && std::string_view(argv[1]) == "--display-tiles-tests") {
		try { DisplayTests::Run(); return 0; }
		catch (const std::exception& ex) { std::cerr << ex.what() << "\n"; return 1; }
	}
	if (argc == 2 && (std::string_view(argv[1]) == "--display-preview" || std::string_view(argv[1]) == "--display-tiles-preview")) {
		try {
			// Static UI fixture: no simulation runs or output files are written.
			const GroFile molecule{ AutomatedTestsDir() / "T4Lysozyme" / "molecule" / "conf.gro" };
			Display display;
			const bool tilePreview = std::string_view(argv[1]) == "--display-tiles-preview";
			display.allowUserInputs = !tilePreview;
			for (int i = 0; i < (tilePreview ? 9 : 4); ++i) {
				auto task = std::make_unique<Rendering::AtomRenderTask>(molecule, false);
				task->simStatus.step = 24000 + i * 1000;
				task->simStatus.temperature = 300.12f + i;
				task->simStatus.maxForce = 1.23e3f;
				task->simStatus.expectedTimeToFinish = std::chrono::duration<double>{ 154. };
				task->simStatus.avgStepTime = .842f;
				task->simStatus.simulationPerformance = 205.23f;
				display.Submit(i, std::move(task), false, nullptr, "Preview " + std::to_string(i + 1));
			}
			while (!display.DisplaySelfTerminated())
				std::this_thread::sleep_for(std::chrono::milliseconds(50));
			if (display.displayThreadException)
				std::rethrow_exception(display.displayThreadException);
			return 0;
		}
		catch (const std::exception& ex) { std::cerr << ex.what() << "\n"; return 1; }
	}
	if (argc == 2 && std::string_view(argv[1]) == "--environment-batch-tests") {
		try { BatchingTests::RunSchedulerTests(Environment::Get()); return 0; }
		catch (const std::exception& ex) { std::cerr << ex.what() << "\n"; return 1; }
	}
	if (argc == 2 && std::string_view(argv[1]) == "--engine-batch-tests") {
		try { EngineBatchTests::RunAll(); return 0; }
		catch (const std::exception& ex) { std::cerr << ex.what() << "\n"; return 1; }
	}
	if (argc == 2 && std::string_view(argv[1]) == "--showcase-multisim") {
		try { ShowcaseMultisim(); return 0; }
		catch (const std::exception& ex) { std::cerr << ex.what() << "\n"; return 1; }
	}
	try {
		constexpr auto envmode = EnvMode::Full;
		Environment& env = Environment::Get();
		//TestDisplayT4();
		//ProgramsTests::TestToGmx_ciffile(envmode);
		//LiveEditTest();
		//BuildCellTest();
		//Benchmarks::ToGmxLargeCif(envmode);

		ShowcaseMultisim();
		//TestDisplayT4Batch();
		//Display::RenderGrofile(TestUtils::AutomatedTestsDir() / "BuildMembraneSmall" / "molecule" / "membrane.gro");
		return 0;
		//TestFourT4BatchMatchReference(env, envmode).RunToCompletion();
		//return 0;
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
		//LoadEnergyMinAndRunBasicSimulation(env, envmode, "T4Lysozyme", "T4Lysozyme");
		//loadAndRunBasicSimulation("T4Lysozyme", envmode,"T4Lysozyme");
		//LoadAndRunBasicSimulation(env, envmode, "Solvents", "SolventBenchmark");

		//const fs::path workDir = simulations_dir / "test";
		//Lipids::Selection lipids;
		//lipids.emplace_back(Lipids::Select{ "DPPE", workDir, 30.5 });
		//lipids.emplace_back(Lipids::Select{ "DMPG", workDir, 39.5 });
		//lipids.emplace_back(Lipids::Select{ "Cholesterol", workDir, 10 });
		//lipids.emplace_back(Lipids::Select{ "SM18", workDir, 20 });
		//GroFile grofile;
		//grofile.box_size = Float3{ 20.f };
		//TopologyFile topfile;
		//topfile.SetSystem("Membrane");
		//SimulationBuilder::CreateMembrane(grofile, topfile, lipids, 5.f);
		//SimulationBuilder::CreateMembrane(grofile, topfile, lipids, 15.f);
		//grofile->printToFile(workDir / "membrane.gro");
		//topfile->printToFile(workDir / "membrane.top");

		//TestBuildmembraneWithCustomlipidAndCustomForcefield(envmode);
		//TestBuildmembraneSmall(envmode, false);
		//TestAllStockholmlipids(env, envmode);

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

		//	grofile.printToFile("membranesolvated_em.gro");
		//	topfile.printToFile("membranesolvated_em.top");
		//}

		/*GroFile grofile{ R"(C:\Users\Daniel\git_repo\LIMA_data\T4Lysozyme\molecule\out.gro)" };
		Display::RenderGrofile(grofile, true);*/
		

		/*GroFile grofile{ R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\stmv\em.gro)" };
		TopologyFile topfile{ R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\stmv\topol.top)" };
		Environment& env = Environment::Get();
		env.CreateSimulation(grofile, topfile, SimParams{});
		env.run();*/
		
		//ForceComparisons::DoAllForceComparisons(envmode);

		//KernelAlgorithms::WarpSort64_Unittest(envmode);
		//Benchmarks::Psome(envmode);
//Benchmarks::ManyT4(envmode);
//Benchmarks::PrepareSimulation_stmv(envmode);
		//TestBuildmembraneSmall(envmode, false);
		// 
		//Benchmarks::STMV(env, envmode, 200, 3).RunToCompletion();
		// 
		RunAllUnitTests();
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


// Every test receives the shared Environment and mode. Extra arguments are only
// written for tests that genuinely need them.
#define ADD_TEST(description, test_function, ...) \
	testman.AddTest(description, [&] { \
		return test_function(environment, envmode __VA_OPT__(,) __VA_ARGS__); \
	})

#define ADD_SERIAL_TEST(description, test_function, ...) \
	ADD_TEST(description, test_function __VA_OPT__(,) __VA_ARGS__)

// Runs all unit tests with the fastest/crucial ones first
void RunAllUnitTests() {
	TimeIt timer("RunAllUnitTests", true);
	Environment& environment = Environment::Get();
	LimaUnittestManager testman;
	constexpr auto envmode = EnvMode::Headless;

	if (!ALL_PHYSICS_ENABLED) {
		TestUtils::setConsoleTextColorRed();
		std::cout << "WARNING: Not all physics modules are enabled, expect tests to fail!" << std::endl;
		TestUtils::setConsoleTextColorDefault();
	}
#ifdef _DEBUG
	TestUtils::setConsoleTextColorRed();
	std::cout << "WARNING: Running tests in debug mode may result in incorrect VC results due to missing floating point math optimizations" << std::endl;
	TestUtils::setConsoleTextColorDefault();
#endif



	// Isolated forces sanity checks
	ADD_TEST("SinglebondForceAndPotentialSanityCheck", SinglebondForceAndPotentialSanityCheck);
	ADD_TEST("SinglebondOscillationTest", SinglebondOscillationTest);
	ADD_TEST("UreyBradleyForceAndPotentialSanityCheck", UreyBradleyForceAndPotentialSanityCheck);
	ADD_TEST("PairbondForceAndPotentialSanityCheck", PairbondForceAndPotentialSanityCheck);
	ADD_TEST("TestIntegration", TestIntegration);

	// Stability tests
	ADD_TEST("doPoolBenchmark", doPoolBenchmark);
	ADD_TEST("doPoolCompSolBenchmark", doPoolCompSolBenchmark);
	ADD_TEST("doSinglebondBenchmark", doSinglebondBenchmark);
	ADD_TEST("doAnglebondBenchmark", doAnglebondBenchmark);
	ADD_TEST("doDihedralbondBenchmark", doDihedralbondBenchmark);
	ADD_TEST("doImproperDihedralBenchmark", doImproperDihedralBenchmark);

	// Smaller compound tests
	ADD_TEST("doMethionineBenchmark", LoadAndRunBasicSimulation, "Met", "doMethionineBenchmark");
	ADD_TEST("doEightResiduesNoSolvent", LoadAndRunBasicSimulation, "8ResNoSol", "doEightResiduesNoSolvent");

	// Larger tests
	ADD_TEST("SolventBenchmark", LoadAndRunBasicSimulation, "Solvents", "SolventBenchmark");
	ADD_TEST("T4Lysozyme", LoadEnergyMinAndRunBasicSimulation, "T4Lysozyme", "T4Lysozyme");
	ADD_TEST("Deterministic Simulations", TestDeterministic);
	ADD_TEST("Four batched T4 simulations match reference", TestFourT4BatchMatchReference);


	// Electrostatics
	ADD_TEST("CoulombForceSanityCheck", CoulombForceSanityCheck);
	ADD_TEST("TestLongrangeEsNoLJTwoParticles", TestLongrangeEsNoLJTwoParticles);
	ADD_TEST("TestLongrangeEsNoLJManyParticles", TestLongrangeEsNoLJManyParticles);
	//ADD_TEST("TestElectrostaticsManyParticles", TestElectrostaticsManyParticles(envmode));
	ADD_TEST("TestChargedParticlesVelocityInUniformElectricField", TestChargedParticlesVelocityInUniformElectricField);

	// Test Forcefield and compoundbuilder
	ADD_TEST("TestLimaChosesSameBondparametersAsGromacs", TestLimaChosesSameBondparametersAsGromacs);

	// Test Setup
	ADD_TEST("TestBoxIsSavedCorrectlyBetweenSimulations", TestBoxIsSavedCorrectlyBetweenSimulations);

	// Programs test
	ADD_TEST("ToGmx PDB matches GROMACS", ProgramsTests::TestToGmx_pdbfile);
	ADD_TEST("ToGmx handles multiple chains", ProgramsTests::TestToGmx_multichain);
	ADD_TEST("ToGmx CIF matches GROMACS", ProgramsTests::TestToGmx_ciffile);
	ADD_TEST("BuildSmallMembrane", TestBuildmembraneSmall, false);
	ADD_TEST("BuildSphericalMembrane", TestSphericalMembraneBuilder);
	ADD_TEST("TestBuildmembraneWithCustomlipidAndCustomForcefield", TestBuildmembraneWithCustomlipidAndCustomForcefield);
	ADD_TEST("TestAllStockholmlipids", TestAllStockholmlipids);

	// Gromacs correctness
	ADD_TEST("ForceComparisons", ForceComparisons::DoAllForceComparisons);

	//ADD_TEST("InsertMoleculesAndDoStaticbodyEM", TestMinorPrograms::InsertMoleculesAndDoStaticbodyEM(envmode));

	//ADD_TEST("ReorderMoleculeParticles", testReorderMoleculeParticles(envmode));
	//ADD_SERIAL_TEST("TestFilesAreCachedAsBinaries", FileTests::TestFilesAreCachedAsBinaries(envmode)); too slow to run...

	// Performance test
	ADD_TEST("ToGmx large CIF benchmark", Benchmarks::ToGmxLargeCif);
	ADD_TEST("T4", Benchmarks::T4, 200, Benchmarks::automatedTestRuns);
	ADD_TEST("stmv sim performance", Benchmarks::STMV, 200, Benchmarks::automatedTestRuns);

	// Meta tests
	//doPool50x(EnvMode::Headless);


	//ADD_SERIAL_TEST("TestBuildmembranesInterface", UserinterfaceTests::TestBuildmembranesInterface(envmode));

	// Total test status will print as testman is destructed
}
