#include "Tests.h"
#include "Display.h"
#include "Workflow.h"
#include "LimitTesting.h"


using namespace TestUtils;
using namespace ForceCorrectness;
using namespace TestMDStability;
using namespace TestMembraneBuilder;
using namespace TestMinorPrograms;
using namespace ElectrostaticsTests;
using namespace VerletintegrationTesting;
using namespace BatchingTests;

int RunAllUnitTests();

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


	env.QueueLiveEditCommand(LiveEdit::InsertMolecule{ "t4/conf.gro", "t4/topol_T4.itp", Float3{8, 10, 10 } });
	env.QueueLiveEditCommand(LiveEdit::InsertMolecule{ "t4/conf.gro", "t4/topol_T4.itp", Float3{16, 10, 10 } });
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


// Runs the STMV benchmark system with the display, to visually inspect the simulation
void ShowcaseSTMV(int nSteps = 50000) {
	const fs::path workDir = HeavyTestsDir() / "benchmarking" / "stmv";
	SimulationJob job;
	job.workDir = workDir;
	job.grofile.emplace(workDir / "conf.gro");
	job.topfile.emplace(workDir / "topol.top");
	job.simParams.emplace(workDir / "sim_params.txt");
	job.mode = EnvMode::Full;
	job.mustRunAlone = true;
	job.preprocess = [nSteps](GroFile&, TopologyFile&, SimParams& params) {
		params.data_logging_interval = 20;
		params.enable_electrostatics = true;
		params.n_steps = nSteps;
	};
	auto result = Environment::Get().Submit(std::move(job)).Get();
	if (!result.simulation || result.simulation->getStep() != nSteps)
		throw std::runtime_error("STMV showcase did not run fully");
	std::cout << "STMV showcase completed: " << nSteps << " steps\n";
}

// Demonstrates LIMA's parallel simulation workflow.
// Builds two membrane compositions using three independent seeds each, then energy-minimizes all six systems.
// Each system is subsequently simulated at 300 K and 340 K, yielding 12 production simulations for comparing membrane stability across composition and temperature.
void ShowcaseMultisim() {
	const fs::path workDir = TestUtils::HeavyTestsDir() / "etc" / "showcase_multisim";
	std::vector<Lipids::Selection> lipidSelections{
		{Lipids::Select{ "DPPC", workDir, 70. }, Lipids::Select{ "DOPC", workDir, 30. }},
		{Lipids::Select{ "DPPC", workDir, 40. }, Lipids::Select{ "DOPC", workDir, 60. }}
	};

	Programs::SimulationWorkflow workflow{ workDir, EnvMode::Full };
	workflow.AddInputs(Programs::MakeMembraneInputs(lipidSelections, { 101, 202, 303 },
		Float3{ 12.f }, MembraneGeometry::Plane{ 4.f }, true));
	SimParams minimization = SimParams::BasicEMSimParams(800.f);
	minimization.n_steps = 5000;
	workflow.AddStage({ "minimize", minimization, {},
		{ OutputSelect::InitialCoordinates, OutputSelect::FinalCoordinates, OutputSelect::Topology } });

	SimParams production;
	production.n_steps = 1000;
	production.apply_thermostat = true;
	production.save_energy = true;
	workflow.AddStage({ "production", production, {
		{ "300K", { { "temperature", "300" } }, [](SimParams& params) { params.ref_t = 300.f; } },
		{ "340K", { { "temperature", "340" } }, [](SimParams& params) { params.ref_t = 340.f; } }
	}, { OutputSelect::FinalCoordinates, OutputSelect::DensityProfile }, true });
	workflow.CompareDensityProfiles("composition", "temperature");
	workflow.Run();

	std::cout << "Multisim showcase completed: 6 minimized membranes and 12 production simulations\n";
}



int main(int argc, char** argv) {
	try {
		// Dispatch before creating Environment: fatal GPU probes belong only to the child process.
		if (argc == 3 && std::string_view(argv[1]) == "--limit-case") {
			const int result = LimitTesting::RunChild(argv[2]);
			// RunChild has destroyed every fixture and explicitly tested Engine cleanup. Skip global
			// CUDA/CRT teardown, which can hang after a deliberately poisoned device context.
			std::cout.flush();
			std::cerr.flush();
			std::_Exit(result);
		}
		if (argc >= 2 && std::string_view(argv[1]) == "--limit-tests") {
			if (argc > 3) throw std::invalid_argument("Usage: limatest --limit-tests [NAME-SUBSTRING]");
			LimaUnittestManager testman(false);
			LimitTesting::AddTests(testman, argc == 3 ? argv[2] : "");
			return testman.Finish() == 0 ? 0 : 1;
		}
		constexpr auto envmode = EnvMode::Full;
		Environment& env = Environment::Get();
		//TestDisplayT4();
		//ProgramsTests::TestToGmx_ciffile(envmode);
		//LiveEditTest();
		//BuildCellTest();
		//Benchmarks::ToGmxLargeCif(envmode);

		
		//ShowcaseSTMV();

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
		//Benchmarks::Benchmark(env, envmode, "stmv", std::nullopt, 1000);
		//Benchmarks::STMV(env, envmode, 1000);
		
		//ShowcaseMultisim();
		//Benchmarks::Load3J3Q(env, envmode);
		//ShowcaseSTMV();
		return RunAllUnitTests() == 0 ? 0 : 1;
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

	return 1;
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
int RunAllUnitTests() {
	TimeIt timer("RunAllUnitTests", true);
	Environment& environment = Environment::Get();
	LimaUnittestManager testman;
	constexpr auto envmode = EnvMode::Headless;

	// Run before enqueueing ordinary simulations so isolated GPU probes do not compete with them.
	//LimitTesting::AddTests(testman);

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
	ADD_TEST("TestTopologyPreprocessor", FileTests::TestTopologyPreprocessor);

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
	ADD_TEST("3j3q load benchmark", Benchmarks::Load3J3Q);
	ADD_TEST("LoadT4", Benchmarks::LoadT4);
	ADD_TEST("T4", Benchmarks::T4, 200, Benchmarks::automatedTestRuns);
	ADD_TEST("stmv sim performance", Benchmarks::STMV, 200, Benchmarks::automatedTestRuns);

	// Meta tests
	//doPool50x(EnvMode::Headless);



	return testman.Finish();
}
