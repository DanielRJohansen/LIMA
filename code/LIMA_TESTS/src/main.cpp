#include "ForceCorrectness.h"
#include "MDStability.cuh"
#include "MembraneBuilder.h"
#include "MinorPrograms.h"
#include "ElectrostaticsTests.h"
#include "Benchmarks.h"
#include "FileTests.h"
#include "ForcefieldTests.h"
#include "SetupTests.h"
#include "Userinterface.h"
#include "Display.h"

using namespace TestUtils;
using namespace ForceCorrectness;
using namespace TestMDStability;
using namespace StressTesting;
using namespace TestMembraneBuilder;
using namespace TestMinorPrograms;
using namespace ElectrostaticsTests;
using namespace VerletintegrationTesting;

void RunAllUnitTests();

uint8_t Dir(int x, int y, int z)  {
	return static_cast<uint8_t>(((x + 1) << 4) | ((y + 1) << 2) | (z + 1));
}


struct Direction3 {
	uint8_t data; // store 6 bits: top 2 for x+1, next 2 for y+1, bottom 2 for z+1

	 constexpr Direction3() : data(0) {}

	 constexpr Direction3(int x, int y, int z)
		: data(static_cast<uint8_t>(((x + 1) << 4) | ((y + 1) << 2) | (z + 1)))
	{}

	 inline int x() const { return ((data >> 4) & 0x3) - 1; }
	 inline int y() const { return ((data >> 2) & 0x3) - 1; }
	 inline int z() const { return (data & 0x3) - 1; }

	 bool operator==(const Direction3& o) const { return data == o.data; }
	 bool operator!=(const Direction3& o) const { return data != o.data; }
};

int main() {
	try {
		constexpr auto envmode = EnvMode::Full;

		//PlotPmePotAsFactorOfDistance(envmode);
		//TestConsistentEnergyWhenGoingFromLresToSres(envmode);
		//TestLongrangeEsNoLJTwoParticles(envmode);
		//TestLongrangeEsNoLJManyParticles(envmode);
		//Lipids::_MakeLipids(true, false);
		//PairbondForceAndPotentialSanityCheck(envmode);
		//loadAndRunBasicSimulation("DisplayTest", envmode);
		//Display::TestDisplay();
		//doPoolBenchmark(envmode);			// Two 1-particle molecules colliding
		//loadAndRunBasicSimulation("PoolElectrostatic", envmode);
		//doPoolCompSolBenchmark(envmode);	// One 1-particle molecule colliding with 1 solvent
		//SinglebondForceAndPotentialSanityCheck(envmode);
		//UreyBradleyForceAndPotentialSanityCheck(envmode);
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

		
		//loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 1.8e-3, 2e-5);
		//loadAndRunBasicSimulation("T4Lysozyme", envmode, 1.15e-4, 2.e-6);
		//loadAndRunBasicSimulation("T4Lysozyme", envmode, 1.15e-4, 2.e-6);

		//const fs::path work_dir = simulations_dir / "test";
		//Lipids::Selection lipids;
		//lipids.emplace_back(Lipids::Select{ "DPPE", work_dir, 30.5 });
		//lipids.emplace_back(Lipids::Select{ "DMPG", work_dir, 39.5 });
		//lipids.emplace_back(Lipids::Select{ "Cholesterol", work_dir, 10 });
		//lipids.emplace_back(Lipids::Select{ "SM18", work_dir, 20 });
		//auto [grofile, topfile] = SimulationBuilder::CreateMembrane(lipids, Float3{ 20.f }, 5.f);
		//SimulationBuilder::CreateMembrane(*grofile, *topfile, lipids, 15.f);
		//Programs::EnergyMinimize(*grofile, *topfile, true, work_dir, envmode, true);
		//grofile->printToFile(work_dir / "membrane.gro");
		//topfile->printToFile(work_dir / "membrane.top");

		//TestBuildmembraneWithCustomlipidAndCustomForcefield(envmode);
		//TestBuildmembraneSmall(envmode, false);
		//TestAllStockholmlipids(envmode);




		//TestMinorPrograms::InsertMoleculesAndDoStaticbodyEM(envmode);
		


		//fs::path dir = R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\EAG1-channel_strong-scaling\inputs\)";
		//fs::path dir = R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\stmv\)";
		//GroFile grofile(dir/"conf.gro");
		//grofile.box_size = Float3(std::ceil(std::max(std::max(grofile.box_size.x, grofile.box_size.y), grofile.box_size.z)));
		////Display::RenderGrofile(grofile, false);
		//TopologyFile topfile(dir / "topol.top");
		//Programs::EnergyMinimize(grofile, topfile, true, dir, envmode, false);
		//grofile.printToFile(std::string{ "em.gro" });

		//SimParams simparams(R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\sim_params.txt)");
		//simparams.dt = 100.f;
		//////auto sim = Programs::EnergyMinimize(grofile, topfile, true, R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\EAG1-channel_strong-scaling)", Full, false);
		//Environment env(R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\EAG1-channel_strong-scaling\inputs)", Full);
		//env.CreateSimulation(grofile, topfile, simparams);
		//env.run(false);

		/*GroFile grofile(R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\membrane20\membrane.gro)");
		TopologyFile topfile(R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\membrane20\membrane.top)");
		Programs::EnergyMinimize(grofile, topfile, true, R"(C:\Users\Daniel\git_repo\LIMA_data\benchmarking\membrane20)", envmode, false);
		grofile.printToFile(std::string("membrane_em.gro"));*/

		//TestForces1To1(envmode);

		//Benchmarks::Benchmark({ "t4", "membrane20", "manyt4" });		
		//Benchmarks::Benchmark({ "t4", "manyt4" });
		//Benchmarks::Benchmark("membrane20"); 
		//Benchmarks::Benchmark("manyt4"); 
		RunAllUnitTests();
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
void RunAllUnitTests() {
	LimaUnittestManager testman;
	constexpr auto envmode = EnvMode::Headless;

	
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
	ADD_TEST("doMethionineBenchmark", TestUtils::loadAndRunBasicSimulation("Met", envmode, 6.3e-4, 2e-6));
	ADD_TEST("doEightResiduesNoSolvent", doEightResiduesNoSolvent(envmode));

	// Larger tests
	ADD_TEST("SolventBenchmark", loadAndRunBasicSimulation("Solvents", envmode, 5.85e-6f, 1.1e-7));
	ADD_TEST("T4Lysozyme", loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 1.224e-3, 2e-5));

	// Electrostatics
	ADD_TEST("CoulombForceSanityCheck", CoulombForceSanityCheck(envmode));
	ADD_TEST("TestLongrangeEsNoLJTwoParticles", TestLongrangeEsNoLJTwoParticles(envmode));
	ADD_TEST("TestLongrangeEsNoLJManyParticles", TestLongrangeEsNoLJManyParticles(envmode));
	ADD_TEST("TestElectrostaticsManyParticles", TestElectrostaticsManyParticles(envmode));
	ADD_TEST("TestChargedParticlesVelocityInUniformElectricField", TestChargedParticlesVelocityInUniformElectricField(envmode));
	
	// Test Forcefield and compoundbuilder
	ADD_TEST("TestLimaChosesSameBondparametersAsGromacs", TestLimaChosesSameBondparametersAsGromacs(envmode));

	// Test Setup
	ADD_TEST("TestBoxIsSavedCorrectlyBetweenSimulations", TestBoxIsSavedCorrectlyBetweenSimulations(envmode));

	// Programs test
	ADD_TEST("BuildSmallMembrane", TestBuildmembraneSmall(envmode, false));
	ADD_TEST("TestBuildmembraneWithCustomlipidAndCustomForcefield", TestBuildmembraneWithCustomlipidAndCustomForcefield(envmode));
	ADD_TEST("TestAllStockholmlipids", TestAllStockholmlipids(envmode));

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
