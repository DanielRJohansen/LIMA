

#include "LIMA_TESTS/src/ForceCorrectness.h"
#include "LIMA_TESTS/src/MDStability.cuh"
#include "LIMA_TESTS/src/EnergyMinimization.cuh"




using namespace TestUtils;
using namespace ForceCorrectness;
using namespace TestMDStability;
using namespace StressTesting;

void runAllUnitTests();

int main() {
	constexpr auto envmode = EnvMode::Full;

	
	//doPoolBenchmark(envmode);			// Two 1-particle molecules colliding
	//doPoolCompSolBenchmark(envmode);	// One 1-particle molecule colliding with 1 solvent
	
	//doSinglebondBenchmark(envmode);
	//doAnglebondBenchmark(envmode);
	//doDihedralbondBenchmark(envmode);
	//TestUtils::loadAndRunBasicSimulation("torsion2", envmode, 0.0002);
	//doImproperDihedralBenchmark(envmode);
//	TestUtils::loadAndRunBasicSimulation("improper", envmode, 5e-5);

	//doMethionineBenchmark(envmode);
	//doPhenylalanineBenchmark(envmode);
	//TestUtils::loadAndRunBasicSimulation("TenSolvents", envmode, 0.0004, 1.2e-6);

	//doEightResiduesNoSolvent(envmode);
	//loadAndEMAndRunBasicSimulation("SolventBenchmark", envmode, 2e-6);

	loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 0.023, 2e-5);


	//doPool50x(EnvMode::Headless);

	//runAllUnitTests();


	return 0;
}

// Runs all unit tests with the fastest/crucial ones first
void runAllUnitTests() {
	LimaUnittestManager testman;
	constexpr auto envmode = EnvMode::Headless;


	// Singled out forces test
	testman.addTest(doPoolBenchmark(envmode));			// Two 1-particle molecules colliding
	testman.addTest(doPoolCompSolBenchmark(envmode));	// One 1-particle molecule colliding with 1 solvent
	testman.addTest(doSinglebondBenchmark(envmode));
	testman.addTest(doAnglebondBenchmark(envmode));
	testman.addTest(doDihedralbondBenchmark(envmode));
	testman.addTest(TestUtils::loadAndRunBasicSimulation("torsion2", envmode, 0.0002));
	testman.addTest(doImproperDihedralBenchmark(envmode));
	TestUtils::loadAndRunBasicSimulation("improper", envmode, 5e-5);

	// Smaller compound tests
	testman.addTest(doMethionineBenchmark(envmode));
	testman.addTest(doPhenylalanineBenchmark(envmode));
	testman.addTest(TestUtils::loadAndRunBasicSimulation("TenSolvents", envmode, 0.0004, 1.2e-6));
	testman.addTest(doEightResiduesNoSolvent(envmode));

	// Larger tests
	testman.addTest(loadAndEMAndRunBasicSimulation("SolventBenchmark", envmode, 2e-6f));
	testman.addTest(loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 0.023, 2e-5));

	// Meta tests
	//doPool50x(EnvMode::Headless);

	// Total test status will print as testman is destructed
}