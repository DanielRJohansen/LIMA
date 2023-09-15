

#include "ForceCorrectness.h"
#include "MDStability.cuh"
#include "EnergyMinimization.cuh"




using namespace TestUtils;
using namespace ForceCorrectness;
using namespace TestMDStability;
using namespace StressTesting;

void runAllUnitTests();


int main() {
	try {
		constexpr auto envmode = EnvMode::Full;

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
		//loadAndEMAndRunBasicSimulation("SolventBenchmark", envmode, 2e-6);


		//loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 8e-5, 2e-5);
		//loadAndRunBasicSimulation("T4Lysozyme", envmode, 8e-5);

		//loadAndEMAndRunBasicSimulation("manyt4", envmode, 2e-4);

		//doPool50x(EnvMode::Headless);

		runAllUnitTests();
	}
	catch (std::runtime_error ex) {
		std::cerr << "Caught exception: " << ex.what() << std::endl;
	}
	catch (const std::string& ex) {
		std::cerr << "Caught string: " << ex.c_str() << std::endl;
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
	ADD_TEST(testman, "Dihedral_exaggerated", TestUtils::loadAndRunBasicSimulation("torsion2", envmode, 2e-4, 1.4e-7));
	ADD_TEST(testman, "doImproperDihedralBenchmark", doImproperDihedralBenchmark(envmode));
	ADD_TEST(testman, "Improper_exaggerated_scaled-up", TestUtils::loadAndRunBasicSimulation("improper", envmode, 7e-5, 2.3e-7));


	// Smaller compound tests
	ADD_TEST(testman, "doMethionineBenchmark", doMethionineBenchmark(envmode));
	ADD_TEST(testman, "doPhenylalanineBenchmark", doPhenylalanineBenchmark(envmode));
	ADD_TEST(testman, "TenSolvents", TestUtils::loadAndRunBasicSimulation("TenSolvents", envmode, 5e-7, 1.2e-6));
	ADD_TEST(testman, "doEightResiduesNoSolvent", doEightResiduesNoSolvent(envmode));


	// Larger tests
	ADD_TEST(testman, "SolventBenchmark", loadAndEMAndRunBasicSimulation("SolventBenchmark", envmode, 2.2e-6f));
	ADD_TEST(testman, "T4Lysozyme", loadAndEMAndRunBasicSimulation("T4Lysozyme", envmode, 3e-5, 2e-5));


	// Meta tests
	//doPool50x(EnvMode::Headless);

	// Total test status will print as testman is destructed
}
